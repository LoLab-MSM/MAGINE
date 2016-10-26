# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pygraphviz as pyg
from bioservices import KEGG, UniProt

uniprot = UniProt()
uniprot.TIMEOUT = 100
kegg = KEGG()
kegg.TIMEOUT = 100


def add_attribute_to_network(graph, list_to_add_attribute, attribute,
                             true_term, false_term='false'):
    tmp_g = graph.copy()
    nodes1 = tmp_g.nodes()
    for i in nodes1:
        if i in list_to_add_attribute:
            tmp_g.node[i][attribute] = true_term
        else:
            tmp_g.node[i][attribute] = false_term

    return tmp_g


def append_attribute_to_network(graph, list_to_add_attribute, attribute,
                                true_term):
    tmp_g = graph.copy()
    nodes1 = tmp_g.nodes()
    for i in nodes1:
        if i in list_to_add_attribute:
            if attribute in tmp_g.node[i].keys():
                tmp_g.node[i][attribute] += true_term
            else:
                tmp_g.node[i][attribute] = true_term

    return tmp_g


def paint_network(graph, list_to_paint, color):
    """
    Paints a graph given a list of nodes and a color for that list
    :param graph: pygraphvix.AGraph
    :param list_to_paint: list
    :param color: string
    :return:
    """
    tmp_g = graph.copy()
    nodes1 = tmp_g.nodes()
    for i in list_to_paint:
        if i in nodes1:
            n = tmp_g.get_node(i)
            n.attr['measured'] = 'True'
            n.attr['color'] = 'black'
            n.attr['fillcolor'] = color
            n.attr['style'] = 'filled'
    return tmp_g


def subtract_network_from_network(net1, net2):
    copy_graph1 = net1.copy()
    nodes1 = net1.nodes()
    nodes2 = net2.nodes()
    for i in nodes2:
        if i in nodes1:
            copy_graph1.remove_node(i)
    return copy_graph1


def create_list_of_species(graph, filename, all_prot):
    output_measured = ''
    output_not_measured = ''
    list_only = ''
    f = open(filename, 'w')
    for i in graph.nodes():
        list_only += i + ','
        if np.array(all_prot[all_prot['Gene  '] == i]).shape == (1, 5):
            output_measured += str(
                    np.array(all_prot[all_prot['Gene  '] == i])[0]).strip(
                    '[').strip(']') + '\n'
        else:
            output_not_measured += i.replace(" ", "_") + ' nan nan nan nan\n'
    f.write(output_measured.replace("'", ""))
    f.write(output_not_measured)
    f.close()
    f = open('list_' + filename, 'w')
    f.write(list_only)
    f.close()


def return_gml(graph):
    newGraph = nx.DiGraph()
    for node in graph.nodes():
        newGraph.add_node(str(node))
    for each in graph.edges():
        newGraph.add_edge(each[0], each[1], attr_dict=each.attr)
    return newGraph


def compress_edges(graph, networkx=True):
    if not networkx:
        graph = nx.from_agraph(graph)
    g = graph.copy()
    nodes = g.nodes()
    neighbor_dict = {}
    for n in nodes:
        neigh = tuple(g.neighbors(n))
        if tuple(neigh) in neighbor_dict:
            neighbor_dict[tuple(neigh)].append(n)
        else:
            neighbor_dict[tuple(neigh)] = []
            neighbor_dict[tuple(neigh)].append(n)

    for i in neighbor_dict:
        neigh = neighbor_dict[i]
        # print i,neigh
        if len(i) != 1:
            continue
        if len(neigh) == 1:
            continue
        interaction_types = {}
        for j in neighbor_dict[i]:
            if g.has_edge(i[0], j):
                direction = 'type1'
                edge = g.get_edge(i[0], j)
            elif g.has_edge(j, i[0]):
                direction = 'type2'
                edge = g.get_edge(j, i[0])
            edge_type = edge.attr['arrowhead']
            if edge_type + direction in interaction_types:
                interaction_types[edge_type + direction].append(j)
            else:
                interaction_types[edge_type + direction] = []
                interaction_types[edge_type + direction].append(j)
        for k in interaction_types:
            if len(interaction_types[k]) > 1:
                # print i[0],'->',
                to_join = []
                for each in interaction_types[k]:
                    if len(g.neighbors(each)) == 1:
                        to_join.append(each)
                # print to_join,k
                if len(to_join) > 1:
                    label = "{"
                    for node in to_join:
                        g.remove_node(node)
                        label += ' %s |' % str(node)
                    label = label[:-1]
                    label += "}"
                    # print label
                    g.add_node(label, shape='record', label=label)
                    if k.endswith('type2'):
                        g.add_edge(label, i[0], dir='both', arrowhead=k[:-5],
                                   arrowtail="none")
                    elif k.endswith('type1'):
                        g.add_edge(i[0], label, dir='both', arrowhead=k[:-5],
                                   arrowtail="none")
    return g


def test_compress():
    g = pyg.AGraph()
    g.add_node('A')
    g.add_node('B')
    g.add_node('C')
    g.add_node('D')
    g.add_node('E')
    g.add_edge('A', 'D', dir='both', arrowhead='dot', arrowtail="none", )
    g.add_edge('B', 'D', dir='both', arrowhead='dot', arrowtail="none", )
    g.add_edge('C', 'D', dir='both', arrowhead='dot', arrowtail="none", )
    g.add_edge('E', 'A', dir='both', arrowhead='dot', arrowtail="none", )
    g.add_edge('E', 'B', dir='both', arrowhead='dot', arrowtail="none", )
    g.add_edge('E', 'C', dir='both', arrowhead='odot', arrowtail="none", )
    # g.draw('before.pdf',prog='dot')
    # g = compress_edges(g)
    # g.draw('after.pdf',prog='dot')
    g = nx.from_agraph(g, create_using=nx.DiGraph())
    nodelist = ['A', 'B', 'C', 'D', 'E']
    adj = nx.to_numpy_recarray(g, nodelist=['A', 'B', 'C', 'D', 'E'],
                               dtype=[('arrowhead', 'S50')])

    for n, i in enumerate(range((np.shape(adj)[1]))):
        to_concate = []
        for k, j in enumerate(range(n, (np.shape(adj)[1]))):
            if i == j:
                continue
            if (adj[:, i] == adj[:, j]).all():
                to_concate.append(j)
                print "here", nodelist[i], nodelist[j]
        print to_concate
    print adj
    test = nx.from_numpy_matrix(adj, create_using=nx.DiGraph())
    nx.draw_graphviz(test, prog='dot')
    plt.savefig('test.png')


def get_uniprot_info(name):
    """

    columns are search terms from UniProt, can be choosen from
    http://www.uniprot.org/help/uniprotkb_column_names

    """

    d = uniprot.search('%s+AND+organism:9606+reviewed:yes' % name, \
                       frmt='tab', \
                       columns="entry name,\
             genes(PREFERRED),\
             comment(FUNCTION),\
             go(biological process),\
             go(molecular function),\
             go(cellular component),\
             comment(PTM)"
                       , limit=1)
    output_line = ''
    try:
        d.split('\n')
    except:
        return ''
    for i in d.split('\n'):
        if i.startswith('Entry'):
            continue
        elif i.startswith('\n'):
            continue
        else:
            output_line += i
    output_line = output_line.replace("; ", ";")
    return output_line + '\n'


def generate_curated_subgraphs(network, i):
    header = 'NetworkName\tEntry name\tGene_names_primary\tFunction\tGO_biological_process\tGO_molecular_function\tGO_cellular_component\tPost-translational_modification'
    print "generating curated list for %s" % i
    size = len(network.nodes())
    file_to_write = open(
            'List_of_species_subgroup_%s_size_%s.txt' % (i + 1, size), 'w')
    output = header + '\n'
    for j in network.nodes():
        tmp = get_uniprot_info(j)
        output += str(j) + '\t' + tmp
    file_to_write.write(output)
    file_to_write.close()
    print "Created curated list for %s", i


def trim_nodes(network, list_of_nodes):
    tmp1 = trim(network, list_of_nodes)
    tmp2 = trim(network, list_of_nodes)
    while tmp1 != tmp2:
        tmp2 = tmp1
        tmp1 = trim(network, list_of_nodes)


def trim(network, list_of_nodes):
    list_of_nodes = set(list_of_nodes)
    found, not_found = 0., 0.
    nodes = set(network.nodes())
    in_dict = network.in_degree(nbunch=nodes)
    out_dict = network.out_degree(nbunch=nodes)
    for i in nodes:
        if i in list_of_nodes:
            found += 1
        else:
            if in_dict[i] == 1 and out_dict[i] == 0:
                network.remove_node(i)
            elif in_dict[i] == 0 and out_dict[i] == 1:
                network.remove_node(i)
            else:
                not_found += 1

    len_nodes = len(network.nodes())
    print("{} found, {} not found".format(found, not_found))
    print("{}% of {} nodes".format(found / len_nodes * 100, len_nodes))
    return len_nodes


def create_lists_of_subgraphs(network, save_name, exp_data):
    G = network.to_undirected()
    sorted_graphs = sorted(nx.connected_component_subgraphs(G), key=len,
                           reverse=True)
    counter = 0
    data = []
    cnt = 0
    subgraph_species = []
    for i in sorted_graphs:

        size = len(i.nodes())
        with open('%s_%s_size_%s.txt' % (save_name, str(cnt), str(size)),
                  'w') as f:
            for j in i.nodes():
                f.write('%s,' % j)

        cnt += 1
        if size == 1:
            counter += 1
        else:
            data.append(size)
        if size == 1:
            continue
        measured = ''
        measured_species = []
        for j in i.nodes():
            if j in exp_data:
                measured_species.append(j)
                measured += '%s,' % str(j)
        if measured == '':
            continue
        else:
            subgraph_species.append(measured_species)
            print(measured)
    data.remove(max(data))
    data = np.asarray(data)
    print(data, 10)
    plt.hist(data)
    plt.title("Distribution of subgraphs with canonical removed")
    # plt.xlim(1,100)
    # plt.ylim(0,20)
    plt.xlabel("Number of nodes")
    plt.ylabel("Count")
    plt.savefig("histogram_mega_minus_canonical_subgraphs.png", dpi=200)
    plt.show()
    print "Number of subgraphs with 1 node = %s" % counter
    return subgraph_species


def add_reactome(network):
    combined_network = network.copy()
    existing_nodes = set(network.nodes())
    ddn2 = nx.read_gml(
            '/home/pinojc/PycharmProjects/Magine/magine/reactome_expansion/reactome_network.gml')
    added_nodes = 0
    for i in ddn2.nodes(data=True):
        if i[0] in existing_nodes:
            continue
        else:
            network.add_node(i[0], i[1])
            added_nodes += 1
    print("Added {} nodes from REACTOME".format(added_nodes))

    added_edges = 0
    kegg_edges = set(combined_network.edges())
    reactome_edges = ddn2.edges(data=True)
    for i, j, k in reactome_edges:
        if (i, j) in kegg_edges:
            continue
        else:
            combined_network.add_edge(i, j, k)
            added_edges += 1
    print("Added {} edges from REACTOME".format(added_edges))
    return combined_network
if __name__ == '__main__':
    test_compress()
