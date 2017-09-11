# -*- coding: utf-8 -*-
import itertools
import os

import networkx as nx
import pathos.multiprocessing as mp

try:
    import igraph as ig

    NO_IGRAPH = False
except ImportError:
    NO_IGRAPH = True


def delete_disconnected_network(full_graph, verbose=False):
    """
    Delete disconnected parts of a provided network

    Parameters
    ----------
    full_graph : networkx.DiGraph
    verbose : bool
        Prints all nodes that are removed

    """

    tmp_g = full_graph.to_undirected()
    sorted_graphs = sorted(nx.connected_component_subgraphs(tmp_g), key=len,
                           reverse=True)
    for i in range(1, len(sorted_graphs)):
        for node in sorted_graphs[i].nodes():
            if verbose:
                print("Removing disconnected node %s", node)
            full_graph.remove_node(node)


def paint_network_overtime(graph, list_of_lists, color_list, save_name,
                           labels=None):
    """
    Adds color attribute to network over time.
    
    Parameters
    ----------
    graph : pygraphviz.AGraph
        Network
    list_of_lists : list_like
        List of lists, where the inner list contains the node to add the
        color 
    color_list : list_like
        list of colors for each time point
    save_name : str
        prefix for images to be saved
    labels: list
        list of labels to add to graph per sample
    
    Returns
    -------

    """
    if len(list_of_lists) != len(color_list):
        print('Length of list of data must equal len of color list')
        return
    if labels is not None:
        if len(labels) != len(list_of_lists):
            print('Length of labels must be equal to len of data')
            return
    string = 'convert -delay 100 '
    tmp_graph = graph.copy()
    tmp_graph = _format_to_directions(tmp_graph)
    if isinstance(tmp_graph, nx.DiGraph):
        try:
            tmp_graph = nx.nx_agraph.to_agraph(tmp_graph)
        except ImportError:
            print("Please install pygraphviz")
            return
    for n, i in enumerate(list_of_lists):
        graph2 = paint_network(tmp_graph, i, color_list[n])

        if labels is not None:
            graph2.graph_attr.update(label=labels[n], ranksep='0.2',
                                     fontsize=13)
        graph2.draw('%s_%04i.png' % (save_name, n), prog='dot')
        string += '%s_%04i.png ' % (save_name, n)
    string1 = string + '  %s.gif' % save_name
    string2 = string + '  %s.pdf' % save_name

    os.system(string1)
    os.system(string2)


def paint_network_overtime_up_down(graph, list_up, list_down, save_name,
                                   color_up='red', color_down='blue',
                                   labels=None):
    """
    Adds color attribute to network over time.

    Parameters
    ----------
    graph : pygraphviz.AGraph
        Network
    list_up : list_like
        List of lists, where the inner list contains the node to add the
        color 
    list_down : list_like
        list of colors for each time point
    color_up : str
        color for first list of species
    color_down : str
        color of second list of species
    save_name : str
        prefix for images to be saved
    labels: list
        list of labels to add to graph per sample

    Returns
    -------

    """
    if len(list_up) != len(list_down):
        print('Length of list of data must equal len of color list')
        return
    if labels is not None:
        if len(labels) != len(list_down):
            print('Length of labels must be equal to len of data')
            return
    string = 'convert -delay 100 '
    tmp_graph = graph.copy()
    tmp_graph = _format_to_directions(tmp_graph)

    if isinstance(tmp_graph, nx.DiGraph):
        try:
            tmp_graph = nx.nx_agraph.to_agraph(tmp_graph)
        except ImportError:
            print("Please install pygraphviz")
            return
    for n, (up, down) in enumerate(zip(list_up, list_down)):
        graph2 = paint_network(tmp_graph, up, color_up)
        graph2 = paint_network(graph2, down, color_down)
        up_s = set(up)
        down_s = set(down)
        both = up_s.intersection(down_s)
        graph2 = paint_network(graph2, both, 'yellow')

        if labels is not None:
            graph2.graph_attr.update(label=labels[n], ranksep='0.2',
                                     fontsize=13)
        graph2.draw('%s_%04i.png' % (save_name, n), prog='dot')
        string += '%s_%04i.png ' % (save_name, n)
    string1 = string + '  %s.gif' % save_name
    string2 = string + '  %s.pdf' % save_name

    os.system(string1)
    os.system(string2)


def paint_network(graph, list_to_paint, color):
    """
    
    Parameters
    ----------
    graph: pygraphvix.AGraph
    list_to_paint : list
        
    color : str

    Returns
    -------

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


def _format_to_directions(network):
    if isinstance(network, nx.DiGraph):
        try:
            network = nx.nx_agraph.to_agraph(network)
        except ImportError:
            print("Need to install pygraphviz")
            return network
    activators = ['activation', 'expression', 'phosphorylation']
    inhibitors = ['inhibition', 'repression', 'dephosphorylation',
                  'ubiquitination']
    physical_contact = ['binding/association', 'dissociation', 'state change',
                        'compound', 'glycosylation']
    indirect_types = ['missing interaction', 'indirect effect']
    for i in network.edges():
        n = network.get_edge(i[0], i[1])
        edge_type = str(i.attr['interactionType'])

        def _find_edge_type():
            for j in activators:
                if edge_type.startswith(j):
                    n.attr['arrowhead'] = 'normal'
                    return
            for j in inhibitors:
                if edge_type.startswith(j):
                    n.attr['arrowhead'] = 'tee'
                    return
            for j in physical_contact:
                if edge_type.startswith(j):
                    n.attr['dir'] = 'both'
                    n.attr['arrowtail'] = 'diamond'
                    n.attr['arrowhead'] = 'diamond'
                    return
            for j in indirect_types:
                if edge_type.startswith(j):
                    n.attr['arrowhead'] = 'diamond'
                    n.attr['style'] = 'dashed'
                    return
                    # print(n, n.attr)

        _find_edge_type()
    network = nx.nx_agraph.from_agraph(network)
    return network


def create_legend(graph):
    """
    adds a legend to a graph
    :param graph:
    :return:
    """
    dict_of_types = {
        'activation': 'onormal',
        'indirect effect': 'odiamondodiamond',
        'expression': 'normal',
        'inhibition': 'tee',
        'binding/association': 'curve',
        'phosphorylation': 'dot',
        'missing interaction': 'odiamond',
        'compound': 'dotodot',
        'dissociation': 'diamond',
        'ubiquitination': 'oldiamond',
        'state change': 'teetee',
        'dephosphorylation': 'onormal',
        'repression': 'obox',
        'glycosylation': 'dot'
    }
    subgraph = []
    len_dic = len(dict_of_types)
    for n, i in enumerate(dict_of_types):
        subgraph.append(n)
        subgraph.append(n + len_dic)
        graph.add_node(n, label="")
        graph.add_node(n + len_dic, label="")
        graph.add_edge(n, n + len_dic, dir='both', arrowhead=dict_of_types[i],
                       arrowtail="none", label=i)
    graph.add_subgraph(subgraph, name='cluster_legend', rank="LR")


def export_to_dot(graph, save_name, image_format='png', engine='dot',
                  dpi=200, concentrate=False):
    """
    Converts networkx graph to graphviz dot

    Parameters
    ----------
    graph : networkx.DiGraph
    save_name : str
        name of file to save
    image_format : str
        format of output( pdf, png, svg)
    engine : str
        graphviz engine
            dot, twopi,
    dpi: int
        resolution of figure

    Returns
    -------

    """
    try:
        py_dot = nx.nx_agraph.to_agraph(graph)
        py_dot.write('{}.dot'.format(save_name))
        arg = '-Gdpi={} -Gconcentrate={}'.format(
            dpi, 'true' if concentrate else 'false'
        )
        py_dot.draw('{}.{}'.format(save_name, image_format), prog=engine,
                    args=arg)
    except ImportError:
        print("No pygraphivz installed")


def add_attribute_to_network(graph, list_to_add_attribute, attribute,
                             true_term, false_term='false'):
    """

    Parameters
    ----------
    graph : networkx graph
    list_to_add_attribute : list
        list of nodes in graph to add attribute to
    attribute : str
        attribute to add to graph
    true_term : str
        value to add for the attribute provided True
    false_term : str
        value to add if attribute is false

    Returns
    -------
    out : networkx graph


    >>> from networkx import DiGraph
    >>> from magine.networks.network_tools import add_attribute_to_network
    >>> g = DiGraph()
    >>> g.add_nodes_from(['a', 'b', 'c'])
    >>> new_g = add_attribute_to_network(g, ['a','b'], 'inTest', 'true')
    >>> new_g.node['a']
    {'inTest': 'true'}
    >>> new_g.node['c']
    {'inTest': 'false'}

    """
    tmp_g = graph.copy()
    nodes = set(tmp_g.nodes())
    set_of_positive = set(list_to_add_attribute)
    for i in nodes:
        if i in set_of_positive:
            tmp_g.node[i][attribute] = true_term
        else:
            tmp_g.node[i][attribute] = false_term

    return tmp_g


def append_attribute_to_network(graph, list_to_add_attribute, attribute,
                                true_term, delimiter=''):
    """

    Parameters
    ----------
    graph : networkx graph
    list_to_add_attribute : list
        list of nodes in graph to add attribute to
    attribute : str
        attribute to add to graph
    true_term : str
        value to add for the attribute provided True

    Returns
    -------
    out : networkx graph


    >>> from networkx import DiGraph
    >>> from magine.networks.network_tools import append_attribute_to_network
    >>> g = DiGraph()
    >>> g.add_nodes_from(['a'])
    >>> new_g = append_attribute_to_network(g, ['a'], 'attribute', 'one')
    >>> new_g.node['a']
    {'attribute': 'one'}
    >>> new_g = append_attribute_to_network(new_g, ['a'], 'attribute', 'two', delimiter=',')
    >>> new_g.node['a']
    {'attribute': 'one,two'}

    """
    tmp_g = graph.copy()
    nodes = set(tmp_g.nodes())
    set_of_positive = set(list_to_add_attribute)
    for i in nodes:
        if i in set_of_positive:
            if attribute in tmp_g.node[i].keys():
                new_attr = tmp_g.node[i][attribute] + delimiter + true_term
                tmp_g.node[i][attribute] = new_attr
            else:
                tmp_g.node[i][attribute] = true_term

    return tmp_g


def trim_sink_source_nodes(network, list_of_nodes):
    """
    Trim graph by removing nodes that are not in provided list if source/sink


    Parameters
    ----------
    network : networkx.DiGraph

    list_of_nodes : list_like
        list of species that are important if sink/source

    Returns
    -------

    """
    tmp_network = network.copy()
    edges = set(tmp_network.edges())
    for i, j in edges:
        if i == j:
            tmp_network.remove_edge(i, j)
            print("removed {} {}".format(i, j))
    tmp1 = _trim(tmp_network, list_of_nodes)
    tmp2 = _trim(tmp_network, list_of_nodes)
    while tmp1 != tmp2:
        tmp2 = tmp1
        tmp1 = _trim(tmp_network, list_of_nodes)
    return tmp_network


def _trim(network, list_of_nodes):
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


def networkx_to_igraph(network):
    if NO_IGRAPH:
        print("igraph not installed")
        return False
    igraph_network = ig.Graph(directed=True)
    for i, data in network.nodes(data=True):
        igraph_network.add_vertex(name=i, **data)
    for edge in network.edges(data=True):
        e1 = edge[0]
        e2 = edge[1]
        data = edge[2]
        igraph_network.add_edge(e1, e2, kwds=data)
    return igraph_network


def subtract_network_from_network(net1, net2):
    """
    subtract one network from another


    Parameters
    ----------
    net1 : networkx.DiGraph
    net2 : networkx.DiGraph


    Returns
    -------
    networkx.DiGraph

    """
    copy_graph1 = net1.copy()
    nodes1 = net1.nodes()
    nodes2 = net2.nodes()
    for i in nodes2:
        if i in nodes1:
            copy_graph1.remove_node(i)
    return copy_graph1


def compress_edges(graph):
    """
    compress edges of networks by finding common paths

    Parameters
    ----------
    graph: networkx.DiGraph

    Returns
    -------

    """

    g = graph.copy()
    nodes = set(g.nodes())
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
        # print(i,neigh)
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
                # print(i[0],'->',)
                to_join = []
                for each in interaction_types[k]:
                    if len(g.neighbors(each)) == 1:
                        to_join.append(each)
                # print(to_join,k)
                if len(to_join) > 1:
                    label = "{"
                    for node in to_join:
                        g.remove_node(node)
                        label += ' %s |' % str(node)
                    label = label[:-1]
                    label += "}"
                    # print(label)
                    g.add_node(label, shape='record', label=label)
                    if k.endswith('type2'):
                        g.add_edge(label, i[0], dir='both', arrowhead=k[:-5],
                                   arrowtail="none")
                    elif k.endswith('type1'):
                        g.add_edge(i[0], label, dir='both', arrowhead=k[:-5],
                                   arrowtail="none")
    return g


def merge_nodes(graph):
    """ merges nodes into single node if same neighbors

    Parameters
    ----------
    graph: nx.DiGraph

    Returns
    -------
    nx.DiGraph
    """
    neighbors2node = dict()

    nodes = set(graph.nodes())
    edges = set(graph.edges())

    new_g = graph.copy()
    for i in nodes:
        down = '|'.join(sorted(graph.successors(i)))
        up = '|'.join(sorted(graph.predecessors(i)))

        if (down, up) in neighbors2node:
            neighbors2node[(down, up)].add(i)
        else:
            neighbors2node[(down, up)] = {i}

    for node, neigh in neighbors2node.items():
        # print(node, neigh)
        if len(neigh) > 1:
            print(node, neigh)
            new_name = '|'.join(sorted(neigh))
            for x in node[0].split('|'):
                if x != '':
                    for n in neigh:
                        if (x, n) in edges:
                            new_g.remove_edge(x, n)
                            new_g.add_edge(x, new_name, **graph.edge[x][n])
                            edges.remove((x, n))

            for x in node[1].split('|'):
                if x != '':
                    for n in neigh:
                        if (n, x) in edges:
                            new_g.remove_edge(n, x)
                            new_g.add_edge(new_name, x, **graph.edge[n][x])
                            edges.remove((n, x))
    print(
    "{} nodes and {} edges".format(len(graph.nodes()), len(graph.edges())))
    print(
    "{} nodes and {} edges".format(len(new_g.nodes()), len(new_g.edges())))
    for n in new_g.nodes():
        if len(new_g.successors(n)) == 0 and len(new_g.predecessors(n)) == 0:
            new_g.remove_node(n)

    print(
    "{} nodes and {} edges".format(len(new_g.nodes()), len(new_g.edges())))
    return new_g


def remove_unmeasured_nodes(graph, measured):
    new_g = graph.copy()
    edge_info_dict = dict()
    for i, j, data in graph.edges(data=True):
        edge_info_dict[(i, j)] = data
    nodes = set(graph.nodes())
    include = set(measured)
    include.intersection_update(nodes)

    # for i, j in :
    def find(d):
        i, j = d
        paths = []
        if nx.has_path(graph, i, j):
            for p in nx.all_shortest_paths(graph, i, j):
                path = []
                label = []
                for n in p:
                    if n in include:
                        path.append(n)
                    else:
                        label.append(n)
                if len(path) == 2:
                    paths.append((path, '|'.join(l for l in label)))
        return paths

    x = mp.Pool(4)
    paths = x.map(find, itertools.combinations(include, 2))
    to_remove = set()
    for p in paths:
        if len(p) != 2:
            continue
        for path, label in p:
            print(path, label)
            for n in label.split('|'):
                to_remove.add(n)
            new_g.add_edge(path[0], path[1], label=label)
    for n in to_remove:
        new_g.remove_node(n)
    return new_g


def _nx_to_dot(network):
    if isinstance(network, nx.DiGraph):
        try:
            network = nx.nx_agraph.to_agraph(network)
        except ImportError:
            print("Need to install pygraphviz")
            return network


def print_network_stats(network, exp_data):
    """

    Parameters
    ----------
    network : nx.DiGraph
    exp_data : magine.data.datatypes.ExperimentalData

    Returns
    -------

    """
    nodes = set(network.nodes())
    g_nodes = set()
    m_nodes = set()
    for i, d in network.nodes(data=True):

        if 'speciesType' not in d:
            print("Doesn't have species type... {} = {}".format(i, d))
            continue
        if d['speciesType'] == 'gene':
            g_nodes.add(i)
        elif d['speciesType'] == 'compound' or d[
            'speciesType'] == 'metabolite':
            m_nodes.add(i)
        else:
            g_nodes.add(i)
            print(i, d)

    total_species = set(exp_data.list_species)
    sig_total_species = set(exp_data.list_sig_species)
    sig_genes = set(exp_data.list_sig_proteins)
    genes = set(exp_data.list_proteins)
    sig_met = set(exp_data.list_sig_metabolites)
    met = set(exp_data.list_metabolites)
    n_measured = len(nodes.intersection(total_species))
    n_genes = len(nodes.intersection(sig_genes))
    n_meta = len(nodes.intersection(sig_met))
    n_sig_measured = len(nodes.intersection(sig_total_species))
    n_nodes = len(nodes)
    st = ''
    st += 'Number of nodes = {}\n'.format(n_nodes)
    st += 'Number of edges = {}\n'.format(len(network.edges()))
    st += "Number of total nodes measured = {}\n".format(n_measured)
    st += "Fraction of total nodes measured = {}\n".format(
        100. * n_measured / n_nodes)
    st += "Number of total nodes sig. changed = {}\n".format(n_sig_measured)
    st += "Fraction of total nodes sig. changed = {}\n".format(
        100. * n_sig_measured / n_nodes)
    st += "Number of protein nodes sig. changed = {}\n".format(n_genes)
    st += "Fraction of protein nodes measured = {}\n".format(
        100. * n_genes / len(g_nodes))
    st += "Number of metabolite nodes sig. changed = {}\n".format(n_meta)
    st += "Fraction of total nodes measured = {}\n".format(
        100. * n_meta / len(m_nodes))

    print(st)


def add_pvalue_and_fold_change():
    # -log10pvalue,
    # # log2fold_change,
    return


def compose(G, H, name=None):
    """Return a new graph of G composed with H.

    Composition is the simple union of the node sets and edge sets.
    The node sets of G and H do not need to be disjoint.

    Parameters
    ----------
    G,H : nx.DiGraph
       A NetworkX graph

    name : str
       Specify name for new graph

    Returns
    -------
    C: A new graph  with the same type as G

    """

    if name is None:
        name = "compose( %s, %s )" % (G.name, H.name)
    new_g = nx.DiGraph()
    new_g.name = name
    new_g.add_nodes_from(H.nodes())
    new_g.add_nodes_from(G.nodes())
    for i, j, data in G.edges_iter(data=True):
        new_g.add_edge(i, j, **data)
        # print(i, j, data)
    edges = set(new_g.edges())
    for i, j, data in H.edges_iter(data=True):
        if (i, j) in edges:
            existing_info = new_g.edge[i][j]
            for n, d in data.items():
                if n in existing_info:
                    if isinstance(existing_info[n], float):
                        print(n, d, i, j)
                    current = set(existing_info[n].split('|'))
                    additions = set(d.split('|'))
                    additions.update(current)
                    new = '|'.join(sorted(additions))
                    new_g[i][j][n] = new
                else:
                    new_g[i][j][n] = d
        else:
            new_g.add_edge(i, j, **data)
    nodes = set(new_g.nodes())
    for i, data in G.nodes_iter(data=True):
        if i in nodes:
            existing_info = new_g.node[i]
            for n, d in data.items():
                if n in existing_info:
                    current = set(existing_info[n].split('|'))
                    additions = set(d.split('|'))
                    additions.update(current)
                    new = '|'.join(sorted(additions))
                    new_g.node[i][n] = new
                else:
                    new_g.node[i][n] = d
        else:
            new_g.add_node(i, **data)
    nodes2 = set(new_g.nodes())
    for i, data in H.nodes_iter(data=True):
        if i in nodes2:
            existing_info = new_g.node[i]
            for n, d in data.items():
                if n in existing_info:
                    current = existing_info[n]
                    new = '|'.join({str(current), str(d)})
                    new_g.node[i][n] = new
                else:
                    new_g.node[i][n] = d
        else:
            new_g.add_node(i, **data)
    return new_g


def compose_all(graphs, name=None):
    """Return the composition of all graphs.

    Composition is the simple union of the node sets and edge sets.
    The node sets of the supplied graphs need not be disjoint.

    Parameters
    ----------
    graphs : list
       List of NetworkX graphs

    name : str
       Specify name for new graph

    Returns
    -------
    C : A graph with the same type as the first graph in list

    """
    graphs = iter(graphs)
    C = next(graphs)
    for H in graphs:
        C = compose(C, H, name=name)
    return C


'''
# deprecated
def get_uniprot_info(name):
    """

    columns are search terms from UniProt, can be choosen from
    http://www.uniprot.org/help/uniprotkb_column_names

    """

    d = uniprot.search('%s+AND+organism:9606+reviewed:yes' % name, frmt='tab',
                       columns="entry name,\
                                genes(PREFERRED),\
                                comment(FUNCTION),\
                                go(biological process),\
                                go(molecular function),\
                                go(cellular component),\
                                comment(PTM)",
                       limit=1)

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


# deprecated
def generate_curated_subgraphs(network, i):
    header = 'NetworkName\tEntry name\tGene_names_primary\tFunction\tGO_biological_process\tGO_molecular_function\tGO_cellular_component\tPost-translational_modification'
    print("generating curated list for %s" % i)
    size = len(network.nodes())
    file_to_write = open(
        'List_of_species_subgroup_%s_size_%s.txt' % (i + 1, size), 'w')
    output = header + '\n'
    for j in network.nodes():
        tmp = get_uniprot_info(j)
        output += str(j) + '\t' + tmp
    file_to_write.write(output)
    file_to_write.close()
    print("Created curated list for %s", i)


# deprecated
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
    print("Number of subgraphs with 1 node = %s" % counter)
    return subgraph_species
'''

if __name__ == '__main__':
    g = nx.DiGraph()
    g.add_edge('A', 'B')
    g.add_edge('B', 'C')
    g.add_edge('B', 'E')
    g.add_edge('C', 'D')
    test_g = remove_unmeasured_nodes(g, ['A', 'D', 'E'])
    export_to_dot(test_g, 'merged_node')
    # new_g = remove_unmeasured_nodes(g, ['A', 'C', 'D'])
    # export_to_dot(new_g, 'merged_node2')
