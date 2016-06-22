import itertools
from textwrap import wrap

import networkx as nx
import numpy as np
import pygraphviz as pyg
from orangecontrib.bio import go


ontology = go.Ontology()
annotations = go.Annotations('hsa', ontology=ontology)


def calculate_terms_between_a_b(nodes, ddn, term_a, term_b):
    path_a_to_b = 0
    path_b_to_a = 0
    in_term_1 = 0
    in_term_2 = 0
    nodes = set(nodes)
    term_a = set(term_a)
    term_b = set(term_b)

    for i in term_a:
        if i in nodes:
            in_term_1 += 1
            neigh = set(ddn.neighbors(i))
            for j in neigh:
                if j in term_b:
                    if ddn.has_edge(i, j):
                        path_a_to_b += 1

    for i in term_b:
        if i in nodes:
            in_term_2 += 1
            neigh = set(ddn.neighbors(i))
            for j in neigh:
                if j in term_a:
                    if ddn.has_edge(i, j):
                        path_b_to_a += 1

    return path_a_to_b, path_b_to_a

def create_go_network(ddn, term1, term2, term3, save_name):

    nodes = ddn.nodes()
    term_1 = annotations.get_all_genes(term1)
    term_2 = annotations.get_all_genes(term2)
    term_3 = annotations.get_all_genes(term3)
    label_1 = ontology[term1].name
    label_2 = ontology[term2].name
    label_3 = ontology[term3].name

    # create_venn3(term_1, term_2, term_3, label_1, label_2, label_3,
    #             '%s_venn_diagram' % save_name)

    g = pyg.AGraph(directed=True, rankdir='LR')

    calculate_terms_between_a_b(nodes, ddn, g, term_1, term_2, label_1, label_2)
    calculate_terms_between_a_b(nodes, ddn, g, term_1, term_3, label_1, label_3)
    calculate_terms_between_a_b(nodes, ddn, g, term_2, term_3, label_2, label_3)
    #g.draw('%s.png' % save_name, prog='dot')
    return g


def create_network_from_list(network, list_of_go_terms, save_name):
    list_of_go_terms = np.array(np.unique(list_of_go_terms))
    graph = nx.DiGraph()
    nodes = network.nodes()
    gene_annotations_dict = dict()
    go_names = dict()
    for i in list_of_go_terms:
        gene_annotations_dict[i] = annotations.get_all_genes(i)
        go_names[i] = "\n".join(wrap(ontology[i].name, 20))
    for n, i in enumerate(itertools.combinations(list_of_go_terms, 2)):
        term1 = i[0]
        term2 = i[1]
        term_1 = gene_annotations_dict[term1]
        term_2 = gene_annotations_dict[term2]
        label_1 = go_names[term1]
        label_2 = go_names[term2]
        a_to_b, b_to_a = calculate_terms_between_a_b(nodes, network, term_1, term_2)
        if (a_to_b and b_to_a) != 0:
            print('both', a_to_b, b_to_a)
            graph.add_node(label_1, go=str(term1).replace(':', ''))
            graph.add_node(label_2, go=str(term2).replace(':', ''))
            graph.add_edge(label_1, label_2, label=str(a_to_b) + '/' + str(b_to_a), dir='both')
        elif a_to_b != 0:
            graph.add_node(label_1, go=str(term1).replace(':', ''))
            graph.add_node(label_2, go=str(term2).replace(':', ''))
            graph.add_edge(label_1, label_2, label=str(a_to_b))

        elif b_to_a != 0:
            graph.add_node(label_1, go=str(term1).replace(':', ''))
            graph.add_node(label_2, go=str(term2).replace(':', ''))
            graph.add_edge(label_2, label_1, label=str(b_to_a))

    nx.nx.write_dot(graph, '{0}.dot'.format(save_name))

    g = pyg.AGraph()
    g.read('{0}.dot'.format(save_name))
    g.draw('{0}.pdf'.format(save_name), prog='dot')
    return graph


if __name__ == '__main__':
    ddn = nx.read_gml('/home/pinojc/git/Network_projects/Cisplatin_project/Network_files/ddn3.gml')
    # test_list = ['GO:0008219', 'GO:0006281', 'GO:0008283', 'GO:0043066', 'GO:0043065', 'GO:1902175','GO:0006805', 'GO:0006766']
    test_list = ['GO:1902175', 'GO:0006805', 'GO:0006766', 'GO:0015893', 'GO:0006936']
    create_network_from_list(ddn, test_list, 'xeno')
    # create_go_network,save_name='death_dnaRepair_proliferation')

    # create_go_network('GO:0043066', 'GO:0043065', 'GO:1902175', 'apoptosis')
