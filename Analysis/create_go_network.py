import networkx as nx
import pygraphviz as pyg
from orangecontrib.bio import go

from Analysis.venn_diagram_maker import create_venn3

ontology = go.Ontology()
annotations = go.Annotations('hsa', ontology=ontology)


def calculate_edges(nodes, ddn, graph, term1, term2, t1, t2):
    path_1_2 = 0
    path_1_1 = 0
    path_2_1 = 0
    path_2_2 = 0
    for i in nodes:
        if i in term1:
            neigh = ddn.neighbors(i)
            for j in neigh:
                if j in term2:
                    if nx.has_path(ddn, i, j):
                        path_1_2 += 1
                if j in term1:
                    if nx.has_path(ddn, i, j):
                        path_1_1 += 1
                    if nx.has_path(ddn, j, i):
                        path_1_1 += 1

        if i in term2:
            neigh = ddn.neighbors(i)
            for j in neigh:
                if j in term1:
                    if nx.has_path(ddn, i, j):
                        path_2_1 += 1
                if j in term2:
                    if nx.has_path(ddn, i, j):
                        path_2_2 += 1
                    if nx.has_path(ddn, j, i):
                        path_2_2 += 1
    graph.add_edge(t1, t2, label=path_1_2)
    graph.add_edge(t2, t1, label=path_2_1)
    graph.add_edge(t1, t1, label=path_1_1)
    graph.add_edge(t2, t2, label=path_2_2)


# cell death GO:0008219
# DNA repair GO:0006281
# proliferation GO:0008283

def create_go_network(term1, term2, term3, save_name):
    ddn = nx.read_gml('Network_files/ddn3.gml')

    nodes = ddn.nodes()
    term_1 = annotations.get_all_genes(term1)
    term_2 = annotations.get_all_genes(term2)
    term_3 = annotations.get_all_genes(term3)
    label_1 = ontology[term1].name
    label_2 = ontology[term2].name
    label_3 = ontology[term3].name
    ddn = nx.read_gml('Network_files/ddn3.gml')

    create_venn3(term_1, term_2, term_3, 'Cell Death', 'DNA Repair', 'Proliferation',
                 '%s_venn_diagram' % save_name)

    g = pyg.AGraph(directed=True, rankdir='LR')

    calculate_edges(nodes, ddn, g, term_1, term_2, label_1, label_2)
    calculate_edges(nodes, ddn, g, term_1, term_3, label_1, label_3)
    calculate_edges(nodes, ddn, g, term_2, term_3, label_2, label_3)
    g.draw('%s.png' % save_name, prog='dot')


# create_go_network('GO:0008219','GO:0006281','GO:0008283',save_name='death_dnaRepair_proliferation')

create_go_network('GO:0043066', 'GO:0043065', 'GO:1902175', 'apoptosis')
