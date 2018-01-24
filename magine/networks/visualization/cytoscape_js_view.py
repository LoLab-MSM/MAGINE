from magine.networks.ontology_network import OntologyNetworkGenerator
from magine.networks.visualization.cytoscapejs_tools import viewer as cyjs
from magine.networks.visualization.util_networkx import from_networkx
import magine.ontology.enrichment_tools as et
from magine.networks.network_subgraphs import NetworkSubgraphs
from magine.networks.visualization.igraph_tools import create_igraph_figure, create_figure
import numpy as np
from IPython.display import SVG, display
import matplotlib.pyplot as plt


def create_subnetwork(terms, df, network, save_name=None, draw_png=False,
                      cytoscape_js=False):
    df = df[df['term_name'].isin(terms)].copy()
    df['combined_score'] = np.abs(df['combined_score'])
    df['combined_score'] = np.log2(df['combined_score'])
    df.loc[df['combined_score'] > 150, 'combined_score'] = 150

    term_dict = dict()
    label_dict = dict()
    for i in terms:
        genes = set(et.term_to_genes(df, i))
        term_dict[i] = genes
        label_dict[i] = i
        print(i, len(genes))
    all_genes = set()
    for i, j in term_dict.items():
        all_genes.update(j)
    ong = OntologyNetworkGenerator(molecular_network=network)
    print("Looking for direct edges")
    term_g, molecular_g = ong.create_network_from_list(
        terms, term_dict, label_dict, save_name=save_name, draw=draw_png)

    if cytoscape_js:
        display_graph(term_g)
        display_graph(molecular_g)
    else:
        return term_g, molecular_g


def display_graph(graph, add_parent=False, display_format='cytoscape'):
    g_copy = graph.copy()
    if display_format == 'cytoscape':
        if add_parent:
            new_nodes = set()
            for i, data in graph.nodes(data=True):
                if 'termName' in data:
                    g_copy.node[i]['parent'] = data['termName']
                    new_nodes.add(data['termName'])
            for each in new_nodes:
                g_copy.add_node(each, )
        g_cyjs = from_networkx(g_copy)
        cyjs.render(g_cyjs, style='Directed', layout_algorithm='cose-bilkent')
    elif display_format == 'igraph':
        display(SVG(create_figure(g_copy)))





def shortest_paths(graph, node_1, node_2, bidirectional):
    ns = NetworkSubgraphs(network=graph)
    display_graph(
        ns.shortest_paths_between_two_proteins(node_1, node_2,
                                               bidirectional=bidirectional)
    )
