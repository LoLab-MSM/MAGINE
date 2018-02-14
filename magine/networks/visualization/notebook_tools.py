import numpy as np

import magine.ontology.enrichment_tools as et
from magine.networks.network_subgraphs import NetworkSubgraphs
from magine.networks.ontology_network import OntologyNetworkGenerator
from magine.networks.visualization.notebooks.view import display_graph, \
    render_graph


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
    return term_g, molecular_g


def find_neighbors(g, start, up_stream=True, down_stream=True, max_dist=1,
                   render=False):
    subgraph_gen = NetworkSubgraphs(g)
    sg = subgraph_gen.neighbors(start, up_stream, down_stream, max_dist)
    if render:
        render_graph(sg)
    else:
        return sg


def find_neighbors_of_list(g, list_start, up_stream=True, down_stream=True,
                           max_dist=1, render=False):
    subgraph_gen = NetworkSubgraphs(g)

    new_g = subgraph_gen.neighbors_of_list(list_start, up_stream, down_stream,
                                           max_dist)

    if render:
        render_graph(new_g)
    else:
        return new_g


def subgraph_from_list(g, list_of_nodes, render=False):
    subgraph_gen = NetworkSubgraphs(g)
    new_g = subgraph_gen.shortest_paths_between_lists(list_of_nodes)
    if render:
        render_graph(new_g)
    else:
        return new_g


if __name__ == '__main__':
    import networkx as nx

    g = nx.DiGraph()
    g.add_edge('BAX', 'CASP3')
    g.add_edge('BCL2', 'BAX')
    g.add_edge('BAX', 'CASP3')
    find_neighbors(g, 'BAX', True, True)
