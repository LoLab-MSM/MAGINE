import numpy as np
import pandas as pd

import magine.ontology.enrichment_tools as et
from magine.networks.network_subgraphs import NetworkSubgraphs
from magine.networks.ontology_network import OntologyNetworkGenerator

pd.set_option('display.precision', 5)
pd.set_option('display.max_colwidth', 100)


def create_subnetwork(terms, df, network, save_name=None, draw_png=False,
                      threshold=0):
    """

    Parameters
    ----------
    terms : list
    df : pd.DataFrame
    network : nx.DiGraph
    save_name : str
    draw_png : bool
    threshold

    Returns
    -------

    """
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
        terms, term_dict, label_dict, save_name=save_name, draw=draw_png,
        threshold=threshold
    )

    return term_g, molecular_g


def find_neighbors(g, start, up_stream=True, down_stream=True, max_dist=1):
    subgraph_gen = NetworkSubgraphs(g)
    sg = subgraph_gen.neighbors(start, up_stream, down_stream, max_dist)
    return sg


def find_neighbors_of_list(g, list_start, up_stream=True, down_stream=True,
                           max_dist=1):
    subgraph_gen = NetworkSubgraphs(g)

    new_g = subgraph_gen.neighbors_of_list(list_start, up_stream, down_stream,
                                           max_dist)
    return new_g


def subgraph_from_list(g, list_of_nodes):
    subgraph_gen = NetworkSubgraphs(g)
    new_g = subgraph_gen.shortest_paths_between_lists(list_of_nodes)

    return new_g


if __name__ == '__main__':
    import networkx as nx

    g = nx.DiGraph()
    g.add_edge('BAX', 'CASP3')
    g.add_edge('BCL2', 'BAX')
    g.add_edge('BAX', 'CASP3')
    find_neighbors(g, 'BAX', True, True)
