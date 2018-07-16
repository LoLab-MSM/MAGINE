import networkx as nx

from exp_data import exp_data
from magine.networks.network_generator import build_network
from magine.networks.subgraphs import Subgraph
from magine.networks.utils import add_data_to_graph, add_attribute_to_network

genes = [

    # dna damage recognition
    'CHEK1', 'CHEK2', 'ATR', 'TP53', 'MDM2', 'PIK3CD', 'AKT1', 'ABL1', 'MAPK1',

    # cell cycle
    'CDKN1A', 'TP73', 'GADD45A',

    # apoptosis
    'FAS', 'FASLG', 'CASP8', 'CASP3', 'CASP9', 'CYCS', 'BAX',
]


def shortest_paths(network, save_name):
    network = add_attribute_to_network(network, genes, 'seed', 'true')
    network = add_data_to_graph(network, exp_data)
    sg = Subgraph(network)
    shortest = sg.paths_between_list(genes)
    print("{} has {} nodes, {} edges".format(save_name,
                                             len(shortest.nodes),
                                             len(shortest.edges)))

    nx.write_gml(network, 'Networks/{}.gml'.format(save_name))
    nx.write_gml(shortest, 'Networks/{}_trimmed.gml'.format(save_name))
    del network


def iterative_build(seed_list, background_list, save_name):
    kegg_only_canonical = build_network(
        seed_list,
        all_measured_list=background_list,
        use_hmdb=False,
        use_biogrid=False,
        use_reactome=False,
        use_signor=False
    )
    shortest_paths(kegg_only_canonical, save_name + '_kegg_only')

    kegg_hmdb_canonical = build_network(
        seed_list,
        all_measured_list=background_list,
        use_hmdb=True,
        use_biogrid=False,
        use_reactome=False,
        use_signor=False
    )

    shortest_paths(kegg_hmdb_canonical, save_name + '_kegg_hmdb')

    kegg_hmdb_biogrid_canonical = build_network(
        seed_list,
        all_measured_list=background_list,
        use_hmdb=True,
        use_biogrid=True,
        use_reactome=False,
        use_signor=False
    )

    shortest_paths(kegg_hmdb_biogrid_canonical,
                   save_name + '_kegg_hmdb_biogrid')

    kegg_hmdb_biogrid_reactome_canonical = build_network(
        seed_list,
        all_measured_list=background_list,
        use_hmdb=True,
        use_biogrid=True,
        use_reactome=True,
        use_signor=False
    )

    shortest_paths(kegg_hmdb_biogrid_reactome_canonical,
                   save_name + '_kegg_hmdb_biogrid_reactome')

    kegg_hmdb_biogrid_reactome_signor_canonical = build_network(
        seed_list,
        all_measured_list=background_list,
        use_hmdb=True,
        use_biogrid=True,
        use_reactome=True,
        use_signor=True
    )

    shortest_paths(kegg_hmdb_biogrid_reactome_signor_canonical, 'ecn')


if __name__ == '__main__':
    iterative_build(genes, None, 'canonical')
