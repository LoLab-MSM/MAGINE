import networkx as nx
from nose.tools import ok_

import magine.networks.databases.kegg_kgml as kegg_tools
from magine.mappings.maps import convert_all

digraph = nx.classes.digraph.DiGraph


def test_construction():
    network = kegg_tools.pathway_id_to_network('hsa00040')
    ok_(type(network) == digraph)
    network = kegg_tools.pathway_id_to_network('hsa04978')
    ok_(type(network) == digraph)
    network = kegg_tools.pathway_id_to_network('hsa04071')
    ok_(type(network) == digraph)
    network = kegg_tools.pathway_id_to_network('hsa04973')
    ok_(type(network) == digraph)


def test_convert_all():
    network = kegg_tools.pathway_id_to_network('hsa00040')
    network = convert_all(network, verbose=True)

    ans = {'chemName': 'D-Ribulose 5-phosphate', 'keggName': 'C00199',
           'hmdbNames': 'HMDB0000618', 'databaseSource': 'KEGG',
           'speciesType': 'compound'}
    # print(network.node['HMDB0000618'])
    # print(network.node['DHDH'])
