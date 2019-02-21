import networkx as nx
from nose.tools import ok_

import magine.networks.databases.kegg_kgml as kegg_tools
from magine.mappings.maps import convert_all

digraph = nx.classes.digraph.DiGraph


def test_construction():
    network = kegg_tools.pathway_id_to_network('hsa00040')
    ok_(isinstance(network, digraph))
    network = kegg_tools.pathway_id_to_network('hsa04978')
    ok_(isinstance(network, digraph))
    network = kegg_tools.pathway_id_to_network('hsa04071')
    ok_(isinstance(network, digraph))
    network = kegg_tools.pathway_id_to_network('hsa04973')
    ok_(isinstance(network, digraph))


def test_convert_all():
    network = kegg_tools.pathway_id_to_network('hsa00040')
    convert_all(network, verbose=True)
