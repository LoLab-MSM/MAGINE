import networkx as nx

import magine.networks.databases.kegg_kgml as kegg_tools
from magine.mappings.maps import convert_all


def test_construction():
    network = kegg_tools.pathway_id_to_network('hsa00040')
    assert (type(network) == nx.classes.digraph.DiGraph)
    network = kegg_tools.pathway_id_to_network('hsa04978')
    assert (type(network) == nx.classes.digraph.DiGraph)
    network = kegg_tools.pathway_id_to_network('hsa04071')
    assert (type(network) == nx.classes.digraph.DiGraph)
    network = kegg_tools.pathway_id_to_network('hsa04973')
    assert (type(network) == nx.classes.digraph.DiGraph)
    network = convert_all(network, verbose=True)
    print(network.nodes)
