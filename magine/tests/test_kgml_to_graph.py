import networkx

from magine.networks.kgml_to_networkx_parser import kgml_to_graph
from magine.networks.network_generator import download_kegg_pathway


def test_construction():
    download_kegg_pathway("hsa04071")
    pathway, pathways_to_add, compounds = kgml_to_graph("hsa04071.xml")

    assert (type(pathway) == networkx.classes.digraph.DiGraph)

