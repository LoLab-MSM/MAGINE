from magine.kgml_to_networkx_parser import kgml_to_graph
from magine.network_generator import download_kegg_pathway


def test_construction():
    """

    :return:
    """
    download_kegg_pathway("hsa04071")
    pathway, pathways_to_add, compounds = kgml_to_graph("hsa04071.xml")
    for i in pathway.nodes():
        print(i)

test_construction()