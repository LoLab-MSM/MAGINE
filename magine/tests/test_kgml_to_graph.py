import networkx

from magine.networks.databases.kegg_kgml import download_kegg_pathway, \
    find_kegg_pathways, kgml_to_graph


def test_find_kegg_pathways():
    """
    tests to see if we can create a list of pathways from a list of proteins
    """
    proteins = ['tp53', 'bax']
    pathways = find_kegg_pathways(proteins, 1, download=False)
    for i in pathways:
        print(i)


def test_construction():
    download_kegg_pathway("hsa04071")
    pathway, pathways_to_add = kgml_to_graph("hsa04071.xml", output_dir='KEGG')

    assert (type(pathway) == networkx.classes.digraph.DiGraph)

