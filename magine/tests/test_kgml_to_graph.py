import networkx as nx

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

    verbose = False

    if verbose:
        for i, data in pathway.nodes(data=True):
            print(i, data)

        for i, j, data in pathway.edges(data=True):
            print(i, j, data)
        print(len(pathway.nodes()))
        print(len(pathway.edges()))

    assert (type(pathway) == nx.classes.digraph.DiGraph)
