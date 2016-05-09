from magine.network_generator import find_kegg_pathways


def test_find_kegg_pathways():
    """
    tests to see if we can create a list of pathways from a list of proteins
    :return:
    """
    proteins = ['tp53','bax']
    pathways = find_kegg_pathways(proteins,1)
    for i in pathways:
        print(i)


test_find_kegg_pathways()