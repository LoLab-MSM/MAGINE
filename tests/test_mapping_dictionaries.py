# -*- coding: utf-8 -*-
import networkx as nx

from magine.mappings.maps import create_gene_dictionaries, \
    create_compound_dictionary


def test_kegg_to_uniprot():
    """
    tests kegg gene to gene name
    :return:
    """
    g = nx.DiGraph()
    g.add_edge('hsa:217','hsa:219')
    g.add_edge('hsa:217', 'hsa:223')
    g.add_edge('hsa:501', 'hsa:219')
    g.add_edge('hsa:224', 'hsa:219')
    g.add_node('hsa:857')
    dic = create_gene_dictionaries(g, species='hsa')
    g = nx.relabel_nodes(g,dic)
    for i in g.nodes():
        print(i)
    assert (g.node['ALDH3A2']['keggName'] == 'hsa:224')


def test_kegg_to_hmdb():
    """
    tests kegg compound to hmdb
    :return:
    """
    g = nx.DiGraph()
    g.add_edge('hsa:224', 'hsa:219')
    g.add_edge('hsa:219', 'cpd:C00197')
    g.add_edge('cpd:C00197','cpd:C00197')
    g.add_edge('cpd:C15972', 'cpd:C00197')
    g.add_edge('cpd:C15972', 'cpd:C00469')
    dic = create_compound_dictionary(g)
    g = nx.relabel_nodes(g,dic)
    nx.write_gml(g, 'test.gml')
    assert (g.node['HMDB60180']['chemName'] == '(2R)-2-Hydroxy-3-(phosphonatooxy)propanoate')
    assert (g.node['HMDB60180']['keggName'] == 'cpd:C00197')


# def test_hugo():
#     """
#
#     :return:
#     """
#     g = nx.DiGraph()
#     g.add_node('hsa:406915')
#     g.add_node('hsa:406913')
#     g.add_node('hsa:406912')
#     g.add_node('hsa:406911')
#     g.add_node('hsa:406910')
#     dic = hugo_mapper(g)
#     g = nx.relabel_nodes(g, dic)
#     for i in g.nodes():
#         print(i)

if __name__ =='__main__':
    #test_kegg_to_uniprot()
    test_kegg_to_hmdb()
    #test_hugo()
