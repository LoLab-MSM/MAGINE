"""
import textwrap

import networkx as nx

import magine.enrichment.ontology_tools as ot
from magine.networks.network_tools import export_to_dot

example_cisplatin = ['GO:0097062', 'GO:0033523', 'GO:0031441', 'GO:1900364',
                     'GO:1900363', 'GO:0022027', 'GO:0090557', 'GO:1901077',
                     'GO:1903299', 'GO:0050686', 'GO:0031468', 'GO:0010666',
                     'GO:0010663', 'GO:0007084', 'GO:0048025', 'GO:0002190',
                     'GO:0046931', 'GO:1903044', 'GO:0035376', 'GO:0006999',
                     'GO:0022027', 'GO:0042776', 'GO:0051292', 'GO:0032594']
topolc_example = ['GO:0044260', 'GO:0006139']


def test_path_to_root():
    go_example = 'GO:0044260'
    g = ot.path_to_root(go_example)
    nodes = set(g.nodes())
    root_nodes = {"GO:0008152",
                  "GO:0043170",
                  "GO:0044260",
                  "GO:0008150",
                  "GO:0071704",
                  "GO:0009987",
                  "GO:0044237"}
    assert nodes == root_nodes


def test_find_common_ancestor():
    topolc_example = ['GO:0044260', 'GO:0006139']
    common_graph = ot.find_disjunction_common_ancestor(topolc_example)
    n_nodes = len(common_graph.nodes())
    print(n_nodes)
    # assert n_nodes == 12
    pretty_format(common_graph, 'go_find_common')


def test_find_common_ancestor_3_terms():
    # This is for the case where one of the terms is in the path to root of
    # one of the other terms
    three_term = ['GO:0006893', 'GO:0031577', 'GO:0033523', 'GO:0010923',
                  'GO:0032786', 'GO:1904667', 'GO:0007094', 'GO:0010828']
    # common_graph = ot.find_disjunction_common_ancestor(topolc_example)
    common_graph = ot.find_disjunction_common_ancestor(three_term)
    # n_nodes = len(common_graph.nodes())
    export_to_dot(common_graph, 'cisplatin_example')
    # assert n_nodes == 12


def test_create_graph_from_list():
    x = ot.create_graph_to_root_from_list_terms(example_cisplatin)
    export_to_dot(x, 'cisplatin_example')

def test_only_one_common():
    # This is for the case where one of the terms is in the path to root of
    # one of the other terms
    common_graph = ot.find_disjunction_common_ancestor(topolc_example)
    # n_nodes = len(common_graph.nodes())
    export_to_dot(common_graph, 'single_ancestor')
    # assert n_nodes == 12


def test_hard_example():
    # Test path that doesn't overlap
    list_of_terms = ['GO:0042776', 'GO:0022027', 'GO:0006999', 'GO:1903044',
                     'GO:0051292']
    list_of_terms = ot.check_term_list(topolc_example)

    common_graph = ot.find_disjunction_common_ancestor(list_of_terms)
    export_to_dot(common_graph, 'common_graph_test')
    ic_dict = nx.get_node_attributes(common_graph, 'ic')
    GOname = nx.get_node_attributes(common_graph, 'GOname')
    new_label = dict()
    for i in ic_dict:
        print(i, ic_dict[i], GOname[i])
        new_label[i] = textwrap.fill("{}\n{}".format(GOname[i], ic_dict[i]),
                                     20)
    nx.set_node_attributes(common_graph, 'label', new_label)
    export_to_dot(common_graph, 'common_graph')


def test_get_common_ancestors():
    list_of_terms = ['GO:0042776', 'GO:0022027', 'GO:0006999', 'GO:1903044',
                     'GO:0051292']
    x = ot.common_parent_go_ids(topolc_example)
    for i in x:
        print(i)


def cisplatin_example():
    common_graph = ot.create_graph_to_root_from_list_terms(example_cisplatin)
    print("Number of nodes = {}".format(len(common_graph.nodes())))
    export_to_dot(common_graph, 'cisplatin_example')

    pretty_format(common_graph, 'cisplatin_all_to_top')

    list_of_terms = ot.check_term_list(example_cisplatin)
    common_graph = ot.find_disjunction_common_ancestor(list_of_terms)
    print("Number of nodes = {}".format(len(common_graph.nodes())))

    export_to_dot(common_graph, 'common_graph_test')
    pretty_format(common_graph, 'cisplatin_term_list_check')


def pretty_format(graph, save_name):
    ic = nx.get_node_attributes(graph, 'ic')
    go_name = nx.get_node_attributes(graph, 'GOname')
    new_label = dict()
    for i in ic:
        new_label[i] = textwrap.fill("{}\n{}".format(go_name[i], ic[i]), 20)

    nx.set_node_attributes(graph, 'label', new_label)
    export_to_dot(graph, save_name)


def slimmed_cisplatin_example():
    example_cisplatin = [
        'GO:0031441',
        'GO:0048025',
        'GO:0050686',
        'GO:1900363',
    ]

    common_graph = ot.create_graph_to_root_from_list_terms(example_cisplatin)
    pretty_format(common_graph, 'cisplatin_rna_example_path_to_top')

    list_of_terms = ot.check_term_list(example_cisplatin)
    common_graph = ot.create_graph_to_root_from_list_terms(example_cisplatin)
    pretty_format(common_graph, 'cisplatin_rna_example_checked_list')

    common_graph = ot.find_disjunction_common_ancestor(list_of_terms)
    pretty_format(common_graph, 'cisplatin_rna_example_checked_list_compress')


if __name__ == "__main__":
    # test_path_to_root()
    # test_find_common_ancestor()
    # test_get_common_ancestors()
    # test_hard_example()
    cisplatin_example()
    # slimmed_cisplatin_example()
    # test_create_graph_from_list()
    # test_find_common_ancestor_3_terms()
    # cisplatin_example()
"""
