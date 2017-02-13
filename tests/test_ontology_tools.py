import magine.ontology.ontology_tools as ot
from magine.network_tools import export_to_dot
def test_path_to_root():
    go_example = 'GO:0044260'
    g = ot.path_to_root(go_example)
    nodes = set(g.nodes())
    # note that GO ids are returned as str within ""
    root_nodes = {'"GO:0008152"',
                  '"GO:0043170"',
                  '"GO:0044260"',
                  '"GO:0008150"',
                  '"GO:0071704"',
                  '"GO:0009987"',
                  '"GO:0044237"'}
    assert nodes == root_nodes


def test_find_common_ancestor():
    topolc_example = ['GO:0044260', 'GO:0006139']
    common_graph = ot.find_disjunction_common_ancestor(topolc_example)
    n_nodes = len(common_graph.nodes())
    assert n_nodes == 12


def test_find_common_ancestor_3_terms():
    # This is for the case where one of the terms is in the path to root of
    # one of the other terms
    topolc_example = ['GO:0061591', 'GO:0061588', 'GO:0061590']
    example_cisplatin = ['GO:0006893', 'GO:0031577', 'GO:0033523',
                         'GO:0010923', 'GO:0032786', 'GO:1904667',
                         'GO:0007094', 'GO:0010828']
    # common_graph = ot.find_disjunction_common_ancestor(topolc_example)
    common_graph = ot.find_disjunction_common_ancestor(example_cisplatin)
    # n_nodes = len(common_graph.nodes())
    export_to_dot(common_graph, 'cisplatin_example')
    # assert n_nodes == 12


def test_create_graph_from_list():
    topolc_example = ['GO:0061591', 'GO:0061588', 'GO:0061590']
    example_cisplatin = ['GO:0006893', 'GO:0031577', 'GO:0033523',
                         'GO:0010923', 'GO:0032786', 'GO:1904667',
                         'GO:0007094', 'GO:0010828']
    x = ot.create_graph_to_root_from_list_terms(example_cisplatin)
    export_to_dot(x, 'cisplatin_example')

def test_only_one_common():
    # This is for the case where one of the terms is in the path to root of
    # one of the other terms
    topolc_example = ['GO:0071840', 'GO:0009987']
    common_graph = ot.find_disjunction_common_ancestor(topolc_example)
    # n_nodes = len(common_graph.nodes())
    export_to_dot(common_graph, 'single_ancestor')
    # assert n_nodes == 12


def test_hard_example():
    # Test path that doesn't overlap
    list_of_terms = ['GO:0042776', 'GO:0022027', 'GO:0006999', 'GO:1903044',
                     'GO:0051292']
    common_graph = ot.find_disjunction_common_ancestor(list_of_terms)


def cisplatin_example():
    example_cisplatin = ['GO:0097062', 'GO:0033523', 'GO:0031441',
        'GO:1900364', 'GO:1900363', 'GO:0022027', 'GO:0090557', 'GO:1901077',
        'GO:1903299', 'GO:0050686', 'GO:0031468', 'GO:0010666', 'GO:0010663',
        'GO:0007084', 'GO:0048025', 'GO:0002190', 'GO:0046931', 'GO:1903044',
        'GO:0035376', 'GO:0006999', 'GO:0022027', 'GO:0042776', 'GO:0051292',
        'GO:0032594']
    x = ot.create_graph_to_root_from_list_terms(example_cisplatin)
    export_to_dot(x, 'cisplatin_example')
    quit()
    common_graph = ot.find_disjunction_common_ancestor(example_cisplatin)
    export_to_dot(common_graph, 'cisplatin_example_painted')
    example_cisplatin = ot.check_term_list(example_cisplatin)
    common_graph = ot.find_disjunction_common_ancestor(example_cisplatin)
    export_to_dot(common_graph, 'cisplatin_example_painted2')


# test_hard_example()
# test_create_graph_from_list()
# test_find_common_ancestor_3_terms()
cisplatin_example()
