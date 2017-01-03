from magine.ontology.ontology_tools import path_to_root, \
    find_disjunction_common_ancestor, export_to_dot


def test_path_to_root():
    go_example = 'GO:0044260'
    g = path_to_root(go_example)
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
    common_graph = find_disjunction_common_ancestor(topolc_example)
    n_nodes = len(common_graph.nodes())
    assert n_nodes == 12


def test_find_common_ancestor_3_terms():
    # This is for the case where one of the terms is in the path to root of
    # one of the other terms
    topolc_example = ['GO:0061591', 'GO:0061588', 'GO:0061590']
    common_graph = find_disjunction_common_ancestor(topolc_example)
    # n_nodes = len(common_graph.nodes())
    export_to_dot(common_graph, '3_terms')
    # assert n_nodes == 12


def test_only_one_common():
    # This is for the case where one of the terms is in the path to root of
    # one of the other terms
    topolc_example = ['GO:0071840', 'GO:0009987']
    common_graph = find_disjunction_common_ancestor(topolc_example)
    # n_nodes = len(common_graph.nodes())
    export_to_dot(common_graph, 'single_ancestor')
    # assert n_nodes == 12


test_only_one_common()
