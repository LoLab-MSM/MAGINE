import itertools
import os

import networkx as nx
import pandas as pd
from goatools import obo_parser
from goatools.semantic import TermCounts, ic, resnik_sim, semantic_similarity
from magine.enrichment.ontology_analysis import MagineGO

from magine.data.storage import id_mapping_dir

obo_file = os.path.join(id_mapping_dir, 'go.obo')

if not os.path.exists(obo_file):
    print("Using ontology for first time")
    print("Downloading files")
    from magine.enrichment.databases.gene_ontology import \
        download_and_process_go
    download_and_process_go()
    assert os.path.exists(obo_file)

go = obo_parser.GODag(obo_file)

mg = MagineGO()
print("Loading termcounts")
associations = mg.gene_to_go
termcounts = TermCounts(go, associations)
print("Loaded termcounts")


def path_to_root(go_term):
    """
    Creates networkx graph from provided term to root term

    Parameters
    ----------
    go_term : str
        source GO term

    Returns
    -------
    graph : nx.DiGraph
        Graph containing all paths from provided GO term to root

    """

    graph = nx.DiGraph()
    paths = go.paths_to_top(go_term)
    all_terms = set()
    for i in paths:
        for j in i:
            all_terms.add(j.id)
    for i in paths:
        for j in i:
            add_children(j.id, graph, all_terms)
    return graph


def print_path_to_root(go_term):
    paths = go.paths_to_top(go_term)
    all_terms = set()
    for i in paths:
        print('\n')
        for j in i:
            all_terms.add(j.id)
            print("GO id = {}, GO name = {}".format(j.id, j.name))


def print_children(go_term):
    print("({}) {} has children of :".format(go[go_term].id, go[go_term].name))
    for i in go[go_term].children:
        print("\t({}) {}".format(i.id, i.name))


def common_parent_go_ids(terms):
    """ Find common parents of list of GO terms

    Parameters
    ----------
    terms: list
    """
    # Find candidates from first
    terms = list(terms)
    path_parents = go[terms[0]]
    candidates = path_parents.get_all_parents()
    candidates.update({terms[0]})

    # Find intersection with second to nth term
    for term in terms[1:]:
        path_parents = go[term]
        parents = path_parents.get_all_parents()
        parents.update({term})

        # Find the intersection with the candidates, and update.
        candidates.intersection_update(parents)

    return candidates


def deepest_common_ancestor(terms):
    """
        This function gets the nearest common ancestor
        using the above function.
        Only returns single most specific - assumes unique exists.
    """
    # Take the element at maximum depth.
    return max(common_parent_go_ids(terms), key=lambda t: go[t].depth)


def min_branch_length(go_id1, go_id2):
    """
        Finds the minimum branch length between two terms in the GO DAG.
    """
    # First get the deepest common ancestor
    dca = deepest_common_ancestor([go_id1, go_id2])

    # Then get the distance from the DCA to each term
    dca_depth = go[dca].depth
    d1 = go[go_id1].depth - dca_depth
    d2 = go[go_id2].depth - dca_depth
    # print(d1, d2)
    # Return the distance from one term to the deepest common ancestor and
    # back to the other term.
    return d1 + d2


def add_children(go_term, graph, gene_set_of_interest):
    """
    add children to graph of provided GO term

    Parameters
    ----------
    go_term : str
        parent GO term
    graph : networkx.DiGraph
        graph to add children
    gene_set_of_interest : set
        all terms to consider

    Returns
    -------

    """
    nodes_in_graph = set(graph.nodes())
    graph.add_node(go[go_term].id,
                   depth=go[go_term].depth,
                   level=go[go_term].level,
                   GOname=go[go_term].name,
                   ic=ic(go[go_term].id, termcounts))
    for i in go[go_term].children:
        if i.id in gene_set_of_interest:
            if i.id not in nodes_in_graph:
                graph.add_node(i.id,
                               depth=i.depth,
                               level=i.level,
                               GOname=i.name,
                               ic=ic(i.id, termcounts))
            graph.add_edge(go[go_term].id, i.id)
            add_children(i.id, graph, gene_set_of_interest)


def check_if_children(list_of_terms, verbose=False):
    """
    Checks to see if any term is a child of another


    Parameters
    ----------
    list_of_terms : list
        list of all terms to consider
    verbose : bool

    Returns
    -------

    """
    list_of_terms = set(list_of_terms)
    to_remove = set()
    for i in list_of_terms:
        g = path_to_root(i)
        nodes = set(g.nodes())
        for j in list_of_terms:
            if not i == j:
                if j in nodes:
                    if verbose:
                        print("{} is a parent of {}, "
                              "removing from list".format(i, j))
                    to_remove.add(i)
    return list_of_terms.difference(to_remove)


def check_depth_level(list_of_terms, verbose=False):
    """
    Checks to see if any terms are too deep in graph

    Parameters
    ----------
    list_of_terms : list
        list of GO terms
    verbose : bool


    Returns
    -------

    """
    list_of_terms = set(list_of_terms)
    to_remove = set()

    for i in list_of_terms:
        if go[i].depth > 11:
            if verbose:
                print("Depth of {} is greater than 11.".format(i),
                      "Removing from list for now as path to root is VAST")
            to_remove.add(i)
    return list_of_terms.difference(to_remove)


def check_depth_and_children(list_of_terms):
    set_of_terms = set(list_of_terms)
    return check_if_children(check_depth_level(set_of_terms))


def create_graph_to_root_from_list_terms(list_of_terms, check_child=False,
                                         check_depth=True):
    """ calculate a path to root for a list of terms

    Parameters
    ----------
    list_of_terms: list
        list of GO terms
    check_child: bool
        if you want to check if there are terms that are children of other
        terms within your list
    check_depth: bool
        check for maximum depth value. Depth larger than 10 usually
        takes a long time to find all paths.


    Returns
    -------
    networkx.DiGraph

    """
    if check_child:
        list_of_terms = check_if_children(list_of_terms)
    if check_depth:
        list_of_terms = check_if_children(list_of_terms)
    combined_graph = nx.DiGraph()
    for i in list_of_terms:
        g = path_to_root(i)
        comb_nodes = set(combined_graph.nodes())
        combined_graph.add_nodes_from(g.nodes(data=True))
        combined_graph.add_edges_from(g.edges(data=True))
        for n in g.nodes():
            if n in comb_nodes:
                combined_graph.node[n]['counter'] += 1
            else:
                combined_graph.node[n]['counter'] = 1
    for each in list_of_terms:
        combined_graph.add_node(each, style='filled', fillcolor='green',
                                **combined_graph[each])

    return combined_graph


def check_term_list(list_of_terms, verbose=False):
    """

    Parameters
    ----------
    list_of_terms : list
    verbose : bool

    Returns
    -------

    """
    list_of_terms = check_depth_and_children(list_of_terms)
    combos = itertools.combinations(list_of_terms, 2)
    to_remove = set()
    for i, j in combos:

        sim2 = semantic_similarity(i, j, go)
        dca = deepest_common_ancestor([i, j])

        sim = resnik_sim(i, j, go, termcounts)
        ic_i = ic(i, termcounts)
        ic_j = ic(j, termcounts)
        ic_dca = ic(dca, termcounts)

        if verbose:
            print("\nGO 1\t\t GO 2\t\t GO3")
            print("Min branch length = {}".format(min_branch_length(i, j)))
            print("{}\t{}\t{}".format(go[i].name, go[j].name, go[dca].name))
            print("{}\t{}\t{}".format(i, j, dca))
            print("{}\t{}\t{}".format(go[i].depth, go[j].depth, go[dca].depth))
            print("Resnik {}\tSemantic {}".format(sim, sim2))
            # print("{}\t{}\n".format(sim, sim2))

            print("IC\t\t{}\t\t\t{}\t\t\t{}".format(ic_i, ic_j, ic_dca))

        # remove if two terms of similar by distance in graph
        if sim2 > .5:
            if go[i].depth > go[j].depth:
                to_remove.add(i)
            else:
                to_remove.add(j)
        # remove if deepest common ancestor has at least 4 IC
        elif sim > 9:
            if go[dca].depth < 3:
                continue
            list_of_terms.add(dca)
            to_remove.add(i)
            to_remove.add(j)
            if dca in to_remove or dca in list_of_terms:
                list_of_terms.add(dca)
                to_remove.add(i)
                to_remove.add(j)
    if verbose:
        print(to_remove)
        print('Number of term before = {}'.format(len(list_of_terms)))

    list_of_terms = list_of_terms.difference(to_remove)

    if verbose:
        print('Number of term after = {}'.format(len(list_of_terms)))

    return list_of_terms


def find_disjunction_common_ancestor(list_of_terms):
    """

    Parameters
    ----------
    list_of_terms : list_like
        list of GO terms

    Returns
    -------
    graph : networkx.DiGraph
        containing sources to root, disjunction common ancestor show in red
    """
    list_of_terms = check_depth_and_children(list_of_terms)
    combined_graph = create_graph_to_root_from_list_terms(list_of_terms)

    # dca = deepest_common_ancestor(list(list_of_terms))

    # find the nodes that are in all terms
    n_terms = len(list_of_terms)
    # print("N terms", n_terms)
    term_counter = nx.get_node_attributes(combined_graph, 'counter')
    common_nodes = set()
    for i in term_counter:
        # print(i, term_counter[i], go[i.replace('"', '')].depth)
        if term_counter[i] == n_terms:
            common_nodes.add(i)

    if len(common_nodes) == 1:
        node_0 = combined_graph.nodes()[0]
        level = combined_graph.node[node_0]['level']
        depth = combined_graph.node[node_0]['depth']
        name = combined_graph.node[node_0]['GOname']
        print("Warning : the most common ancestor is {} ({}) with a depth of "
              "{} and level of {}".format(name, node_0, depth, level))

    # Check to see if any of the provided terms are in all graphs
    # This means it was an ancestor node and should be removed from list
    # to_remove = []
    # for each in list_of_terms:
    #     if each in common_nodes:
    #         print("Warning : {} term is already in graph! Thus will treat as "
    #               "upstream term and remove from list of terms"
    #               " provided".format(each))
    #         to_remove.append(each)
    # for each in to_remove:
    #     list_of_terms.remove(each)
    print(common_nodes)
    # return subgraph containing only nodes that are in all
    common_graph = combined_graph.subgraph(common_nodes)

    # find the common dijunction nodes and shortest paths between
    disjunction_nodes = set()
    for i in common_graph.nodes():
        if common_graph.out_degree(i) == 0:
            disjunction_nodes.add(i)

    # color disjunction nodes and find all shortest paths between source nodes
    # and each disjunction node
    for disjun_node in disjunction_nodes:
        common_graph.node[disjun_node]['style'] = 'filled'
        common_graph.node[disjun_node]['fillcolor'] = 'red'
        for each in list_of_terms:
            if nx.has_path(combined_graph, disjun_node, each):
                print("From {} to {}".format(disjun_node, each))
                sh_paths = [p for p in
                            nx.all_shortest_paths(combined_graph, disjun_node,
                                                  each)]

                # ic_paths = []
                # for path in sh_paths:
                #     ic_path = 0
                #     for g in path:
                #         ic_local = combined_graph.node[g]['ic']
                #         ic_path += ic_local
                #     ic_paths.append(ic_path)
                # min_ic = max(ic_paths)
                for n, path in enumerate(sh_paths):
                    # if ic_paths[n] == min_ic:
                    #     for g in path:
                    #         common_graph.add_node(g, **combined_graph.node[g])
                    #     common_graph.add_path(path)
                    for g in path:
                        common_graph.add_node(g, **combined_graph.node[g])
                    common_graph.add_path(path)

    # color source nodes
    for each in list_of_terms:
        common_graph.add_node(each, **combined_graph.node[each])

    return common_graph


def filter_ontology_df(data, n_hits_per_time=None, go_aspects=None,
                       trim_nodes=False, additional_ids_to_include=None,
                       pvalue=0.05, min_depth=3, max_depth=10, min_ref=5,
                       max_ref=500, variable_of_interest='enrichment_score'):
    """

    Parameters
    ----------
    data :  pandas.DataFrame 
        output from magine.ontology_analysis.calculate_enrichment
    trim_nodes : bool, optional, default=False
        remove GO terms that are similar
    n_hits_per_time : int, optional
        number of terms for each sample
    go_aspects : list, optional
        aspects to plots, options are
        {"biological_process", "cellular_compartment", "molecular_function"}
    additional_ids_to_include: list
        list of additional GO ids to keep
    pvalue : float
        upper pvalue limit to filter terms
    min_depth : int
        minimum depth of GO ontology
    max_depth : int
        maximum depth of GO ontology
    min_ref : int
        minimum number of references for GO ontology term
    max_ref : int
        maximum number of references for GO ontology term
    variable_of_interest : str, {"enrichment_score", "pvalue}



    Returns
    -------

    """

    if isinstance(data, str):
        data = pd.read_csv(data)

    # make sure all columns we need are defined
    for i in ['GO_id', 'pvalue', 'enrichment_score', 'sample_index', 'aspect']:
        if i not in data.columns:
            raise AssertionError("Must have {} in df".format(i))

    data['genes'] = data['genes'].astype(set)

    # removes aspects of GO that are not wanted
    if go_aspects is None:
        go_aspects = ['biological_process']
    if isinstance(go_aspects, str):
        data = data[data['aspect'] == go_aspects]
    elif isinstance(go_aspects, list):
        data = data[data['aspect'].isin(go_aspects)]
    else:
        print("go_aspect must be a list! \n"
              "biological_process, cellular_component, molecular_function")

    # filter terms based on reference and depth
    if max_ref is not None:
        data = data[data['ref'] <= max_ref]
    if min_ref is not None:
        data = data[data['ref'] >= min_ref]
    if max_depth is not None:
        data = data[data['depth'] <= max_depth]
    if min_depth is not None:
        data = data[data['depth'] >= min_depth]

    # filter out non-signficant terms
    tmp = data[data['pvalue'] <= pvalue]

    # create sample labels
    labels = tmp['sample_index'].unique()

    enrichment_list = [(variable_of_interest, i) for i in labels]
    ascend = False
    if variable_of_interest == 'pvalue':
        ascend = True

    list_all_go = tmp['GO_id'].unique()
    if n_hits_per_time is not None:
        tmp = pd.pivot_table(tmp, index=['GO_id', ], columns='sample_index')

        def find_n_go_terms(terms, updated_index):
            n_needed = n_hits_per_time - len(terms)
            terms.update(set(list(tmp.index)[:updated_index + n_needed]))
            return check_term_list(terms)

        list_all_go = set()
        for i in enrichment_list:
            tmp = tmp.sort_values(by=i, ascending=ascend)
            list_of_go = set(tmp.head(n_hits_per_time).index)
            if trim_nodes:
                list_of_go = check_term_list(list_of_go)
                count = n_hits_per_time
                while len(list_of_go) < n_hits_per_time:
                    list_of_go = find_n_go_terms(list_of_go, count)
                    count += 1

            list_all_go.update(list_of_go)

    if trim_nodes:
        list_all_go = check_term_list(list_all_go)

    if additional_ids_to_include is not None:
        assert isinstance(additional_ids_to_include, list)
        additional_ids_to_include = set(additional_ids_to_include)
        if n_hits_per_time is not None:
            additional_ids_to_include.update(list_all_go)
        data = data[data['GO_id'].isin(additional_ids_to_include)]

    elif list_all_go is not None:
        data = data[data['GO_id'].isin(list_all_go)]

    return data


def slim_ontology(data, pvalue=0.05, n_top_hits=None, go_aspects=None,
                  trim_nodes=False, additional_ids_to_include=None,
                  min_depth=3, max_depth=10, min_ref=5, max_ref=500,
                  variable_of_interest='enrichment_score'):
    """

    Parameters
    ----------
    data :  pandas.DataFrame 
        output from magine.ontology_analysis.calculate_enrichment
    pvalue : float
        upper pvalue limit to filter terms
    trim_nodes : bool, optional, default=False
        remove GO terms that are similar
    n_top_hits : int, optional
        number of terms for each sample
    go_aspects : list, optional
        aspects to plots, options are
        {"biological_process", "cellular_compartment", "molecular_function"}
    additional_ids_to_include: list
        list of additional GO ids to keep
    min_depth : int
        minimum depth of GO ontology
    max_depth : int
        maximum depth of GO ontology
    min_ref : int
        minimum number of references for GO ontology term
    max_ref : int
        maximum number of references for GO ontology term
    variable_of_interest : str, {"enrichment_score", "pvalue"}
    Returns
    -------

    """

    if isinstance(data, str):
        data = pd.read_csv(data)

    # make sure all columns we need are defined
    for i in ['GO_id', 'pvalue', 'enrichment_score', 'aspect']:
        assert i in data.columns, "Must have {} in df".format(i)

    data['genes'] = data['genes'].astype(set)

    # removes aspects of GO that are not wanted
    if go_aspects is None:
        go_aspects = ['biological_process']
    if isinstance(go_aspects, str):
        data = data[data['aspect'] == go_aspects]
    elif isinstance(go_aspects, list):
        data = data[data['aspect'].isin(go_aspects)]
    else:
        print("go_aspect must be a list! \n"
              "biological_process, cellular_component, molecular_function")

    # remove terms with reference smaller than 5
    if max_ref is not None:
        data = data[data['ref'] <= max_ref]
    if min_ref is not None:
        data = data[data['ref'] >= min_ref]
    if max_depth is not None:
        data = data[data['depth'] <= max_depth]
    if min_depth is not None:
        data = data[data['depth'] >= min_depth]

    # normalize maximum enrichment to 20
    # done for simplicity, may not be desired
    # data.loc[data['enrichment_score'] > 20, 'enrichment_score'] = 20

    # filter out non-significant terms
    tmp = data[data['pvalue'] <= pvalue].copy()
    assert variable_of_interest in ('enrichment_score', 'pvalue'), \
        "variable_of_interest must be 'enrichment_score' or 'pvalue'"

    enrichment_list = [variable_of_interest]

    list_all_go = tmp['GO_id'].unique()
    if n_top_hits is not None:
        def find_n_go_terms(terms, updated_index):
            n_needed = n_top_hits - len(terms)
            terms.update(set(list(tmp.index)[:updated_index + n_needed]))
            return check_term_list(terms)

        list_all_go = set()
        tmp.set_index('GO_id', inplace=True)
        tmp.sort_values(by=enrichment_list, ascending=False, inplace=True)
        list_of_go = set(tmp.head(n_top_hits).index)
        if trim_nodes:
            list_of_go = check_term_list(list_of_go)
            print("start {}", format(list_of_go))
            count = n_top_hits
            while len(list_of_go) < n_top_hits:
                list_of_go = find_n_go_terms(list_of_go, count)
                count += 1
        list_all_go.update(list_of_go)

    if trim_nodes:
        list_all_go = check_term_list(list_all_go)

    if additional_ids_to_include is not None:
        assert isinstance(additional_ids_to_include, list)
        additional_ids_to_include = set(additional_ids_to_include)
        if n_top_hits is not None:
            additional_ids_to_include.update(list_all_go)
        data = data[data['GO_id'].isin(additional_ids_to_include)]

    elif list_all_go is not None:
        data = data[data['GO_id'].isin(list_all_go)]

    return data
