import itertools
from sys import modules

import networkx as nx
import pandas as pd

from magine.ontology.go_from_goatools import go
from goatools.semantic import deepest_common_ancestor, resnik_sim, \
    semantic_similarity, ic

try:
    termcounts = modules['termcounts']
except:
    from magine.ontology.go_from_goatools import load_termcount
    print("Loading termcounts")
    termcounts = load_termcount()
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


# number of edges between each pair of nodes
def common_parent_go_ids(terms):
    """ Find common parents of list of GO terms

    Parameters
    ----------
    terms: list
    """
    # Find candidates from first
    terms = list(terms)
    rec = go[terms[0]]
    candidates = rec.get_all_parents()
    candidates.update({terms[0]})

    # Find intersection with second to nth term
    for term in terms[1:]:
        rec = go[term]
        parents = rec.get_all_parents()
        parents.update({term})

        # Find the intersection with the candidates, and update.
        candidates.intersection_update(parents)

    return candidates


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


def check_if_children(list_of_terms):
    to_remove = set()
    for i in list_of_terms:
        name = '"{}"'.format(i)
        g = path_to_root(i)
        nodes = set(g.nodes())
        for j in list_of_terms:
            name2 = '"{}"'.format(j)
            if not name == name2:
                if name2 in nodes:
                    print(
                        "{} is a parent of {}, removing from list".format(i,
                                                                          j))
                    to_remove.add(i)
    return list_of_terms.difference(to_remove)


def check_depth(list_of_terms):
    list_of_terms = set(list_of_terms)
    to_remove = set()
    for i in list_of_terms:
        # print(go[i].depth)
        if go[i].depth > 10:
            print("Depth of {} is greater than 10.".format(i),
                  "Removing from list for now as path to root is VAST")
            to_remove.add(i)
    return list_of_terms.difference(to_remove)


def check_depth_and_children(list_of_terms):
    set_of_terms = set(list_of_terms)
    return check_if_children(check_depth(set_of_terms))


def check_term_list(list_of_terms):
    list_of_terms = check_depth_and_children(list_of_terms)
    combos = itertools.combinations(list_of_terms, 2)
    to_remove = set()
    for i, j in combos:
        sim2 = semantic_similarity(i, j, go)
        dca = deepest_common_ancestor([i, j], go)

        sim = resnik_sim(i, j, go, termcounts)
        # this is the IC for biological_process,
        # subtracting it to normalize to 0
        sim -= 4.175976034663472

        print("\nGO 1\t\t GO 2\t\t GO3")
        print("{}\t{}\t{}".format(i, j, dca))
        print("{}\t{}\t{}".format(go[i].name, go[j].name, go[dca].name))
        print("{}\t{}\t{}".format(go[i].depth, go[j].depth, go[dca].depth))
        print("Resnik\tSemantic")
        print("{}\t{}\n".format(sim, sim2))
        ic_i = ic(i, termcounts)
        ic_j = ic(j, termcounts)
        ic_dca = ic(dca, termcounts)
        print("IC\t\t{}\t\t\t{}\t\t\t{}".format(ic_i, ic_j, ic_dca))

        # remove if two terms of similar by distance in graph
        if sim2 > .5:
            if go[i].depth > go[j].depth:
                to_remove.add(i)
            else:
                to_remove.add(j)
        # remove if deepest common ancestor has at least 4 IC
        elif sim > 4:
            if go[dca].depth < 3:
                continue
            if dca in to_remove or dca in list_of_terms:
                list_of_terms.add(dca)
                to_remove.add(i)
                to_remove.add(j)
            print('Removing')
            print('-' * 20)
            print(i, j, dca, go[i].depth, go[j].depth, go[dca].depth)
            print(go[i].name, go[j].name, go[dca].name, sim, sim2)
            print('-' * 20)
                # else:
            # I was thinking I could compare IC for regulators,
            # but it doesn't really provide much
                # ic_i = ic(i, termcounts)
                # ic_j = ic(j, termcounts)
                # print("GO\t\t{}\t\t{}".format(go[i].name, go[j].name))
            #
                # print("Score = {}\n".format(sim2))
            # if go[dca].depth > 3:
            #     print("Depth\t{}\t\t\t{}".format(go[i].depth, go[j].depth))
            #     print("DCA\tID\tDepth\tIC")
            #     print("{}\t{}\t{}\t{}\t".format(go[dca].name, dca, go[dca].depth, sim))

    print(to_remove)
    print('Number of term before = {}'.format(len(list_of_terms)))
    list_of_terms = list_of_terms.difference(to_remove)
    print('Number of term after = {}'.format(len(list_of_terms)))
    return list_of_terms


def create_graph_to_root_from_list_terms(list_of_terms):
    """ calculate a path to root for a list of terms

    Parameters
    ----------
    list_of_terms: list

    Returns
    -------

    networkx.DiGraph
    """

    list_of_terms = check_depth_and_children(list_of_terms)
    combined_graph = nx.DiGraph()
    for i in list_of_terms:
        g = path_to_root(i)
        comb_nodes = set(combined_graph.nodes())
        combined_graph.add_nodes_from(g.nodes(data=True))
        combined_graph.add_edges_from(g.edges(data=True))
        for i in g.nodes():
            if i in comb_nodes:
                if 'counter' in combined_graph.node[i]:
                    combined_graph.node[i]['counter'] += 1
            else:
                combined_graph.node[i]['counter'] = 1
    for each in list_of_terms:
        combined_graph.add_node(each, style='filled', fillcolor='green',
                                **combined_graph[each])

    return combined_graph


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

    # find the nodes that are in all terms
    n_terms = len(list_of_terms)
    term_counter = nx.get_node_attributes(combined_graph, 'counter')
    common_nodes = set()
    for i in term_counter:
        print(i, term_counter[i], go[i.replace('"', '')].depth)
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
    to_remove = []
    for each in list_of_terms:
        if each in common_nodes:
            print("Warning : {} term is already in graph! Thus will treat as "
                  "upstream term and remove from list of terms"
                  " provided".format(each))
            to_remove.append(each)
    for each in to_remove:
        list_of_terms.remove(each)

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
                sh_paths = [p for p in
                            nx.all_shortest_paths(combined_graph, disjun_node,
                                                  each)]
                for path in sh_paths:
                    for g in path:
                        common_graph.add_node(g, **combined_graph.node[g])
                    common_graph.add_path(path)

    # color source nodes
    for each in list_of_terms:
        common_graph.add_node(each, style='filled',
                              fillcolor='green', **combined_graph.node[each])
    gonames = nx.get_node_attributes(combined_graph, 'GOname')
    nodes_in_graph = common_graph.nodes()
    to_remove = set()
    for i in gonames:
        if i not in nodes_in_graph:
            to_remove.add(i)
    for i in to_remove:
        gonames.pop(i)
    # nx.set_node_attributes(common_graph, 'label', gonames)
    # export_to_dot(common_graph, 'common', view=False)
    return common_graph


def filter_ontology_df(data, n_hits_per_time=None, go_aspects=None,
                       trim_nodes=False, additional_ids_to_include=None):
    """

    Parameters
    ----------
    data :  pandas.DataFrame or str
        output from magine.ontology_analysis.create_enrichment_array
    trim_nodes : bool, optional, default=False
        remove GO terms that are similar
    n_hits_per_time : int, optional
        number of terms for each sample
    go_aspects : list, optional
        aspects to plots, options are
        {"biological_process", "cellular_compartment", "molecular_function"}
    additional_ids_to_include: list
        list of additional GO ids to keep

    Returns
    -------

    """

    if isinstance(data, str):
        data = pd.read_csv(data)

    # make sure all columns we need are defined
    for i in ['GO_id', 'pvalue', 'enrichment_score', 'sample_index', 'aspect']:
        assert i in data.columns, "Must have {} in df".format(i)

    data['genes'] = data['genes'].astype(set)

    # removes aspects of GO that are not wanted
    if go_aspects is None:
        go_aspects = ['biological_process']
    data = data[data['aspect'].isin(go_aspects)]

    # remove terms with reference smaller than 5
    data = data[data['ref'] >= 10]
    data = data[data['ref'] <= 300]

    # normalize maximum enrichment to 20
    # done for simplicity, may not be desired
    data.loc[data['enrichment_score'] > 20, 'enrichment_score'] = 20

    # filter out non-signficant terms
    tmp = data[data['pvalue'] < 0.05]

    # create sample labels
    labels = tmp['sample_index'].unique()

    def find_n_go_terms(n_terms):
        terms = set(tmp.head(n_terms).index)
        terms = check_term_list(terms)
        return terms

    enrichment_list = [('enrichment_score', i) for i in labels]
    list_all_go = tmp['GO_id'].unique()
    if n_hits_per_time is not None:
        tmp = pd.pivot_table(tmp, index=['GO_id', ], columns='sample_index')
        list_all_go = set()
        for i in enrichment_list:
            tmp = tmp.sort_values(by=i, ascending=False)
            list_of_go = find_n_go_terms(n_hits_per_time)
            count = 1
            while len(list_of_go) < n_hits_per_time:
                list_of_go = find_n_go_terms(n_hits_per_time + count)
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
