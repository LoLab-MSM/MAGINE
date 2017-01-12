import itertools
import os

import networkx as nx
from orangecontrib.bio.utils import serverfiles

from goatools import obo_parser
from goatools.associations import read_gaf
from goatools.semantic import TermCounts, ic, resnik_sim
from goatools.semantic import semantic_similarity, deepest_common_ancestor
from magine.network_tools import export_to_dot

# sys.setrecursionlimit(100000)
os.environ[
    "PATH"] += os.pathsep + "C:\Users\James Pino\Miniconda2\envs\MAGINE\Library\\bin\graphviz"

default_database_path = os.path.join(serverfiles.localpath(), "GO")

short_path = filename = os.path.join(default_database_path,
                                     "gene_ontology_edit.obo.tar.gz",
                                     "gene_ontology_edit.obo")
if not os.path.exists(short_path):
    from orangecontrib.bio.go import Ontology

    print("Using ontology for first time")
    print("Downloading files via Orange.bio")
    ontology = Ontology()
    assert os.path.exists(short_path)

go = obo_parser.GODag(short_path)
associations = read_gaf(
    "http://geneontology.org/gene-associations/goa_human.gaf.gz")

to_remove = set()
for gene, terms in associations.items():
    terms_copy = terms.copy()
    for go_id in terms:
        if go_id not in go:
            terms_copy.remove(go_id)
    associations[gene] = terms_copy
associations = {key: value for key, value in associations.items()
                if value is not to_remove}
termcounts = TermCounts(go, associations)


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


# def dendrogram(graph):
# number of edges between each pair of nodes


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
    for i in go[go_term].children:
        if i.id in gene_set_of_interest:
            if i.id not in nodes_in_graph:
                graph.add_node('"{}"'.format(i.id), depth=i.depth,
                               level=i.level, GOname=i.name)
            # if go[go_term].id not in nodes_in_graph:
            graph.add_node('"{}"'.format(go[go_term].id),
                           depth=go[go_term].depth,
                           level=go[go_term].level,
                           GOname=go[go_term].name)
            edge = ('"{}"'.format(go[go_term].id), '"{}"'.format(i.id))
            graph.add_edge(*edge)
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


def check_term_list(list_of_terms):
    list_of_terms = set(list_of_terms)

    to_remove = set()
    for i in list_of_terms:
        # print(go[i].depth)
        if go[i].depth > 10:
            print("Depth of {} is greater than 10.".format(i),
                  "Removing from list for now as path to root is VAST")
            to_remove.add(i)
    list_of_terms = list_of_terms.difference(to_remove)
    list_of_terms = check_if_children(list_of_terms)
    combos = itertools.combinations(list_of_terms, 2)
    for i, j in combos:
        sim2 = semantic_similarity(i, j, go)
        dca = deepest_common_ancestor([i, j], go)

        sim = resnik_sim(i, j, go, termcounts)
        # this is the IC for biological_process,
        # subtracting it to normalize to 0
        sim -= 4.175976034663472
        # print(i, j, dca, go[i].depth, go[j].depth, go[dca].depth)
        # print(go[i].name, go[j].name, go[dca].name, sim, sim2)

        # remove if two terms of similar by distance in graph
        if sim2 > .5:
            if go[i].depth > go[j].depth:
                to_remove.add(i)
            else:
                to_remove.add(j)
        # remove if deepest common ancestor has at least 4 IC
        elif sim > 4:
            # if go[dca].depth < 3:
            #     continue
            if dca in to_remove or dca in list_of_terms:
                continue
            list_of_terms.add(dca)
            to_remove.add(i)
            to_remove.add(j)
            # print('Removing')
            # print('-' * 20)
            # print(i, j, dca, go[i].depth, go[j].depth, go[dca].depth)
            # print(go[i].name, go[j].name, go[dca].name, sim, sim2)
            # print('-' * 20)
        else:
            # I was thinking I could compare IC for regulators,
            # but it doesn't really provide much
            ic_i = ic(i, termcounts)
            ic_j = ic(j, termcounts)
            print("GO\t\t{}\t\t{}".format(go[i].name, go[j].name))
            print("IC\t\t{}\t\t\t{}".format(ic_i, ic_j))
            print("Score = {}\n".format(sim2))
            # if go[dca].depth > 3:
            #     print("Depth\t{}\t\t\t{}".format(go[i].depth, go[j].depth))
            #     print("DCA\tID\tDepth\tIC")
            #     print("{}\t{}\t{}\t{}\t".format(go[dca].name, dca, go[dca].depth, sim))

    # quit()
    print(to_remove)
    print('Number of term before = {}'.format(len(list_of_terms)))
    list_of_terms = list_of_terms.difference(to_remove)
    # list_of_terms = [i for i in list_of_terms if i not in to_remove]


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

    combined_graph = nx.DiGraph()

    node_list = []
    # list_of_terms = check_term_list(list_of_terms)
    for i in list_of_terms:
        g = path_to_root(i)
        node_list.append(set(g.nodes()))
        combined_graph.add_nodes_from(g.nodes(data=True))
        combined_graph.add_edges_from(g.edges())
    # count the number of times a node is in the list to find nodes that are
    # common to all terms
    term_counter = {}
    for i in node_list:
        for j in i:
            if j in term_counter:
                term_counter[j] += 1
            else:
                term_counter[j] = 1

    # find the nodes that are in all terms
    n_terms = len(list_of_terms)

    # For graphiv to recongnize ":" within a name, we must put it in a string
    # TODO : find a more clever way to go about this
    formatted_list_of_terms = []
    for i in list_of_terms:
        formatted_list_of_terms.append('"{}"'.format(i))
    common_nodes = set()
    for i in term_counter:
        print(i, term_counter[i], go[i.replace('"', '')].depth)
        if term_counter[i] == n_terms:
            common_nodes.add(i)

    if len(common_nodes) == 1:
        for i in common_nodes:
            level = combined_graph.node[i]['level']
            depth = combined_graph.node[i]['depth']
            name = combined_graph.node[i]['GOname']
        print("Warning : the most common ancestor is {} ({}) with a depth of "
              "{} and level of {}".format(name, i, depth, level))

    # Check to see if any of the provided terms are in all graphs
    # This means it was an ancestor node and should be removed from list
    to_remove = []
    for each in formatted_list_of_terms:
        if each in common_nodes:
            print("Warning : {} term is already in graph! Thus will treat as "
                  "upstream term and remove from list of terms"
                  " provided".format(each))
            to_remove.append(each)
    for each in to_remove:
        formatted_list_of_terms.remove(each)

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
        for each in formatted_list_of_terms:
            if nx.has_path(combined_graph, disjun_node, each):
                sh_paths = [p for p in
                            nx.all_shortest_paths(combined_graph, disjun_node,
                                                  each)]
                for path in sh_paths:
                    common_graph.add_path(path)

    # color source nodes
    for each in formatted_list_of_terms:
        common_graph.add_node(each, style='filled',
                              fillcolor='green', **combined_graph[each])
    gonames = nx.get_node_attributes(common_graph, 'GOname')
    print(gonames)
    nx.set_node_attributes(common_graph, 'label', gonames)
    export_to_dot(common_graph, 'common', view=False)
    return common_graph
