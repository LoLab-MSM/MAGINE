import os

import graphviz as gv
import networkx as nx
from goatools import obo_parser

os.environ[
    "PATH"] += os.pathsep + "C:\Users\James Pino\Miniconda2\envs\MAGINE\Library\\bin\graphviz"

# must download this from Gene ontology
# TODO point to the ontology and annotations from orange.biocontrib
short_name = 'gene_ontology_edit.obo'

short_path = os.path.join(os.path.dirname(__file__),
                          'ontology_files',
                          'go-basic.obo')
go = obo_parser.GODag(short_path)


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
            if i.id not in nodes_in_graph:
                graph.add_node('"{}"'.format(go[go_term].id),
                               depth=go[go_term].depth,
                               level=go[go_term].level,
                               GOname=go[go_term].name)
            edge = ('"{}"'.format(go[go_term].id), '"{}"'.format(i.id))
            graph.add_edge(*edge)
            add_children(i.id, graph, gene_set_of_interest)


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
                              fillcolor='green')

    export_to_dot(common_graph, 'common', view=False)
    return common_graph


def export_to_dot(graph, save_name, format='png', engine='dot', view=False):
    """
    Converts networkx graph to graphviz dot

    Parameters
    ----------
    graph : networkx.DiGraph
    save_name : str
        name of file to save
    format : str
        format of output( pdf, png, svg)
    engine : str
        graphviz engine
            dot, twopi,
    view : bool
        open up the rendered image

    Returns
    -------

    """
    dot_string = nx.nx_pydot.to_pydot(graph).to_string()
    dot_graph = gv.Source(dot_string, format=format, engine=engine)
    dot_graph.render('{}'.format(save_name), view=view, cleanup=True)
    dot_graph.save('{}.dot'.format(save_name))
