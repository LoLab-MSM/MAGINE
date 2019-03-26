import json
import sys
import tempfile

import networkx as nx
import pydot


def nx_to_igraph(network):
    try:
        import igraph
    except ImportError:
        raise ImportError('requires igraph ')
    file_descriptor, file_path = tempfile.mkstemp(suffix='.gml')

    nx.write_gml(network, file_path)
    igraph_network = igraph.read(file_path)
    return igraph_network


def nx_to_json(network):
    """
    Convert from networkx Digraph to cytoscape json format
    Parameters
    ----------
    network : nx.DiGraph

    Returns
    -------
    dict
    """
    nodes, edges = [], []

    for node_id in network.nodes:
        node = network.node[node_id]
        new_node = {}
        node_columns = node.keys()
        data = {col: node[col] for col in node_columns if col != 0}
        if 'position' in node.keys():
            new_node['position'] = node['position']
        data['id'] = str(node_id)
        data['name'] = data['id']
        new_node['data'] = data
        nodes.append(new_node)

    for source, target, data in network.edges(data=True):
        data['source'] = str(source)
        data['target'] = str(target)
        edges.append({'data': data})

    return {'nodes': json.dumps(nodes),
            'edges': json.dumps(edges)}


def nx_to_dot(graph):
    graph = format_to_directions(graph)
    # set Graphviz graph type
    if graph.is_directed():
        graph_type = 'digraph'
    else:
        graph_type = 'graph'
    strict = graph.number_of_selfloops() == 0 and not graph.is_multigraph()

    name = graph.name
    if name is not '':
        name = make_str(name)
    args = {make_str(i): make_str(graph.graph[i]) for i in graph.graph}
    g = pydot.Dot(name, graph_type=graph_type, strict=strict, **args)
    try:
        g.set_node_defaults(**{make_str(i): make_str(graph.graph[i]) for i in
                               graph.graph['node']})
    except KeyError:
        pass
    try:
        g.set_edge_defaults(**{make_str(i): make_str(graph.graph[i]) for i in
                               graph.graph['edge']})
    except KeyError:
        pass

    for n, data in graph.nodes(data=True):
        str_data = dict((make_str(k), make_str(v)) for k, v in data.items())
        g.add_node(pydot.Node(make_str(n), **str_data))

    for u, v, data in graph.edges(data=True):
        str_data = dict((make_str(k), make_str(v)) for k, v in data.items()
                        if k != 'key')
        g.add_edge(pydot.Edge(make_str(u), make_str(v), **str_data))
    return g


def export_to_dot(graph, save_name, image_format='png', engine='dot'):
    """
    Converts networkx graph to graphviz dot

    Parameters
    ----------
    graph : networkx.DiGraph
    save_name : str
        name of file to save
    image_format : str
        format of output( pdf, png, svg)
    engine : str
        graphviz engine
            dot, twopi,

    Returns
    -------

    """

    py_dot = nx_to_dot(graph)
    if save_name.endswith('.png'):
        basename = save_name[:-4]
    else:
        basename = save_name
    py_dot.write('{}.{}'.format(basename, image_format), prog=engine,
                 format=image_format)


def check_graphviz(network):
    if isinstance(network, nx.DiGraph):
        network = format_to_directions(network)

    return nx_to_dot(network)


def format_to_directions(network):
    activators = ['activate', 'expression', 'phosphorylate']
    inhibitors = [
        'inhibit', 'repression', 'dephosphorylate', 'deubiquitinate',
        'ubiquitinate'
    ]
    physical_contact = ['binding', 'dissociation', 'stateChange',
                        'compound', 'glycosylation']
    indirect_types = ['indirect']
    for source, target, data in network.edges(data=True):
        if 'interactionType' in data:
            edge_type = data['interactionType']
            for j in activators:
                if j in edge_type:
                    network[source][target]['arrowhead'] = 'normal'
            for j in inhibitors:
                if j in edge_type:
                    network[source][target]['arrowhead'] = 'tee'

            for j in physical_contact:
                if j in edge_type:
                    network[source][target]['dir'] = 'both'
                    network[source][target]['arrowtail'] = 'diamond'
                    network[source][target]['arrowhead'] = 'diamond'

            for j in indirect_types:
                if j in edge_type:
                    network[source][target]['arrowhead'] = 'diamond'
                    network[source][target]['style'] = 'dashed'
    return network


if sys.version_info[0] == 2:
    def make_str(x):
        """Return the string representation of t."""
        if isinstance(x, unicode):
            return x
        else:
            return unicode(str(x), 'unicode-escape')
else:
    def make_str(x):
        """Return the string representation of t."""
        return str(x)
