import json
import tempfile

import networkx as nx


def nx_to_igraph(network):
    try:
        import igraph
    except ImportError:
        raise ImportError('requires igraph ',
                          'http://pygraphviz.github.io/')
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

    for node_id in network.nodes():
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

    return {'nodes': json.dumps(nodes), 'edges': json.dumps(edges)}


def nx_to_dot(graph):
    try:
        import pygraphviz
    except ImportError:
        raise ImportError('requires pygraphviz ',
                          'http://pygraphviz.github.io/')

    graph = format_to_directions(graph)

    directed = graph.is_directed()
    strict = graph.number_of_selfloops() == 0 and not graph.is_multigraph()
    new_g = pygraphviz.AGraph(name=graph.name,
                              strict=strict,
                              directed=directed,
                              encoding='utf8')
    # default graph attributes
    new_g.graph_attr.update(graph.graph.get('graph', {}))
    new_g.node_attr.update(graph.graph.get('node', {}))
    new_g.edge_attr.update(graph.graph.get('edge', {}))

    # add nodes
    for n, data in graph.nodes(data=True):
        new_g.add_node(n, **data)

    # loop over edges
    for u, v, data in graph.edges(data=True):
        new_g.add_edge(u, v, **dict((k, str(v)) for k, v in data.items()))
    return new_g


def export_to_dot(graph, save_name, image_format='png', engine='dot',
                  dpi=100, concentrate=False):
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
    dpi: int
        resolution of figure
    concentrate : bool
        Concentrate bi-edges into one edge

    Returns
    -------

    """

    py_dot = nx_to_dot(graph)
    if save_name.endswith('.png'):
        basename = save_name[:-4]
    else:
        basename = save_name
    py_dot.write('{}.dot'.format(basename))
    arg = '-Gdpi={} -Gconcentrate={}'.format(
        dpi, 'true' if concentrate else 'false'
    )
    py_dot.draw('{}.{}'.format(basename, image_format), prog=engine, args=arg)


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
