import json
try:
    import igraph as ig
except ImportError:
    ig = None
import networkx as nx


def networkx_to_igraph(network):
    try:
        import igraph as ig
    except ImportError:
        raise ImportError('requires igraph ',
                          'http://pygraphviz.github.io/')

    nx.write_gml(network, 'test.gml')
    igraph_network = ig.read('test.gml')
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

    for node_id in network.nodes_iter():
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

    for source, target, data in network.edges_iter(data=True):
        data['source'] = str(source)
        data['target'] = str(target)
        edges.append({'data': data})

    return {'nodes': json.dumps(nodes), 'edges': json.dumps(edges)}


def export_to_dot(graph, save_name, image_format='png', engine='dot',
                  dpi=200, concentrate=False):
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
    try:
        py_dot = nx.nx_agraph.to_agraph(graph)
        py_dot.write('{}.dot'.format(save_name))
        arg = '-Gdpi={} -Gconcentrate={}'.format(
            dpi, 'true' if concentrate else 'false'
        )
        py_dot.draw('{}.{}'.format(save_name, image_format), prog=engine,
                    args=arg)
    except ImportError:
        print("No pygraphivz installed")


def nx_to_dot(graph):

    try:
        import pygraphviz
    except ImportError:
        raise ImportError('requires pygraphviz ',
                          'http://pygraphviz.github.io/')
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
    for u, v, data in graph.edges_iter(data=True):
        new_g.add_edge(u, v, **dict((k, str(v)) for k, v in data.items()))
    return new_g