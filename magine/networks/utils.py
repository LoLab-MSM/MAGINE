import json

import igraph as ig
import networkx as nx


def networkx_to_igraph(network):
    nx.write_gml(network, 'test.gml')
    igraph_network = ig.read('test.gml')
    return igraph_network


def from_networkx(g):
    """

    Parameters
    ----------
    g : nx.DiGraph

    Returns
    -------
    dict
    """
    graph = dict(data={}, elements=dict(nodes=[], edges=[]))

    graph['data'] = {col: g.graph[col] for col in g.graph.keys() if col != 0}
    for node_id in g.nodes_iter():
        node = g.node[node_id]
        new_node = {}
        node_columns = node.keys()
        data = {col: node[col] for col in node_columns if col != 0}

        if 'position' in node.keys():
            new_node['position'] = node['position']
        data['id'] = str(node_id)
        data['name'] = str(node_id)
        new_node['data'] = data
        graph['elements']['nodes'].append(new_node)

    for source, target, data in g.edges_iter(data=True):
        data['source'] = str(source)
        data['target'] = str(target)
        graph['elements']['edges'].append(data)

    return graph


def nx_to_json(network):
    graph = from_networkx(network)
    nodes = graph['elements']['nodes']
    edges = graph['elements']['edges']
    data = {
        'nodes': json.dumps(nodes),
        'edges': json.dumps(edges),
    }
    return data