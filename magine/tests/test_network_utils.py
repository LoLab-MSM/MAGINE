import networkx as nx

import magine.networks.utils as utils

"""
def test_paint_over_time():
    graph = nx.DiGraph()
    graph.add_nodes_from(['a', 'b', 'c'])
    labels = ['1', '2', '3']
    colors = ['red', 'green', 'blue']
    nodes = [['a'], ['c', 'b'], ['a', 'b', 'c']]
    utils.paint_network_overtime(graph, nodes, colors,
                                 labels=labels,
                                 save_name='ok',
                                 create_gif=False)


def test_paint_over_time_up_down():
    graph = nx.DiGraph()
    graph.add_nodes_from(['a', 'b', 'c'])
    labels = ['1', '2', '3']
    up_nodes = [['a'], ['c', 'b'], ['a', 'b', 'c']]
    down_nodes = [['b'], ['a'], []]
    utils.paint_network_overtime_up_down(graph, list_up=up_nodes,
                                         list_down=down_nodes,
                                         color_down='blue',
                                         color_up='red',
                                         labels=labels,
                                         save_name='ok',
                                         create_gif=False)

"""


def test_delete_disconnected():
    g = nx.DiGraph()
    g.add_nodes_from(['a', 'b', 'c'])
    g.add_edge('a', 'b')
    large_g_connected = utils.delete_disconnected_network(g)
    large_g_connected.nodes()
    assert set(large_g_connected.nodes()) == {'b', 'a'}
    assert set(large_g_connected.edges()) == {('a', 'b')}


def test_add_attribute_to_network():

    g = nx.DiGraph()
    g.add_nodes_from(['a', 'b', 'c'])
    new_g = utils.add_attribute_to_network(g, ['a', 'b'], 'inTest', 'true')
    assert new_g.node['a'] == {'inTest': 'true'}

    assert new_g.node['c'] == {'inTest': 'false'}
