import networkx as nx
from nose.tools import ok_

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
    ok_(set(large_g_connected.nodes()) == {'b', 'a'})
    ok_(set(large_g_connected.edges()) == {('a', 'b')})


def test_add_attribute_to_network():
    g = nx.DiGraph()
    g.add_nodes_from(['a', 'b', 'c'])
    new_g = utils.add_attribute_to_network(g, ['a', 'b'], 'inTest', 'true')
    ok_(new_g.node['a'] == {'inTest': 'true'})
    ok_(new_g.node['c'] == {'inTest': 'false'})


def test_trim_source_sink():
    graph = nx.DiGraph()
    graph.add_edge('A', 'B')
    graph.add_edge('B', 'C')
    graph.add_edge('B', 'B')
    graph.add_edge('D', 'C')
    graph.add_edge('C', 'N')
    trimmed = utils.trim_sink_source_nodes(graph, ['A', 'B', 'D'],
                                           remove_self_edge=True)
    ok_(set(trimmed.edges()) == {('B', 'C'), ('A', 'B'), ('D', 'C')})
    return
    import matplotlib.pyplot as plt
    plt.subplot(121)
    nx.draw(graph,
            pos=nx.nx_pydot.graphviz_layout(graph, 'dot'),
            with_labels=True)
    plt.subplot(122)
    nx.draw(trimmed,
            pos=nx.nx_pydot.graphviz_layout(graph, 'dot'),
            with_labels=True)
    plt.show()

# test_trim_source_sink()
