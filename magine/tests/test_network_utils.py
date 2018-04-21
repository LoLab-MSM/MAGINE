import networkx as nx

import magine.networks.utils as utils


def test_paint_over_time():
    graph = nx.DiGraph()
    graph.add_nodes_from(['a', 'b', 'c'])
    labels = ['1', '2', '3']
    colors = ['red', 'green', 'blue']
    nodes = [['a'], ['c', 'b'], ['a', 'b', 'c']]
    utils.paint_network_overtime(graph, nodes, colors,
                                 labels=labels,
                                 save_name='ok',
                                 create_gif=True)


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
                                         create_gif=True)
