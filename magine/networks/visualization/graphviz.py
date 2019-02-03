import networkx as nx
from IPython.display import Image, display

from magine.networks.exporters import nx_to_dot


def draw_graphviz(network, layout='dot', width=100, save_name=None):
    """

    Parameters
    ----------
    network : nx.DiGraph
    layout : str
        Which graphviz engine to use to layout graph

    Returns
    -------

    """
    net_copy = network.copy()
    net_copy.graph.update({'K': '1'})
    net_copy.graph.update({'repulsiveforce ': '1'})
    net_copy.graph.update({'overlap ': 'false'})
    net_copy.graph.update({'splines ': 'true'})
    if save_name is not None:
        nx_to_dot(net_copy).write(save_name + '.png', format='png',
                                  prog=layout)
        display(Image(save_name + '.png', width=width))
    else:
        img = nx_to_dot(net_copy).create(format='png', prog=layout)
        display(Image(img, width=width))
