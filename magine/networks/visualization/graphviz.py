from IPython.display import Image, display

from magine.networks.exporters import nx_to_dot


def draw_graphviz(network, layout='dot'):
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
    display(Image(nx_to_dot(net_copy).create(format='png', prog=layout)))
