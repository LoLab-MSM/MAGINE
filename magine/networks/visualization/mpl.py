import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pydot

from magine.networks.exporters import nx_to_dot


def pydot_layout(G, prog='neato'):
    """Create node positions using :mod:`pydot` and Graphviz.

    Parameters
    --------
    G : Graph
        NetworkX graph to be laid out.
    prog : optional[str]
        Basename of the GraphViz command with which to layout this graph.
        Defaults to `neato`, the default GraphViz command for undirected graphs.

    Returns
    --------
    dict
        Dictionary of positions keyed by node.

    """

    net = nx_to_dot(G)
    dot_file = pydot.graph_from_dot_data(net.create(prog=prog, format='dot'))
    net = dot_file[0]
    node_pos = {}
    for n in G.nodes():
        pydot_node = pydot.Node(n).get_name()
        node = net.get_node(pydot_node)

        if isinstance(node, list):
            node = node[0]
        pos = node.get_pos()[1:-1]  # strip leading and trailing double quotes
        if pos is not None:
            xx, yy = pos.split(",")
            node_pos[n] = (float(xx), float(yy))

    return node_pos


def render_mpl(network, layout='dot'):
    """ Draw network using networkx and matplotlib

    Parameters
    ----------
    network : nx.DiGraph
    layout : str

    """
    for (n, d) in network.nodes(data=True):
        if 'keggName' in d:
            del d["keggName"]

    if layout in ['dot', 'neato', 'fdp']:
        pos = pydot_layout(network, prog=layout, )
    else:
        assert layout in ['circular_layout', 'random_layout', 'shell_layout',
                          'spring_layout', 'spectral_layout',
                          'fruchterman_reingold_layout']
        if layout == 'spring_layout':
            pos = nx.drawing.layout.spring_layout(network)
        elif layout == 'random_layout':
            pos = nx.drawing.layout.random_layout(network)
        elif layout == 'shell_layout':
            pos = nx.drawing.layout.shell_layout(network)
        elif layout == 'spectral_layout':
            pos = nx.drawing.layout.spectral_layout(network)
        elif layout == 'fruchterman_reingold_layout':
            pos = nx.drawing.layout.fruchterman_reingold_layout(network)
        elif layout == 'circular_layout':
            pos = nx.drawing.layout.circular_layout(network)

    # some layout algorithms (graphviz ones) can provide large position values
    # normalize the positions to be from 0 to 1
    # this allows us to add space between pode and label
    positions = np.array([p for _, p in pos.items()])
    x_min = positions[:, 0].min()
    x_max = positions[:, 0].max()
    y_min = positions[:, 0].min()
    y_max = positions[:, 0].max()

    positions[:, 0] = (positions[:, 0] - x_min) / (x_max - x_min)
    positions[:, 1] = (positions[:, 1] - y_min) / (y_max - y_min)

    new_pos, label_pos = {}, {}
    for num, n in enumerate(pos):
        x, y = positions[num]
        new_pos[n] = [x, y]
        # raise text positions
        label_pos[n] = [x, y + .05]
    try:
        node_colors = [n['color'] for _, n in network.nodes(data=True)]
    except:
        node_colors = ['red' for _ in network.nodes()]

    nx.draw_networkx(network, pos=new_pos,
                     with_labels=False,
                     node_color=node_colors,
                     node_size=500,
                     alpha=.4,
                     )
    nx.draw_networkx_labels(network, label_pos, font_size=16)
    plt.axis('equal')
    plt.xticks([])
    plt.yticks([])
    plt.tight_layout()


if __name__ == '__main__':
    x = nx.DiGraph()
    x.add_edge('A', 'B')
    x = pydot_layout(x)
