from locale import getpreferredencoding

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pydot
from networkx.drawing.layout import rescale_layout

from magine.networks.exporters import nx_to_dot

try:
    basestring
except NameError:
    basestring = str
    unicode = str


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
    tmp_g = G.copy()
    tmp_g.graph.update({'K': '1'})
    tmp_g.graph.update({'repulsiveforce ': '1'})
    tmp_g.graph.update({'overlap ': 'false'})
    tmp_g.graph.update({'splines ': 'true'})
    dot_fmt = nx_to_dot(tmp_g).create(prog=prog, format='dot', encoding='utf8')
    dot_fmt = unicode(dot_fmt, encoding=getpreferredencoding())
    net = pydot.graph_from_dot_data(dot_fmt)[0]
    node_pos = {}

    for n in tmp_g.nodes():
        pydot_node = pydot.Node(n).get_name()
        node = net.get_node(pydot_node)

        if isinstance(node, list):
            node = node[0]
        pos = node.get_pos()[1:-1]  # strip leading and trailing double quotes
        if pos is not None:
            xx, yy = pos.split(",")
            node_pos[n] = (float(xx), float(yy))

    return node_pos


def draw_mpl(network, layout='dot', scale=1, node_size=500, font_size=16):
    """ Draw network using networkx and matplotlib

    Parameters
    ----------
    network : nx.DiGraph
    layout : str

    """
    for (n, d) in network.nodes(data=True):
        if 'keggName' in d:
            del d["keggName"]
    mpl_layouts = ['circular_layout', 'random_layout', 'shell_layout',
                   'spring_layout', 'spectral_layout', 'kamada_kawai_layout',
                   'fruchterman_reingold_layout']
    dot_layouts = ['dot', 'neato', 'fdp', 'twopi', 'circo', 'sfdp']
    if layout in dot_layouts:
        pos = pydot_layout(network, prog=layout, )
    else:
        if layout not in mpl_layouts:
            raise AssertionError("Please provide valid layout option:{}"
                                 "".format(mpl_layouts + dot_layouts))
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
        elif layout == 'kamada_kawai_layout':
            pos = nx.drawing.layout.kamada_kawai_layout(network)

    # some layout algorithms (graphviz ones) can provide large position values
    # normalize the positions to be from 0 to 1
    # this allows us to add space between pode and label
    positions = np.array([p for _, p in pos.items()])

    positions = rescale_layout(positions, scale=1)
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
    fig = plt.figure(figsize=(4 * scale, 4 * scale))
    ax = fig.add_subplot()
    nx.draw_networkx(network, pos=new_pos,
                     with_labels=False,
                     node_color=node_colors,
                     node_size=node_size,
                     alpha=.4, ax=ax
                     )
    nx.draw_networkx_labels(network, label_pos, font_size=font_size, ax=ax)
    # plt.axis('equal')
    plt.axis('off')
    # plt.xticks([])
    # plt.yticks([])
    # plt.tight_layout()
    return fig


if __name__ == '__main__':
    x = nx.DiGraph()
    x.add_edge('A', 'B')
    x.add_edge('A', 'C')
    x.add_edge('C', 'D')
    x.add_edge('C', 'E')
    x.add_edge('E', 'A')
    # x = pydot_layout(x)
    draw_mpl(x, 'circo', scale=1)
    draw_mpl(x, 'circo', scale=2)
    plt.show()
