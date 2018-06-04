import matplotlib.pyplot as plt
import networkx as nx
import numpy as np


def render2(network):
    for (n, d) in network.nodes(data=True):
        if 'keggName' in d:
            del d["keggName"]
    node_colors = [n['color'] for _, n in network.nodes(data=True)]

    pos = nx.nx_pydot.graphviz_layout(network, prog='fdp')

    nx.draw_networkx(network, pos=pos,
                     with_labels=False,
                     node_color=node_colors,
                     node_size=500,
                     alpha=.4,
                     )
    label_pos = {}
    for n, p in pos.items():  # raise text positions
        label_pos[n] = [p[0], p[1] + 10]

    nx.draw_networkx_labels(network, label_pos, font_size=16)


def render_mpl(network, layout='dot'):
    for (n, d) in network.nodes(data=True):
        if 'keggName' in d:
            del d["keggName"]

    if layout in ['dot', 'neato', 'fdp']:
        pos = nx.nx_pydot.graphviz_layout(network, prog=layout, )
    else:
        assert layout in ['circular_layout', 'random_layout', 'shell_layout',
                          'spring_layout', 'spectral_layout',
                          'fruchterman_reingold_layout']

        pos = nx.drawing.layout.spring_layout(network)

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
