import json
import uuid

import networkx as nx
import numpy as np
from IPython.display import display, HTML

from magine.html_templates.html_tools import env
from magine.networks.exporters import nx_to_json
from .cyjs_options import layouts, styles


def draw_cyjs(graph, add_parent=False, layout='cose-bilkent',
              bg_color='white', default_node_color='white', **layout_args):
    """ Renders a graph using cytoscape.js

    Parameters
    ----------
    graph : nx.DiGraph
    add_parent : bool
        Group together nodes that share 'termName'.
    layout : str
        Layout
    bg_color
    default_node_color : str
        Will only be used if color is not already set for node
    layout_args : dict
        Any additional arguments for use in cytoscape.js. Examples can be found
        in cyjs_options.py

    Returns
    -------

    """
    # Dont think users need these
    height = 700
    width = 100
    g_copy = graph.copy()
    if add_parent:
        g_copy = _add_parent_term(g_copy)

    _scale_edges(g_copy)
    _set_node_color(g_copy, default_node_color)

    d = nx_to_json(g_copy)
    d['background'] = bg_color
    d['uuid'] = "cy" + str(uuid.uuid4())
    d['widget_width'] = str(width)
    d['widget_height'] = str(height)

    if layout not in layouts:
        raise Exception(
            "layout {} is not in {}".format(layout, layouts.keys())
        )
    else:
        layout_opts = layouts[layout].copy()
        layout_opts.update(layout_args)

    d['layout_json'] = json.dumps(layout_opts)
    d['style_json'] = json.dumps(styles['default'])
    d['fitbutton'] = "fit{}".format(uuid.uuid4())

    template = env.get_template('subgraph_2.html')
    display(HTML(template.render(d)))
 

def render_graph(graph, add_parent=False, default_color='white'):
    g_copy = graph.copy()

    _scale_edges(g_copy)

    if add_parent:
        g_copy = _add_parent_term(g_copy)
    _set_node_color(g_copy, default_color)
    d = nx_to_json(g_copy)
    u_name = "cy{}".format(uuid.uuid4())
    d['uuid'] = u_name
    d['style_json'] = json.dumps(styles['default'])

    edge_types = set()
    for i, j, data in graph.edges(data=True):
        if 'interactionType' in data:
            for e in data['interactionType'].split('|'):
                edge_types.add(e)

    d['edge_list'] = list(edge_types)

    fname_temp = '{}.html'.format(u_name)
    subgraph_html = env.get_template('subgraph.html')
    template = env.get_template('main_view.html')

    with open(fname_temp, 'w') as f:
        f.write(subgraph_html.render(d))

    display(HTML(template.render(name=u_name, filename=fname_temp)))


def _set_node_color(net, color):
    for i in net.nodes:
        if 'color' not in net.node[i]:
            net.node[i]['color'] = color


def _scale_edges(net):
    # store in a weird format. An EdgeView should maintain same ordering when
    # enumerating, but just in case going to keep index pairs for placing width
    vals = [[[i, j], val] for i, j, val in net.edges.data('weight', default=5)]
    pairs = [pair for pair, _ in vals]
    widths = np.array([val for _, val in vals])

    # normalize between 1, 100
    edge_width = np.interp(widths, (widths.min(), widths.max()), (3, 10))
    for width, (i, j) in zip(edge_width, pairs):
        net.edges[i, j]['width'] = width


def _add_parent_term(graph):
    g_copy = graph.copy()
    new_nodes = set()
    for i, data in graph.nodes(data=True):
        if 'termName' in data:
            term = data['termName']
            g_copy.node[i]['parent'] = term
            new_nodes.add(term)
    for each in new_nodes:
        g_copy.add_node(each)
    return g_copy


if __name__ == '__main__':
    g = nx.DiGraph()
    g.add_edge('a', 'b', weight=10)
    g.add_edge('b', 'c', weight=2)
    g.add_edge('b', 'a')
    g.add_edge('b', 'p', weight=100)
    _scale_edges(g)

    # render_graph(g)
