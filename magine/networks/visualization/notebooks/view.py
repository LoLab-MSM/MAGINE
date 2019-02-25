import json
import uuid

import networkx as nx
from IPython.display import display, HTML

from magine.html_templates.html_tools import env
from magine.networks.exporters import nx_to_json
from magine.networks.visualization.notebooks.cy_stypes import styles, layouts


def display_graph(graph, add_parent=False, layout='cose-bilkent',
                  background='#FFFFFF', height=700, width=100,
                  default_color='white', **layout_args):
    g_copy = graph.copy()
    if add_parent:
        g_copy = _add_parent_term(g_copy)

    _set_node_color(g_copy, default_color)

    d = nx_to_json(g_copy)
    d['background'] = background
    d['uuid'] = "cy" + str(uuid.uuid4())
    d['widget_width'] = str(width)
    d['widget_height'] = str(height)

    if layout not in layout:
        layout_opts = dict()
    else:
        layout_opts = layouts[layout].copy()
    layout_opts.update(layout_args)

    d['layout_json'] = json.dumps(layout_opts)
    d['style_json'] = json.dumps(styles['default'])
    d['fitbutton'] = "fit" + str(uuid.uuid4())

    template = env.get_template('subgraph_2.html')
    display(HTML(template.render(d)))


def render_graph(graph, add_parent=False, default_color='white'):
    g_copy = graph.copy()
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


def render_graph2(graph, add_parent=False):
    g_copy = graph.copy()
    if add_parent:
        g_copy = _add_parent_term(g_copy)

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
    subgraph_html = env.get_template('subgraph_3.html')

    with open(fname_temp, 'w') as f:
        f.write(subgraph_html.render(d))

    display(HTML(subgraph_html.render(d)))


def _set_node_color(net, color):
    for i in net.nodes:
        if 'color' not in net.node[i]:
            net.node[i]['color'] = color


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
    g.add_edge('a', 'b')

    render_graph(g)
