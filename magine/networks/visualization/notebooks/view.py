import json
import uuid

import networkx as nx
from IPython.display import display, HTML, Javascript

from magine.html_templates.cy_stypes import styles
from magine.html_templates.html_tools import env
from magine.networks.exporters import nx_to_json


def init():
    _urls = [
        "https://ajax.googleapis.com/ajax/libs/jquery/3.2.1/jquery.min.js",
        "https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js",
        "https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.1.4/cytoscape.js",
        "https://cdn.rawgit.com/cytoscape/cytoscape.js-cose-bilkent/1.6.5/cytoscape-cose-bilkent.js",
        "https://cdn.rawgit.com/cpettitt/dagre/v0.7.4/dist/dagre.min.js",
        "https://cdn.rawgit.com/cytoscape/cytoscape.js-dagre/1.5.0/cytoscape-dagre.js",
        "https://cdnjs.cloudflare.com/ajax/libs/jquery/2.2.4/jquery.min"
    ]

    # JS_LOADER_FILE = "loader.js"
    # path = os.path.abspath(os.path.dirname(__file__)) + "/" + JS_LOADER_FILE
    # js_loader = open(path).read()
    Javascript(lib=_urls)


# init()

layouts = {
    'breadthfirst': {
        'name': 'breadthfirst',
        'directed': 'true',
        'spacingFactor': .5
    },
    'cose': {
        'name': 'cose',
    },
    'cose-bilkent': {
        'name': 'cose-bilkent',
    }
}


def display_graph(graph, add_parent=False, layout='cose-bilkent',
                  background='#FFFFFF', height=700, width=100,
                  default_color='white', **layout_args):
    g_copy = graph.copy()
    if add_parent:
        g_copy = _add_parent_term(g_copy)
    for i in g_copy.nodes:
        if 'color' not in g_copy.node[i]:
            g_copy.node[i]['color'] = default_color
    d = nx_to_json(g_copy)
    d['background'] = background
    d['uuid'] = "cy" + str(uuid.uuid4())
    d['widget_width'] = str(width)
    d['widget_height'] = str(height)

    layout_opts = layouts[layout].copy()
    layout_opts.update(layout_args)

    d['layout_json'] = json.dumps(layout_opts)
    d['style_json'] = json.dumps(styles['default'])
    d['fitbutton'] = "fit" + str(uuid.uuid4())

    template = env.get_template('subgraph_2.html')
    display(HTML(template.render(d)))


def render_graph(graph, add_parent=False):
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
