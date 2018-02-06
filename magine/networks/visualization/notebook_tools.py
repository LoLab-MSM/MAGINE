import json
from IPython.display import HTML, Javascript
# from gui.network_functions import neighbors
from magine.networks.visualization.util_networkx import from_networkx
from magine.html_templates.html_tools import env
from magine.networks.network_subgraphs import NetworkSubgraphs


Javascript("https://ajax.googleapis.com/ajax/libs/jquery/3.2.1/jquery.min.js")
Javascript(
    "https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js")
Javascript(
    "https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.1.4/cytoscape.js")
Javascript(
    "https://cdn.rawgit.com/cytoscape/cytoscape.js-cose-bilkent/1.6.5/cytoscape-cose-bilkent.js")

Javascript("https://cdn.rawgit.com/cpettitt/dagre/v0.7.4/dist/dagre.min.js")
Javascript("https://cdn.rawgit.com/cytoscape/cytoscape.js-dagre/1.5.0/cytoscape-dagre.js")


def find_neighbors(g, start, up_stream=True, down_stream=True, max_dist=1):

    subgraph_gen = NetworkSubgraphs(g)
    sg = subgraph_gen.neighbors(start, up_stream, down_stream, max_dist)
    return render_graph(sg)


def render_graph(graph):
    graph = from_networkx(graph)
    d = _json_graph(graph)
    # """

    prefix = 'style_main'
    fname_temp = '{}.html'.format(prefix)

    with open('{}.html'.format(prefix), 'w') as f:
        subgraph_html = env.get_template('subgraph.html')
        f.write(subgraph_html.render(d))

    template = env.get_template('main_view.html')
    outfile = template.render(name=prefix, filename=fname_temp)

    # with open('test.html', 'w') as f:
    #     f.write(outfile)

    return HTML(outfile)


def _json_graph(graph):
    nodes = graph['elements']['nodes']
    edges = graph['elements']['edges']
    data = {
        'nodes': json.dumps(nodes),
        'edges': json.dumps(edges),
    }
    return data


if __name__ == '__main__':
    find_neighbors('BAX', True, True)
