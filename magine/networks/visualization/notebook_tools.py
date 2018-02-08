import json

import numpy as np
from IPython.display import HTML, Javascript
from IPython.display import SVG, display

import magine.ontology.enrichment_tools as et
from magine.html_templates.html_tools import env
from magine.networks.network_subgraphs import NetworkSubgraphs
from magine.networks.ontology_network import OntologyNetworkGenerator
from magine.networks.visualization.cytoscapejs_tools import viewer as cyjs
from magine.networks.visualization.igraph_tools import create_figure
from magine.networks.visualization.util_networkx import from_networkx

Javascript("https://ajax.googleapis.com/ajax/libs/jquery/3.2.1/jquery.min.js")
Javascript(
    "https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js")
Javascript(
    "https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.1.4/cytoscape.js")
Javascript(
    "https://cdn.rawgit.com/cytoscape/cytoscape.js-cose-bilkent/1.6.5/cytoscape-cose-bilkent.js")

Javascript("https://cdn.rawgit.com/cpettitt/dagre/v0.7.4/dist/dagre.min.js")
Javascript("https://cdn.rawgit.com/cytoscape/cytoscape.js-dagre/1.5.0/cytoscape-dagre.js")


def create_subnetwork(terms, df, network, save_name=None, draw_png=False,
                      cytoscape_js=False):
    df = df[df['term_name'].isin(terms)].copy()
    df['combined_score'] = np.abs(df['combined_score'])
    df['combined_score'] = np.log2(df['combined_score'])
    df.loc[df['combined_score'] > 150, 'combined_score'] = 150

    term_dict = dict()
    label_dict = dict()
    for i in terms:
        genes = set(et.term_to_genes(df, i))
        term_dict[i] = genes
        label_dict[i] = i
        print(i, len(genes))
    all_genes = set()
    for i, j in term_dict.items():
        all_genes.update(j)
    ong = OntologyNetworkGenerator(molecular_network=network)
    print("Looking for direct edges")
    term_g, molecular_g = ong.create_network_from_list(
        terms, term_dict, label_dict, save_name=save_name, draw=draw_png)

    if cytoscape_js:
        display_graph(term_g)
        display_graph(molecular_g)
    return term_g, molecular_g


def display_graph(graph, add_parent=False, display_format='cytoscape'):
    g_copy = graph.copy()
    if display_format == 'cytoscape':
        if add_parent:
            new_nodes = set()
            for i, data in graph.nodes(data=True):
                if 'termName' in data:
                    term = data['termName']
                    g_copy.node[i]['parent'] = term
                    new_nodes.add(term)
            for each in new_nodes:
                g_copy.add_node(each, )
        g_cyjs = from_networkx(g_copy)
        cyjs.render(g_cyjs, style='Directed', layout_algorithm='cose-bilkent')
    elif display_format == 'igraph':
        display(SVG(create_figure(g_copy)))


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
    import networkx as nx

    g = nx.DiGraph()
    g.add_edge('BAX', 'CASP3')
    g.add_edge('BCL2', 'BAX')
    g.add_edge('BAX', 'CASP3')
    find_neighbors(g, 'BAX', True, True)
