import igraph
import networkx as nx
import pydotplus

import magine.networks.exporters as exporters


def test_nx_to_igraph():
    g = nx.DiGraph()
    g.add_edge('a', 'b')
    i_g = exporters.nx_to_igraph(g)
    assert isinstance(i_g, igraph.Graph)


def test_nx_to_jsonh():
    g = nx.DiGraph()
    g.add_edge('a', 'b')
    json_g = exporters.nx_to_json(g)
    answer = {'edges': '[{"data": {"source": "a", "target": "b"}}]',
              'nodes': '[{"data": {"id": "a", "name": "a"}}, {"data": {"id": "b", "name": "b"}}]'}

    assert json_g['edges'] == answer['edges']
    assert json_g['nodes'] == answer['nodes']


def test_nx_to_dot():
    g = nx.DiGraph()
    g.add_edge('a', 'b')
    dot_g = exporters.nx_to_dot(g)
    assert isinstance(dot_g, pydotplus.Dot)
