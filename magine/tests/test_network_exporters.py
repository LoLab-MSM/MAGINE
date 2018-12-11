import json

import igraph
import networkx as nx
import pydot

import magine.networks.exporters as exporters


def test_nx_to_igraph():
    g = nx.DiGraph()
    g.add_edge('a', 'b')
    i_g = exporters.nx_to_igraph(g)
    assert isinstance(i_g, igraph.Graph)


def ordered(obj):
    if isinstance(obj, dict):
        return sorted((k, ordered(v)) for k, v in obj.items())
    if isinstance(obj, list):
        return sorted(ordered(x) for x in obj)
    else:
        return obj


def test_nx_to_json():
    g = nx.DiGraph()
    g.add_edge('a', 'b')
    json_g = exporters.nx_to_json(g)
    answer = {'edges': '[{"data": {"source": "a", "target": "b"}}]',
              'nodes': '[{"data": {"id": "a", "name": "a"}},'
                       ' {"data": {"id": "b", "name": "b"}}]'}
    a_edges = json.loads(json_g['edges'])
    b_edges = json.loads(answer['edges'])
    a_nodes = json.loads(json_g['nodes'])
    b_nodes = json.loads(answer['nodes'])

    assert ordered(a_edges) == ordered(b_edges)
    assert ordered(a_nodes) == ordered(b_nodes)


def test_nx_to_dot():
    g = nx.DiGraph()
    g.add_edge('a', 'b')
    dot_g = exporters.nx_to_dot(g)
    assert isinstance(dot_g, pydot.Dot)
