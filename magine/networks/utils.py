import igraph as ig
import networkx as nx

ID = 'id'
NAME = 'name'
DATA = 'data'
ELEMENTS = 'elements'
NODES = 'nodes'
EDGES = 'edges'
SOURCE = 'source'
TARGET = 'target'
DEF_SCALE = 100


def networkx_to_igraph(network):
    nx.write_gml(network, 'test.gml')
    igraph_network = ig.read('test.gml')
    return igraph_network


def __map_table_data(columns, graph_obj):
    # data = {}
    # data[col] = graph_obj[col]
    # data =
    # for col in columns:
    #     if col == 0:
    #         break
    return {col: graph_obj[col] for col in columns if col != 0}


def __create_node(node, node_id):
    new_node = {}
    node_columns = node.keys()
    data = __map_table_data(node_columns, node)
    # Override special keys
    data[ID] = str(node_id)
    data[NAME] = str(node_id)

    if 'position' in node.keys():
        new_node['position'] = node['position']

    new_node[DATA] = data
    return new_node


def __build_multi_edge(edge_tuple, g):
    source = edge_tuple[0]
    target = edge_tuple[1]
    key = edge_tuple[2]
    data = edge_tuple[3]

    data['source'] = str(source)
    data['target'] = str(target)
    data['interaction'] = str(key)
    return {DATA: data}


def __build_edge(edge_tuple, g):
    source = edge_tuple[0]
    target = edge_tuple[1]
    data = edge_tuple[2]

    data['source'] = str(source)
    data['target'] = str(target)
    return {DATA: data}


def __build_empty_graph():
    return {
        DATA: {},
        ELEMENTS: {
            NODES: [],
            EDGES: []
        }
    }


def from_networkx(g):
    """

    Parameters
    ----------
    g : nx.DiGraph

    Returns
    -------
    dict
    """
    cygraph = __build_empty_graph()

    # Map network table data
    cygraph[DATA] = __map_table_data(g.graph.keys(), g.graph)

    for node_id in g.nodes_iter():
        new_node = __create_node(g.node[node_id], node_id)
        cygraph['elements']['nodes'].append(new_node)

    for edge in g.edges_iter(data=True):
        cygraph['elements']['edges'].append(__build_edge(edge, g))

    return cygraph
