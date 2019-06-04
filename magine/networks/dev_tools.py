import itertools

import networkx as nx
import pathos.multiprocessing as mp


def remove_unmeasured_nodes(graph, measured):
    new_g = graph.copy()
    edge_info_dict = dict()
    for i, j, data in graph.edges(data=True):
        edge_info_dict[(i, j)] = data
    nodes = set(graph.nodes())
    include = set(measured)
    include.intersection_update(nodes)

    # for i, j in :
    def find(d):
        i, j = d
        paths = []
        if nx.has_path(graph, i, j):
            for p in nx.all_shortest_paths(graph, i, j):
                path = []
                label = []
                for n in p:
                    if n in include:
                        path.append(n)
                    else:
                        label.append(n)
                if len(path) == 2:
                    paths.append((path, '|'.join(l for l in label)))
        return paths

    x = mp.Pool(4)
    paths = x.map(find, itertools.combinations(include, 2))
    to_remove = set()
    for p in paths:
        if len(p) != 2:
            continue
        for path, label in p:
            print(path, label)
            for n in label.split('|'):
                to_remove.add(n)
            new_g.add_edge(path[0], path[1], label=label)
    for n in to_remove:
        new_g.remove_node(n)
    return new_g


def merge_nodes(graph):
    """ merges nodes into single node if same neighbors

    Parameters
    ----------
    graph: nx.DiGraph

    Returns
    -------
    nx.DiGraph
    """
    neighbors2node = dict()

    nodes = set(graph.nodes())
    edges = set(graph.edges())

    new_g = graph.copy()
    for i in nodes:
        down = '|'.join(sorted(graph.successors(i)))
        up = '|'.join(sorted(graph.predecessors(i)))

        if (down, up) in neighbors2node:
            neighbors2node[(down, up)].add(i)
        else:
            neighbors2node[(down, up)] = {i}

    for node, neigh in neighbors2node.items():
        # print(node, neigh)
        if len(neigh) > 1:
            print(node, neigh)
            new_name = '|'.join(sorted(neigh))
            for x in node[0].split('|'):
                if x != '':
                    for n in neigh:
                        if (x, n) in edges:
                            new_g.remove_edge(x, n)
                            new_g.add_edge(x, new_name, **graph.edge[x][n])
                            edges.remove((x, n))

            for x in node[1].split('|'):
                if x != '':
                    for n in neigh:
                        if (n, x) in edges:
                            new_g.remove_edge(n, x)
                            new_g.add_edge(new_name, x, **graph.edge[n][x])
                            edges.remove((n, x))
    print("{} nodes and {} edges".format(len(graph.nodes), len(graph.edges)))
    print("{} nodes and {} edges".format(len(new_g.nodes), len(new_g.edges)))
    for n in new_g.nodes():
        if len(new_g.successors(n)) == 0 and len(new_g.predecessors(n)) == 0:
            new_g.remove_node(n)

    print("{} nodes and {} edges".format(len(new_g.nodes), len(new_g.edges)))
    return new_g


def compress_edges(graph):
    """
    compress edges of networks by finding common paths

    Parameters
    ----------
    graph: networkx.DiGraph

    Returns
    -------

    """

    g = graph.copy()
    nodes = set(g.nodes())
    neighbor_dict = {}
    for n in nodes:
        neigh = tuple(g.neighbors(n))
        if tuple(neigh) in neighbor_dict:
            neighbor_dict[tuple(neigh)].append(n)
        else:
            neighbor_dict[tuple(neigh)] = []
            neighbor_dict[tuple(neigh)].append(n)

    for i in neighbor_dict:
        neigh = neighbor_dict[i]
        # print(i,neigh)
        if len(i) != 1 or len(neigh) == 1:
            continue
        interaction_types = {}
        for j in neighbor_dict[i]:
            if g.has_edge(i[0], j):
                direction = 'type1'
                edge = g.get_edge(i[0], j)
            elif g.has_edge(j, i[0]):
                direction = 'type2'
                edge = g.get_edge(j, i[0])
            edge_type = edge.attr['arrowhead']
            if edge_type + direction in interaction_types:
                interaction_types[edge_type + direction].append(j)
            else:
                interaction_types[edge_type + direction] = []
                interaction_types[edge_type + direction].append(j)
        for k in interaction_types:
            if len(interaction_types[k]) > 1:
                # print(i[0],'->',)
                to_join = []
                for each in interaction_types[k]:
                    if len(g.neighbors(each)) == 1:
                        to_join.append(each)
                # print(to_join,k)
                if len(to_join) > 1:
                    label = "{"
                    for node in to_join:
                        g.remove_node(node)
                        label += ' %s |' % str(node)
                    label = label[:-1]
                    label += "}"
                    # print(label)
                    g.add_node(label, shape='record', label=label)
                    if k.endswith('type2'):
                        g.add_edge(label, i[0], dir='both', arrowhead=k[:-5],
                                   arrowtail="none")
                    elif k.endswith('type1'):
                        g.add_edge(i[0], label, dir='both', arrowhead=k[:-5],
                                   arrowtail="none")
    return g
