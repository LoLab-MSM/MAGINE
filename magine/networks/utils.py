import os

import networkx as nx

import magine.networks.exporters as exporters

try:
    from IPython.display import Image, display
    IPYTHON = True
except RuntimeError:
    IPYTHON = False


def delete_disconnected_network(full_graph):
    """
    Delete disconnected parts of a provided network

    Parameters
    ----------
    full_graph : networkx.DiGraph

    Examples
    --------
    >>> from magine.networks.utils import delete_disconnected_network
    >>> from networkx import DiGraph
    >>> g = DiGraph()
    >>> g.add_nodes_from(['a', 'b', 'c'])
    >>> g.add_edge('a', 'b')
    >>> large_g_connected = delete_disconnected_network(g)
    >>> large_g_connected.nodes()
    ['b', 'a']
    >>> large_g_connected.edges()
    [('a', 'b')]
    """
    copy_network = full_graph.copy()
    tmp_g = copy_network.to_undirected()
    sorted_graphs = sorted(nx.connected_component_subgraphs(tmp_g), key=len,
                           reverse=True)
    for i in range(1, len(sorted_graphs)):
        copy_network.remove_nodes_from(sorted_graphs[i].nodes())
    return copy_network


def add_attribute_to_network(graph, list_to_add_attribute, attribute,
                             true_term, false_term='false'):
    """

    Parameters
    ----------
    graph : networkx graph
    list_to_add_attribute : list
        list of nodes in graph to add attribute to
    attribute : str
        attribute to add to graph
    true_term : str
        value to add for the attribute provided True
    false_term : str
        value to add if attribute is false

    Returns
    -------
    nx.DiGraph

    Examples
    --------
    >>> from magine.networks.utils import add_attribute_to_network
    >>> from networkx import DiGraph
    >>> g = DiGraph()
    >>> g.add_nodes_from(['a', 'b', 'c'])
    >>> new_g = add_attribute_to_network(g, ['a','b'], 'inTest', 'true')
    >>> new_g.node['a']
    {'inTest': 'true'}
    >>> new_g.node['c']
    {'inTest': 'false'}

    """
    tmp_g = graph.copy()
    nodes = set(tmp_g.nodes())
    set_of_positive = set(list_to_add_attribute)
    for i in nodes:
        if i in set_of_positive:
            tmp_g.node[i][attribute] = true_term
        else:
            tmp_g.node[i][attribute] = false_term

    return tmp_g


def add_color_graphviz_fmt(graph, list_to_paint, color):
    """

    Parameters
    ----------
    graph: networkx.DiGraph
    list_to_paint : list_like

    color : str


    Examples
    --------
    >>> from magine.networks.utils import add_color_graphviz_fmt
    >>> from networkx import DiGraph
    >>> g = DiGraph()
    >>> g.add_nodes_from(['a', 'b', 'c'])
    >>> new_g = add_color_graphviz_fmt(g, ['a','b'], 'red')
    >>> new_g = add_color_graphviz_fmt(new_g, ['c'], 'blue')
    >>> new_g.node['a']
    {'style': 'filled', 'color': 'black', 'measured': True, 'fillcolor': 'red'}
    >>> new_g.node['c']
    {'style': 'filled', 'color': 'black', 'measured': True, 'fillcolor': 'blue'}


    Returns
    -------

    """
    tmp_g = graph.copy()
    nodes1 = set(tmp_g.nodes())
    for i in list_to_paint:
        if i in nodes1:
            tmp_g.node[i]['measured'] = True
            tmp_g.node[i]['color'] = 'black'
            tmp_g.node[i]['fillcolor'] = color
            tmp_g.node[i]['style'] = 'filled'

    return tmp_g


def paint_network_overtime_up_down(graph, list_up, list_down, save_name,
                                   color_up='red', color_down='blue',
                                   labels=None, create_gif=False):
    """
    Adds color attribute to network over time.

    Parameters
    ----------
    graph : nx.DiGraph
        Network
    list_up : list_like
        List of lists, where the inner list contains the node to add the
        color
    list_down : list_like
        list of colors for each time point
    color_up : str
        color for first list of species
    color_down : str
        color of second list of species
    save_name : str
        prefix for images to be saved
    labels: list_like
        list of labels to add to graph per sample
    create_gif : bool
        Create a gif from series of images
    Returns
    -------

    """

    if len(list_up) != len(list_down):
        print('Length of list of data must equal len of color list')
        return
    if labels is not None:
        if len(labels) != len(list_down):
            print('Length of labels must be equal to len of data')
            return
    string = 'convert -delay 100 '
    tmp_graph = graph.copy()

    for n, (up, down) in enumerate(zip(list_up, list_down)):
        graph2 = add_color_graphviz_fmt(tmp_graph, up, color_up)
        graph2 = add_color_graphviz_fmt(graph2, down, color_down)
        both = set(up).intersection(set(down))
        graph2 = add_color_graphviz_fmt(graph2, both, 'yellow')

        if labels is not None:
            graph2.graph['label'] = labels[n]
            graph2.graph['fontsize'] = 13

        s_name = '%s_%04i.png' % (save_name, n)

        exporters.export_to_dot(graph2, s_name, 'png', 'dot')
        string += s_name + ' '
        if IPYTHON and run_from_ipython():
            display(Image(s_name))

    string1 = string + '  %s.gif' % save_name
    string2 = string + '  %s.pdf' % save_name
    if create_gif:
        os.system(string1)
        os.system(string2)


def paint_network_overtime(graph, list_of_lists, color_list, save_name,
                           labels=None, create_gif=False):
    """
    Adds color attribute to network over time.

    Parameters
    ----------
    graph : nx.DiGraph
        Network
    list_of_lists : list_like
        List of lists, where the inner list contains the node to add the
        color
    color_list : list_like
        list of colors for each time point
    save_name : str
        prefix for images to be saved
    labels: list_like
        list of labels to add to graph per sample

    """

    if len(list_of_lists) != len(color_list):
        print('Length of list of data must equal len of color list')
        return
    if labels is not None:
        if len(labels) != len(list_of_lists):
            print('Length of labels must be equal to len of data')
            return

    string = 'convert -delay 100 '
    tmp_graph = graph.copy()

    for n, i in enumerate(list_of_lists):
        graph2 = add_color_graphviz_fmt(tmp_graph, i, color_list[n])

        if labels is not None:
            graph2.graph['label'] = labels[n]
            graph2.graph['fontsize'] = 13

        s_name = '%s_%04i.png' % (save_name, n)
        exporters.export_to_dot(graph2, s_name, 'png', 'dot')
        if IPYTHON and run_from_ipython():
            display(Image(s_name))

        string += s_name + ' '
    string1 = string + '  %s.gif' % save_name
    string2 = string + '  %s.pdf' % save_name
    if create_gif:
        os.system(string1)
        os.system(string2)


def compose(g, g_1):
    """Return a new graph of G composed with H.

    Composition is the simple union of the node sets and edge sets.
    The node sets of G and H do not need to be disjoint.

    Parameters
    ----------
    g, : nx.DiGraph
    g_1 : nx.DiGraph
       A NetworkX graph

    Returns
    -------
    C: A new graph  with the same type as G

    """

    new_g = nx.DiGraph()

    _add_nodes(g, new_g)
    _add_nodes(g_1, new_g)

    _add_edges(g, new_g)
    _add_edges(g_1, new_g)

    return new_g


def compose_all(graphs):
    """Return the composition of all graphs.

    Composition is the simple union of the node sets and edge sets.
    The node sets of the supplied graphs need not be disjoint.

    Parameters
    ----------
    graphs : list
       List of NetworkX graphs

    Returns
    -------
    C : A graph with the same type as the first graph in list

    """
    graphs = iter(graphs)
    g = next(graphs)
    for h in graphs:
        g = compose(g, h)
    return g


def _add_nodes(old_network, new_network):
    new_nodes = set(new_network.nodes())
    for i, data in old_network.nodes_iter(data=True):
        if i not in new_nodes:
            new_network.add_node(i, **data)
        else:
            existing_info = new_network.node[i]
            for n, d in data.items():
                if n not in existing_info:
                    new_network.node[i][n] = d
                else:
                    additions = set(d.split('|'))
                    if isinstance(existing_info[n], list):
                        old = set(existing_info[n][0].split('|'))
                    else:
                        old = set(existing_info[n].split('|'))
                    additions.update(old)
                    new_network.node[i][n] = '|'.join(sorted(additions))


def _add_edges(current_network, new_network):
    edges = set(new_network.edges())
    for i, j, data in current_network.edges_iter(data=True):
        if (i, j) not in edges:
            new_network.add_edge(i, j, **data)
        else:
            existing_info = new_network.edge[i][j]
            for n, d in data.items():
                if n not in existing_info:
                    new_network[i][j][n] = d
                else:
                    additions = set(d.split('|'))
                    additions.update(set(existing_info[n].split('|')))
                    new_network[i][j][n] = '|'.join(sorted(additions))


_maps = {
    'activation': 'activate',
    'activator': 'activate',
    'potentiator': 'activate',

    'inducer': 'expression',
    'stimulator': 'expression',
    'suppressor': 'repression',

    'blocker': 'inhibit',
    'inhibitor': 'inhibit',
    'inhibition': 'inhibit',
    'inhibitor, competitive': 'inhibit',

    'proteolytic processing': 'cleavage',

    # binding
    'binding/association': 'binding',
    'binder': 'binding',
    'complex': 'binding',
    'dissociation': 'binding',

    # indirect/missing
    'indirect effect': 'indirect',
    'missing interaction': 'indirect',

    'state change': 'stateChange',

    'ubiquitination': 'ubiquitinate',
    'methylation': 'methylate',
    'glycosylation': 'glycosylate',
    'sumoylation': 'sumoylate',
    'ribosylation': 'ribosylate',
    'neddylation': 'neddylate',
    'desumoylation': 'desumoylate',
    'deneddylation': 'deneddylate',
    'demethylation': 'demethylate',
    'deacetylation': 'deacetylate',
    'desensitize the target': 'inhibit',
    'deubiquitination': 'deubiquitinate',
    'nedd(rub1)ylation': 'neddy(rub1)late',

    'dephosphorylation': 'dephosphorylate',
    'phosphorylation': 'phosphorylate',

    'negative modulator': 'inhibit',
    'inhibitory allosteric modulator': 'allosteric|inhibit',
    'allosteric modulator': 'allosteric|modulate',
    'positive allosteric modulator': 'activate|allosteric',
    'positive modulator': 'activate',
    'partial agonist': 'activate|chemical',
    'inverse agonist': 'activate|chemical',
    'agonist': 'activate|chemical',

    'antagonist': 'inhibit|chemical',
    'partial antagonist': 'inhibit|chemical',

    # chemical related
    'compound': 'chemical',
    'product of': 'chemical',
    'ligand': 'chemical',
    'cofactor': 'chemical',
    'multitarget': 'chemical',
}


def standardize_edge_types(network):
    to_remove = set()
    for source, target, data in network.edges_iter(data=True):
        if 'interactionType' in data:
            edge_type = data['interactionType']
            edge_type = set(i for i in edge_type.split('|'))

            for k, v in _maps.items():
                if k in edge_type:
                    edge_type.remove(k)
                    edge_type.add(v)
            for i in ['', ' ']:
                if i in edge_type:
                    edge_type.remove(i)

            if 'reaction' in edge_type:
                if len(edge_type) != 1:
                    edge_type.remove('reaction')
            if 'catalyze' in edge_type:
                if len(edge_type) != 1:
                    edge_type.remove('catalyze')
            edge_type = '|'.join(sorted(edge_type))
            if edge_type == '':
                # network.remove_edge(source, target)
                to_remove.add((source, target))
            else:
                network[source][target]['interactionType'] = edge_type
    for source, target in to_remove:
        network.remove_edge(source, target)


def append_attribute_to_network(graph, list_to_add_attribute, attribute,
                                true_term, delimiter=''):
    """

    Parameters
    ----------
    graph : networkx graph
    list_to_add_attribute : list
        list of nodes in graph to add attribute to
    attribute : str
        attribute to add to graph
    true_term : str
        value to add for the attribute provided True
    delimiter : str


    Returns
    -------
    out : networkx graph


    >>> from networkx import DiGraph
    >>> from magine.networks.utils import append_attribute_to_network
    >>> g = DiGraph()
    >>> g.add_nodes_from(['a'])
    >>> new_g = append_attribute_to_network(g, ['a'], 'attribute', 'one')
    >>> new_g.node['a']
    {'attribute': 'one'}
    >>> new_g = append_attribute_to_network(new_g, ['a'], 'attribute', 'two', delimiter=',')
    >>> new_g.node['a']
    {'attribute': 'one,two'}

    """
    tmp_g = graph.copy()
    nodes = set(tmp_g.nodes())
    set_of_positive = set(list_to_add_attribute)
    for i in nodes:
        if i in set_of_positive:
            if attribute in tmp_g.node[i].keys():
                new_attr = tmp_g.node[i][attribute] + delimiter + true_term
                tmp_g.node[i][attribute] = new_attr
            else:
                tmp_g.node[i][attribute] = true_term

    return tmp_g


def trim_sink_source_nodes(network, list_of_nodes, remove_self_edge=False):
    """
    Trim graph by removing nodes that are not in provided list if source/sink


    Parameters
    ----------
    network : networkx.DiGraph
    list_of_nodes : list_like
        list of species that are important if sink/source
    remove_self_edge : bool
        Remove self edges
    Returns
    -------

    """
    tmp_network = network.copy()
    if remove_self_edge:
        tmp_network = remove_self_edges(tmp_network)
    tmp1 = _trim(tmp_network, list_of_nodes)
    tmp2 = _trim(tmp_network, list_of_nodes)
    while tmp1 != tmp2:
        tmp2 = tmp1
        tmp1 = _trim(tmp_network, list_of_nodes)
    return tmp_network


def remove_self_edges(network):
    tmp_network = network.copy()
    tmp_network.remove_edges_from(tmp_network.selfloop_edges())
    return tmp_network


def _trim(network, list_of_nodes):
    list_of_nodes = set(list_of_nodes)
    found, not_found = 0., 0.
    nodes = set(network.nodes())
    in_dict = network.in_degree(nbunch=nodes)
    out_dict = network.out_degree(nbunch=nodes)
    for i in nodes:
        if i in list_of_nodes:
            found += 1
        else:
            if in_dict[i] == 1 and out_dict[i] == 0:
                network.remove_node(i)
            elif in_dict[i] == 0 and out_dict[i] == 1:
                network.remove_node(i)
            else:
                not_found += 1

    len_nodes = len(network.nodes())
    print("{} found, {} not found".format(found, not_found))
    print("{}% of {} nodes".format(found / len_nodes * 100, len_nodes))
    return len_nodes


def subtract_network_from_network(net1, net2):
    """
    subtract one network from another

    Parameters
    ----------
    net1 : networkx.DiGraph
    net2 : networkx.DiGraph

    Returns
    -------
    networkx.DiGraph

    """
    copy_graph1 = net1.copy()
    nodes1 = set(net1.nodes())
    nodes2 = set(net2.nodes())
    overlap = nodes2.intersection(nodes1)
    copy_graph1.remove_edges_from(overlap)

    return copy_graph1


def run_from_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False


if __name__ == '__main__':
    g = nx.DiGraph()
    g.add_edge('A', 'B')
    g.add_edge('B', 'C')
    g.add_edge('B', 'E')
    g.add_edge('C', 'D')
    # test_g = remove_unmeasured_nodes(g, ['A', 'D', 'E'])
    # export_to_dot(test_g, 'merged_node')
    # new_g = remove_unmeasured_nodes(g, ['A', 'C', 'D'])
    # export_to_dot(new_g, 'merged_node2')
    g = nx.DiGraph()
    g.add_node('A', color='red', intType='ugly')
    g.add_edge('A', 'B', iType='no')
    gg = nx.DiGraph()
    gg.add_node('A', color='green', intType='ugly')
    gg.add_edge('A', 'B', iType='yes')
    fg = compose(g, gg)
    for i in fg.nodes(data=True):
        print(i)
    for i in fg.edges(data=True):
        print(i)
