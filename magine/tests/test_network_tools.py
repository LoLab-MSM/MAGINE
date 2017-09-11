"""
import magine.networks.network_tools as nt
import networkx as nx



g = nx.read_gpickle('../networks/background_network.p')

g = nt.merge_nodes(g)
nx.write_gml(g, 'merged_nodes_test.gml')




def test_compress():
    g = nx.DiGraph()
    g.add_node('A')
    g.add_node('B')
    g.add_node('C')
    g.add_node('D')
    g.add_node('E')
    g.add_edge('A', 'D', dir='both', arrowhead='dot', arrowtail="none", )
    g.add_edge('B', 'D', dir='both', arrowhead='dot', arrowtail="none", )
    g.add_edge('C', 'D', dir='both', arrowhead='dot', arrowtail="none", )
    g.add_edge('E', 'A', dir='both', arrowhead='dot', arrowtail="none", )
    g.add_edge('E', 'B', dir='both', arrowhead='dot', arrowtail="none", )
    g.add_edge('E', 'C', dir='both', arrowhead='odot', arrowtail="none", )
    nodelist = ['A', 'B', 'C', 'D', 'E']
    adj = nx.to_numpy_recarray(g, nodelist=['A', 'B', 'C', 'D', 'E'],
                               dtype=[('arrowhead', 'S50')])

    for n, i in enumerate(range((np.shape(adj)[1]))):
        to_concate = []
        for k, j in enumerate(range(n, (np.shape(adj)[1]))):
            if i == j:
                continue
            if (adj[:, i] == adj[:, j]).all():
                to_concate.append(j)
                print("here", nodelist[i], nodelist[j])
        print(to_concate)
    print(adj)
    # test = nx.from_numpy_matrix(adj, create_using=nx.DiGraph())
    # nx.draw_graphviz(test, prog='dot')
    # plt.savefig('test.png')


def test_trim_sink_source():
    g = nx.DiGraph()

    sink_nodes = ['sink_1', 'sink_2', 'sink_3', 'sink_4']

    g.add_path(['source_1', 'b', 'c', 'sink_1'])
    g.add_path(['source_2', 'e', 'f', 'sink_2'])
    g.add_path(['source_2', 'f', 'sink_3'])
    g.add_path(['source_3', 'f', 'l', 'sink_4'])
    # source_1 and source_3 now have two edges
    g.add_edge('source_1', 'b')
    g.add_edge('source_3', 'b')
    g.add_edge('b', 'c')
    g.add_edge('a', 'c')
    export_to_dot(g, 'test_before')
    test_set = ['source_1']
    # test should result in only
    # ['a', 'c', 'b', 'source_3', 'source_1']
    test_1 = trim_sink_source_nodes(g, test_set)
    test_1_nodes = test_1.nodes()
    print(test_1_nodes)

    # export_to_dot(test_1, 'test_after')
    for i in sink_nodes:
        assert i not in test_1_nodes
"""
