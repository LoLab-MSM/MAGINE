import matplotlib.pyplot as plt
import networkx as nx

from magine.networks.visualization import render_igraph, render_mpl


def te_igraph():
    g = nx.DiGraph()
    g.add_node('A', termName='one')
    g.add_node('B', termName='one')
    g.add_node('C', termName='two')
    g.add_node('D', termName='two')
    g.add_edge('A', 'B')
    g.add_edge('C', 'B')
    g.add_edge('C', 'D')

    plot, pos = render_igraph(g, )
    plot.show()
    plot, pos = render_igraph(g, cluster=True)
    plot.show()


def test_mpl():
    g = nx.DiGraph()
    g.add_node('A', termName='one')
    g.add_node('B', termName='one')
    g.add_node('C', termName='two')
    g.add_node('D', termName='two')
    g.add_edge('A', 'B')
    g.add_edge('C', 'B')
    g.add_edge('C', 'D')

    render_mpl(g, )
    plt.show()


if __name__ == '__main__':
    test_mpl()
