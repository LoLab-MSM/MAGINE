import matplotlib.pyplot as plt
import networkx as nx

from magine.networks.visualization import render_igraph, render_mpl


def test_igraph():
    g = nx.DiGraph()
    g.add_node('A', termName='one')
    g.add_node('B', termName='one')
    g.add_node('C', termName='two')
    g.add_node('D', termName='two')
    g.add_edge('A', 'B')
    g.add_edge('C', 'B')
    g.add_edge('C', 'D')

    plot, pos = render_igraph(g, )
    plot, pos = render_igraph(g, cluster=True)
    plot, pos = render_igraph(g, cluster=True, positions=pos)


def test_mpl():
    g = nx.DiGraph()
    g.add_node('A', termName='one')
    g.add_node('B', termName='one')
    g.add_node('C', termName='two')
    g.add_node('D', termName='two')
    g.add_edge('A', 'B')
    g.add_edge('C', 'B')
    g.add_edge('C', 'D')
    for i in ['circular_layout', 'random_layout', 'shell_layout',
              'spring_layout', 'spectral_layout',  # 'dot', 'neato', 'fdp',
              'fruchterman_reingold_layout']:
        fig = render_mpl(g, i)
        plt.close()


if __name__ == '__main__':
    test_igraph()
