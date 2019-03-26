import matplotlib.pyplot as plt
import networkx as nx

from magine.networks.visualization import draw_igraph, draw_mpl


def test_igraph():
    g = nx.DiGraph()
    g.add_node('A', termName='one')
    g.add_node('B', termName='one')
    g.add_node('C', termName='two')
    g.add_node('D', termName='two')
    g.add_edge('A', 'B')
    g.add_edge('C', 'B')
    g.add_edge('C', 'D')

    plot, pos = draw_igraph(g, )
    plot, pos = draw_igraph(g, cluster=True)
    plot, pos = draw_igraph(g, cluster=True, positions=pos)


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
        fig = draw_mpl(g, i)
        plt.close()
