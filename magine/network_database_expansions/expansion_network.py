import networkx as nx
import os


def load_reactome_fi():
    path = os.path.join(os.path.dirname(__file__),
                        'reactome_expansion',
                        'reactome_fi.gml')
    return nx.read_gml(path)
