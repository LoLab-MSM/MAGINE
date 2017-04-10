import networkx as nx
import os
from reactome_expansion.download_reactome_functional_interaction import ReactomeFunctionalInteraction


def load_reactome_fi():
    path = os.path.join(os.path.dirname(__file__),
                        'reactome_expansion',
                        'reactome_fi.gml')
    if not os.path.exists(path):
        print("Downloading Reactome Functional interaction network!")
        rfi = ReactomeFunctionalInteraction()
        rfi.setup()
        rfi.parse_network()
    if not os.path.exists(path):
        print('Failed!')
        quit()
    return nx.read_gml(path)
