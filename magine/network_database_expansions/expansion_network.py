import os

import networkx as nx

from magine.network_database_expansions.reactome.reactome_functional_interaction \
    import ReactomeFunctionalInteraction


def load_reactome_fi():
    """
    Load reactome functional interaction network
    Returns
    -------

    """
    path = os.path.join(os.path.dirname(__file__),
                        'reactome',
                        'reactome_fi.gml')
    if not os.path.exists(path):
        print("Downloading Reactome Functional interaction network!")
        rfi = ReactomeFunctionalInteraction()
        rfi.parse_network()
    if not os.path.exists(path):
        print('Failed!')
        quit()
    return nx.read_gml(path)
