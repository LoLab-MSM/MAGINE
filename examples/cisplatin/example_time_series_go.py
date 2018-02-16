import networkx as nx
import pandas as pd

from magine.networks.visualization.notebook_tools import create_subnetwork

if __name__ == '__main__':

    df = pd.read_csv('Data/cisplatin_enrichment.csv.gz')
    network = nx.read_gpickle('Data/cisplatin_based_network.p')

    shorten_names = {
        'Cell Cycle_Homo sapiens_R-HSA-1640170': 'Cell Cycle',
        'DNA Repair_Homo sapiens_R-HSA-73894': 'DNA Repair',
        'MAP2K and MAPK activation_Homo sapiens_R-HSA-5674135': 'MAPK signals',
        'Interleukin-2 signaling_Homo sapiens_R-HSA-451927': 'IL2',
        'positive regulation of apoptotic process ': 'apoptosis'
    }
    renamed = df.copy()
    renamed['term_name'] = renamed['term_name'].replace(shorten_names)
    term_net, mol_net = create_subnetwork(shorten_names.values(), renamed,
                                          network, 'subgraphs_from_one_hour')




