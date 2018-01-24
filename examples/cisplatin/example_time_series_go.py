import pandas as pd
import networkx as nx
from magine.data.datatypes import ExperimentalData
from magine.networks.visualization.cytoscape_js_view import create_subnetwork, display_graph


if __name__ == '__main__':
    data = pd.read_csv('Data/norris_et_al_2017_cisplatin_data.csv.gz',
                       low_memory=False)
    exp_data = ExperimentalData(data)

    df = pd.read_csv('Data/all_cisplatin_out.csv.gz')
    df.sort_values(by='pvalue', ascending=True, inplace=True)

    # TODO update with network painted by timepoint, modality
    network = nx.read_gpickle('Data/cisplatin_network.p')
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
                                          network, '1hr', cytoscape_js=False)




