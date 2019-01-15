import os

import networkx as nx

import magine.networks.utils as utils
from magine.examples.cisplatin.exp_data import exp_data
from magine.networks.network_generator import build_network

if __name__ == '__main__':

    network = build_network(
        seed_species=exp_data.species.sig.id_list,  # genes seed species
        all_measured_list=exp_data.species.id_list,  # all data measured
        use_biogrid=True,  # expand with biogrid
        use_hmdb=True,  # expand with hmdb
        use_reactome=True,  # expand with reactome
        use_signor=True,  # expand with signor
        trim_source_sink=True,  # remove all source and sink nodes not measured
        save_name='Data/cisplatin_based_network_new'
    )

    # Load the network, note that it is returned above but for future use
    # we will use load in
    network = nx.read_gpickle('Data/cisplatin_based_network.p')

    utils.add_data_to_graph(network, exp_data)
    print("Saving network")
    # write to GML for cytoscape or other program
    nx.write_gml(

        network,
        os.path.join(os.path.dirname(__file__), 'Data',
                     'cisplatin_network_w_attributes.gml')
    )

    # write to gpickle for fast loading in python
    nx.write_gpickle(
        network,
        os.path.join(os.path.dirname(__file__), 'Data',
                     'cisplatin_based_network.p'),
    )
