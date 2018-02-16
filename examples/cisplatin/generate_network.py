import os

import networkx as nx

import magine.networks.network_tools as nt
from exp_data import exp_data
from magine.networks.network_generator import build_network

if __name__ == '__main__':


    network = build_network(
        gene_list=exp_data.list_sig_proteins,  # genes seed species
        metabolite_list=exp_data.list_metabolites,  # metabolites seed species
        all_measured_list=exp_data.list_species,  # all data measured
        use_biogrid=True,  # expand with biogrid
        use_hmdb=True,  # expand with hmdb
        use_reactome=True,  # expand with reactome
        trim_source_sink=True,  # remove all source and sink nodes not measured
        save_name='Data/cisplatin_based_network'
    )

    # Load the network, note that it is returned above but for future use
    # we will use load in
    network = nx.read_gpickle('Data/cisplatin_based_network.p')

    # label all nodes that are measure and significantly changed in the data
    measured = set(exp_data.list_species)
    sig_measured = set(exp_data.list_sig_species)

    network = nt.add_attribute_to_network(network, sig_measured, 'sigMeasured',
                                          'red', 'blue')
    network = nt.add_attribute_to_network(network, measured, 'measured', 'red',
                                          'blue')
    # add labels for if node is measured in any of our data_types
    m, sig_m = exp_data.get_measured_by_datatype()

    for exp_type, spec in m.items():
        attr_name = exp_type.replace('_', '')
        attr_name = attr_name.replace('-', '')
        network = nt.add_attribute_to_network(network, spec, attr_name, 'red',
                                              'blue')
    # add labels for if node is measured in any of our time points
    for time, spec in exp_data.sig_species_over_time.items():
        time = 'timepoint{}'.format(time)
        network = nt.add_attribute_to_network(network, spec, time, 'red', 'blue')

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
