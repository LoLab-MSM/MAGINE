from magine.networks.network_generator import  build_network
from magine.data.datatypes import ExperimentalData
import magine.networks.network_tools as nt
import networkx as nx
import os


data_file = os.path.join(
    os.path.dirname(__file__),
    'Data',
    'norris_et_al_2017_cisplatin_data.csv.gz'
)

exp_data = ExperimentalData(data_file)

network = build_network(gene_list = exp_data.list_sig_proteins,
                        metabolite_list=exp_data.list_metabolites,
                        all_measured_list=exp_data.list_species,
                        use_biogrid=True,
                        use_hmdb=True,
                        use_reactome=True
                        )

print("node {} has attributes {}".format('EPB41', network.node['EPB41']))

nx.write_gpickle(
    network,
    os.path.join(os.path.dirname(__file__), 'Data', 'cisplatin_based_network.p'),
)


network = nx.read_gpickle('Data/cisplatin_based_network.p')

measured = set(exp_data.list_species)
sig_measured = set(exp_data.list_sig_species)
network = nt.add_attribute_to_network(network, sig_measured, 'sigMeasured', 'red', 'blue')
network = nt.add_attribute_to_network(network, measured, 'measured', 'red', 'blue')


# nx.write_gml(network, os.path.join(os.path.dirname(__file__), 'Data',
#                                    'cisplatin_network_w_attributes.gml')
#              )


m, sig_m = exp_data.get_measured_by_datatype()
nodes = set(network.nodes())

network = nt.trim_sink_source_nodes(network, exp_data.list_species)


for exp_type, spec in m.items():
    attr_name = exp_type.replace('_', '')
    attr_name = attr_name.replace('-', '')
    network = nt.add_attribute_to_network(network, spec, attr_name, 'red', 'blue')
    # network = nt.add_attribute_to_network(network, spec, 'color', 'red', 'lightblue')


# for exp_type, spec in sig_m.items():
#     attr_name = exp_type.replace('_', '')
#     attr_name = attr_name.replace('-', '')
#     true = 'sig'
#     network = nt.add_attribute_to_network(network, spec, attr_name, true)
#     network = nt.add_attribute_to_network(network, spec, 'color', 'red',
#                                           'lightblue')


for time, spec in exp_data.sig_species_over_time.items():
    time = 'timepoint{}'.format(time)
    network = nt.add_attribute_to_network(network, spec, time, 'red', 'blue')



print("Saving network")
nx.write_gml(network, os.path.join(os.path.dirname(__file__), 'Data',
                                   'cisplatin_network_w_attributes.gml')
             )
