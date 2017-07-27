from magine.data.datatypes import ExperimentalData
from magine.plotting.species_plotting import plot_dataframe
from magine.networks.network_generator import build_network
from magine.ontology.enrichr import Enrichr


data = ExperimentalData('Data/norris_et_al_2017_cisplatin_data.csv.gz')

# want data organized
# data
# proteomics
# metabolomics
plot_dataframe(data.metabolites, 'metabolites', 'metabolites',
               type_of_species='metabolites')

plot_dataframe(data.proteomics, 'genes', 'protein',
               type_of_species='protein')


e = Enrichr()
e.run_key_dbs(data.proteomics_over_time, data.proteomics_time_points,
              save_name='proteomics_changed')

e.run_key_dbs(data.rna_over_time, data.rna_time_points,
              save_name='rnaseq_changed')

build_network(data.list_sig_proteins,
              save_name='example_network',
              all_measured_list=data.list_species,
              metabolite_list=data.list_metabolites,
              use_reactome=True, use_hmdb=True
              )

