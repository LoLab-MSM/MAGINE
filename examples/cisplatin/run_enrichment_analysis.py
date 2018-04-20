from exp_data import exp_data
from magine.enrichment.enrichr import run_enrichment_for_project

# exp_data.create_table_of_data(sig=True, save_name='test', plot=True)

run_enrichment_for_project(exp_data, 'cisplatin_enrichment')
