from magine.enrichment.enrichr import run_enrichment_for_project
from magine.examples.cisplatin.exp_data import exp_data

run_enrichment_for_project(exp_data, 'cisplatin_enrichment')
