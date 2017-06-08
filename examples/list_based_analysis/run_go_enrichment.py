import pandas as pd

from magine.ontology.ontology_analysis import GoAnalysis
from magine.ontology.ontology_tools import slim_ontology

pd.set_option('display.width', 10000)
pd.set_option('max_colwidth', 100)

list_of_proteins = ['CASP3', 'CASP6', 'FAS', 'FADD', 'CASP8', 'CFLAR',
                    'BFAR', 'BAD', 'BID', 'PMAIP1', 'MCL1', 'BCL2', 'BCL2L1',
                    'BAX', 'BAK1', 'DIABLO', 'CYCS', 'PARP1', 'APAF1', 'XIAP']

# go analyzer class
go = GoAnalysis(species='hsa', output_directory='output',
                save_name='example_earm_proteins',
                )

results = go.calculate_enrichment(list_of_proteins)

# just want to look at biological processes
results = results[results['aspect'] == 'biological_process']
# want to look at these columns
cols = ['GO_name', 'GO_id', 'pvalue', 'enrichment_score']

results.sort_values(by='enrichment_score', ascending=False, inplace=True)
results.to_csv('enrichment_results.csv', index=False)
print(results[cols].head(10))

results.sort_values(by='pvalue', ascending=True, inplace=True)
print(results[cols].head(10))

sig_results = slim_ontology(results,
                            pvalue=0.05,
                            trim_nodes=True,
                            n_top_hits=10)
sig_results.sort_values(by='enrichment_score', ascending=False, inplace=True)
print(sig_results[cols].head(10))

sig_results.sort_values(by='pvalue', ascending=True, inplace=True)
print(sig_results[cols].head(10))
sig_results.to_csv('enrichment_results_slimmed.csv', index=False)

"""

CASP3, CASP6, FAS, FADD, CASP8, CFLAR, BFAR, BAD, BID, PMAIP1, MCL1, BCL2, 
BCL2L1, BAX, BAK1, DIABLO, CYCS, PARP1, APAF1, XIAP

CASP3
CASP6
FAS
FADD
CASP8
CFLAR
BFAR
BAD
BID
PMAIP1
MCL1
BCL2
BCL2L1
BAX
BAK1
DIABLO
CYCS
PARP1
APAF1
XIAP

"""
