from magine.ontology.ontology_analysis import GoAnalysis
from magine.ontology.ontology_tools import slim_ontology

pro_death = ['BID', 'BAD', 'PMAIP1',  # NOXA
             'BAX', 'BAK1',  # effector proteins
             ]
anti_death = ['MCL1', 'BCL2', 'BCL2L1', ]  # anti-apoptotic
caspases = ['CASP3', 'CASP6', 'CASP8', ]
downstream_mito = ['DIABLO', 'CYCS',
                   'PARP1', 'APAF1', 'XIAP', ]
upstream_mito = ['FAS', 'FADD',  # DISC
                    'CFLAR',  # FLIP
                    'BFAR',  # BAR
                 ]
list_of_proteins = [pro_death, anti_death, caspases, downstream_mito,
                    upstream_mito]

go = GoAnalysis(species='hsa', output_directory='output',
                save_name='multi_sample_earm_proteins',
                )

results = go.calculate_enrichment(list_of_proteins)

sig_results = slim_ontology(results,
                            pvalue=0.05,
                            go_aspects='biological_process',
                            trim_nodes=True,
                            n_top_hits=10)
sig_results.sort_values(by='enrichment_score', ascending=False, inplace=True)
print(sig_results.head(10))

