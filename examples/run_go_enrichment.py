from magine.ontology.ontology_analysis import GoAnalysis
from magine.ontology.ontology_tools import slim_ontology

list_of_proteins = ['CASP3', 'CASP6',

                    'FAS', 'FADD', 'CASP8',  # DISC
                    'CFLAR',  # FLIP
                    'BFAR',  # BAR
                    'BAD',
                    # pro-apoptotic
                    'BID',
                    'PMAIP1',  # NOXA
                    'MCL1', 'BCL2', 'BCL2L1',  # anti-apoptotic
                    'BAX', 'BAK1',  # effector proteins
                    'DIABLO', 'CYCS',
                    'PARP1', 'APAF1', 'XIAP',
                    ]

go = GoAnalysis(species='hsa', output_directory='output',
                save_name='example_earm_proteins',
                )

results = go.enrichment_analysis_of_single_sample(list_of_proteins)

sig_results = slim_ontology(results,
                            pvalue=0.05,
                            go_aspects='biological_process',
                            trim_nodes=True,
                            n_top_hits=10)
sig_results.sort_values(by='enrichment_score', ascending=False, inplace=True)
print(sig_results.head(10))
