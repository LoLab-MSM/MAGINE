import matplotlib.pyplot as plt

from magine.ontology.enrichr import Enrichr

pro_death = ['BID', 'BAD', 'PMAIP1', 'BAX', 'BAK1']
anti_death = ['MCL1', 'BCL2', 'BCL2L1']
caspases = ['CASP3', 'CASP6', 'CASP8']
downstream_mito = ['DIABLO', 'CYCS', 'PARP1', 'APAF1', 'XIAP']
upstream_mito = ['FAS', 'FADD', 'CFLAR', 'BFAR']

list_of_proteins = [pro_death, anti_death, caspases, downstream_mito,
                    upstream_mito]
sample_labels = ['pro_death', 'anti_death', 'caspases', 'down_mito', 'up_mito']

enrichr = Enrichr()

results = enrichr.run_samples(sample_lists=list_of_proteins,
                              sample_ids=sample_labels,
                              gene_set_lib='GO_Biological_Process_2017')

enrichment = results['combined_score'].copy()[sample_labels]
enrichment.sort_values(by=sample_labels, ascending=False, inplace=True)
print(enrichment.head(10))

plt.imshow(enrichment.as_matrix(), aspect='auto')
plt.xticks(range(len(sample_labels)), sample_labels, rotation=45)
plt.colorbar()
plt.show()
