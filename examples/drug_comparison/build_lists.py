import pandas as pd

from magine.ontology.enrichr import get_background_list, Enrichr

info = get_background_list('Drug_Perturbations_from_GEO_up')

e = Enrichr()

terms = [
    'chlorambucil db00291 human gse8832 sample 2691',
    'cisplatin db00515 human gse47856 sample 3146',
    'carboplatin db00958 human gse49577 sample 3177',
    'methotrexate db00563 human gse11440 sample 3072',
    'doxorubicin db00997 human gse12972 sample 2720',
]

df = pd.DataFrame(info, columns=['term', 'gene_list', 'n_genes'])

df.to_csv('drug_pertubations_from_geo_up.csv')

df = df[df['term'].isin(terms)]

terms, gene_sets = [], []

for t in df.to_dict('records'):
    terms.append(t['term'].split(' ')[0])
    gene_sets.append(t['gene_list'])

enrichment = e.run_sample_set_of_dbs(gene_sets, terms, pivot=False)

enrichment.to_csv('compare_drug_dbs.csv', encoding='utf8')
