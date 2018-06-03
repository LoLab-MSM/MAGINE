import ast

import pandas as pd

from magine.enrichment.enrichr import get_background_list, Enrichr

e = Enrichr()


def download_database():
    info = get_background_list('Drug_Perturbations_from_GEO_up')

    df = pd.DataFrame(info, columns=['term', 'gene_list', 'n_genes'])

    def split(row):
        names = row['term'].split(' ')
        for i in names:
            if i.startswith('gse'):
                return i

    def get_sample(row):
        return row['term'].split('sample')[1]

    def only_drug(row):
        return row['term'].split(' ')[0]

    df['gse'] = df.apply(split, axis=1)
    df['sample'] = df.apply(get_sample, axis=1)
    df['drug'] = df.apply(only_drug, axis=1)
    df = df[['term', 'drug', 'gse', 'sample', 'gene_list', 'n_genes']]
    df.to_csv('drug_pertubations_from_geo_up.csv')


df = pd.read_csv('drug_pertubations_from_geo_up.csv', index_col=0,
                 converters={"gene_list": lambda x: ast.literal_eval(x)})


def get_dna_damage_drugs():
    terms = [
        'chlorambucil db00291 human gse8832 sample 2691',
        'cisplatin db00515 human gse47856 sample 3146',
        'carboplatin db00958 human gse49577 sample 3177',
        'methotrexate db00563 human gse11440 sample 3072',
        'doxorubicin db00997 human gse12972 sample 2720',
    ]

    data = df[df['term'].isin(terms)]

    terms, gene_sets = [], []

    for t in data._dict('records'):
        terms.append(t['term'].split(' ')[0])
        gene_sets.append(t['gene_list'])

    enrichment = e.run_sample_set_of_dbs(gene_sets, terms, pivot=False)

    enrichment.to_csv('compare_drug_dbs.csv', encoding='utf8')


def run_gse_6907():
    data = df[df['gse'] == 'gse6907'].copy()

    terms, gene_sets = [], []

    for t in data._dict('records'):
        terms.append(t['term'].split(' ')[0])
        gene_sets.append(t['gene_list'])

    enrichment = e.run_sample_set_of_dbs(gene_sets, terms, pivot=False)
    enrichment.to_csv('gse6907_enrichment.csv', encoding='utf8')


for i, g in df.groupby(['gse']):
    drugs = g['drug'].unique()
    n_drugs = len(drugs)
    if n_drugs > 2:
        print("{} : {}".format(i, ','.join(drugs)))
    # print(i, n_drugs)
quit()


def print_containing():
    def p(row):
        if 'cisplatin' in row['term']:
            print(row['term'])


df.apply(p, axis=1)
