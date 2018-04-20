import os

import pandas as pd

import magine.enrichment.tools as et

data_dir = os.path.dirname(__file__)

df = pd.read_csv(os.path.join(data_dir, 'Data', 'enrichr_test_enrichr.csv'))

print(df.head(10))


def test_filter_row():
    terms = ['apoptotic process',
             'regulation of mitochondrial membrane potential'],

    # checks if single entry
    slimmed = et.filter_rows(df, 'term_name', terms[0])
    assert slimmed.shape[0] == 4

    # checks if list
    slimmed = et.filter_rows(df, 'term_name', terms)
    assert slimmed.shape[0] == 111


test_filter_row()
