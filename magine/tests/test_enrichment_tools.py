import os

import pandas as pd

import magine.enrichment.tools as et

data_dir = os.path.dirname(__file__)

df = pd.read_csv(os.path.join(data_dir, 'Data', 'enrichr_test_enrichr.csv'))


def test_filter_row():
    terms = ['apoptotic process',
             'regulation of mitochondrial membrane potential'],

    # checks if single entry
    slimmed = et.filter_rows(df, 'term_name', terms[0])
    assert slimmed.shape[0] == 4

    # checks if list
    slimmed = et.filter_rows(df, 'term_name', terms)
    assert slimmed.shape[0] == 111


def test_jaccard_index():
    term1 = ['BAX', 'BCL2', 'MCL1', 'CASP3']
    term2 = ['BAX', 'BCL2', 'MCL1', 'TP53']
    score = et.jaccard_index(term1, term2)

    assert score == 0.6


def test_filter_sim_terms():
    slimmed = et.filter_similar_terms(df)
    assert slimmed.shape == (3, 10)

    sim2 = et.remove_redundant(df)
    assert sim2.shape == (0, 10)


def test_genes_from_df():
    genes = et.all_genes_from_df(df)
    assert len(genes) == 5


def test_find_similar_terms():
    sim = et.find_similar_terms('apoptotic process', df)


def test_dist():
    import matplotlib.figure
    dist = et.dist_matrix(df)
    assert isinstance(dist, matplotlib.figure.Figure)


def test_dist2():
    import matplotlib.figure
    dist = et.dist_matrix2(df)
    assert isinstance(dist, matplotlib.figure.Figure)
