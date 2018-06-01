import os

import pandas as pd

import magine.enrichment.tools as et

data_dir = os.path.dirname(__file__)

df = pd.read_csv(os.path.join(data_dir, 'Data', 'enrichr_test_enrichr.csv'))


class TestEnrichmentResult(object):
    def setUp(self):
        d_name = os.path.join(data_dir, 'Data', 'enrichr_test_enrichr.csv')
        self.data = et.load_enrichment_csv(d_name)

    def test_filter_row(self):
        terms = ['apoptotic process',
                 'regulation of mitochondrial membrane potential'],

        # checks if single entry
        slimmed = self.data.filter_rows('term_name', terms[0])
        assert slimmed.shape[0] == 4

        # checks if list
        slimmed = self.data.filter_rows('term_name', terms)
        assert slimmed.shape[0] == 111

    def test_filter_multi(self):
        slimmed = self.data.filter_multi(p_value=0.05, combined_score=20)
        assert slimmed.shape == (20, 10)

    def test_term_to_gene(self):
        genes = self.data.term_to_genes('apoptotic process')
        assert genes == {'CASP8', 'CASP10', 'BCL2', 'BAX', 'CASP3'}

    def test_filter_based_on_word(self):
        slimmed = self.data.filter_based_on_words('apoptotic')
        assert slimmed.shape == (40, 10)

        slimmed = self.data.filter_based_on_words('mitochondrial')
        assert slimmed.shape == (9, 10)

    def test_all_genes(self):
        all_g = self.data.all_genes_from_df()
        assert all_g == {'CASP8', 'CASP10', 'BCL2', 'BAX', 'CASP3'}


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
    slimmed = et.remove_redundant(df, level='all')
    assert slimmed.shape == (3, 10)

    sim2 = et.remove_redundant(df, level='sample')
    assert sim2.shape == (0, 10)

    sim2 = et.remove_redundant(df, level='dataframe', verbose=True)
    assert sim2.shape == (19, 10)


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
