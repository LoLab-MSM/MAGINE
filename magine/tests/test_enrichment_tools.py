import os

import matplotlib.figure
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

    def test_filter_sim_terms(self):

        sim2 = self.data.remove_redundant(level='sample')
        assert sim2.shape == (18, 10)

        slimmed = self.data.remove_redundant(level='all')
        assert slimmed.shape == (18, 10)

        sim2 = self.data.remove_redundant(level='dataframe', verbose=True)
        assert sim2.shape == (29, 10)

    def test_dist(self):

        dist = self.data.dist_matrix()

        assert isinstance(dist, matplotlib.figure.Figure)

        dist = self.data.dist_matrix(fig_size=(3, 3), level='each')
        assert isinstance(dist, matplotlib.figure.Figure)

    def test_find_similar_terms(self):
        sim = self.data.find_similar_terms('apoptotic process')
        print(sim)

    def test_jaccard_index(self):
        term1 = ['BAX', 'BCL2', 'MCL1', 'CASP3']
        term2 = ['BAX', 'BCL2', 'MCL1', 'TP53']
        score = self.data.jaccard_index(term1, term2)

        assert score == 0.6
