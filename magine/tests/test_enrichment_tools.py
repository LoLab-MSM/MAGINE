import os

import matplotlib.figure
import matplotlib.pyplot as plt
from nose.tools import ok_, raises

import magine.enrichment.enrichment_result as et
from magine.data.experimental_data import load_data
from magine.plotting.heatmaps import heatmap_by_terms

data_dir = os.path.dirname(__file__)


class TestEnrichmentResult(object):
    def setUp(self):
        d_name = os.path.join(data_dir, 'Data', 'enrichr_test_enrichr.csv')
        self.data = et.load_enrichment_csv(d_name)

    def test_filter_row(self):
        terms = ['p53 signaling pathway_hsa_hsa04115',
                 'apoptosis_hsa_hsa04210'],

        # checks if single entry
        slimmed = self.data.filter_rows('term_name', terms[0])
        ok_(slimmed.shape[0] == 6)

        # checks if list
        slimmed = self.data.filter_rows('term_name', terms)
        ok_(slimmed.shape[0] == 128)

    def test_filter_multi(self):
        slimmed = self.data.filter_multi(p_value=0.05, combined_score=20)
        ok_(slimmed.shape[0] == 8)

    def test_term_to_gene(self):
        genes = self.data.term_to_genes('apoptosis_hsa_hsa04210')
        ok_(genes == {'CASP8', 'CASP10', 'BCL2', 'BAX', 'CASP3'})

    def test_filter_based_on_word(self):
        slimmed = self.data.filter_based_on_words('p53')
        ok_(slimmed.shape[0] == 5)

        slimmed = self.data.filter_based_on_words(['apop', 'p53'])
        ok_(slimmed.shape[0] == 11)

    def test_all_genes(self):
        all_g = self.data.all_genes_from_df()
        ok_(all_g == {'CASP8', 'CASP10', 'BCL2', 'BAX', 'CASP3'})

    def test_filter_sim_terms(self):
        sim2 = self.data.remove_redundant(level='sample', sort_by='rank')
        ok_(sim2.shape[0] == 7)

        sim2 = self.data.remove_redundant(level='sample')
        ok_(sim2.shape[0] == 7)

        sim2 = self.data.remove_redundant(level='dataframe', verbose=True)
        ok_(sim2.shape[0] == 15)

        copy_data = self.data.copy()
        copy_data.remove_redundant(level='sample', verbose=True, inplace=True)
        ok_(copy_data.shape[0] == 7)

    def test_dist(self):
        # dist = self.data.dist_matrix()
        # assert isinstance(dist, matplotlib.figure.Figure)
        copy_data = self.data.copy()
        copy_data.remove_redundant(level='sample', verbose=False, inplace=True)
        copy_data.dist_matrix()

        copy_data.heatmap(convert_to_log=True, index='term_name',
                          columns='sample_id', figsize=(6, 14),
                          annotate_sig=True, linewidths=.01, cluster_row=True)
        plt.close()
        dist = self.data.dist_matrix(figsize=(3, 3), level='each')
        ok_(isinstance(dist, matplotlib.figure.Figure))
        plt.close()

        df = load_data(os.path.join(os.path.dirname(__file__), 'Data',
                                    'example_apoptosis.csv'))
        terms = [{'BAX'}, {'PARP4'}]
        heatmap_by_terms(df.species, ['1', '2'], terms)
        plt.close()

    @raises(DeprecationWarning)
    def test_deprecated_warning(self):
        # dist = self.data.dist_matrix()
        # assert isinstance(dist, matplotlib.figure.Figure)
        copy_data = self.data.copy()
        copy_data.remove_redundant(level='sample', verbose=False, inplace=True)

        copy_data.heatmap(convert_to_log=True, index='term_name',
                          columns='sample_id', figsize=(6, 14),
                          rank_index=True,
                          )
        plt.close()

    def test_find_similar_terms(self):
        sim = self.data.find_similar_terms('apoptotic process')
        print(sim)

    def test_jaccard_index(self):
        term1 = {'BAX', 'BCL2', 'MCL1', 'CASP3'}
        term2 = {'BAX', 'BCL2', 'MCL1', 'TP53'}
        score = et.jaccard_index(term1, term2)

        ok_(score == 0.6)

    def test_term_to_dict(self):
        g = self.data.term_to_genes_dict(
            ['apoptosis_hsa_hsa04210', 'p53 signaling pathway_hsa_hsa04115']
        )

        ok_(g['apoptosis_hsa_hsa04210'] == {'BCL2', 'CASP10'})
        ok_(g['apoptosis_hsa_hsa04210,p53 signaling pathway_hsa_hsa04115'] == {
            'CASP3', 'CASP8', 'BAX'})
