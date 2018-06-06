import os
import shutil
import tempfile

import matplotlib.pyplot as plt
import pandas as pd

from magine.data.experimental_data import ExperimentalData, load_data_csv


class TestExpData(object):
    def setUp(self):
        self._dir = os.path.join(os.path.dirname(__file__), 'Data')
        self.exp_data = load_data_csv(os.path.join(self._dir,
                                                   'example_apoptosis.csv'))
        self.out_dir = tempfile.mkdtemp()

    def tearDown(self):
        self.exp_data = None
        shutil.rmtree(self.out_dir)

    def test_load_from_df(self):
        df = pd.read_csv(os.path.join(self._dir, 'example_apoptosis.csv'))
        exp_data = ExperimentalData(df)

    def test_protein(self):
        assert self.exp_data.proteins.id_list == {'AHR', 'PARP4', 'ADORA1',
                                                  'BAX',
                                               'TP53', 'PARP1', 'ADRA1A',
                                               'ADORA2A', 'CASP3', 'AGTR2'}
        assert self.exp_data.proteins.sig.id_list == {'ADRA1A', 'BAX', 'TP53',
                                                   'AGTR2', 'PARP4', 'PARP1',
                                                   'CASP3'}

        assert self.exp_data.proteins.up.id_list == {'PARP4', 'BAX', 'PARP1',
                                                  'CASP3', 'TP53', 'ADRA1A'}

        assert self.exp_data.proteins.down.id_list == {'AGTR2', 'BAX'}

    def test_gene(self):
        assert self.exp_data.genes.id_list == {'ADORA1', 'PARP4', 'AKT1',
                                               'CASP3',
                                            'ADRA1A', 'AIF1', 'PARP1', 'AGTR2',
                                            'BAX', 'AKT2', 'ADORA2A', 'TP53',
                                            'AHR'}

        assert self.exp_data.genes.sig.id_list == {'TP53', 'AIF1', 'AKT1',
                                                   'BAX',
                                                'CASP3', 'PARP4', 'PARP1',
                                                'AGTR2', 'ADRA1A'}

    def test_rna(self):
        assert self.exp_data.rna.id_list == {'AIF1', 'AKT1', 'AKT2'}
        assert self.exp_data.rna.sig.id_list == {'AIF1', 'AKT1'}

    def test_compounds(self):
        compounds = {'HMDB0009901', 'HMDB0000001'}
        assert self.exp_data.compounds.sig.id_list == compounds
        assert self.exp_data.compounds.id_list == compounds

    def test_species(self):
        assert self.exp_data.species.id_list == {'BID', 'TP53', 'PARP4',
                                                 'HMDB0009901', 'BAX', 'AKT1',
                                                 'AKT2', 'PARP1', 'AGTR2',
                                                 'AHR',
                                                 'AIF1', 'ADORA1', 'CASP3',
                                                 'HMDB0000001', 'ADORA2A',
                                                 'ADRA1A'}

        assert self.exp_data.species.sig.id_list == {'BID', 'PARP4', 'BAX',
                                                     'AKT1', 'PARP1', 'AGTR2',
                                                     'AIF1', 'CASP3', 'TP53',
                                                     'HMDB0000001', 'ADRA1A',
                                                     'HMDB0009901'}

    def test_plot_list(self):
        self.exp_data.plot_species(self.exp_data.rna.sig.id_list,
                                   'del_test', self.out_dir)
        plt.close()

    def test_html_output(self):
        self.exp_data.plot_all_proteins(html_file_name='del',
                                        out_dir=self.out_dir)
        plt.close()

    def test_volcano(self):
        self.exp_data.volcano_analysis(out_dir=self.out_dir)
        plt.close()

    def test_time_series_volcano(self):
        self.exp_data.time_series_volcano('label_free', 'test_label_free',
                                          out_dir=self.out_dir, bh_critera=True)
        plt.close()

    def test_plot_all_metabolites(self):
        self.exp_data.plot_all_metabolites('metab', out_dir=self.out_dir,
                                           plot_type='matplotlib')
        plt.close()

    def test_list_metabolites(self):
        l = ['HMDB0000001', 'HMDB0009901']
        self.exp_data.plot_species(l, save_name='metab', out_dir=self.out_dir,
                                   plot_type='plotly')
        plt.close()

    def test_histogram(self):
        self.exp_data.create_histogram_measurements('label_free',
                                                    save_name='lf_hist',
                                                    y_range=[0, 100],
                                                    out_dir=self.out_dir)
        plt.close()

    def test_table(self):
        self.exp_data.create_table_of_data()
        self.exp_data.create_table_of_data(write_latex=True, save_name='latex')
        self.exp_data.create_table_of_data(sig=True)
        self.exp_data.create_table_of_data(sig=True, unique=True)

    def test_log2(self):
        x = self.exp_data.rna.log2_normalize_df('fold_change')
        assert x.to_dict() == \
               {'sample_id': {8: 'Time_3', 9: 'Time_3', 10: 'Time_3'},
                'source': {8: 'rna_seq', 9: 'rna_seq', 10: 'rna_seq'},
                'significant': {8: True, 9: True, 10: False},
                'fold_change': {8: -0.5849625007211562,
                                9: -1.8073549220576042,
                                10: -0.2630344058337938},
                'p_value': {8: 0.01, 9: 0.06, 10: 0.22},
                'identifier': {8: 'AIF1', 9: 'AKT1', 10: 'AKT2'},
                'species_type': {8: 'protein', 9: 'protein', 10: 'protein'},
                'label': {8: 'AIF1', 9: 'AKT1', 10: 'AKT2'}}

    def test_pivotor(self):
        x = self.exp_data.rna.pivoter(True, columns='sample_id',
                                      values='fold_change', fill_value='-',
                                      min_sig=1)
        assert x.to_dict('index') == {'fold_change':
                                          {'Time_3': -0.8851172762041847}}
