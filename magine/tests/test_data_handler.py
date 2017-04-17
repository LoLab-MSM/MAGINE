import os

from magine.data.datatypes import ExperimentalData

data_directory = os.path.join(os.path.dirname(__file__), 'Data')

exp_data = ExperimentalData('large_example.csv', data_directory)


def test_table():
    exp_data.create_table_of_data(sig=True, save_name='sig_table')


def test_time_series_volcano():
    exp_data.time_series_volcano('LF', 'test_label_free',
                                 bh_critera=True)


def test_plot_list():
    # x = list(exp_data.list_proteins)
    x = ['CAV1', 'SRSF5', 'H2AFV', 'TPR', 'HIST1H2AA', 'GPRC5A']
    exp_data.plot_list_of_genes(x, 'del_test', 'DELETE', 'test')

def test_html_output():
    exp_data.plot_all_proteins()



if __name__ == '__main__':
    # test_plot_list()
    # test_table()
    test_time_series_volcano()
    # test_plot_list_of_species()
