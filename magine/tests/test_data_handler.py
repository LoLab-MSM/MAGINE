import os

import pandas as pd

from magine.data.datatypes import ExperimentalData

data_directory = os.path.join(os.path.dirname(__file__), 'Data')

exp_data = ExperimentalData('example_apoptosis.csv', data_directory)


def test_table():
    df = pd.read_csv(os.path.join(data_directory, 'test_data.csv'))
    df['compound_id'] = df['compound']
    exp_data = ExperimentalData(df)
    # exp_data.create_table_of_data(sig=True, save_name='sig_table')
    # exp_data.create_table_of_data(sig=False, save_name='all_table')
    exp_data.create_table_of_data(unique=True, save_name='unique_table')
    exp_data.create_table_of_data(unique=True, sig=True,
                                  save_name='sig_unique_table')


def test_time_series_volcano():
    exp_data.volcano_plot('LF', 'test_label_free',
                          bh_critera=True)


def test_plot_list():
    x = list(exp_data.list_proteins)
    exp_data.plot_list_of_genes(x, 'del_test', 'DELETE', 'test')


def test_html_output():
    exp_data.plot_all_proteins(html_file_name='del', out_dir='DEL')


if __name__ == '__main__':
    test_table()
    # test_time_series_volcano()
    # test_plot_list()
