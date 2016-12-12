import os

from magine.data_handler import ExperimentalData

data_directory = os.path.join(os.path.dirname(__file__), 'Data')
exp_data = ExperimentalData('example_apoptosis.csv', data_directory)


def test_table():
    exp_data.create_table_of_data(sig=True, save_name='sig_table')


def test_time_series_volcano():
    exp_data.time_series_volcano('label_free', 'test_label_free',
                                 bh_critera=True)


test_time_series_volcano()
