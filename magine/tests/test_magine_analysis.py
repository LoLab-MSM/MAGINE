import os

from magine.data.datatypes import ExperimentalData
from magine.magine_analysis import Analyzer

data_directory = os.path.join(os.path.dirname(__file__), 'Data')
exp_data = ExperimentalData('example_apoptosis.csv', data_directory)


def test_html_array():
    go_analyzer = Analyzer(exp_data, network=None, species='hsa',
                           output_directory='test_html',
                           save_name='test')

    go_analyzer.run_proteomics_go()


def test_run_go():
    go_analyzer = Analyzer(exp_data, network=None, species='hsa',
                           output_directory='test_html',
                           save_name='test')
    data_dict = go_analyzer.run_proteomics_go()


if __name__ == '__main__':

    test_run_go()

    # test_html_array()
