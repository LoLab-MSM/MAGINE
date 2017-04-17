import os

from magine.data.datatypes import ExperimentalData
from magine.magine_analysis import Analyzer

data_directory = os.path.join(os.path.dirname(__file__), 'Data')


def test_html_array():
    exp_data = ExperimentalData('example_apoptosis.csv', data_directory)
    go_analyzer = Analyzer(exp_data, network=None, species='hsa',
                           output_directory='test_html',
                           save_name='test')

    go_analyzer.run_proteomics_go()



if __name__ == '__main__':
    test_html_array()
