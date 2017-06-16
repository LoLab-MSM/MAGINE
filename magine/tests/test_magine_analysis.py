import os

from magine.data.datatypes import ExperimentalData
from magine.magine_analysis import Analyzer


def test_html_array():
    data_directory = os.path.join(os.path.dirname(__file__), 'Data')
    exp_data = ExperimentalData('example_apoptosis.csv', data_directory)
    go_analyzer = Analyzer(exp_data, network=None, species='hsa',
                           output_directory='test_html',
                           save_name='test')

    go_analyzer.generate_network('test_2')
    go_analyzer.run_go_and_create_html('test_del', create_figure=False)


if __name__ == '__main__':
    test_html_array()
