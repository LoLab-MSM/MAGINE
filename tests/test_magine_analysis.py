
import os

from magine.data_handler import ExperimentalData
from magine.magine_analysis import Analyzer

data_directory = os.path.join(os.path.dirname(__file__), 'Data')
exp_data = ExperimentalData('example_apoptosis.csv', data_directory)

if __name__ == '__main__':
    def test_html_array():
        go_analyzer = Analyzer(exp_data, network=None, species='hsa',
                               output_directory='test_html',
                               save_name='test')
        go_analyzer.run_go_and_create_html('html_df_test')


    test_html_array()
