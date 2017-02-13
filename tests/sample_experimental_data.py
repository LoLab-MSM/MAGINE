import os

from magine.data_handler import ExperimentalData

data_dir = os.path.dirname(__file__)
exp_data = ExperimentalData(proteomics_file='example_apoptosis.csv',
                            data_directory=os.path.join(data_dir, 'Data'))
