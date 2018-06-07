import os

from magine.data.experimental_data import load_data_csv

data_dir = os.path.dirname(__file__)
exp_data = load_data_csv(os.path.join(os.path.join(data_dir, 'Data'),
                                      'example_apoptosis.csv'))
