import os

from magine import load_data_csv

file_path = os.path.join(os.path.dirname(__file__), 'Data',
                         'norris_et_al_2017_cisplatin_data.csv.gz')

exp_data = load_data_csv(file_path, low_memory=False)
