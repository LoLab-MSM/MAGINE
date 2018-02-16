import os

import pandas as pd

from magine.data.datatypes import ExperimentalData

_file = os.path.join(os.path.dirname(__file__), 'Data',
                     'norris_et_al_2017_cisplatin_exp_data.csv.gz')
df = pd.read_csv(_file, low_memory=False)

exp_data = ExperimentalData(df)
