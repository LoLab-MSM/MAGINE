import pandas as pd

from magine.data.datatypes import ExperimentalData
from magine.magine_analysis import Analyzer

if __name__ == '__main__':
    data = pd.read_csv('Data/norris_et_al_2017_cisplatin_data.csv.gz',
                       low_memory=False)
    exp_data = ExperimentalData(data)

    magine = Analyzer(experimental_data=exp_data,
                      save_name='cisplatin_go_analysis')

    # Runs go analysis
    magine.run_go_and_create_html('cisplatin_main')
