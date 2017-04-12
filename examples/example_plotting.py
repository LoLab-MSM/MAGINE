import pandas as pd
from magine import ExperimentalData

data = pd.read_csv('Data/norris_et_al_2017_cisplatin_data.csv.gz',
                   low_memory=False)

exp_data = ExperimentalData(data)

print(exp_data.proteomics[exp_data.proteomics['protein'].str.contains('BAX')])
from magine.plotting.species_plotting import plot_list_of_genes2

"""
# create a plotly interactive plot
plot_list_of_genes2(data, ['BAX', 'BID'],
                    save_name='cisplatin_example',
                    out_dir='plots',
                    plot_type='plotly')
"""
# create a matplotlib plot
plot_list_of_genes2(data, ['BAX', 'BID'],
                    save_name='cisplatin_example',
                    out_dir='plots',
                    plot_type='matplotlib',
                    image_format='png')
# """
# create a volcano plot

exp_data.time_series_volcano(exp_data_type='LF',
                             save_name='LF_volcano',
                             p_value=0.05,
                             fold_change_cutoff=1.5,
                             out_dir='plots',

                             )
