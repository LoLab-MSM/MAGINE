from magine.plotting.species_plotting import plot_list_of_genes2
import pandas as pd
import os

d_name = os.path.join(os.path.dirname(__file__), 'Data', 'large_example.csv')
# d_name = os.path.join(os.path.dirname(__file__), 'Data', 'example_apoptosis.csv')

d = pd.read_csv(d_name)

# list_species = ['CAV1', 'GPRC5A', 'PLEC']
list_species = ['PARP1', 'TP53']
plot_list_of_genes2(d, list_of_genes=list_species,
                    out_dir='DELETE', save_name='test_me',
                    plot_type='matplotlib')
