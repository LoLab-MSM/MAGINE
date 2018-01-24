from magine.data.datatypes import ExperimentalData
from magine.plotting.venn_diagram_maker import create_venn3

exp_data = ExperimentalData('norris_et_al_2017_cisplatin_data.csv.gz',
                            data_directory='Data')

exp_data.time_series_volcano(exp_data_type='LF',
                             save_name='LF_volcano',
                             p_value=0.05,
                             fold_change_cutoff=1.5,
                             out_dir='plots')

exp_data.create_table_of_data(save_name='cisplatin_data_counts',
                              sig=True)

create_venn3(exp_data.proteomics_over_time[0],
             exp_data.proteomics_over_time[1],
             exp_data.proteomics_over_time[2],
             '1hr', '24hr', '48hr',
             save_name='protein_over_time'
             )
