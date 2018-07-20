from exp_data import exp_data

from magine.plotting.heatmaps import heatmap_from_array

# create a plotly interactive plot
exp_data.plot_species(['BAX', 'BID'],
                      save_name='cisplatin_example',
                      out_dir='plots',
                      plot_type='plotly')

# create a matplotlib plot
exp_data.plot_species(['BAX', 'BID'],
                      save_name='cisplatin_example',
                      out_dir='plots',
                      plot_type='matplotlib',
                      image_format='png')

# create a volcano plot

exp_data.time_series_volcano(exp_data_type='LF',
                             save_name='LF_volcano',
                             p_value=0.05,
                             fold_change_cutoff=1.5,
                             out_dir='plots',
                             )

# create clustered heatmap of significantly changed proteins only.
sig_data = exp_data.data[exp_data.data['significant_flag']]

heatmap_from_array(sig_data, cluster_row=False, convert_to_log=True,
                   index='protein', values='treated_control_fold_change',
                   columns='time_points', div_colors=True)
