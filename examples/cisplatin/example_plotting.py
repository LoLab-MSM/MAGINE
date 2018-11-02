from exp_data import exp_data

# create a plotly interactive plot
exp_data.species.plot_species(['BAX', 'BID'],
                              save_name='cisplatin_example',
                              out_dir='plots',
                              plot_type='plotly')

# create a matplotlib plot
exp_data.species.plot_species(['BAX', 'BID'],
                              save_name='cisplatin_example',
                              out_dir='plots',
                              plot_type='matplotlib',
                              image_format='png')

# create a volcano plot
exp_data.label_free.volcano_by_sample(
    save_name='LF_volcano',
    p_value=0.05,
    fold_change_cutoff=1.5,
    out_dir='plots',
)

# create clustered heatmap of significantly changed proteins only.
sig_data = exp_data.species.sig.copy()

sig_data.heatmap(cluster_row=False, convert_to_log=True,
                 columns='sample_id', div_colors=True)
