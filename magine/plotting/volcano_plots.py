import os
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np

from magine.data.formatter import log2_normalize_df

fold_change = 'treated_control_fold_change'
flag = 'significant_flag'
exp_method = 'data_type'
p_val = 'p_value_group_1_and_group_2'
rna = 'rna_seq'
gene = 'gene'
protein = 'protein'
metabolites = 'metabolites'
species_type = 'species_type'
sample_id = 'time'


def filter_data(data, use_sig=True, p_value=0.1, fold_change_cutoff=1.5):
    tmp = data.loc[:, (p_val, fold_change, flag)].copy()
    # convert to log10 scale
    tmp[p_val] = np.log10(data[p_val]) * -1
    # convert to log2 space
    tmp = log2_normalize_df(tmp, fold_change=fold_change)

    # Visual example of volcano plot
    # section 0 are significant criteria

    #    0    #    1   #      0     #
    #         #        #            #
    #################################
    #         #        #            #
    #    2    #    2   #      2     #
    #         #        #            #
    #################################

    if use_sig:
        sec_0 = tmp[tmp[flag]]
        sec_2 = tmp[~tmp[flag]]
        sec_1 = None
    else:
        fc = np.log2(fold_change_cutoff)
        p_value = -1 * np.log10(p_value)
        criteria_1 = tmp[p_val] >= p_value
        sec_0 = tmp[criteria_1 & (np.abs(tmp['log2fc']) >= fc)]
        sec_1 = tmp[criteria_1 & (np.abs(tmp['log2fc']) < fc)]
        sec_2 = tmp[(tmp[p_val] < p_value)]
    return sec_0, sec_1, sec_2


def add_volcano_plot(fig_axis, section_0, section_1, section_2):
    fig_axis.scatter(section_0['log2fc'], section_0[p_val], marker='.',
                     color='blue')
    if section_1 is not None:
        fig_axis.scatter(section_1['log2fc'], section_1[p_val],
                         marker='.', color='red')
    fig_axis.scatter(section_2['log2fc'], section_2[p_val], s=1,
                     marker=',', color='gray')
    fig_axis.set_ylabel('-log$_{10}$ p-value', fontsize=16)
    fig_axis.set_xlabel('log$_2$ Fold Change', fontsize=16)


def volcano_plot(data, save_name, out_dir=None,
                 bh_criteria=False, p_value=0.1, fold_change_cutoff=1.5,
                 x_range=None, y_range=None):
    """ Create a volcano plot of data
    Creates a volcano plot of data type provided

    Parameters
    ----------
    data : pandas.Dataframe
        data to create volcano plots from
    save_name: str
        name to save figure
    out_dir: str, directory
        Location to save figure
    bh_criteria: bool, optional
        If to use significant flags of data
    p_value: float, optional
        Criteria for significant
    fold_change_cutoff: float, optional
        Criteria for significant
    y_range: array_like
        upper and lower bounds of plot in y direction
    x_range: array_like
        upper and lower bounds of plot in x direction

    Returns
    -------


    """

    data = data.dropna(subset=[p_val])
    data = data[np.isfinite(data[fold_change])]
    filtered_data = filter_data(data, bh_criteria, p_value,
                                fold_change_cutoff)
    sec_0, sec_1, sec_2 = filtered_data
    fig = plt.figure()
    ax = fig.add_subplot(111)
    add_volcano_plot(ax, sec_0, sec_1, sec_2)

    if not bh_criteria:
        fc = np.log2(fold_change_cutoff)
        log_p_val = -1 * np.log10(p_value)
        ax.axvline(x=fc, linestyle='--')
        ax.axvline(x=-1 * fc, linestyle='--')
        ax.axhline(y=log_p_val, linestyle='--')
    if y_range is not None:
        ax.set_ylim(y_range[0], y_range[1])
    if x_range is not None:
        ax.set_xlim(x_range[0], x_range[1])
        fig.tight_layout()
    save_plot(fig, save_name=save_name, out_dir=out_dir)


def save_plot(fig, save_name, out_dir=None, image_type='png'):
    """
    
    Parameters
    ----------
    fig : matplotlib.pyplot.figure
        Figure to be saved
    save_name : str
        output file name
    out_dir : str, optional
        output path
    image_type : str, optional
        output type of file, {"png", "pdf", etc..}

    Returns
    -------

    """
    fig.tight_layout()
    save_name = '{}.{}'.format(save_name, image_type)
    if out_dir is not None:
        save_name = os.path.join(out_dir, save_name)
    fig.savefig(save_name, bbox_inches='tight')
    plt.close()
