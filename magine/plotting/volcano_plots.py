import matplotlib.pyplot as plt
import numpy as np
import os

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


def _filter_data(data, use_sig=True, p_value=0.1, fold_change_cutoff=1.5):
    tmp = data.loc[:, (p_val, fold_change, flag)].copy()
    # convert to log10 scale
    tmp[p_val] = np.log10(data[p_val]) * -1
    # convert to log2 space
    greater_than = tmp[fold_change] > 0
    less_than = tmp[fold_change] < 0
    tmp.loc[greater_than, fold_change] = \
        np.log2(tmp[greater_than][fold_change])
    tmp.loc[less_than, fold_change] = -np.log2(-tmp[less_than][fold_change])
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
        sec_0 = tmp[criteria_1 & (np.abs(tmp[fold_change]) >= fc)]
        sec_1 = tmp[criteria_1 & (np.abs(tmp[fold_change]) < fc)]
        sec_2 = tmp[(tmp[p_val] < p_value)]
    return sec_0, sec_1, sec_2


def _add_volcano_plot(fig_axis, section_0, section_1, section_2):
    fig_axis.scatter(section_0[fold_change], section_0[p_val], marker='.',
                     color='blue')
    if section_1 is not None:
        fig_axis.scatter(section_1[fold_change], section_1[p_val],
                         marker='.', color='red')
    fig_axis.scatter(section_2[fold_change], section_2[p_val], s=1,
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
    filtered_data = _filter_data(data, bh_criteria, p_value,
                                 fold_change_cutoff)
    sec_0, sec_1, sec_2 = filtered_data
    fig = plt.figure()
    ax = fig.add_subplot(111)
    _add_volcano_plot(ax, sec_0, sec_1, sec_2)

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
    tmp_save_name_pdf = '{0}.pdf'.format(save_name)
    tmp_save_name_png = '{0}.png'.format(save_name)
    if out_dir is not None:
        tmp_save_name_pdf = os.path.join(out_dir, tmp_save_name_pdf)
        tmp_save_name_png = os.path.join(out_dir, tmp_save_name_png)
    fig.savefig(tmp_save_name_pdf, bbox_inches='tight')
    fig.savefig(tmp_save_name_png, bbox_inches='tight')
    plt.close()
