import os
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pandas
from pandas.plotting import table

import magine.plotting.volcano_plots as v_plot
from magine.data.formatter import log2_normalize_df
from magine.plotting.species_plotting import plot_dataframe, \
    plot_list_of_genes, \
    plot_list_of_metabolites

# pandas.set_option('display.max_colwidth', -1)
# column definitions
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
valid_cols = [fold_change, flag, p_val, protein, gene, species_type, sample_id]


# sample_id = 'sample_id'


class ExperimentalData(object):
    """
    Manages all experimental data

    Attributes
    ----------
    data : pd.DataFrame
        Original dataframe provided to class

    proteins : pd.DataFrame
        pandas dataframe containing only species=='proteins'

    exp_methods : list
        List of all experimental methods in data
    list_metabolites
    metabolites
    list_sig_metabolites
    sample_ids
    rna_seq
    proteins_non_rna
    proteins_sig
    rna_seq_sig
    proteins_non_rna_sig
    list_proteins
    list_rna
    list_proteins_non_rna
    list_species
    list_sig_proteins
    list_sig_rna
    list_sig_proteins_non_rna
    list_sig_species
    n_sig_proteins
    n_sig_rna_seq
    n_sig_proteins_non_rna
    n_sig_metabolites
    n_sig_species
    proteins_up
    proteins_down
    sign_changed_proteomics
    proteomics_up
    proteomics_down
    proteomics_sign_changed
    proteomics_sample_ids
    rna_sample_ids
    proteomics_by_sample_id
    proteomics_up_by_sample_id
    proteomics_down_by_sample_id
    sig_species_over_time
    sig_species_up_over_time
    sig_species_down_over_time
    genes_over_time
    genes_up_over_time
    genes_down_over_time
    rna_up
    rna_down
    rna_sig_changed
    rna_over_time
    rna_down_over_time
    rna_up_over_time
    rna_down_over_time


    """

    def __init__(self, data_file, data_directory=os.getcwd(),
                 file_object=None):
        """

        Parameters
        ----------
        data_file : str, pandas.DataFrame
            Name of file, generally csv.
            If provided a str, the file will be read in as a pandas.DataFrame
        data_directory : str
            Directory of csv file
        file_object : object
            File object used for web interface




        """
        if file_object:
            df = pandas.read_csv(
                file_object,
                parse_dates=False, low_memory=False
            )
        elif isinstance(data_file, pandas.DataFrame):
            df = data_file.copy()
        else:
            df = pandas.read_csv(
                os.path.join(data_directory, data_file),
                parse_dates=False, low_memory=False
            )
        for i in valid_cols:
            if i not in df.dtypes:
                print("{} not in columns.")
        df = df[df[fold_change].notnull()]
        df.loc[df[fold_change] == np.inf, fold_change] = 1000
        df.loc[df[fold_change] == -np.inf, fold_change] = -1000
        self.data = df
        self.exp_methods = list(self.data[exp_method].unique())

        self.metabolites = []
        self.list_metabolites = []
        self.list_sig_metabolites = []
        self.metabolite_sign = []
        self.exp_methods_metabolite = []

        if metabolites in self.data[species_type].unique():
            self._set_up_metabolites()
        else:
            self.data['compound_id'] = None
        # return
        self.exp_methods = self.data[exp_method].unique()
        self.sample_ids = sorted(list(self.data[sample_id].unique()))

        self.proteins = df[df[species_type] == protein].dropna(subset=[gene])

        self.proteins_sig = self.proteins[self.proteins[flag]]
        self.rna_seq = self.proteins[self.proteins[exp_method] == rna]
        self.rna_seq_sig = self.rna_seq[self.rna_seq[flag]]

        self.proteins_non_rna = self.proteins[self.proteins[exp_method] != rna]
        self.proteins_non_rna_sig = self.proteins_non_rna[
            self.proteins_non_rna[flag]]

        self.list_proteins = list(self.proteins[gene].astype(str).unique())
        self.list_rna = list(self.rna_seq[gene].unique())
        self.list_proteins_non_rna = list(self.proteins_non_rna[gene].unique())
        self.list_species = list(self.list_proteins + self.list_metabolites)

        self.list_sig_proteins = list(self.proteins_sig[gene].unique())
        self.list_sig_rna = self.return_rna(significant=True)
        self.list_sig_proteins_non_rna = self.return_proteomics(
            significant=True)

        self.list_sig_species = list(self.list_sig_proteins +
                                     self.list_sig_metabolites)

        self.n_sig_proteins = len(self.list_sig_proteins)
        self.n_sig_rna_seq = len(self.list_sig_rna)
        self.n_sig_proteins_non_rna = len(self.list_sig_proteins_non_rna)
        self.n_sig_metabolites = len(self.list_sig_metabolites)
        self.n_sig_species = len(self.list_sig_species)

        self.proteins_up = self.return_proteomics(significant=True,
                                                  fold_change_value=1)
        self.proteins_down = self.return_proteomics(significant=True,
                                                    fold_change_value=-1)

        self.sign_changed_proteomics = [
            set(self.proteins_down + self.proteins_up)]

        self.proteomics_up = dict()
        self.proteomics_down = dict()
        self.proteomics_sign_changed = dict()

        self.proteomics_sample_ids = np.sort(
            self.proteins_non_rna[sample_id].unique())
        self.rna_sample_ids = np.sort(self.rna_seq[sample_id].unique())

        self.proteomics_by_sample_id = []
        self.proteomics_up_by_sample_id = []
        self.proteomics_down_by_sample_id = []

        self.sig_species_over_time = dict()
        self.sig_species_up_over_time = dict()
        self.sig_species_down_over_time = dict()

        self.genes_over_time = []
        self.genes_up_over_time = []
        self.genes_down_over_time = []

        self.rna_up = {}
        self.rna_down = {}
        self.rna_sig_changed = {}
        self.rna_over_time = []
        self.rna_down_over_time = []
        self.rna_up_over_time = []
        self.rna_down_over_time = []
        for i in self.proteomics_sample_ids:
            self.proteomics_up[i] = self.return_proteomics(sample_id_name=i,
                                                           significant=True,
                                                           fold_change_value=1)
            self.proteomics_down[i] = self.return_proteomics(
                sample_id_name=i, significant=True, fold_change_value=-1.
            )
            self.proteomics_sign_changed[i] = list(
                set(self.proteomics_up[i] + self.proteomics_down[i]))

            self.proteomics_by_sample_id.append(
                self.proteomics_sign_changed[i])
            self.proteomics_up_by_sample_id.append(self.proteomics_up[i])
            self.proteomics_down_by_sample_id.append(self.proteomics_down[i])

            self.genes_over_time.append(self.filter(i, True, mol_type='genes'))
            self.genes_up_over_time.append(self.filter(i, True, 'up', 'genes'))
            self.genes_down_over_time.append(
                self.filter(i, True, 'down', 'genes'))

            self.sig_species_over_time[i] = self.filter(i, True)
            self.sig_species_up_over_time[i] = self.filter(i, True, 'up')
            self.sig_species_down_over_time[i] = self.filter(i, True, 'down')

        for i in self.rna_sample_ids:
            self.rna_up[i] = self.return_rna(sample_id_name=i,
                                             significant=True,
                                             fold_change_value=1)
            self.rna_down[i] = self.return_rna(sample_id_name=i,
                                               significant=True,
                                               fold_change_value=-1.)
            self.rna_sig_changed[i] = list(
                set(self.rna_up[i] + self.rna_down[i]))
            self.rna_over_time.append(self.rna_sig_changed[i])
            self.rna_up_over_time.append(self.rna_up[i])
            self.rna_down_over_time.append(self.rna_down[i])

    def _set_up_metabolites(self):
        """
        sets up internal attributes of metabolites
        Returns
        -------

        """
        self.metabolites = self.data[self.data[species_type] == metabolites]
        self.metabolites = self.metabolites.dropna(subset=['compound'])
        self.exp_methods_metabolite = self.metabolites[exp_method].unique()

        if 'compound_id' not in self.metabolites.columns.values:
            self.metabolites['compound_id'] = self.metabolites['compound']
            self.data['compound_id'] = self.data['compound']
            self.data['compound_id'] = self.data['compound_id'].astype(str)

        self.metabolites.loc[:, 'compound_id'] = \
            self.metabolites['compound_id'].astype(str)

        self.metabolite_sign = self.metabolites[self.metabolites[flag]]
        self.list_metabolites = list(self.metabolites['compound_id'].unique())
        self.list_sig_metabolites = \
            list(self.metabolite_sign['compound_id'].astype(str).unique())

        self.metabolites_time_points = \
            np.sort(self.metabolites[sample_id].unique())

    def return_proteomics(self, sample_id_name=0.0, significant=False,
                          fold_change_value=None):
        """
        Returns list of proteins species according to criteria

        Parameters
        ----------
        sample_id_name: str or float
            The sample_id which to filter data
        significant: bool
            Where you want to return significant data or not
        fold_change_value: float
            If you want the positive or negative fold change

        Returns
        -------
        proteins species: list
            List of all proteins species that are of provided criteria

        """
        if sample_id_name == 0.0:
            sample_id_name = False
        tmp = self.proteins.dropna(subset=[gene]).copy()
        tmp = tmp[tmp[exp_method] != rna]
        if significant:
            tmp = tmp[tmp[flag]]
        if fold_change_value is not None:
            if fold_change_value > 0.:
                tmp = tmp[tmp[fold_change] >= fold_change_value]
            else:
                tmp = tmp[tmp[fold_change] <= fold_change_value]
        if sample_id_name:
            return list(tmp[tmp[sample_id] == sample_id_name][gene].unique())
        else:
            return list(tmp[gene].unique())

    def return_rna(self, sample_id_name=None, significant=False,
                   fold_change_value=None):
        """
        Returns list of rna species according to criteria

        Parameters
        ----------
        sample_id_name: str or float
            The sample_id which to filter data
        significant: bool
            Where you want to return significant data or not
        fold_change_value: float
            If you want the positive or negative fold change

        Returns
        -------
        rna species: list
            List of all rna species that are of provided criteria

        """
        if sample_id_name is None:
            sample_id_name = False
        tmp = self.rna_seq.dropna(subset=[gene]).copy()
        if significant:
            tmp = tmp[tmp[flag]]
        if fold_change_value is not None:
            if fold_change_value > 0.:
                tmp = tmp[tmp[fold_change] >= fold_change_value]
            else:
                tmp = tmp[tmp[fold_change] <= fold_change_value]
        if sample_id_name:
            return list(tmp[tmp[sample_id] == sample_id_name][gene].unique())
        else:
            return list(tmp[gene].unique())

    def filter(self, sample_id_name=0.0, significant=False,
               fold_change_value=None, mol_type=None):
        """
        Returns list of proteins species according to criteria

        Parameters
        ----------
        sample_id_name: str or float
            The sample_id which to filter data
        significant: bool
            Where you want to return significant data or not
        fold_change_value: str
            If you want the positive or negative fold change
        mol_type : str
            can be 'gene', 'metabolite', or None for both
        Returns
        -------
        gene list : list
            List of all gene species that match criteria

        """
        if sample_id_name == 0.0:
            sample_id_name = False
        tmp = self.data.copy()

        if sample_id_name:
            tmp = tmp[tmp[sample_id] == sample_id_name]

        if significant:
            tmp = tmp[tmp[flag]]

        if fold_change_value == 'up':
            tmp = tmp[tmp[fold_change] >= 0]
        elif fold_change_value == 'down':
            tmp = tmp[tmp[fold_change] <= 0]
        if mol_type == 'genes':
            return list(tmp[gene].unique())
        elif mol_type == 'metabolite':
            return list(tmp[metabolites].unique())
        else:
            if 'compound_id' not in tmp.dtypes:
                return list(tmp[gene].unique())
            else:
                return list(tmp[gene].unique()) + list(
                    tmp['compound_id'].unique())

    def get_measured_by_datatype(self):
        """
        Returns dict of species per data type

        Returns
        -------
        dict

        """
        return get_measured_by_datatype(self.data)

    def create_table_of_data(self, sig=False, unique=False, save_name=None,
                             plot=False, write_latex=False):
        """
        Creates a summary table of data.


        Parameters
        ----------
        sig: bool
            Flag to summarize significant species only
        save_name: str
            Name to save csv and .tex file
        unique: bool
            If you want to only consider unique species
            ie count gene species rather than PTMs
        plot: bool
            If you want to create a plot of the table
        write_latex: bool
            Create latex file of table


        Returns
        -------
        pandas.DataFrame

        """
        return create_table_of_data(self.data, sig=sig, unique=unique,
                                    save_name=save_name, plot=plot,
                                    write_latex=write_latex)

    def plot_list_of_genes(self, list_of_genes, save_name, out_dir=None,
                           title=None, plot_type='plotly', image_format='png'):
        """
        Creates an HTML table of plots provided a list gene names
        Parameters
        ----------
    
        list_of_genes: list
            List of genes to be plotter
        save_name: str
            Filename to be saved as
        out_dir: str
            Path for output to be saved
        title: str
            Title of plot, useful when list of genes corresponds to a GO term
        plot_type : str
            Use plotly to generate html output or matplotlib to generate pdf
        image_format : str
            pdf or png, only used if plot_type="matplotlib"

        Returns
        -------

        """
        plot_list_of_genes(
            self.proteins, genes=list_of_genes,
            save_name=save_name,
            out_dir=out_dir, title=title, plot_type=plot_type,
            image_format=image_format
        )

    def plot_list_of_metabolites(self, list_of_metabolites, save_name,
                                 out_dir=None, title=None, plot_type='plotly',
                                 image_format='png'):
        """
        Creates an HTML table of plots provided a list metabolites
        
        Parameters
        ----------
        list_of_metabolites : list
            list of compounds 
        save_name : str
            Name of html output file
        out_dir : str
            Location to place plots
        title : str
            Title for HTML page
        plot_type : str
            Type of plot outputs, can be "plotly" or "matplotlib"
        image_format : str
            pdf or png, only used if plot_type="matplotlib"
        Returns
        -------

        """
        if not isinstance(self.metabolites, list):
            plot_list_of_metabolites(
                self.metabolites, list_of_metab=list_of_metabolites,
                save_name=save_name, out_dir=out_dir, title=title,
                plot_type=plot_type, image_format=image_format
            )

    def plot_all_proteins(self, html_file_name, out_dir='proteins',
                          plot_type='plotly', run_parallel=False):
        """
        Creates a plot of all proteins

        Parameters
        ----------
        html_file_name : str
            filename to save html of all plots
        out_dir: str, path
            Directory that will contain all proteins
        plot_type : str
            plotly or matplotlib output
        run_parallel : bool
            Create the plots in parallel
        Returns
        -------

        """

        plot_dataframe(self.data, html_filename=html_file_name,
                       out_dir=out_dir, plot_type=plot_type,
                       type_of_species='protein',
                       run_parallel=run_parallel)

    def plot_all_metabolites(self, html_file_name, out_dir='metabolites',
                             plot_type='plotly', run_parallel=False):
        """
        Creates a plot of all metabolites

        Parameters
        ----------
        html_file_name : str
            filename to save html of all plots
        out_dir: str, path
            Directory that will contain all proteins
        plot_type : str
            plotly or matplotlib output
        run_parallel : bool
            Create the plots in parallel
        Returns
        -------

        """

        plot_dataframe(self.data, html_filename=html_file_name,
                       out_dir=out_dir, plot_type=plot_type,
                       type_of_species='metabolites',
                       run_parallel=run_parallel)

    def volcano_analysis(self, out_dir, use_sig_flag=True,
                         p_value=0.1, fold_change_cutoff=1.5):
        """
        Creates a volcano plot for each experimental method

        Parameters
        ----------
        out_dir: str, path
            Path to where the output figures will be saved
        use_sig_flag: bool
            Use significant flag of data
        p_value: float, optional
            p value criteria for significant
            Will not be used if use_sig_flag
        fold_change_cutoff: float, optional
            fold change criteria for significant
            Will not be used if use_sig_flag

        Returns
        -------

        """
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        for i in self.exp_methods:
            self.volcano_plot(i, i, out_dir=out_dir, bh_critera=use_sig_flag,
                              p_value=p_value,
                              fold_change_cutoff=fold_change_cutoff)

    def time_series_volcano(self, exp_data_type, save_name, p_value=0.1,
                            out_dir=None, fold_change_cutoff=1.5,
                            y_range=None, x_range=None, bh_critera=False):
        """
        Creates a figure of subplots of provided experimental method

        Parameters
        ----------
        exp_data_type: str
            Type of experimental method for plot
        save_name: str
            name to save figure
        out_dir: str, directory
            Location to save figure
        bh_critera: bool, optional
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
        if not self._check_experiment_type_existence(exp_type=exp_data_type):
            return

        data = self.data[self.data[exp_method] == exp_data_type].copy()
        n_sample = np.sort(data[sample_id].unique())

        if len(n_sample) > 8:
            n_cols = 3
        else:
            n_cols = 2
        n_rows = np.rint(np.rint(len(n_sample) / float(n_cols)))
        if n_cols * n_rows < len(n_sample):
            if n_cols >= n_rows:
                n_rows += 1
            else:
                n_cols += 1

        fig = plt.figure(figsize=(4 * n_rows, 3 * n_cols))
        for n, i in enumerate(n_sample):
            sample = data[data[sample_id] == i].copy()

            sample = sample.dropna(subset=[p_val])
            sample = sample[np.isfinite(sample[fold_change])]
            sample = sample.dropna(subset=[fold_change])
            sec_0, sec_1, sec_2 = v_plot.create_mask(sample, bh_critera,
                                                     p_value,
                                                     fold_change_cutoff)

            ax = fig.add_subplot(n_rows, n_cols, n + 1)
            ax.set_title(i)
            v_plot.add_volcano_plot(ax, sec_0, sec_1, sec_2)
            if not bh_critera:
                fc = np.log2(fold_change_cutoff)
                log_p_val = -1 * np.log10(p_value)
                ax.axvline(x=fc, linestyle='--')
                ax.axvline(x=-1 * fc, linestyle='--')
                ax.axhline(y=log_p_val, linestyle='--')
            if y_range is not None:
                ax.set_ylim(y_range[0], y_range[1])
            if x_range is not None:
                ax.set_xlim(x_range[0], x_range[1])

        v_plot.save_plot(fig, save_name=save_name, out_dir=out_dir)

    def volcano_plot(self, exp_data_type, save_name, out_dir=None,
                     bh_critera=False, p_value=0.1, fold_change_cutoff=1.5,
                     x_range=None, y_range=None):
        """ Create a volcano plot of data
        Creates a volcano plot of data type provided

        Parameters
        ----------
        exp_data_type: str
            Type of experimental method for plot
        save_name: str
            name to save figure
        out_dir: str, directory
            Location to save figure
        bh_critera: bool, optional
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

        if not self._check_experiment_type_existence(exp_type=exp_data_type):
            return
        data = self.data[self.data[exp_method] == exp_data_type].copy()
        fig = v_plot.volcano_plot(data, save_name=save_name, out_dir=out_dir,
                                  bh_criteria=bh_critera, p_value=p_value,
                                  fold_change_cutoff=fold_change_cutoff,
                                  x_range=x_range, y_range=y_range)
        return fig

    def create_histogram_measurements(self, exp_data_type, save_name=None,
                                      y_range=None, out_dir=None):
        """
        Plots a histogram of data

        Parameters
        ----------
        exp_data_type: str
            Which data to plot
        save_name: str
            Name of figure
        out_dir: str, path
            Path to location to save figure
        y_range: array_like
            range of data


        Returns
        -------

        """

        if not self._check_experiment_type_existence(exp_type=exp_data_type):
            return
        data = self.data[self.data[exp_method] == exp_data_type]
        data = data.dropna(subset=[p_val])
        data = data[np.isfinite(data[fold_change])]
        data = data.dropna(subset=[fold_change])

        tmp = np.array(log2_normalize_df(data, fold_change)[fold_change])

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(tmp, 50, color='gray')
        if y_range is not None:
            plt.xlim(y_range[0], y_range[1])

        ax.set_yscale('log', basey=10)
        ax.set_xlabel('log$_2$ Fold Change', fontsize=16)
        ax.set_ylabel('Count', fontsize=16)
        if save_name is not None:
            v_plot.save_plot(fig, save_name, out_dir)
        return fig

    def _check_experiment_type_existence(self, exp_type):
        if exp_type not in self.exp_methods:
            print("Must provide experimental method for volcano plot")
            # raise Warning
            return False
        return True


template = r'''
\documentclass[12pt, letterpaper]{{article}}
\usepackage{{booktabs}}
\usepackage{{geometry}}
\usepackage{{pdflscape}}
%\usepackage{{nopageno}}
\geometry{{ papersize={{4.444in,12.681in}},total={{2.8in,6.8in}} }}
\begin{{document}}
\begin{{landscape}}
\begin{{table}}
{}
\end{{table}}
\end{{landscape}}
\end{{document}}
'''


def get_measured_by_datatype(data):
    d = data.copy()
    measured = dict()
    sig_measured = dict()
    for i, r in d.groupby('data_type'):
        if i in ['HILIC', 'C18']:
            key = 'compound_id'
            if key not in r:
                key = 'compound'
        else:
            key = 'gene'
        r = r.dropna(subset=[key])
        measured[i] = set(r[key].unique())
        sig_measured[i] = set(r[r['significant_flag']][key].unique())

    return measured, sig_measured


def create_table_of_data(data, sig=False, unique=False, save_name=None,
                         plot=False, write_latex=False):
    """
    Creates a summary table of data.


    Parameters
    ----------
    data : pandas.DataFrame
    sig: bool
        Flag to summarize significant species only
    save_name: None, str
        Name to save csv and .tex file
    unique: bool
        If you want to only consider unique species
        ie count gene species rather than PTMs
    plot: bool
        If you want to create a plot of the table
    write_latex: bool
        Create latex file of table


    Returns
    -------
    pandas.DataFrame

    """
    tmp_data = data.copy()

    if sig:
        tmp_data = tmp_data[tmp_data[flag]]
    meta_d = tmp_data[tmp_data[species_type] == metabolites].copy()
    exp_methods_metabolite = meta_d[exp_method].unique()
    do_metab = True
    if len(meta_d) == 0:
        do_metab = False

    gene_d = tmp_data[tmp_data[species_type] == protein].copy()

    def _calculate_table(d, value, index):
        d_return = d.pivot_table(values=value, index=index, columns='time',
                                 aggfunc=lambda x: int(x.dropna().nunique()))
        return d_return

    if unique:
        gene_index = 'gene'
        compound_index = 'compound'
    else:
        gene_index = 'protein'
        compound_index = 'compound_id'

    tmp_data2 = _calculate_table(gene_d, gene_index, exp_method)

    if do_metab:
        tmp_data1 = _calculate_table(meta_d, compound_index, exp_method)
        t = pandas.concat([tmp_data1, tmp_data2]).fillna('-')
    else:
        t = tmp_data2.fillna('-')

    unique_col = {}
    for i in t.index:
        loc = tmp_data[tmp_data[exp_method] == i]
        if i in exp_methods_metabolite:
            n = len(loc[compound_index].dropna().unique())
        else:
            n = len(loc[gene_index].dropna().unique())
        unique_col[i] = int(n)

    t['Total Unique Across'] = pandas.Series(unique_col, index=t.index)

    if plot:
        ax = plt.subplot(111, frame_on=False)

        table(ax, t, loc='center')
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        plt.tight_layout()
        if save_name is not None:
            plt.savefig('{}.png'.format(save_name), dpi=300,
                        bbox_inches='tight')

    if save_name is not None:
        t.to_csv('{0}.csv'.format(save_name))
    if write_latex and save_name is not None:
        _write_to_latex(pd_table=t, save_name=save_name)
    return t


def _write_to_latex(pd_table, save_name):
    filename = '{0}.tex'.format(save_name)

    with open(filename, 'wt') as f:
        st = pd_table.to_latex(
            column_format='*{{{}}}{{c}}'.format(str(pd_table.shape[0] + 2)))
        f.write(template.format(st))

    if _which('pdflatex'):
        print('Compiling table')
        with open(os.devnull, "w") as fnull:
            subprocess.call(['pdflatex', filename], stderr=subprocess.STDOUT,
                            stdout=fnull)

            if _which('pdfcrop'):
                pdffile = '{0}.pdf'.format(save_name)
                subprocess.call(['pdfcrop', pdffile, pdffile],
                                stderr=subprocess.STDOUT, stdout=fnull)
                if _which('convert'):
                    tmp_png_name = '{0}.png'.format(save_name)
                    subprocess.call(['convert', '-density', '300', pdffile,
                                     '-quality', '90', tmp_png_name],
                                    stderr=subprocess.STDOUT, stdout=fnull)
    else:
        print('Install pdflatex to compile to pdf or png\n'
              'You can use the csv file for use in outside tools')


def _which(program):
    def _is_exe(filepath):
        return os.path.isfile(filepath) and os.access(filepath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if _is_exe(program):
            return program
        elif _is_exe(program + '.exe'):
            return program + '.exe'
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if _is_exe(exe_file):
                return exe_file

    return None
