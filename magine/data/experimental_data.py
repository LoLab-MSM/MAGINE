import os
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.plotting import table

import magine.plotting.volcano_plots as v_plot
from magine.data import Data
from magine.data.formatter import log2_normalize_df
from magine.plotting.species_plotting import plot_dataframe, plot_species

# pandas.set_option('display.max_colwidth', -1)
# column definitions

fold_change = 'fold_change'
flag = 'significant'
exp_method = 'source'
p_val = 'p_value'
rna = 'rna_seq'
gene = 'gene'
protein = 'protein'
metabolites = 'metabolites'
species_type = 'species_type'
sample_id = 'sample_id'
identifier = 'identifier'
label = 'label'
valid_cols = [fold_change, flag, p_val, species_type, sample_id]


def load_data_csv(file_name, **kwargs):
    """ Load data into EnrichmentResult data class

    Parameters
    ----------
    file_name : str
    kwargs :
        Flags to pass to pandas.

    Returns
    -------
    df : EnrichmentResult

    """
    df = pd.read_csv(file_name, **kwargs)
    df = df[df[fold_change].notnull()]
    return ExperimentalData(df)


class Sample(Data):
    def __init__(self, *args, **kwargs):
        super(Sample, self).__init__(*args, **kwargs)
        self._identifier = identifier
        self._label = label
        self._up = None
        self._down = None
        self._sig = None

    @property
    def _constructor(self):
        return Sample

    @property
    def sample_ids(self):
        return sorted(set(self[sample_id].values))

    @property
    def up(self):
        # return up regulated species
        return self.loc[self[flag] & (self[fold_change] > 0)]

    @property
    def down(self):
        # return down regulated species
        return self.loc[self[flag] & (self[fold_change] < 0)]

    @property
    def sig(self):
        return self.loc[self[flag]]

    @property
    def id_list(self):
        return set(self[self._identifier].values)

    @property
    def label_list(self):
        return set(self[self._label].values)

    @property
    def up_by_sample(self):
        over_time = []
        for i in self.sample_ids:
            cur_slice = self.copy()
            cur_slice = cur_slice.loc[cur_slice[sample_id] == i]
            over_time.append(cur_slice.up.id_list)
        return over_time

    @property
    def down_by_sample(self):
        over_time = []
        for i in self.sample_ids:
            cur_slice = self.copy()
            cur_slice = cur_slice.loc[cur_slice[sample_id] == i]
            over_time.append(cur_slice.down.id_list)
        return over_time

    @property
    def by_sample(self):
        over_time = []
        for i in self.sample_ids:
            cur_slice = self.copy()
            cur_slice = cur_slice.loc[cur_slice[sample_id] == i]
            over_time.append(cur_slice.down.id_list)
        return over_time

    def volcano_plot(self, save_name, out_dir=None, bh_critera=False,
                     p_value=0.1, fold_change_cutoff=1.5, x_range=None,
                     y_range=None):
        """ Create a volcano plot of data
        Creates a volcano plot of data type provided

        Parameters
        ----------
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
        matplotlib.Figure

        """
        fig = v_plot.volcano_plot(self, save_name=save_name, out_dir=out_dir,
                                  bh_criteria=bh_critera, p_value=p_value,
                                  fold_change_cutoff=fold_change_cutoff,
                                  x_range=x_range, y_range=y_range)
        return fig


class ExperimentalData(object):
    """
    Manages all experimental data

    """

    def __init__(self, data_file):
        """

        Parameters
        ----------
        data_file : str, pandas.DataFrame
            Name of file, generally csv.
            If provided a str, the file will be read in as a pandas.DataFrame


        """
        if isinstance(data_file, pd.DataFrame):
            df = data_file.copy()
        else:
            df = pd.read_csv(data_file, parse_dates=False, low_memory=False)
        for i in valid_cols:
            if i not in df.dtypes:
                print("{} not in columns.".format(i))

        self.data = Data(df)
        self._index = 'identifier'
        self.__proteins = None
        self.__genes = None
        self.__species = None
        self.__rna = None
        self.__compounds = None
        for i in self.exp_methods:
            self.__setattr__(i, Sample(
                self.data.loc[self.data[exp_method] == i]))

    def __setattr__(self, name, value):
        super(ExperimentalData, self).__setattr__(name, value)

    def __getitem__(self, name):
        return super(ExperimentalData, self).__getattribute__(name)

    @property
    def genes(self):
        """ All data tagged with gene

        Includes protein and RNA.

        Returns
        -------

        """
        if self.__genes is None:
            tmp = self.data.copy()
            tmp = tmp.loc[tmp[species_type] == protein]

            self.__genes = Sample(tmp)
        return self.__genes

    @property
    def proteins(self):
        """ Protein level data

        Tagged with "gene" identifier that is not RNA

        Returns
        -------

        """
        if self.__proteins is None:
            tmp = self.data.copy()
            tmp = tmp.loc[(self.data[species_type] == protein) &
                          ~(tmp[exp_method] == rna)]
            self.__proteins = Sample(tmp)
        return self.__proteins

    @property
    def rna(self):
        """ RNA level data

        Tagged with "RNA"

        Returns
        -------

        """
        if self.__rna is None:
            tmp = self.data.copy()
            tmp = tmp.loc[tmp[exp_method] == rna]
            self.__rna = Sample(tmp)
        return self.__rna

    @property
    def compounds(self):
        if self.__compounds is None:
            tmp = self.data.copy()
            tmp = tmp.loc[tmp[species_type] == metabolites]
            self.__compounds = Sample(tmp)
        return self.__compounds

    @property
    def species(self):
        if self.__species is None:
            self.__species = Sample(self.data.copy())
        return self.__species

    @property
    def exp_methods(self):
        return list(self.data[exp_method].unique())

    @property
    def sample_ids(self):
        return sorted(list(self.data[sample_id].unique()))

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
        return create_table_of_data(self, sig=sig, unique=unique,
                                    save_name=save_name, plot=plot,
                                    write_latex=write_latex)

    def plot_species(self, species_list, save_name, out_dir=None, title=None,
                     plot_type='plotly', image_format='png'):
        """
        Creates an HTML table of plots provided a list metabolites

        Parameters
        ----------
        species_list : list
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

        return plot_species(
            self.species, species_list=species_list,
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

        plot_dataframe(self.genes, html_filename=html_file_name,
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

        plot_dataframe(self.compounds, html_filename=html_file_name,
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
        data = self.data[self.data[exp_method] == exp_data_type].copy()
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
    data : ExperimentalData
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

    if sig:
        data_copy = data.species.sig.copy()
    else:
        data_copy = data.species.copy()

    if unique:
        index = identifier
    else:
        index = label

    count_table = data_copy.pivot_table(values=index, index=exp_method,
                                        columns=sample_id, fill_value='-',
                                        aggfunc=lambda x: x.dropna().nunique())

    unique_col = {}
    for i in data.exp_methods:
        if unique:
            if sig:
                unique_col[i] = len(set(data[i].sig.label_list))
            else:
                unique_col[i] = len(set(data[i].label_list))
        else:
            if sig:
                unique_col[i] = len(set(data[i].sig.id_list))
            else:
                unique_col[i] = len(set(data[i].id_list))

    count_table['Total Unique Across'] = pd.Series(unique_col,
                                                   index=count_table.index)
    if plot:
        ax = plt.subplot(111, frame_on=False)

        table(ax, count_table, loc='center')
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        plt.tight_layout()
        if save_name is not None:
            plt.savefig('{}.png'.format(save_name), dpi=300,
                        bbox_inches='tight')

    if save_name is not None:
        count_table.to_csv('{0}.csv'.format(save_name))
    if write_latex and save_name is not None:
        _write_to_latex(pd_table=count_table, save_name=save_name)
    return count_table


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
