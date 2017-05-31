import subprocess

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas
from magine.html_templates.html_tools import write_single_table
from magine.plotting.species_plotting import plot_list_of_genes2

pandas.set_option('display.max_colwidth', -1)
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
    """

    def __init__(self, proteomics_file, data_directory=os.getcwd(),
                 file_object=None):
        if file_object:
            raw_df = pandas.read_csv(
                file_object,
                parse_dates=False, low_memory=False
                )
        elif isinstance(proteomics_file, pandas.DataFrame):
            raw_df = proteomics_file.copy()
        else:
            raw_df = pandas.read_csv(
                os.path.join(data_directory, proteomics_file),
                parse_dates=False, low_memory=False
                )
        for i in valid_cols:
            if i not in raw_df.dtypes:
                print("{} not in columns.")
        raw_df = raw_df[raw_df[fold_change].notnull()]
        raw_df.loc[raw_df[fold_change] == np.inf, fold_change] = 1000
        raw_df.loc[raw_df[fold_change] == -np.inf, fold_change] = -1000
        self.data = raw_df
        self.exp_methods = list(self.data[exp_method].unique())

        self.metabolites = []
        self.list_metabolites = []
        self.list_sig_metabolites = []
        self.metabolite_sign = []

        if metabolites in self.data[species_type].unique():
            self._set_up_metabolites()
        else:
            self.data['compound_id'] = None
        # return
        self.proteomics = raw_df[raw_df[species_type] == protein]
        self.proteomics = self.proteomics.dropna(subset=[gene])

        self.exp_methods = self.data[exp_method].unique()
        self.timepoints = sorted(list(self.data[sample_id].unique()))

        self.rna_seq = self.proteomics[self.proteomics[exp_method] == rna]
        self.proteomics_non_rna = self.proteomics[
            self.proteomics[exp_method] != rna]

        self.proteomics_sign = self.proteomics[self.proteomics[flag]]

        self.rna_seq_sign = self.rna_seq[self.rna_seq[flag]]
        self.proteomics_non_rna_sig = self.proteomics_non_rna[
            self.proteomics_non_rna[flag]]

        self.list_proteins = list(self.proteomics[gene].unique())
        self.list_rna = list(self.rna_seq[gene].unique())
        self.list_proteins_non_rna = list(
                self.proteomics_non_rna[gene].unique())

        self.list_species = list(self.list_proteins + self.list_metabolites)

        self.list_sig_proteins = list(self.proteomics_sign[gene].unique())
        self.list_sig_rna = self.return_rna(significant=True)
        self.list_sig_proteins_non_rna = self.return_rna(significant=True)

        self.list_sig_species = list(self.list_sig_proteins +
                                     self.list_sig_metabolites)

        self.n_sig_proteins = len(self.list_sig_proteins)
        self.n_sig_rnaseq = len(self.list_sig_rna)
        self.n_sig_proteins_non_rna = len(self.list_sig_proteins_non_rna)
        self.n_sig_metabolites = len(self.list_sig_metabolites)
        self.n_sig_species = len(self.list_sig_species)

        self.total_unique_genes_proteomics = len(
                self.proteomics_non_rna[gene].unique())

        self.total_sign_rna = len(
                self.rna_seq[self.rna_seq[flag]][gene].unique())

        self.gene_fold_change = None

        self.proteins_measured = self.return_proteomics()

        self.prot_up = self.return_proteomics(significant=True,
                                              fold_change_value=1)
        self.prot_down = self.return_proteomics(significant=True,
                                                fold_change_value=-1)

        self.sign_changed_proteomics = list(set(self.prot_down + self.prot_up))

        self.proteomics_up = {}
        self.proteomics_down = {}
        self.proteomics_sign_changed = {}
        self.protomics_time_points = \
            np.sort(self.proteomics[sample_id].unique())
        self.rna_time_points = np.sort(self.rna_seq[sample_id].unique())

        self.proteomics_over_time = []
        self.proteomics_up_over_time = []
        self.proteomics_down_over_time = []
        self.sig_species_over_time = []
        self.sig_species_up_over_time = []
        self.sig_species_down_over_time = []
        self.genes_over_time = []
        self.genes_up_over_time = []
        self.genes_down_over_time = []
        self.rna_up = {}
        self.rna_down = {}
        self.rna_sign_changed = {}
        self.rna_over_time = []
        self.rna_down_over_time = []
        self.rna_up_over_time = []
        self.rna_down_over_time = []
        for i in self.protomics_time_points:
            self.proteomics_up[i] = self.return_proteomics(sample_id_name=i,
                                                           significant=True,
                                                           fold_change_value=1)
            self.proteomics_down[i] = self.return_proteomics(
                    sample_id_name=i, significant=True, fold_change_value=-1.)
            self.proteomics_sign_changed[i] = list(
                    set(self.proteomics_up[i] + self.proteomics_down[i]))

            self.proteomics_over_time.append(self.proteomics_sign_changed[i])
            self.proteomics_up_over_time.append(self.proteomics_up[i])
            self.proteomics_down_over_time.append(self.proteomics_down[i])

            self.genes_over_time.append(
                    self.filter_measurements(i, True, mol_type='genes'))
            self.genes_up_over_time.append(
                    self.filter_measurements(i, True, 'up', 'genes'))
            self.genes_down_over_time.append(
                    self.filter_measurements(i, True, 'down', 'genes'))

            self.sig_species_over_time.append(
                    self.filter_measurements(i, True))
            self.sig_species_up_over_time.append(
                    self.filter_measurements(i, True, 'up'))
            self.sig_species_down_over_time.append(
                self.filter_measurements(i, True, 'down'))

        for i in self.rna_time_points:
            self.rna_up[i] = self.return_rna(sample_id_name=i,
                                             significant=True,
                                             fold_change_value=1)
            self.rna_down[i] = self.return_rna(sample_id_name=i,
                                               significant=True,
                                               fold_change_value=-1.)
            self.rna_sign_changed[i] = list(
                set(self.rna_up[i] + self.rna_down[i]))
            self.rna_over_time.append(self.rna_sign_changed[i])
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
        if 'compound_id' not in self.metabolites.dtypes:
            self.metabolites['compound_id'] = self.metabolites['compound']
        self.metabolite_sign = self.metabolites[self.metabolites[flag]]
        self.list_metabolites = list(self.metabolites['compound_id'].unique())
        self.list_sig_metabolites = list(
                self.metabolite_sign['compound_id'].unique())
        self.exp_methods_metabolite = self.metabolites[exp_method].unique()

    def return_proteomics(self, sample_id_name=0.0, significant=False,
                          fold_change_value=None):
        """
        Returns list of proteomics species according to criteria

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
        proteomics species: list
            List of all proteomics species that are of provided criteria

        """
        if sample_id_name == 0.0:
            sample_id_name = False
        tmp = self.proteomics.dropna(subset=[gene]).copy()
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

    def filter_measurements(self, sample_id_name=0.0, significant=False,
                            fold_change_value=None, mol_type=None):
        """
        Returns list of proteomics species according to criteria

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

    def create_table_of_data(self, sig=False, save_name='tmp', unique=False):
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


        Returns
        -------

        """
        timepoints = list(self.timepoints)
        tmp_data = self.data.copy()

        if sig:
            tmp_data = tmp_data[tmp_data[flag]]
        meta_d = tmp_data[tmp_data[species_type] == metabolites].copy()
        gene_d = tmp_data[tmp_data[species_type] == protein].copy()

        if unique:

            tmp_data1 = meta_d.pivot_table(values='compound', index=exp_method,
                                           columns='time',
                                           aggfunc=lambda x:
                                           len(x.dropna().unique()))

            tmp_data2 = gene_d.pivot_table(values='gene', index=exp_method,
                                           columns='time',
                                           aggfunc=lambda x:
                                           len(x.dropna().unique()))
            unique_col = []
            for i in self.exp_methods:
                if i in self.exp_methods_metabolite:
                    n = len(tmp_data[tmp_data[exp_method] == i][
                                'compound'].dropna().unique())
                else:
                    n = len(tmp_data[tmp_data[exp_method] == i][
                                'gene'].dropna().unique())
                unique_col.append(n)
        else:

            tmp_data1 = meta_d.pivot_table(values='compound_id',
                                           index=exp_method,
                                           columns='time',
                                           aggfunc=lambda x:
                                           len(x.dropna().unique()))
            tmp_data2 = gene_d.pivot_table(values='protein', index=exp_method,
                                           columns='time',
                                           aggfunc=lambda x:
                                           len(x.dropna().unique()))
            unique_col = []
            for i in self.exp_methods:
                if i in self.exp_methods_metabolite:
                    n = len(tmp_data[tmp_data[exp_method] == i][
                                'compound_id'].dropna().unique())
                else:
                    n = len(tmp_data[tmp_data[exp_method] == i][
                                'protein'].dropna().unique())
                unique_col.append(n)
        t = pandas.concat([tmp_data1, tmp_data2]).fillna('-')

        t['Total Unique Across '] = pandas.Series(unique_col, index=t.index)

        ax = plt.subplot(111, frame_on=False)
        pandas.tools.plotting.table(ax, t, loc='center')
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        # plt.tight_layout()
        plt.savefig('{}.png'.format(save_name))
        plt.close()
        t.to_csv('{0}.csv'.format(save_name))
        filename = '{0}.tex'.format(save_name)

        with open(filename, 'wt') as f:
            f.write(template.format(t.to_latex(
                    column_format='*{}{{c}}'.format(str(len(timepoints) + 1)))))

        if _which('pdflatex'):
        # t = True
        # if t:
            print('Compiling table')
            with open(os.devnull, "w") as fnull:
                subprocess.call(['pdflatex', filename],
                                stderr=subprocess.STDOUT,
                                stdout=fnull)

                # subprocess.call(['rm',
                #                  '{0}.aux'.format(save_name),
                #                  '{0}.log'.format(save_name)],
                #                 stderr=subprocess.STDOUT)
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

    def plot_list_of_genes(self, list_of_genes, save_name, out_dir=None,
                           title=None):

        plot_list_of_genes2(self.proteomics, list_of_genes=list_of_genes,
                            save_name=save_name, out_dir=out_dir, title=title)

    def plot_all_proteins(self, out_dir='proteins', plot_type='plotly'):
        """
        Creates a plot of all proteins

        Parameters
        ----------
        out_dir: str, path
            Directory that will contain all proteins
        plot_type : str
            plotly or matplotlib output

        Returns
        -------

        """
        from magine.plotting.species_plotting import plot_dataframe

        plot_dataframe(self.data, 'proteins', out_dir, plot_type)
        return
        if os.path.exists(out_dir):
            pass
        else:
            os.mkdir(out_dir)
        html_pages = []
        genes_to_plot = self.list_sig_proteins
        for i in np.sort(genes_to_plot):

            plot_list_of_genes2(self.proteomics, list_of_genes=[i], save_name=i,
                                out_dir=out_dir, title=i, plot_type=plot_type)
            if plot_type == 'plotly':
                html_pages.append(
                        '<a href="{0}/{1}.html">{1}</a>'.format(out_dir, i))
            else:
                html_pages.append(
                    '<a href="{0}/{1}.pdf">{1}</a>'.format(out_dir, i))

        proteins = pandas.DataFrame(html_pages, columns=['Genes'])
        write_single_table(proteins, 'proteins', 'All proteins')

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
        fig = plt.figure(figsize=(3 * n_rows, 3 * n_cols))
        for n, i in enumerate(n_sample):
            sample = data[data[sample_id] == i].copy()
            sample = sample.dropna(subset=[p_val])
            sample = sample[np.isfinite(sample[fold_change])]
            sample = sample.dropna(subset=[fold_change])
            filtered_data = self._filter_data(sample, bh_critera, p_value,
                                              fold_change_cutoff)
            sec_0, sec_1, sec_2 = filtered_data
            ax = fig.add_subplot(n_rows, n_cols, n + 1)
            ax.set_title(i)
            self._add_volcano_plot(ax, sec_0, sec_1, sec_2)
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

        self._save_plot(fig, save_name=save_name, out_dir=out_dir)

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
        data = self.data[self.data[exp_method] == exp_data_type]
        data = data.dropna(subset=[p_val])
        data = data[np.isfinite(data[fold_change])]
        filtered_data = self._filter_data(data, bh_critera, p_value,
                                          fold_change_cutoff)
        sec_0, sec_1, sec_2 = filtered_data
        fig = plt.figure()
        ax = fig.add_subplot(111)
        self._add_volcano_plot(ax, sec_0, sec_1, sec_2)
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
        self._save_plot(fig, save_name=save_name, out_dir=out_dir)

    @staticmethod
    def _save_plot(fig, save_name, out_dir):
        fig.tight_layout()
        tmp_save_name_pdf = '{0}.pdf'.format(save_name)
        tmp_save_name_png = '{0}.png'.format(save_name)
        if out_dir is not None:
            tmp_save_name_pdf = os.path.join(out_dir, tmp_save_name_pdf)
            tmp_save_name_png = os.path.join(out_dir, tmp_save_name_png)
        fig.savefig(tmp_save_name_pdf, bbox_inches='tight')
        fig.savefig(tmp_save_name_png, bbox_inches='tight')
        plt.close()

    @staticmethod
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

    @staticmethod
    def _filter_data(data, use_sig=True, p_value=0.1, fold_change_cutoff=1.5):
        tmp = data.loc[:, (p_val, fold_change, flag)]

        # convert to log10 scale
        tmp[p_val] = np.log10(data[p_val]) * -1

        # convert to log2 space
        greater_than = tmp[fold_change] > 0
        less_than = tmp[fold_change] < 0

        tmp.loc[greater_than, fold_change] = \
            np.log2(tmp[greater_than][fold_change])

        tmp.loc[less_than, fold_change] = \
            -np.log2(-tmp[less_than][fold_change])

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

            sec_0 = tmp[
                (tmp[p_val] >= p_value) & (np.abs(tmp[fold_change]) >= fc)]
            sec_1 = tmp[
                (tmp[p_val] >= p_value) & (np.abs(tmp[fold_change]) < fc)]
            sec_2 = tmp[(tmp[p_val] < p_value)]
        return sec_0, sec_1, sec_2

    def create_histogram_measurements(self, exp_data_type, save_name,
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
        tmp = np.array(data[fold_change])

        # tmp = np.where(tmp > 0, np.log2(tmp), -np.log2(-tmp))

        tmp[tmp > 0] = np.log2(tmp[tmp > 0])
        tmp[tmp < 0] = -np.log2(-tmp[tmp < 0])

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(tmp, 50, color='gray')
        if y_range is not None:
            plt.xlim(y_range[0], y_range[1])

        ax.set_yscale('log', base=10)
        ax.set_ylabel('Count', fontsize=16)
        ax.set_xlabel('log$_2$ Fold Change', fontsize=16)
        self._save_plot(fig, save_name, out_dir)

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
\geometry{{ papersize={{4.444in,8.681in}},total={{3.8in,6.8in}} }}
\begin{{document}}
\begin{{landscape}}
\begin{{table}}
{}
\end{{table}}
\end{{landscape}}
\end{{document}}
'''


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
