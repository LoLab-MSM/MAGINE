import os
import subprocess
from textwrap import wrap

import matplotlib
import numpy as np
import pandas

from magine.html_templates.html_tools import write_single_table

matplotlib.use('Agg')
import matplotlib.pyplot as plt

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


# sample_id = 'sample_id'


class ExperimentalData:
    """
    Manages all experimental data
    """

    def __init__(self, proteomics_file, data_directory):
        raw_df = pandas.read_csv(
                os.path.join(data_directory, proteomics_file),
                parse_dates=False, low_memory=False)

        raw_df = raw_df[raw_df[fold_change].notnull()]
        raw_df.loc[raw_df[fold_change] == np.inf, fold_change] = 1000
        raw_df.loc[raw_df[fold_change] == -np.inf, fold_change] = -1000
        self.data = raw_df
        self.exp_methods = list(self.data[exp_method].unique())
        self.proteomics = raw_df[raw_df[species_type] == protein]
        self.proteomics = self.proteomics.dropna(subset=[gene])
        self.metabolites = []
        self.list_metabolites = []
        self.list_sig_metabolites = []
        self.metabolite_sign = []
        self.exp_methods = self.data[exp_method].unique()
        self.timepoints = sorted(list(self.data[sample_id].unique()))

        if metabolites in self.data[species_type].unique():
            self._set_up_metabolites()

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

        self.gene_fold_change = {}
        self._set_up_gene_fold_change_access()
        self.proteins_measured = self.return_proteomics()

        self.prot_up = self.return_proteomics(significant=True,
                                              fold_change_value=1)
        self.prot_down = self.return_proteomics(significant=True,
                                                fold_change_value=-1)

        self.sign_changed_proteomics = list(set(self.prot_down + self.prot_up))

        self.proteomics_up = {}
        self.proteomics_down = {}
        self.proteomics_sign_changed = {}
        self.protomics_time_points = np.sort(self.proteomics[
                                                 sample_id].unique())
        self.rna_time_points = np.sort(self.rna_seq[sample_id].unique())

        self.proteomics_over_time = []
        self.proteomics_up_over_time = []
        self.proteomics_down_over_time = []
        self.rna_up = {}
        self.rna_down = {}
        self.rna_sign_changed = {}
        self.rna_over_time = []
        self.rna_down_over_time = []
        self.rna_over_time = []
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

        for i in self.rna_time_points:
            self.rna_up[i] = self.return_rna(sample_id_name=i,
                                             significant=True,
                                             fold_change_value=1)
            self.rna_down[i] = self.return_rna(sample_id_name=i,
                                               significant=True,
                                               fold_change_value=-1.)
            self.rna_sign_changed[i] = [set(self.rna_up[i] + self.rna_down[i])]
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
        self.metabolite_sign = self.metabolites[self.metabolites[flag]]
        self.list_metabolites = list(self.metabolites['compound'].unique())
        self.list_sig_metabolites = list(
                self.metabolite_sign['compound'].unique())
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
        tmp = self.proteomics.dropna(subset=[gene])
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
        tmp = self.rna_seq.dropna(subset=[gene])
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

    def _set_up_gene_fold_change_access(self):
        """
        sets up easy access to fold change across time per species
        Returns
        -------

        """
        # tmp = self.proteomics[self.proteomics[flag]]
        tmp = self.proteomics.copy()
        tmp = tmp.dropna(subset=[gene])
        group = tmp.groupby(gene)
        gene_fold_change = {}
        for i, j in group:
            tmp_dict = {}
            group2 = j.groupby(protein)
            for n, m in group2:
                x = np.array(m[sample_id])
                y = np.array(m[fold_change])
                sig_flag = np.array(m[flag])
                tmp_dict[n] = [x, y, sig_flag]
            gene_fold_change[i] = tmp_dict
        self.gene_fold_change = gene_fold_change

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
        measured_table = []
        genes = set()
        meta = set()
        timepoints = list(self.timepoints)
        for i in self.exp_methods:
            tmp = self.data[self.data[exp_method] == i]
            if sig:
                tmp = tmp[tmp[flag]]
            local = set()
            values = []

            for j in timepoints:
                tmp2 = tmp[tmp[sample_id] == j]
                if unique:
                    if i in self.exp_methods_metabolite:
                        tmp2 = tmp2['compound'].unique()
                        for l in list(tmp2):
                            local.add(l)
                            meta.add(l)
                    else:
                        tmp2 = tmp2[gene].unique()
                        for l in list(tmp2):
                            local.add(l)
                            genes.add(l)
                if len(tmp2) == 0:
                    values.append('-')
                else:
                    values.append(len(tmp2))
            if unique:
                values.append(len(local))
            measured_table.append(values)
            if unique:
                print(i, len(local))

        measured_table = np.array(measured_table)

        if unique:
            print('Unique genes = {}'.format(len(genes)))
            print('Unique metabolites = {}'.format(len(meta)))
            timepoints.append('Total Unique')

        t = pandas.DataFrame(data=measured_table,
                             index=self.exp_methods,
                             columns=timepoints)
        from pandas.tools.plotting import table
        ax = plt.subplot(111, frame_on=False)
        table(ax, t, loc='center')
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        # plt.tight_layout()
        plt.savefig('test_table.png')
        t.to_csv('{0}.csv'.format(save_name))
        filename = '{0}.tex'.format(save_name)

        with open(filename, 'wb') as f:
            f.write(template.format(t.to_latex(
                    column_format='*{%s}{c}' % str(len(timepoints) + 1))))

        if _which('pdflatex'):
            print('Compiling table')
            with open(os.devnull, "w") as fnull:
                subprocess.call(['pdflatex', filename],
                                stderr=subprocess.STDOUT,
                                stdout=fnull)
                subprocess.call(['rm',
                                 '{0}.aux'.format(save_name),
                                 '{0}.log'.format(save_name)],
                                stderr=subprocess.STDOUT)
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
            print('Install pdflatex to compile to pdf or png'
                  'You can used the csv file for use in outside tools')

    def plot_all_proteins(self, out_dir='proteins'):
        """
        Creates a plot of all proteins

        Parameters
        ----------
        out_dir: str, path
            Directory that will contain all proteins

        Returns
        -------

        """
        if os.path.exists(out_dir):
            pass
        else:
            os.mkdir(out_dir)
        html_pages = []
        genes_to_plot = self.list_sig_proteins
        for i in np.sort(genes_to_plot):
            self.plot_list_of_genes(list_of_genes=[i], save_name=i,
                                    out_dir=out_dir, title=i, plot_all_x=True)

            html_pages.append(
                    '<a href="{0}/{1}.pdf">{1}</a>'.format(out_dir, i))

        proteins = pandas.DataFrame(html_pages, columns=['Genes'])
        write_single_table(proteins, 'proteins', 'All proteins')

    def plot_list_of_genes(self, list_of_genes, save_name='test', out_dir='.',
                           title=None, plot_all_x=False, log_scale=False):
        """

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
        plot_all_x: bool
            Used if data samples is of time. This ensures all plots have same
            x axis.

        Returns
        -------

        """

        if os.path.exists(out_dir):
            pass
        else:
            os.mkdir(out_dir)
        if type(list_of_genes) in (list, tuple) and save_name == 'test':
            list_of_genes, save_name, out_dir, title, plot_all_x, log_scale = list_of_genes
        x_values = []
        y_values = []
        x_points = np.array(self.timepoints)

        if isinstance(x_points[0], np.float):
            x_point_dict = {i: x_points[n] for n, i
                            in enumerate(self.timepoints)}
        else:
            x_point_dict = {i: 2 ** n + 1 for n, i
                            in enumerate(self.timepoints)}
        sig_flags = []
        labels = []
        for i in sorted(list_of_genes):
            if i in self.gene_fold_change:
                for each in sorted(self.gene_fold_change[i]):
                    x, y, sig = self.gene_fold_change[i][each]
                    if len(x) == 0:
                        continue
                    x_values.append(x)
                    y_values.append(y)
                    sig_flags.append(sig)
                    labels.append(str(each))

        cm = plt.get_cmap('jet')
        num_colors = len(x_values)

        if num_colors != 0:
            ax = plt.subplot(111)
            ax.set_prop_cycle(
                    plt.cycler(
                            'color',
                            [cm(1. * i / num_colors) for i in
                             range(num_colors)]))
            for i in range(num_colors):
                x_index = []
                x_labels = []
                label = "\n".join(wrap(labels[i], 40))
                index = np.argsort(x_values[i])
                x = x_values[i][index]
                y = y_values[i][index]
                s_flag = sig_flags[i][index]

                # y = np.where(y > 0, np.log2(y), -np.log2(-y))

                y[y > 0] = np.log2(y[y > 0])
                y[y < 0] = -np.log2(-y[y < 0])

                for i in x:
                    x_index.append(x_point_dict[i])
                    x_labels.append(i)

                x_index = np.array(x_index)
                p = ax.plot(x_index, y, '.-', label=label)
                if len(s_flag) != 0:
                    color = p[0].get_color()
                    ax.plot(x_index[s_flag], y[s_flag], '^', color=color)

            if title is not None:
                plt.title(title)

            if plot_all_x:
                plt.xlim(min(x_point_dict.values()) - 2,
                         max(x_point_dict.values()) + 2)
                if log_scale:
                    ax.set_xscale('log', basex=2)
                ax.get_xaxis().get_major_formatter().labelOnlyBase = False
                ax.set_xticks(sorted(x_point_dict.values()))
                ax.set_xticklabels(x_points)
                plt.xlabel('Time(hr)')
            else:
                plt.xlim(min(x_point_dict.values()) - 2,
                         max(x_point_dict.values()) + 2)
                if log_scale:
                    ax.set_xscale('log', basex=2)
                    ax.get_xaxis().get_major_formatter().labelOnlyBase = False

                ax.set_xticks(x_index)
                ax.set_xticklabels(x_labels)
                plt.xlabel('Sample index')

            if type(x[0]) == float:
                ax.xaxis.set_major_formatter(
                        matplotlib.ticker.FormatStrFormatter('%.4f'))

            plt.ylabel('log$_2$ Fold Change')
            locs, labels = plt.xticks()
            plt.setp(labels, rotation=90)
            plt.axhline(y=np.log2(1.5), linestyle='--')
            plt.axhline(y=-np.log2(1.5), linestyle='--')
            handles, labels = ax.get_legend_handles_labels()
            lgd = ax.legend(handles, labels, loc='best',
                            bbox_to_anchor=(1.01, 1.0))
            tmp_savename = os.path.join(out_dir, "{}.pdf".format(save_name))
            plt.savefig(tmp_savename, bbox_extra_artists=(lgd,),
                        bbox_inches='tight')
            plt.close()

    def plot_list_of_genes_plotly(self, list_of_genes, save_name='test',
                                  out_dir='.',
                                  title=None, plot_all_x=False,
                                  log_scale=False):

        from plotly.offline import plot

        from matplotlib.pyplot import cm
        if os.path.exists(out_dir):
            pass
        else:
            os.mkdir(out_dir)
        if type(list_of_genes) in (list, tuple) and save_name == 'test':
            list_of_genes, save_name, out_dir, title, plot_all_x, log_scale = list_of_genes
        x_values = []
        y_values = []
        x_points = np.array(self.timepoints)

        if isinstance(x_points[0], np.float):
            x_point_dict = {i: x_points[n] for n, i
                            in enumerate(self.timepoints)}
        else:
            x_point_dict = {i: 2 ** (n + 1) for n, i
                            in enumerate(self.timepoints)}
        sig_flags = []
        labels = []
        for i in sorted(list_of_genes):
            if i in self.gene_fold_change:
                for each in sorted(self.gene_fold_change[i]):
                    x, y, sig = self.gene_fold_change[i][each]
                    if len(x) == 0:
                        continue
                    x_values.append(x)
                    y_values.append(y)
                    sig_flags.append(sig)
                    labels.append(str(each))
        num_colors = len(x_values)
        jet = plt.get_cmap('jet')
        cNorm = colors.Normalize(vmin=0, vmax=num_colors)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

        plotly_list = []
        if num_colors != 0:
            color = cm.rainbow(np.linspace(0, 1, num_colors))
            for i in range(num_colors):
                print(color[i])
                colorVal = scalarMap.to_rgba(i)
                l_color = (
                    'rgba(%4.2f,%4.2f,%4.2f)' % (
                        colorVal[0], colorVal[1], colorVal[2])
                )

                print(l_color)
                x_index = []
                x_labels = []
                label = "\n".join(wrap(labels[i], 40))
                index = np.argsort(x_values[i])
                x = x_values[i][index]
                y = y_values[i][index]
                s_flag = sig_flags[i][index]

                # y = np.where(y > 0, np.log2(y), -np.log2(-y))

                y[y > 0] = np.log2(y[y > 0])
                y[y < 0] = -np.log2(-y[y < 0])

                for k in x:
                    x_index.append(x_point_dict[k])
                    x_labels.append(k)

                x_index = np.array(x_index)

                plotly_list.append(go.Scatter(x=x_index,
                                              y=y,
                                              hoveron='points',
                                              name=label,
                                              mode='lines+markers',
                                              legendgroup='group_{}'.format(i),

                                              marker=dict(symbol='circle',
                                                          color=l_color,
                                                          line=dict(
                                                              color=l_color), ),

                                              )
                                   )
                if len(s_flag) != 0:
                    plotly_list.append(go.Scatter(
                            x=x_index[s_flag],
                            y=y[s_flag],
                            hoveron='points',
                            name=label,
                            legendgroup='group_{}'.format(i),
                            mode='markers',
                            marker=dict(symbol='triangle-up',
                                        size=10,
                                        color=l_color),
                    )
                    )

            if title is not None:
                title = None
            layout = dict(title=title,
                          # yaxis=dict(title='$\\text{log}_2 \\text{Fold Change}$'),
                          xaxis=dict(title='Sample index'))
            fig = dict(data=plotly_list, layout=layout)
            plot(fig, filename='{}.html'.format(save_name), image='png')



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

    def time_series_volcano(self, exp_date_type, save_name, p_value=0.1,
                            out_dir=None, fold_change_cutoff=1.5,
                            y_range=None, x_range=None, bh_critera=False):
        """
        Creates a figure of subplots of provided experimental method

        Parameters
        ----------
        exp_date_type: str
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

        data = self.data[self.data[exp_method] == exp_date_type]
        n_sample = np.sort(data[sample_id].unique())
        if len(n_sample) > 8:
            n_cols = 3
        else:
            n_cols = 2
        n_rows = len(n_sample) / n_cols
        print(n_rows, n_cols)

        fig = plt.figure(figsize=(10, 10))
        for n, i in enumerate(n_sample):
            sample = data[data[sample_id] == i]
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

    def volcano_plot(self, exp_date_type, save_name, out_dir=None,
                     bh_critera=False, p_value=0.1, fold_change_cutoff=1.5,
                     x_range=None, y_range=None):
        """ Create a volcano plot of data
        Creates a volcano plot of data type provided

        Parameters
        ----------
        exp_date_type: str
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

        if exp_date_type not in self.exp_methods:
            print("Must provide experimental method for volcano plot")
            quit()
        data = self.data[self.data[exp_method] == exp_date_type]
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
        tmp[fold_change] = np.where(tmp[fold_change] > 0,
                                    np.log2(tmp[fold_change]),
                                    -np.log2(-tmp[fold_change]))
        # tmp.loc[tmp[fold_change] >= 0, fold_change] = np.log2(tmp[tmp[fold_change] >= 0]][fold_change])
        # tmp.loc[tmp[fold_change] < 0, fold_change] = np.log2(tmp[tmp[fold_change] < 0]][fold_change])
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

    def create_histogram_measurements(self, exp_date_type, save_name,
                                      y_range=None, out_dir=None):
        """
        Plots a histogram of data

        Parameters
        ----------
        exp_date_type: str
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
        if exp_date_type not in self.exp_methods:
            print("Must provide experimental method for volcano plot")
            quit()
        data = self.data[self.data[exp_method] == exp_date_type]
        data = data.dropna(subset=[p_val])
        data = data[np.isfinite(data[fold_change])]
        data = data.dropna(subset=[fold_change])
        tmp = np.array(data[fold_change])

        tmp = np.where(tmp > 0, np.log2(tmp), -np.log2(-tmp))
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(tmp, 50, color='gray')
        if y_range is not None:
            plt.xlim(y_range[0], y_range[1])

        ax.yscale('log', base=10)
        ax.ylabel('Count', fontsize=16)
        ax.xlabel('log$_2$ Fold Change', fontsize=16)
        self._save_plot(fig, save_name, out_dir)


template = r'''
\documentclass[12pt, letterpaper]{{article}}
\usepackage{{booktabs}}
\usepackage{{geometry}}
\usepackage{{pdflscape}}
\usepackage{{nopageno}}
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
    import os

    def _is_exe(filepath):
        return os.path.isfile(filepath) and os.access(filepath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if _is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if _is_exe(exe_file):
                return exe_file

    return None
