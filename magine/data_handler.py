import os

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas
from textwrap import wrap
import subprocess

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


class ExperimentalData:
    def __init__(self, proteomics_file, data_directory, P_VALUE=0.1,
                 FOLD_CHANGE=.1):
        raw_df = pandas.read_csv(
                os.path.join(data_directory, proteomics_file),
                parse_dates=False, low_memory=False)

        raw_df = raw_df[raw_df[fold_change].notnull()]
        raw_df.loc[raw_df[fold_change] == np.inf, fold_change] = 1000
        raw_df.loc[raw_df[fold_change] == -np.inf, fold_change] = -1000
        self.data = raw_df
        self.exp_methods = list(self.data[exp_method].unique())
        self.proteomics = raw_df[raw_df['species_type'] == protein]
        self.proteomics = self.proteomics.dropna(subset=[gene])
        self.metabolites = []
        self.list_metabolites = []
        self.list_sig_metabolites = []

        if metabolites in self.data['species_type'].unique():
            self.set_up_metabolites()

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
        self.set_up_gene_fold_change_access()
        self.proteins_measured = self.return_proteomics()

        self.prot_up = self.return_proteomics(significant=True,
                                              fold_change_value=1)
        self.prot_down = self.return_proteomics(significant=True,
                                                fold_change_value=-1)

        self.sign_changed_proteomics = list(set(self.prot_down + self.prot_up))

        self.proteomics_up = {}
        self.proteomics_down = {}
        self.proteomics_sign_changed = {}
        self.protomics_time_points = np.sort(self.proteomics['time'].unique())
        self.rna_time_points = np.sort(self.rna_seq['time'].unique())

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
            self.proteomics_up[i] = self.return_proteomics(time=i,
                                                           significant=True,
                                                           fold_change_value=1)
            self.proteomics_down[i] = self.return_proteomics(time=i,
                                                             significant=True,
                                                             fold_change_value=-1.)
            self.proteomics_sign_changed[i] = list(
                    set(self.proteomics_up[i] + self.proteomics_down[i]))
            self.proteomics_over_time.append(self.proteomics_sign_changed[i])
            self.proteomics_up_over_time.append(self.proteomics_up[i])
            self.proteomics_down_over_time.append(self.proteomics_down[i])
        for i in self.rna_time_points:
            self.rna_up[i] = self.return_rna(time=i, significant=True,
                                             fold_change_value=1)
            self.rna_down[i] = self.return_rna(time=i, significant=True,
                                               fold_change_value=-1.)
            self.rna_sign_changed[i] = list(
                    set(self.rna_up[i] + self.rna_down[i]))
            self.rna_over_time.append(self.rna_sign_changed[i])
            self.rna_up_over_time.append(self.rna_up[i])
            self.rna_down_over_time.append(self.rna_down[i])

    def set_up_metabolites(self):
        self.metabolites = self.data[self.data['species_type'] == metabolites]
        self.metabolites = self.metabolites.dropna(subset=['compound'])
        self.metabolite_sign = self.metabolites[self.metabolites[flag]]
        self.list_metabolites = list(self.metabolites['compound'].unique())
        self.list_sig_metabolites = list(
            self.metabolite_sign['compound'].unique())

    def return_proteomics(self, time=0.0, significant=False,
                          fold_change_value=None):
        """ return proteomic data

        Can provide time and/or if significant

        :param fold_change_value:
        :param significant:
        :param time:
        :return:
        """
        if time == 0.0:
            time = False
        tmp = self.proteomics.dropna(subset=[gene])
        tmp = tmp[tmp[exp_method] != rna]
        if significant:
            tmp = tmp[tmp[flag]]
        if fold_change_value is not None:
            if fold_change_value > 0.:
                tmp = tmp[tmp[fold_change] >= fold_change_value]
            else:
                tmp = tmp[tmp[fold_change] <= fold_change_value]
        if time:
            return list(tmp[tmp['time'] == time][gene].unique())
        else:
            return list(tmp[gene].unique())

    def return_rna(self, time=0.0, significant=False, fold_change_value=None):
        """ return proteomic data

        Can provide time and/or if significant

        :param fold_change_value:
        :param significant:
        :param time:
        :return:
        """
        if time == 0.0:
            time = False
        tmp = self.rna_seq.dropna(subset=[gene])
        if significant:
            tmp = tmp[tmp[flag]]
        if fold_change_value is not None:
            if fold_change_value > 0.:
                tmp = tmp[tmp[fold_change] >= fold_change_value]
            else:
                tmp = tmp[tmp[fold_change] <= fold_change_value]
        if time:
            return list(tmp[tmp['time'] == time][gene].unique())
        else:
            return list(tmp[gene].unique())

    def print_info(self):
        """ Prints summary information about data

        :return:
        """
        print(
            "\nNumber of unique genes in measured(non rna seq) = %s" % self.n_sig_proteins_non_rna)
        print(
            "\nNumber of unique genes rna seq  = %s" % self.n_sig_rnaseq)

        print(
        "\nNumber of significantly changed unique genes (RNA&proteomics) = %s" % self.n_sig_proteins)

        print("\nNumber of proteins")
        for i in np.sort(self.protomics_time_points):
            print(
                "{0} hour = {1}".format(i,
                                        len(self.proteomics_sign_changed[i])))

            # print('\nStats on upregulated proteins')
            #
            # for i in np.sort(self.protomics_time_points):
            #     print("{0} hour = {1}".format(i, len(self.proteomics_up[i])))
            #
            # print("\nNumber of down regulated proteins")
            #
            # for i in np.sort(self.protomics_time_points):
            #     print("{0} hour = {1}".format(i, len(self.proteomics_down[i])))

    def set_up_gene_fold_change_access(self):

        tmp = self.proteomics[self.proteomics[flag]]
        tmp = tmp.dropna(subset=[gene])
        group = tmp.groupby(gene)
        gene_fold_change = {}
        for i, j in group:
            tmp_dict = {}
            group2 = j.groupby(protein)
            for n, m in group2:
                x = np.array(m['time'])
                y = np.array(m[fold_change])
                sig_flag = np.array(m[flag])
                tmp_dict[n] = [x, y, sig_flag]
            gene_fold_change[i] = tmp_dict
        self.gene_fold_change = gene_fold_change

    def plot_all_proteins(self, out_dir='proteins'):
        if os.path.exists(out_dir):
            pass
        else:
            os.mkdir(out_dir)
        html_pages = []
        # tmp = self.list_sig_proteins
        # tmp2 = self.lisproteomics_non_rna[gene].unique()
        genes_to_plot = self.list_sig_proteins
        for i in np.sort(genes_to_plot):
            self.create_heat_map_of_gene(out_dir, i)
            html_pages.append(
                    '<a href="{0}/{1}.pdf">{1}</a>'.format(out_dir, i))

        proteins = pandas.DataFrame(html_pages, columns=['Genes'])
        header = '<script src="sorttable_tmp.js"></script>\n<html>\n\t<body>\n'
        footer = "\n\t</body>\n</html>"
        # write out the html file
        with open('proteins.html', 'w') as f:
            f.write(header)
            f.write('\n<div style="float: left">\n')
            f.write(proteins.to_html(classes='sortable', escape=False,
                                     justify='left'))
            f.write(footer)

    def plot_list_of_genes(self, list_of_genes, save_name, title=None):
        X = []
        Y = []
        sig_flags = []
        labels = []
        for i in sorted(list_of_genes):
            if i in self.gene_fold_change:
                for each in sorted(self.gene_fold_change[i]):
                    x, y, sig = self.gene_fold_change[i][each]
                    if len(x) == 0:
                        continue
                    X.append(x)
                    Y.append(y)
                    sig_flags.append(sig)
                    labels.append(str(each))

        cm = plt.get_cmap('jet')
        num_colors = len(X)
        if num_colors != 0:
            ax = plt.subplot(111)
            ax.set_color_cycle(
                    [cm(1. * i / num_colors) for i in range(num_colors)])
            for i in range(num_colors):
                label = "\n".join(wrap(labels[i], 40))
                index = np.argsort(X[i])
                x = X[i][index]
                y = Y[i][index]
                s_flag = sig_flags[i][index]
                y = np.where(y > 0, np.log2(y), -np.log2(-y))
                p = ax.plot(x, y, '.-', label=label)
                if len(s_flag) != 0:
                    color = p[0].get_color()
                    ax.plot(x[s_flag], y[s_flag], '^', color=color)
            plt.xlabel('Time(hr)')
            plt.ylabel('log$_2$ Fold Change')
            if title is not None:
                plt.title(title)

            ax.set_xscale('log', basex=2)
            ax.get_xaxis().get_major_formatter().labelOnlyBase = False
            ax.set_xticks(self.protomics_time_points)
            ax.xaxis.set_major_formatter(
                    matplotlib.ticker.FormatStrFormatter('%.4f'))
            plt.xlim(min(self.protomics_time_points),
                     max(self.protomics_time_points))
            locs, labels = plt.xticks()
            plt.setp(labels, rotation=90)
            plt.axhline(y=np.log2(1.5), linestyle='--')
            plt.axhline(y=-np.log2(1.5), linestyle='--')
            handles, labels = ax.get_legend_handles_labels()
            lgd = ax.legend(handles, labels, loc='best',
                            bbox_to_anchor=(1.01, 1.0))
            plt.savefig('%s' % save_name, bbox_extra_artists=(lgd,),
                        bbox_inches='tight')
            plt.close()

    def create_table_of_data(self, sig=False, save_name='tmp', unique=False):
        # if sig:
        #     proteins = self.proteomics[self.proteomics[flag]]
        #     metabolites = self.metabolites[self.metabolites[flag]]
        # else:
        #     proteins = self.proteomics
        #     metabolites = self.metabolites

        exp_methods = self.data[exp_method].unique()
        timepoints = list(self.data['time'].unique())
        measured_table = []
        genes = set()
        meta = set()
        for i in exp_methods:
            tmp = self.data[self.data[exp_method] == i]
            if sig:
                tmp = tmp[tmp[flag]]
            local = set()
            values = []

            for j in timepoints:
                tmp2 = tmp[tmp['time_points'] == j]
                if unique:
                    if i in ['hilic', 'rplc', 'HILIC', 'C18']:
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
            print('genes', len(genes))
            print('meta', len(meta))
            timepoints.append('Total Unique')
        t = pandas.DataFrame(data=measured_table, index=exp_methods,
                             columns=timepoints)
        filename = '{0}.tex'.format(save_name)
        pdffile = '{0}.pdf'.format(save_name)
        outname = '{0}.png'.format(save_name)

        with open(filename, 'wb') as f:
            f.write(template.format(t.to_latex(
                    column_format='*{%s}{c}' % str(len(timepoints) + 1))))

        with open(os.devnull, "w") as fnull:
            subprocess.call(['pdflatex', filename], stderr=subprocess.STDOUT,
                            stdout=fnull)
            subprocess.call(['pdfcrop', pdffile, pdffile],
                            stderr=subprocess.STDOUT, stdout=fnull)
            subprocess.call(
                    ['convert', '-density', '300', pdffile, '-quality', '90',
                     outname],
                    stderr=subprocess.STDOUT, stdout=fnull)

    def create_heat_map_of_gene(self, out_dir, gene_name):
        X = []
        Y = []
        sig_flags = []
        labels = []
        gene_data = self.proteomics[self.proteomics[gene] == gene_name]
        grouped = gene_data.groupby(protein)
        for i, row in grouped:
            fc = np.array(row[fold_change])
            sig_flag = np.array(row[flag])

            time = np.array(row['time'])

            if len(time) == 0:
                continue
            X.append(time)
            Y.append(fc)
            sig_flags.append(sig_flag)
            labels.append(str(i))
        cm = plt.get_cmap('jet')
        num_colors = len(X)
        if num_colors != 0:
            fig = plt.figure()
            ax = plt.subplot(111)
            ax.set_color_cycle(
                    [cm(1. * i / num_colors) for i in range(num_colors)])
            for i in range(num_colors):
                label = "\n".join(wrap(labels[i], 40))
                index = np.argsort(X[i])
                x = X[i][index]
                y = Y[i][index]
                s_flag = sig_flags[i][index]
                y = np.where(y > 0, np.log2(y), -np.log2(-y))
                p = ax.plot(x, y, '.-', label=label)
                if len(s_flag) != 0:
                    color = p[0].get_color()
                    ax.plot(x[s_flag], y[s_flag], '^', color=color)

            plt.xlabel('Time(hr)')
            plt.ylabel('log$_2$ Fold Change')
            plt.title(gene)
            ax.set_xscale('log', basex=2)
            ax.get_xaxis().get_major_formatter().labelOnlyBase = False
            ax.set_xticks(self.protomics_time_points)
            ax.xaxis.set_major_formatter(
                    matplotlib.ticker.FormatStrFormatter('%.4f'))

            plt.xlim(min(self.protomics_time_points),
                     max(self.protomics_time_points))
            locs, labels = plt.xticks()
            plt.setp(labels, rotation=90)
            plt.axhline(y=np.log2(1.5), linestyle='--')
            plt.axhline(y=-np.log2(1.5), linestyle='--')
            handles, labels = ax.get_legend_handles_labels()
            lgd = ax.legend(handles, labels, loc='best',
                            bbox_to_anchor=(1.01, 1.0))
            fig.savefig('%s/%s.pdf' % (out_dir, gene_name),
                        bbox_extra_artists=(lgd,), bbox_inches='tight')
            plt.close()

    def volcano_analysis(self, out_dir):
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        for i in self.exp_methods:
            self.volcano_plot(i, i, out_dir=out_dir)

    def time_series_volcano(self, exp_date_type, save_name, p_value=0.1,
                            out_dir=None,
                            fold_change_cutoff=1.5, def_range=None,
                            bh_critera=False):

        data = self.data[self.data[exp_method] == exp_date_type]

        n_sample = np.sort(data['time'].unique())
        print(n_sample)
        if len(n_sample) > 8:
            n_cols = 3
        else:
            n_cols = 2
        n_rows = len(n_sample) / n_cols
        print(n_rows, n_cols)

        figure = plt.figure(figsize=(10, 10))
        for n, i in enumerate(n_sample):
            sample = data[data['time'] == i]
            sample = sample.dropna(subset=[p_val])
            sample = sample[np.isfinite(sample[fold_change])]
            sample = sample.dropna(subset=[fold_change])
            tmp = sample.loc[:, (p_val, fold_change, flag)]
            # convert to log10 scale
            tmp[p_val] = np.log10(data[p_val]) * -1
            # convert to log2 space
            tmp[fold_change] = np.where(tmp[fold_change] > 0,
                                        np.log2(tmp[fold_change]),
                                        -np.log2(-tmp[fold_change]))
            if bh_critera:
                sec_0 = tmp[tmp[flag]]
                sec_2 = tmp[~tmp[flag]]
                sec_1 = None
                print('{0} : sec0 = {1} , sec1 = {2}'.format(save_name,
                                                             len(sec_0),
                                                             len(sec_2)))
            else:
                fc = np.log2(fold_change_cutoff)
                p_value = -1 * np.log10(p_value)

                sec_0 = tmp[
                    (tmp[p_val] >= p_value) & (np.abs(tmp[fold_change]) >= fc)]
                sec_1 = tmp[
                    (tmp[p_val] >= p_value) & (np.abs(tmp[fold_change]) < fc)]
                sec_2 = tmp[(tmp[p_val] < p_value)]

                print('{0} : sec0 = {1} , sec1 = {2} , '
                      'sec2 = {3}'.format(save_name, len(sec_0), len(sec_1),
                                          len(sec_2)))
            ax = figure.add_subplot(n_rows, n_cols, n + 1)
            ax.set_title(i)
            self._add_volcano_plot(ax, sec_0, sec_1, sec_2)
            if not bh_critera:
                ax.axvline(x=fc, linestyle='--')
                ax.axvline(x=-1 * fc, linestyle='--')
                ax.axhline(y=p_value, linestyle='--')
            if def_range is not None:
                ax.set_xlim(def_range[0], def_range[1])
        plt.tight_layout()
        if out_dir is None:
            plt.savefig('{0}.pdf'.format(save_name), bbox_inches='tight')
            plt.savefig('{0}.png'.format(save_name), bbox_inches='tight')
        else:
            plt.savefig('{0}/{1}.pdf'.format(out_dir, save_name),
                        bbox_inches='tight')
            plt.savefig('{0}/{1}.png'.format(out_dir, save_name),
                        bbox_inches='tight')
        plt.close()

    def volcano_plot(self, exp_date_type, save_name, p_value=0.1, out_dir=None,
                     fold_change_cutoff=1.5, def_range=None, bh_critera=False):
        # Create a volcano plot
        # We declare signficant p-values as less than 0.05
        # We declare fold change as  less than -1 or greater than 1
        """ create a volcano plot of data

        :param exp_date_type:
        :param save_name:
        :param def_range:
        :param p_value:
        :param fold_change:
        :return:
        """

        if exp_date_type not in self.exp_methods:
            print("Must provide experimental method for volcano plot")
            quit()
        data = self.data[self.data[exp_method] == exp_date_type]

        # fold change cutoff
        fc = np.log2(fold_change_cutoff)
        p_value = -1 * np.log10(p_value)
        data = data.dropna(subset=[p_val])
        data = data[np.isfinite(data[fold_change])]
        data = data.dropna(subset=[fold_change])
        tmp = data.loc[:, (p_val, fold_change, flag)]

        # convert to log10 scale
        tmp[p_val] = np.log10(data[p_val]) * -1
        # convert to log2 space
        tmp[fold_change] = np.where(tmp[fold_change] > 0,
                                    np.log2(tmp[fold_change]),
                                    -np.log2(-tmp[fold_change]))

        #    0    #    1   #      0     #
        #         #        #            #
        #################################
        #         #        #            #
        #    2    #    2   #      2     #
        #         #        #            #
        #################################

        sec_0 = tmp[(tmp[p_val] >= p_value) & (np.abs(tmp[fold_change]) >= fc)]
        sec_1 = tmp[(tmp[p_val] >= p_value) & (np.abs(tmp[fold_change]) < fc)]
        sec_2 = tmp[(tmp[p_val] < p_value)]

        if bh_critera:
            sec_0 = tmp[tmp[flag]]
            sec_2 = tmp[~tmp[flag]]
            sec_1 = None
            print('{0} : sec0 = {1} , sec1 = {2}'.format(save_name,
                                                         len(sec_0),
                                                         len(sec_2)))
        else:
            print('{0} : sec0 = {1} , sec1 = {2} , '
                  'sec2 = {3}'.format(save_name, len(sec_0), len(sec_1),
                                      len(sec_2)))

        fig = plt.figure()
        ax = fig.add_subplot(111)
        self._add_volcano_plot(ax, sec_0, sec_1, sec_2)
        if not bh_critera:
            ax.axvline(x=fc, linestyle='--')
            ax.axvline(x=-1 * fc, linestyle='--')
            ax.axhline(y=p_value, linestyle='--')
        if def_range is not None:
            ax.set_xlim(def_range[0], def_range[1])
        if out_dir is None:
            plt.savefig('{0}.pdf'.format(save_name), bbox_inches='tight')
            plt.savefig('{0}.png'.format(save_name), bbox_inches='tight')
        else:
            plt.savefig('{0}/{1}.pdf'.format(out_dir, save_name),
                        bbox_inches='tight')
            plt.savefig('{0}/{1}.png'.format(out_dir, save_name),
                        bbox_inches='tight')
        plt.close()

    def _add_volcano_plot(self, fig_axis, section_0, section_1, section_2):

        fig_axis.scatter(section_0[fold_change], section_0[p_val], marker='.',
                         color='blue')
        if section_1 is not None:
            fig_axis.scatter(section_1[fold_change], section_1[p_val],
                             marker='.', color='red')
        fig_axis.scatter(section_2[fold_change], section_2[p_val], s=1,
                         marker=',', color='gray')
        fig_axis.set_ylabel('-log10(p-value)', fontsize=16)
        fig_axis.set_xlabel('log2(Fold Change)', fontsize=16)

    def create_histogram_measurments(self, exp_date_type, save_name,
                                     def_range=None):
        # Create a volcano plot
        # We declare signficant p-values as less than 0.05
        # We declare fold change as  less than -1 or greater than 1
        if exp_date_type not in self.exp_methods:
            print("Must provide experimental method for volcano plot")
            quit()
        data = self.data[self.data[exp_method] == exp_date_type]
        data = data.dropna(subset=[p_val])
        data = data[np.isfinite(data[fold_change])]
        data = data.dropna(subset=[fold_change])
        tmp = np.array(data[fold_change])

        tmp = np.where(tmp > 0, np.log2(tmp), -np.log2(-tmp))

        plt.hist(tmp, 50, color='gray')
        if def_range is not None:
            plt.xlim(def_range[0], def_range[1])

        plt.yscale('log', base=10)
        plt.ylabel('Count', fontsize=16)
        plt.xlabel('log2(Fold Change)', fontsize=16)
        plt.savefig('{0}.pdf'.format(save_name), bbox_inches='tight')
        plt.savefig('{0}.png'.format(save_name), bbox_inches='tight')
        plt.close()


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
