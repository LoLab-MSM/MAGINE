import os

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas
from textwrap import wrap
import subprocess


class ExperimentalData:
    def __init__(self, proteomics_file, data_directory, P_VALUE=0.1,
                 FOLD_CHANGE=1.5):
        all_prot = pandas.read_csv(
            os.path.join(data_directory, proteomics_file),
            parse_dates=False, low_memory=False)
        # all_prot.dropna(subset=['gene'])
        self.proteomics = all_prot[all_prot['species_type'] == 'protein']
        self.metabolite = all_prot[all_prot['species_type'] == 'metabolite']

        self.proteomics_sign = self.proteomics[
            self.proteomics['significant_flag']]
        self.metabolite_sign = self.metabolite[
            self.metabolite['significant_flag']]
        self.rna_seq = self.proteomics[
            self.proteomics['data_type'] == 'rna_seq']
        self.proteomics_non_rna = self.proteomics[
            self.proteomics['data_type'] != 'rna_seq']

        self.rna_seq_sign = self.proteomics_sign[
            self.proteomics_sign['data_type'] == 'rna_seq']
        self.proteomics_non_rna_sig = self.proteomics_sign[
            self.proteomics_sign['data_type'] != 'rna_seq']

        self.total_unique_genes_proteomics = len(
            self.proteomics_non_rna['gene'].unique())

        self.total_sign_rna = len(
            self.rna_seq[self.rna_seq['significant_flag']]['gene'].unique())
        return
        self.gene_fold_change = {}
        self.set_up_gene_fold_change_access()
        self.proteins_measured = self.return_proteomics()
        self.proteins_sign_measured = self.return_proteomics(significant=True)

        self.prot_up = self.return_proteomics(significant=True,
                                              fold_change_value=FOLD_CHANGE)
        self.prot_down = self.return_proteomics(significant=True,
                                                fold_change_value=-1 * FOLD_CHANGE)
        self.sign_changed_proteomics = list(set(self.prot_down + self.prot_up))
        self.proteomics_up = {}
        self.proteomics_down = {}
        self.proteomics_sign_changed = {}
        self.protomics_time_points = np.sort(self.proteomics['time'].unique())
        self.rna_time_points = np.sort(self.rna_seq['time'].unique())
        print(self.rna_time_points)
        self.proteomics_over_time = []
        self.proteomics_up_over_time = []
        self.proteomics_down_over_time = []
        self.rna_up = {}
        self.rna_down = {}
        self.rna_sign_changed = {}

        self.rna_over_time = []
        self.rna_up_over_time = []
        self.rna_down_over_time = []
        self.rna_over_time = []
        self.rna_up_over_time = []
        self.rna_down_over_time = []
        for i in self.protomics_time_points:
            self.proteomics_up[i] = self.return_proteomics(time=i,
                                                           significant=True,
                                                           fold_change_value=FOLD_CHANGE)
            self.proteomics_down[i] = self.return_proteomics(time=i,
                                                             significant=True,
                                                             fold_change_value=-1. * FOLD_CHANGE)
            self.proteomics_sign_changed[i] = list(
                set(self.proteomics_up[i] + self.proteomics_down[i]))
            self.proteomics_over_time.append(self.proteomics_sign_changed[i])
            self.proteomics_up_over_time.append(self.proteomics_up[i])
            self.proteomics_down_over_time.append(self.proteomics_down[i])
        for i in self.rna_time_points:
            print(i)
            self.rna_up[i] = self.return_rna(time=i, significant=True,
                                             fold_change_value=FOLD_CHANGE)
            self.rna_down[i] = self.return_rna(time=i, significant=True,
                                               fold_change_value=-1. * FOLD_CHANGE)

            self.rna_sign_changed[i] = list(
                set(self.rna_up[i] + self.rna_down[i]))
            self.rna_over_time.append(self.rna_sign_changed[i])
            self.rna_up_over_time.append(self.rna_up[i])
            self.rna_down_over_time.append(self.rna_down[i])

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
        if significant:
            tmp = self.proteomics.dropna(subset=['gene'])
            # tmp = tmp[tmp['data_type'] != 'rna_seq']
            tmp = tmp[tmp['significant_flag']]
            if fold_change_value is not None:
                if fold_change_value > 0.:
                    tmp = tmp[tmp[
                                  'treated_control_fold_change'] >= fold_change_value]
                else:
                    tmp = tmp[tmp[
                                  'treated_control_fold_change'] <= fold_change_value]
            if time:
                return list(tmp[tmp['time'] == time]['gene'].unique())
            else:
                return list(tmp['gene'].unique())
        else:
            tmp = self.proteomics.dropna(subset=['gene'])
            if time:
                return list(tmp[tmp['time'] == time]['gene'].unique())
            else:
                return list(tmp['gene'].unique())

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
        if significant:
            tmp = self.rna_seq.dropna(subset=['gene'])
            tmp = tmp[tmp['significant_flag']]
            if fold_change_value is not None:
                if fold_change_value > 0.:
                    tmp = tmp[tmp[
                                  'treated_control_fold_change'] >= fold_change_value]
                else:
                    tmp = tmp[tmp[
                                  'treated_control_fold_change'] <= fold_change_value]
            if time:
                return list(tmp[tmp['time'] == time]['gene'].unique())
            else:
                return list(tmp['gene'].unique())
        else:
            tmp = self.rna_seq.dropna(subset=['gene'])
            if time:
                return list(tmp[tmp['time'] == time]['gene'].unique())
            else:
                return list(tmp['gene'].unique())

    def print_info(self):
        """ Prints summary information about data

        :return:
        """
        print(
        "\nNumber of unique genes in proteomics, non rnaseq = %s" % self.total_unique_genes_proteomics)
        print(
        "\nNumber of sig , unique rna seq proteomics = %s" % self.total_sign_rna)

        print("\nNumber of significantly change unique proteomics = %s" % len(
            self.sign_changed_proteomics))

        # print("\nNumber of proteins")
        # for i in np.sort(self.protomics_time_points):
        #     print("{0} hour = {1}".format(i, len(self.proteomics_sign_changed[i])))
        #
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

        print(np.shape(self.proteomics))

        tmp = self.proteomics[self.proteomics['significant_flag']]
        tmp = tmp.dropna(subset=['gene'])

        group = tmp.groupby('gene')
        gene_fold_change = {}
        for i, j in group:
            tmp_dict = {}
            group2 = j.groupby('protein')
            for n, m in group2:
                x = np.array(m['time'])
                y = np.array(m['treated_control_fold_change'])
                tmp_dict[n] = [x, y]
            gene_fold_change[i] = tmp_dict
        self.gene_fold_change = gene_fold_change

    # def set_up_gene_fold_change_access(self):
    #
    #     print(np.shape(self.proteomics))
    #
    #     tmp = self.proteomics
    #     tmp = (tmp['significant_flag'] > True) | (tmp['data_type'] > 'rna_seq')
    #     tmp |= self.proteomics['data_type'] != 'rna_seq'
    #     tmp = self.proteomics[tmp]
    #     tmp = tmp.dropna(subset=['gene'])
    #
    #
    #     genes = np.sort(tmp['gene'].unique())
    #     gene_fold_change = {}
    #     for i in genes:
    #         tmp_dict = {}
    #         #group2 = j.groupby('protein')
    #         gene_data = self.proteomics[self.proteomics['gene'] == i]
    #         for n in gene_data['protein'].unique():
    #         #for n, m in group2:
    #             m = gene_data[gene_data['gene'] == n]
    #             x = np.array(m['time'])
    #             y = np.array(m['treated_control_fold_change'])
    #             tmp_dict[n] = [x, y]
    #         gene_fold_change[i] = tmp_dict
    #     self.gene_fold_change = gene_fold_change

    def plot_all_proteins(self, out_dir='proteins'):
        if os.path.exists(out_dir):
            pass
        else:
            os.mkdir(out_dir)
        html_pages = []
        tmp = self.proteomics_sign[
            self.proteomics_sign['data_type'] != 'rna_seq']

        for i in np.sort(tmp['gene'].unique()):
            self.create_heat_map_of_gene('proteins', i)
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
        labels = []
        for i in sorted(list_of_genes):
            if i in self.gene_fold_change:
                for each in sorted(self.gene_fold_change[i]):
                    x, y = self.gene_fold_change[i][each]
                    X.append(x)
                    Y.append(y)
                    labels.append(str(each))

        cm = plt.get_cmap('jet')
        num_colors = len(X)
        if num_colors != 0:
            ax = plt.subplot(111)
            ax.set_color_cycle(
                    [cm(1. * i / num_colors) for i in range(num_colors)])
            for i in range(num_colors):
                label = "\n".join(wrap(labels[i], 40))
                ax.plot(X[i], Y[i], 'o-', label=label)
            plt.xlabel('Time(hr)')
            plt.ylabel('Fold Change')
            plt.title(title)
            ax.set_xscale('log', basex=2)
            ax.get_xaxis().get_major_formatter().labelOnlyBase = False
            ax.set_xticks(self.protomics_time_points)
            ax.xaxis.set_major_formatter(
                matplotlib.ticker.FormatStrFormatter('%.4f'))
            # ax.set_xlim([-2, 60])
            locs, labels = plt.xticks()
            plt.setp(labels, rotation=90)
            plt.axhline(y=1.5, linestyle='--')
            plt.axhline(y=-1.5, linestyle='--')
            handles, labels = ax.get_legend_handles_labels()
            lgd = ax.legend(handles, labels, loc='best',
                            bbox_to_anchor=(1.01, 1.0))
            plt.savefig('%s' % save_name, bbox_extra_artists=(lgd,),
                        bbox_inches='tight')
            plt.close()

    def create_table_of_data(self, sig=False, save_name='tmp', unique=False):
        if sig:
            proteins = self.proteomics[self.proteomics['significant_flag']]
            metabolites = self.metabolite[self.metabolite['significant_flag']]
        else:
            proteins = self.proteomics
            metabolites = self.metabolite

        exp_methods = ['hilic', 'rplc', 'label_free', 'silac', 'ph_silac',
                       'rna_seq']
        timepoints = ['30s', '30min', '1h', '3h', '6h', '12h', '18h', '24h',
                      '36h', '48h', '60h', '72h']
        measured_table = []
        genes = set()
        meta = set()
        for i in exp_methods:
            if i in ['hilic', 'rplc']:
                tmp = metabolites[metabolites['data_type'] == i]
            else:
                tmp = proteins[proteins['data_type'] == i]
            local = set()
            values = []

            for j in timepoints:
                tmp2 = tmp[tmp['time_points'] == j]
                if unique:
                    if i in ['hilic', 'rplc']:
                        tmp2 = tmp2['compound'].unique()
                        for l in list(tmp2):
                            local.add(l)
                            meta.add(l)
                    else:
                        tmp2 = tmp2['gene'].unique()
                        for l in list(tmp2):
                            local.add(l)
                            genes.add(l)
                if len(tmp2) == 0:
                    values.append('-')
                else:
                    values.append(len(tmp2))
            measured_table.append(values)
            print(i, len(local))
        print('genes', len(genes))
        print('meta', len(meta))
        measured_table = np.array(measured_table)
        t = pandas.DataFrame(data=measured_table, index=exp_methods,
                             columns=timepoints)
        filename = '{0}.tex'.format(save_name)
        pdffile = '{0}.pdf'.format(save_name)
        outname = '{0}.png'.format(save_name)

        template = r'''\documentclass[12pt, letterpaper]{{article}}
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

        with open(filename, 'wb') as f:
            # f.write(template.format(t.to_latex(column_format='|*{%s}{c|}|' % str(len(timepoints) + 1))))
            f.write(template.format(t.to_latex(
                column_format='*{%s}{c}' % str(len(timepoints) + 1))))

        with open(os.devnull, "w") as fnull:
            subprocess.call(['pdflatex', filename], stderr=subprocess.STDOUT,
                            stdout=fnull)
            subprocess.call(['pdfcrop', pdffile, pdffile],
                            stderr=subprocess.STDOUT, stdout=fnull)
            # os.system('pdfcrop {0}/go_network_{1}_{2}.pdf {0}/go_network_{1}_{2}_wpr.pdf'.format(directory, prefix, j))
            subprocess.call(
                    ['convert', '-density', '300', pdffile, '-quality', '90',
                     outname],
                    stderr=subprocess.STDOUT, stdout=fnull)

    def create_heat_map_of_gene(self, out_dir, gene):
        X = []
        Y = []
        labels = []

        gene_data = self.proteomics[self.proteomics['gene'] == gene]
        grouped = gene_data.groupby('protein')
        for i, row in grouped:
            fc = np.array(row['treated_control_fold_change'])
            time = np.array(row['time'])
            X.append(time)
            Y.append(fc)
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
                ax.plot(X[i], Y[i], 'o-', label=label)
            plt.xlabel('Time(hr)')
            plt.ylabel('Fold Change')
            plt.title(gene)
            ax.set_xscale('log', basex=2)
            ax.get_xaxis().get_major_formatter().labelOnlyBase = False
            ax.set_xticks(self.protomics_time_points)
            ax.xaxis.set_major_formatter(
                matplotlib.ticker.FormatStrFormatter('%.4f'))
            # ax.set_xlim([-2, 60])
            locs, labels = plt.xticks()
            plt.setp(labels, rotation=90)
            plt.axhline(y=1.5, linestyle='--')
            plt.axhline(y=-1.5, linestyle='--')
            handles, labels = ax.get_legend_handles_labels()
            lgd = ax.legend(handles, labels, loc='best',
                            bbox_to_anchor=(1.01, 1.0))
            fig.savefig('%s/%s.pdf' % (out_dir, gene),
                        bbox_extra_artists=(lgd,), bbox_inches='tight')
            plt.close()


# label_free = [30s,30m,1hr,3hr,6hr,12hr,18hr, 24hr, 36hr,48hr]
# silac = [
# phsilac
# rnaseq

if __name__ == '__main__':
    exp_data = ExperimentalData()
    exp_data.print_info()
    # print(exp_data.proteomics.head(10))

    exp_data.plot_list_of_genes(['SMARCB1', 'APOA5'], 'test.pdf')
