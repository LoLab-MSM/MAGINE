import os

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas
from textwrap import wrap

P_VALUE = 0.1
FOLD_CHANGE = 1.5


class ExperimentalData:
    def __init__(self, proteomics_file, data_directory):
        all_prot = pandas.read_csv(os.path.join(data_directory, proteomics_file),
                                   parse_dates=False)
        all_prot.dropna(subset=['gene'])
        self.proteomics = all_prot
        self.proteomics_sign = all_prot[all_prot['p_value_group_1_and_group_2'] < P_VALUE]
        self.proteomics_sign = self.proteomics_sign[np.abs(
            self.proteomics_sign['treated_control_fold_change']) >= FOLD_CHANGE]
        # self.rna_seq = all_rna
        self.gene_fold_change = {}
        self.set_up_gene_fold_change_access()
        self.proteins_measured = self.return_proteomics()
        self.proteins_sign_measured = self.return_proteomics(significant=True)

        self.prot_up = self.return_proteomics(significant=True, fold_change_value=FOLD_CHANGE)
        self.prot_down = self.return_proteomics(significant=True, fold_change_value=-1 * FOLD_CHANGE)
        self.sign_changed_proteomics = set(self.prot_down + self.prot_up)
        self.proteomics_up = {}
        self.proteomics_down = {}
        self.proteomics_sign_changed = {}
        self.protomics_time_points = self.proteomics['time'].unique()
        self.proteomics_over_time = []
        self.proteomics_up_over_time = []
        self.proteomics_down_over_time = []
        for i in self.protomics_time_points:
            self.proteomics_up[i] = self.return_proteomics(time=i, significant=True, fold_change_value=FOLD_CHANGE)
            self.proteomics_down[i] = self.return_proteomics(time=i, significant=True,
                                                             fold_change_value=-1. * FOLD_CHANGE)
            self.proteomics_sign_changed[i] = set(self.proteomics_up[i] + self.proteomics_down[i])
            self.proteomics_over_time.append(self.proteomics_sign_changed[i])
            self.proteomics_up_over_time.append(self.proteomics_up[i])
            self.proteomics_down_over_time.append(self.proteomics_down[i])

    def return_proteomics(self, time=0.0, significant=False, fold_change_value=None):
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
            tmp = self.proteomics_sign.copy()
            if fold_change_value is not None:
                if fold_change_value > 0.:
                    tmp = tmp[tmp['treated_control_fold_change'] >= fold_change_value]
                else:
                    tmp = tmp[tmp['treated_control_fold_change'] <= fold_change_value]
            if time:
                return list(tmp[tmp['time'] == time]['gene'].unique())
            else:
                return list(tmp['gene'].unique())
        else:
            if time:
                return list(self.proteomics[self.proteomics['time'] == time]['gene'].unique())
            else:
                return list(self.proteomics['gene'].unique())

    def print_info(self):
        """ Prints summary information about data

        :return:
        """
        print("\nNumber of unique proteomics = %s" % len(self.proteins_measured))

        print("\nNumber of significantly change unique proteomics = %s" % len(self.sign_changed_proteomics))

        print("\nNumber of proteins")
        for i in self.protomics_time_points:
            print("{0} hour = {1}".format(i, len(self.proteomics_sign_changed[i])))

        print('\nStats on upregulated proteins')

        for i in self.protomics_time_points:
            print("{0} hour = {1}".format(i, len(self.proteomics_up[i])))

        print("\nNumber of down regulated proteins")

        for i in self.protomics_time_points:
            print("{0} hour = {1}".format(i, len(self.proteomics_down[i])))

    def set_up_gene_fold_change_access(self):
        group = self.proteomics.groupby('gene')
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

    def plot_list_of_genes(self, list_of_genes, save_name, title=None):
        X = []
        Y = []
        labels = []
        for i in sorted(list_of_genes):
            if i in self.gene_fold_change:
                for each in self.gene_fold_change[i]:
                    x, y = self.gene_fold_change[i][each]
                    X.append(x)
                    Y.append(y)
                    labels.append(str(i) + ' : ' + str(each))

        cm = plt.get_cmap('jet')
        num_colors = len(X)
        if num_colors != 0:
            fig = plt.figure()
            ax = plt.subplot(111)
            ax.set_color_cycle([cm(1. * i / num_colors) for i in range(num_colors)])
            for i in range(num_colors):
                label = "\n".join(wrap(labels[i], 40))
                ax.plot(X[i], Y[i], 'o-', label=label)
            plt.xlabel('Time(hr)')
            plt.ylabel('Fold Change')
            plt.title(title)
            ax.set_xscale('log')
            ax.get_xaxis().get_major_formatter().labelOnlyBase = False
            ax.set_xticks([.5, 1, 6, 24, 36, 48], )
            ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            locs, labels = plt.xticks()
            plt.setp(labels, rotation=90)
            plt.axhline(y=1.5, linestyle='--')
            plt.axhline(y=-1.5, linestyle='--')
            handles, labels = ax.get_legend_handles_labels()
            lgd = ax.legend(handles, labels, loc='best', bbox_to_anchor=(1.01, 1.0))
            fig.savefig('%s' % save_name, bbox_extra_artists=(lgd,), bbox_inches='tight')
            plt.close()


if __name__ == '__main__':
    exp_data = ExperimentalData()
    exp_data.print_info()
    # print(exp_data.proteomics.head(10))

    exp_data.plot_list_of_genes(['SMARCB1', 'APOA5'], 'test.pdf')
