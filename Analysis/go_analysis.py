"""
GO analysis function using orange bioinformatics
"""
import numpy as np
from orangecontrib.bio import go
import scipy.cluster.hierarchy as sch
import os
from textwrap import wrap
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt


# noinspection PyUnresolvedReferences
class GoAnalysis:
    """
    Go analysis class.
    Uses the orangecontrib.bio package
    """

    def __init__(self, species='hsa', slim=False, output_directory='tmp', reference=None):
        self.array = None
        self.names = None
        self.ontology = go.Ontology()
        self.slim = False
        self.reference = reference
        if slim:
            self.slim = True
            self.slim_name = slim
            self.ontology.set_slims_subset(slim)  # goslim_pir goslim_generic
        self.annotations = go.Annotations(species, ontology=self.ontology)
        self.global_go = {}
        self.out_dir = output_directory
        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)

    def create_array_per_time(self, data, aspect='F'):
        """

        :param data:
        :param aspect:
        :return:
        """
        res = self.annotations.get_enriched_terms(data, slims_only=self.slim, aspect=aspect, reference=self.reference)
        # Slims goslim_chembl
        sorted_list = np.ones((len(res.items()), 4), dtype='S50')
        sorted_list[:, 2].astype('float')
        n = 0
        for go_id, (genes, p_value, ref) in res.items():
            if self.slim:
                pass
            elif ref < 10.:
                continue
            if p_value > 0.05:
                continue
            else:
                # sorted_list[n, 2] = float(len(genes) / float(ref)) * 100.
                sorted_list[n, 2] = -1. * np.log10(p_value)
            sorted_list[n, 0] = self.ontology[go_id].name
            sorted_list[n, 1] = go_id
            print(go_id, float(len(genes) / float(ref)) * 100, self.ontology[go_id].name, p_value, len(genes), ref,
                  genes, n)
            self.global_go[go_id] = self.ontology[go_id].name
            sorted_list[n, 3] = str(genes)
            n += 1
        if n == 0:
            print("No significant p-values")
        return sorted_list[:n, :]

    def create_data_set(self, list_of_exp, aspect):
        """

        :param list_of_exp:
        :param aspect:
        :return:
        """
        num_of_lists = len(list_of_exp)
        results = []
        for i in list_of_exp:
            tmp_array = self.create_array_per_time(i, aspect)
            results.append(tmp_array)
        all_terms = []
        for each in results:
            all_terms.extend(list(each[:, 1]))
        all_terms = set(all_terms)
        data_set = np.zeros((len(all_terms), num_of_lists + 1), dtype='S50')
        for n, i in enumerate(all_terms):
            data_set[n, 0] = i
            for num in range(1, num_of_lists + 1):
                data_set[n, num] = return_go_number(i, results[num - 1])
        return data_set

    def find_terms(self, data, term):
        """

        :param data:
        :param term:
        :return:
        """
        shape = np.shape(data)[1]
        new_data = np.zeros((len(data), shape), dtype='S50')
        new_names = data[:, 0]
        n = 0
        for i in range(len(new_names)):
            if len(self.global_go[new_names[i]].split(term)) < 2:
                continue
            new_data[n, :] = data[i, :]
            n += 1
        return new_data[:n - 1, :]

    def plot_heatmap(self, data, all_names, start, stop, savename, labels):
        """

        :param labels:
        :param data:
        :param all_names:
        :param start:
        :param stop:
        :param savename:
        :return:
        """
        matrix = data[start:stop, :]
        size_of_data = np.shape(data)[1]
        length_matrix = len(matrix)
        fig = plt.figure(figsize=(6, 12))
        ax1 = fig.add_subplot(111)
        # plt.title(savename)
        im = ax1.imshow(matrix, aspect=.25, interpolation='nearest', extent=(0, size_of_data, 0, length_matrix + 1),
                        origin='lower')
        plt.colorbar(im)
        names_2 = []
        for i in all_names[start:stop]:
            names_2.append(self.global_go[i])
        y_ticks = np.linspace(.5, length_matrix + .5, length_matrix)
        plt.yticks(y_ticks, names_2, fontsize=16)
        x_ticks = np.linspace(.5, size_of_data - .5, size_of_data)
        if labels:
            plt.xticks(x_ticks, labels, fontsize=16, rotation='90')
        else:
            print("Provide labels")
        plt.savefig(os.path.join(self.out_dir, '%s.pdf' % savename), dpi=300, bbox_inches='tight')
        plt.close()

    def analysis_data(self, data, aspect='F', savename='tmp', labels=None, analyze=True):
        """

        :param analyze:
        :param data:
        :param aspect:
        :param savename:
        :param labels:
        :return:
        """
        data = self.create_data_set(data, aspect)
        # self.print_ranked_over_time(data)
        names, array = sort_data(data)
        self.array = array
        self.names = names
        if not analyze:
            return
        if self.slim:
            savename += "_%s" % self.slim_name
        savename += '_%s' % aspect

        self.plot_heatmap(array, names, 0, 50, '%s_top' % savename, labels)
        self.plot_heatmap(array, names, -50, -1, '%s_bottom' % savename, labels)
        tmp_array = array[:, :].copy()
        fig = plt.figure()
        axdendro = fig.add_axes([0.09, 0.1, 0.2, 0.8])
        linkage = sch.linkage(tmp_array, method='centroid')
        dendrogram = sch.dendrogram(linkage, orientation='right')
        axdendro.set_xticks([])
        axdendro.set_yticks([])
        axmatrix = fig.add_axes([0.3, 0.1, 0.6, 0.8])
        index = dendrogram['leaves']
        tmp_array = tmp_array[index, :]
        names_sorted = names[index]
        im = axmatrix.matshow(tmp_array, aspect='auto', origin='lower')
        axmatrix.set_xticks([])
        axmatrix.set_yticks([])
        axcolor = fig.add_axes([0.91, 0.1, 0.02, 0.8])
        plt.colorbar(im, cax=axcolor)
        fig.savefig(os.path.join(self.out_dir, '%s_dendrogram.pdf' % savename))
        self.plot_heatmap(tmp_array, names_sorted, 0, 50, '%s_top_dendrogram' % savename, labels)
        self.plot_heatmap(tmp_array, names_sorted, -50, -1, '%s_bottom_dendrogram' % savename, labels)

    def print_ranked_over_time(self, data):
        """ Prints information about top hits

        :param data:
        :return:
        """
        for i in range(0, np.shape(data)[1] - 2):
            print(i)
            names, tmp = sort_data_by_index(data, i)
            for j in range(len(data), len(data) - 10):
                print(i, self.global_go[names[j]], tmp[j, 1 + i] - tmp[j, i])

    def find_and_plot_subterms(self, term, savename):
        """

        :param term:
        :param savename:
        """
        print(term)
        terms = self.ontology.extract_sub_graph(term)
        for i in terms:
            if i in self.global_go:
                print(i, self.global_go[i])
        self.plot_specific_go(terms, savename)

    def plot_specific_go(self, term, savename):
        """ Plots a scatter plot of selected terms over time

        :param term:
        :param savename:
        :return:
        """
        found_terms = []
        for i in term:
            if i in self.names:
                # print(self.ontology.term_depth(i),i)
                found_terms.append(i)
                # else:
                #    print("%s is not found" % i)
        num_colors = len(found_terms)
        if num_colors == 0:
            print("Did not find any of the terms")
            return
        fig = plt.figure()
        ax = plt.subplot(111)
        cm = plt.get_cmap('jet')

        ax.set_color_cycle([cm(1. * i / num_colors) for i in range(num_colors)])
        x = [1, 6, 24, 48]
        for i in found_terms:
            y_index = np.where(self.names == i)

            y = self.array[y_index][0]
            label = "\n".join(wrap(self.global_go[i], 40))
            ax.plot(x, y, 'o-', label=label)
        plt.ylabel('-log(p-value)', fontsize=16)
        plt.xlabel('Time(hr)', fontsize=16)
        # ax.legend(loc=0, bbox_to_anchor=(0.5, 1.05),
        #          ncol=1, fancybox=True, shadow=True)
        # box = ax.get_position()
        # ax.set_position([box.x0, box.y0 + box.height * 0.1,
        #                  box.width, box.height * 0.9])
        # ax.legend(loc='center', bbox_to_anchor=(0.5, -0.15),
        #          fancybox=True, shadow=True, ncol=5)
        # #plt.legend(loc=0)
        # plt.tight_layout()
        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles, labels, loc='best', bbox_to_anchor=(1.01, 1.0))
        fig.savefig('%s.png' % savename, bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.show()


def return_go_number(go_term, go_array):
    """

    :param go_term:
    :param go_array:
    :return:
    """
    if go_term in list(go_array[:, 1]):
        for i in range(len(go_array[:, 0])):
            if go_term == go_array[i, 1]:
                return go_array[i, 2]
    else:
        return 0


def sort_data(data):
    """

    :param data:
    :return:
    """
    names = data[:, 0].copy()
    array = data[:, 1:].astype(np.float32).copy()
    # step_size = np.average(array,axis=1).argsort()
    # step_size = (array[:, 2] - array[:, 1]).argsort()
    step_size = (array[:, -1]).argsort()
    names = names[step_size]
    array = array[step_size]
    # array_2 = array[:,1:] - array[:,:-1]
    return names, array


def sort_data_by_index(data, index):
    """

    :param index:
    :param data:
    :return:
    """
    names = data[:, 0].copy()
    array = data[:, 1:].astype(np.float32).copy()
    # step_size = np.average(array,axis=1).argsort()
    step_size = (array[:, index + 1] - array[:, index]).argsort()
    names = names[step_size]
    array = array[step_size]
    return names, array
