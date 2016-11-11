"""
GO analysis function using orange bioinformatics
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
from orangecontrib.bio import go
from textwrap import wrap

pd.set_option('display.max_colwidth', -1)

evidence_codes = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC', 'IEA']
# noinspection PyUnresolvedReferences
class GoAnalysis:
    """
    Go analysis class.
    Uses the orangecontrib.bio package
    """

    def __init__(self, species='hsa', slim=False, output_directory='tmp', reference=None, metric='enrichment',
                 verbose=False, experimental_data=None, save_png=False, save_name='tmp'):

        self.array = None
        self.exp_data = experimental_data
        self.data = None
        self.data_2 = None
        self.names = None
        self.savename = save_name
        self.go_terms = None
        self.ontology = go.Ontology()
        self.slim = False
        self.reference = reference
        self.metric = metric
        self.verbose = verbose
        if slim:
            self.slim = True
            self.slim_name = slim
            self.savename += '_slim'
            self.ontology.set_slims_subset(slim)  # goslim_pir goslim_generic goslim_chembl
        self.annotations = go.Annotations(species, ontology=self.ontology)
        self.number_of_total_reference_genes = len(self.annotations.gene_names)
        self.global_go = {}
        self.out_dir = output_directory
        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)
        if not os.path.exists(os.path.join(self.out_dir, 'Figures')):
            os.mkdir(os.path.join(self.out_dir, 'Figures'))
        self.top_hits = []
        self.num_data_sets = 0
        self.save_png = save_png
        self.created_go_pds = set()


    def create_array_per_time(self, data, aspect='F', num=0):
        """

        :param data:
        :param aspect:
        :return:
        """
        genes_present = set()
        for i in data:
            if type(i) != str:
                continue
            elif i in self.annotations.gene_names:
                genes_present.add(i)
            else:
                split_name = i.split(',')
                if len(split_name) > 1:
                    for n in split_name:
                        if n in self.annotations.gene_names:
                            genes_present.add(n)
        number_of_genes = len(genes_present)
        print("Number of genes given = {0}. Number of genes in GO = {1}".format(len(data), number_of_genes))
        res = self.annotations.get_enriched_terms(genes_present, slims_only=self.slim, aspect=aspect,
                                                  # evidence_codes=evidence_codes,
                                                  reference=self.reference)

        sorted_list = np.ones((len(res.items()), 4), dtype='S50')
        sorted_list[:, 2].astype('float')

        n = 0
        sorted_list_2 = []
        for go_id, (genes, p_value, ref) in res.items():
            tmp_entry = []
            if self.slim:
                pass
            # elif ref < 5.:
            #     continue
            if p_value > 0.05:
                continue
            # else:

            if self.metric == 'fraction':
                score = float(len(genes) / float(ref)) * 100.
            if self.metric == 'pvalue':
                score = -1. * np.log10(p_value)
            if self.metric == 'enrichment':
                # expected_value = number_of_genes * num_ref / num_total_reference_genes
                expected_value = number_of_genes * float(
                    ref) / self.number_of_total_reference_genes
                # print(
                # 'number of reference', self.number_of_total_reference_genes)
                # pvalue_2 = scipy.stats.binom_test(len(genes),
                #                                  number_of_genes,
                #                                  float(ref) / self.number_of_total_reference_genes)
                enrichment = float(len(genes)) / expected_value
                # print('their pvalue = {}, my pvalue = {}, enrichment = {}'.format(p_value, pvalue_2, enrichment))
                # score = np.log2(enrichment)
                score = enrichment

            sorted_list[n, 0] = self.ontology[go_id].name
            sorted_list[n, 1] = go_id
            sorted_list[n, 2] = score
            sorted_list[n, 3] = str(genes)
            self.global_go[go_id] = self.ontology[go_id].name

            n += 1
            # This is to return pandas datagrame
            tmp_entry.append(self.ontology[go_id].name)
            tmp_entry.append(go_id)
            tmp_entry.append(score)
            tmp_entry.append(list(np.sort(genes)))
            tmp_entry.append(len(genes))
            tmp_entry.append(ref)
            sorted_list_2.append(tmp_entry)
            # End of pandas upgrade space
            if self.verbose:
                print(go_id, float(len(genes) / float(ref)) * 100, self.ontology[go_id].name, p_value, len(genes), ref,
                      score)

        if n == 0:
            print("No significant p-values")
        cols = ['GO_name', 'GO_id', 'score_{0}'.format(num), 'genes_{0}'.format(num),'num_{0}'.format(num),'ref']
        data = pd.DataFrame(sorted_list_2, columns=cols)
        return sorted_list[:n, :], data

    def create_data_set(self, list_of_exp, aspect):
        """

        :param list_of_exp:
        :param aspect:
        :return:
        """

        results = []
        self.num_data_sets = len(list_of_exp)
        for n, i in enumerate(list_of_exp):
            tmp_array, tmp_array_2 = self.create_array_per_time(i, aspect, n)
            results.append(tmp_array)
            if n == 0:
                data_set_2 = tmp_array_2
            else:
                data_set_2 = pd.merge(data_set_2, tmp_array_2, on=['GO_id', 'GO_name'], how='outer')
        self.data_2 = data_set_2.fillna(0)

        all_terms = []
        for each in results:
            all_terms.extend(list(each[:, 1]))
        all_terms = set(all_terms)
        num_of_lists = len(list_of_exp)
        data_set = np.zeros((len(all_terms), num_of_lists + 1), dtype='S50')
        for n, i in enumerate(all_terms):
            data_set[n, 0] = i
            for num in range(1, num_of_lists + 1):
                data_set[n, num] = return_go_number(i, results[num - 1])
        if np.shape(data_set)[0] == 0:
            print('Warning! No significanance')

        return data_set

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

        if stop is None:
            matrix = data[start:, :]
        else:
            matrix = data[start:stop, :]

        size_of_data = np.shape(data)[1]
        length_matrix = len(matrix)
        fig = plt.figure(figsize=(14, 20))
        ax1 = fig.add_subplot(111)

        im = ax1.imshow(matrix, aspect=.25, interpolation='nearest',
                        extent=(0, size_of_data, 0, length_matrix + 1),
                        origin='lower')

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
        # divider = make_axes_locatable(ax1)
        # cax = divider.append_axes("top", size="5%", pad=0.05)
        plt.colorbar(im, fraction=0.046, pad=0.04)
        # plt.colorbar(im)
        # plt.axes().set_aspect('equal', 'datalim')
        # plt.tight_layout()
        plt.savefig(os.path.join(self.out_dir, 'Figures', '%s.pdf' % savename),
                    dpi=150, bbox_inches='tight')
        if self.save_png:
            plt.savefig(
                os.path.join(self.out_dir, 'Figures', '%s.png' % savename),
                dpi=150, bbox_inches='tight')
        plt.close()

    def analysis_data(self, data, aspect='F', savename='tmp', labels=None, analyze=True, search_term=None):
        """

        :param analyze:
        :param data:
        :param aspect:
        :param savename:
        :param labels:
        :return:
        """
        data = self.create_data_set(data, aspect)
        if data is None:
            print('No data! Error! Returning nothing')
            return
        self.data_2.to_csv('{0}.csv'.format(savename))

        if search_term is not None:
            self.data = self.find_terms(data, search_term)
        else:
            self.data = data

        self.names, self.array = self.sort_by_hierarchy(data)

        self.go_terms = data[:, 0]

        if not analyze:
            return
        savename += '_%s' % aspect
        self.savename = savename

        self.plot_heatmap(self.array, self.names, 0, 100, '%s_top' % savename, labels)
        self.plot_heatmap(self.array, self.names, -100, None, '%s_bottom' % savename, labels)
        tmp_array = self.array[:, :].copy()
        fig = plt.figure()
        axdendro = fig.add_axes([0.09, 0.1, 0.2, 0.8])
        # distance_matrix = dist.pdist(tmp_array, 'euclidean')
        #linkage = sch.linkage(distance_matrix, method='single')
        linkage = sch.linkage(tmp_array, method='centroid')
        dendrogram = sch.dendrogram(linkage, orientation='right')
        axdendro.set_xticks([])
        axdendro.set_yticks([])
        axmatrix = fig.add_axes([0.3, 0.1, 0.6, 0.8])
        index = dendrogram['leaves']
        tmp_array = tmp_array[index, :]
        names_sorted = self.names[index]
        im = axmatrix.matshow(tmp_array, aspect='auto', origin='lower')
        axmatrix.set_xticks([])
        axmatrix.set_yticks([])
        axcolor = fig.add_axes([0.91, 0.1, 0.02, 0.8])
        plt.colorbar(im, cax=axcolor)
        if self.save_png:
            fig.savefig(
                    os.path.join(self.out_dir, 'Figures',
                                 '%s_clustered.png' % savename))
        fig.savefig(os.path.join(self.out_dir, 'Figures',
                                 '%s_clustered.pdf' % savename))
        plt.close()

        self.plot_heatmap(tmp_array, names_sorted, 0, 100,
                          '%s_top_clustered' % savename, labels)
        self.plot_heatmap(tmp_array, names_sorted, -100, None,
                          '%s_bottom_clustered' % savename, labels)

        figures = ['%s_top' % savename,
                   '%s_bottom' % savename,
                   '%s_clustered' % savename,
                   '%s_top_clustered' % savename,
                   '%s_bottom_clustered' % savename, ]
        html_pages = []

        for i in figures:
            html_pages.append(
                    '<a href="Figures/{0}.pdf">{0}</a>'.format(i))

        self.html_pdfs = pd.DataFrame(html_pages, columns=['Clustered output'])
        self.print_ranked_over_time(savename=savename, labels=labels)

    def print_ranked_over_time(self, savename=None, labels=None, number=10,
                               create_plots=True):
        """ Prints information about top hits
        """
        html_pages = []
        self.html_pdfs2 = []
        self.top_hits = []
        n = np.shape(self.data)[1] - 1
        score_names = []

        if labels is None:
            labels = list(range(n))
        for i in range(n):
            score_names.append('score_{0}'.format(i))

        for i in range(0, np.shape(self.data)[1] - 1):
            terms = self.retrieve_top_ranked(i, number)
            tmp = []
            names = []
            for t in terms:
                names.append(t)
                tmp.append(np.array(self.data_2[self.data_2['GO_id'] == t][score_names])[0])

            tmp = np.array(tmp)[::-1]
            names = np.array(names)[::-1]

            if create_plots:
                self.plot_heatmap(tmp, names, -1 * number, None
                                  , 'top_hits_entry_%i_%s' % (i, savename), labels)
                html_pages.append(
                        '<a href="Figures/top_hits_entry_{0}_{1}.pdf">{2}</a>'.format(
                                i, savename, labels[i]))
            terms_dict = {}
            for j in range(number):
                if j > len(tmp) - 1:
                    continue
                terms_dict[names[j]] = tmp[j]
                if tmp[j, i] == 0.0:
                    continue
                else:
                    print("Top hits = {0}, {1}, {2}, {3}".format(
                            i, names[j], self.global_go[names[j]], tmp[j, i]))

            self.top_hits.append(terms_dict)

        if create_plots:
            self.html_pdfs2 = pd.DataFrame(html_pages, columns=['Top hits per time'])

    def retrieve_top_ranked(self, index, number=20):
        """ Prints information about top hits
        """
        names, tmp = sort_data_by_index(self.data, index)

        if len(names) < number:
            number = len(names)

        scores = {}
        for i in range(len(names)):
            scores[names[i]] = tmp[i, index]
        points = []
        counter = 0
        go_terms = list(names[-1 * number:])
        go_terms.reverse()
        #return go_terms

        term_to_add = -1 * number
        while counter < number:
            parents = dict([(term, self.get_parents(term, go_terms)) for term in go_terms])
            top_level_terms = [id for id in parents if not parents[id]]

            terms_to_remove = []
            index_to_add = []
            for term in top_level_terms:
                child = self.get_children(term, go_terms, parents)
                if len(child) == 0:
                    counter += 1
                else:
                    term_to_add -= 1
                    if np.abs(term_to_add) > len(names):
                        continue
                    terms_to_remove.append(term)

                    index_to_add.append(names[term_to_add])
            for t in terms_to_remove:
                go_terms.remove(t)
            for t in index_to_add:
                go_terms.append(t)

        for j in range(number):
            points.append(go_terms[j])
        return points

    def find_and_plot_subterms(self, term, savename, x=[1, 6, 24, 48]):
        """

        :param term: term to plot
        :param savename: name to save the image as
        :param x : labels for x axis
        """
        print(term)
        terms = self.ontology.extract_sub_graph(term)

        self.plot_specific_go(terms, savename, x)

    def plot_specific_go(self, term, savename, x):
        """ Plots a scatter plot of selected terms over time

        :param term:
        :param savename:
        :return:
        """
        if x is None:
            x = [1, 6, 24, 48]
        found_terms = []
        for i in term:
            if i in self.names:
                found_terms.append(i)

        num_colors = len(found_terms)
        if num_colors == 0:
            print("Did not find any of the terms")
            return
        fig = plt.figure()
        ax = plt.subplot(111)
        cm = plt.get_cmap('jet')

        ax.set_color_cycle([cm(1. * i / num_colors) for i in range(num_colors)])

        for i in found_terms:
            y_index = np.where(self.names == i)

            y = self.array[y_index][0]
            label = "\n".join(wrap(self.global_go[i], 40))
            ax.plot(x, y, 'o-', label=label)
        y_label = None
        if self.metric == 'enrichment':
            y_label = 'log2(enrichment)'
        if self.metric == 'pvalue':
            y_label = '-log10(p value)'
        if self.metric == 'fraction':
            y_label = 'fraction(%%)'
        plt.ylabel(y_label, fontsize=16)
        plt.xlabel('Time(hr)', fontsize=16)
        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles, labels, loc='best', bbox_to_anchor=(1.01, 1.0))
        fig.savefig('%s.pdf' % savename, bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.close()

    def export_to_html(self, labels, html_name='tmp'):
        """ exports to html to be viewed by browser

        :param labels:
        :param html_name:
        :param x:
        :return:
        """

        real_names = []
        html_array = self.array.copy()

        # replace GO numbers with GO name
        # create plot of genes over time
        for n, i in enumerate(self.names):

            gene_set = set()
            for k in range(self.num_data_sets):
                genes = list(self.data_2[self.data_2['GO_id'] == i]['genes_{0}'.format(k)])

                for g in genes:
                    if type(g) == list:
                        for e in g:
                            gene_set.add(e)
            save_name = '{0}/Figures/go_{1}_{2}.pdf'.format(self.out_dir, i,
                                                            self.savename).replace(':', '')
            title = "{0} : {1}".format(str(i), self.global_go[i])

            if len(gene_set) > 50:
                real_names.append('<a>{0}</a>'.format(self.global_go[i]))
            else:
                if save_name not in self.created_go_pds:
                    self.created_go_pds.add(save_name)
                    self.exp_data.plot_list_of_genes(list(gene_set), save_name,
                                                     title=title)
                real_names.append(
                    '<a href="Figures/go_{0}_{1}.pdf">{2}</a>'.format(i,
                                                                      self.savename,
                                                                      self.global_go[
                                                                          i]).replace(
                        ':', ''))
        # turn it into a pandas dataframe
        d = pd.DataFrame(data=html_array, index=real_names, columns=labels)
        # duplicate the sort table class so we can sort by column
        with open(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                               'html_additions', 'sorttable.js'), 'r') as f:
            x = f.read()
            with open('{}/Figures/sorttable_tmp.js'.format(self.out_dir),
                      'w') as tmp_file:
                tmp_file.write(x)

        header = '<script src="sorttable_tmp.js"></script>\n<html>\n\t<body>\n\n<div style="float: left">\n'
        footer = "\n\t</body>\n</html>"
        # write out the html file
        with open(os.path.join(self.out_dir, '%s.html' % html_name), 'w') as f:
            f.write(header)
            f.write(self.html_pdfs.to_html(classes='sortable', escape=False, justify='left'))
            f.write('\n</div>\n')
            f.write(self.html_pdfs2.to_html(classes='sortable', escape=False, justify='left'))
            f.write(d.to_html(classes='sortable', float_format='{0:.4}'.format, escape=False, justify='left'))
            f.write(footer)

    def get_parents(self, term, data):
        parents = self.ontology.extract_super_graph([term])
        parents = [id for id in parents if id in data and id != term]
        lst = [set(self.ontology.extract_super_graph([id])) - set([id]) for id in parents]
        c = set.union(*lst) if lst else set()
        parents = [t for t in parents if t not in c]
        return parents

    def get_children(self, term, data, parents):
        return [id for id in data if term in parents[id]]

    def sort_by_hierarchy(self, data):
        """ Sorts according to GO hierarchy

        :param data:
        :return:
        """
        names = data[:, 0].copy()
        array = data[:, 1:].astype(np.float32).copy()
        parents = dict([(term, self.get_parents(term, names)) for term in names])
        top_level_terms = [id for id in parents if not parents[id]]
        term_list = []
        visited = []

        def collect(go_term, parent):
            """  collects all children terms to a hierarchy
            modified from bio.orange.go
            :param go_term:
            :param parent:
            :return:
            """
            term_list.append((go_term, self.global_go[go_term], parent))
            parent = len(term_list) - 1
            for c in self.get_children(go_term, names, parents):
                if c in visited:
                    continue
                else:
                    visited.append(c)
                    collect(c, parent)

        for topTerm in top_level_terms:
            collect(topTerm, None)
        term_list = np.array(term_list)
        list_sorted = list(term_list[:, 0])
        tmp_array = []
        for i in list_sorted:
            ind = np.where(names == i)
            tmp_array.append(ind[0][0])
        array = array[tmp_array]
        names = names[tmp_array]
        return names, array

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


def sort_data_by_index(data, index=0):
    """

    :param index:
    :param data:
    :return:
    """
    names = data[:, 0].copy()
    array = data[:, 1:].astype(np.float32).copy()
    step_size = (array[:, index]).argsort()
    names = names[step_size]
    array = array[step_size]
    return names, array
