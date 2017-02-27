"""
GO analysis function using orange bioinformatics
"""
import os
import time

import numpy as np
import pandas as pd
import pathos.multiprocessing as mp
from orangecontrib.bio import go

from magine.html_templates.html_tools import write_single_table
from magine.plotting.heatmaps import plot_heatmap, plot_heatmap_cluster

from magine.ontology.go_from_goatools import go as goa_tools

pd.set_option('display.max_colwidth', -1)

evidence_codes = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC', 'IEA']


class GoAnalysis:
    """
    Go analysis class.
    Uses the orangecontrib.bio package
    """

    def __init__(self, species='hsa', output_directory='tmp',
                 reference=None, verbose=False, experimental_data=None,
                 save_name='tmp', slim_name='goslim_pir'):

        self.exp_data = experimental_data
        self.data = None
        self.savename = save_name
        self.out_dir = output_directory
        self.ontology = go.Ontology()
        self.annotations = go.Annotations(species, ontology=self.ontology)

        self.reference = reference
        self.verbose = verbose

        # options can be goslim_pir goslim_generic goslim_chembl
        self.slim_name = slim_name

        self.gene_annotations = set(self.annotations.gene_names)
        self.number_of_total_reference_genes = len(self.annotations.gene_names)
        self.slim_set = self.ontology.named_slims_subset(self.slim_name)

        self.global_go = {}

        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)
        if not os.path.exists(os.path.join(self.out_dir, 'Figures')):
            os.mkdir(os.path.join(self.out_dir, 'Figures'))
        self.top_hits = []
        self.num_data_sets = 0
        self.created_go_pds = set()

    def enrichment_analysis_of_single_sample(self, gene_list):
        """
        Performs enrichment analysis of list of genes

        Parameters
        ----------
        gene_list : array like
            List of genes to perform analysis
        """
        # checks first to see if all genes are annotated in GO
        genes_present = set()
        for i in gene_list:
            if type(i) != str:
                continue
            elif i in self.gene_annotations:
                genes_present.add(i)
            else:
                split_name = i.split(',')
                if len(split_name) > 1:
                    for n in split_name:
                        if n in self.gene_annotations:
                            genes_present.add(n)
        number_of_genes = len(genes_present)
        print("Number of genes given = {0}."
              " Number of genes in GO = {1}".format(len(gene_list),
                                                    number_of_genes))

        res = self.annotations.get_enriched_terms(genes_present,
                                                  aspect=['F', 'P', 'C'],
                                                  # evidence_codes=evidence_codes,
                                                  reference=self.reference)

        n = 0
        column_names = dict(GO_id=None, enrichment_score=None, pvalue=None,
                            genes=None, n_genes=None, )

        all_go_rows = []
        for go_id, (genes, p_value, ref) in res.items():
            n += 1
            go_row = column_names.copy()

            expected_value = number_of_genes * float(
                    ref) / self.number_of_total_reference_genes
            enrichment = float(len(genes)) / expected_value
            score = enrichment

            self.global_go[go_id] = self.ontology[go_id].name

            go_row['GO_id'] = go_id
            go_row['enrichment_score'] = score
            go_row['pvalue'] = p_value
            go_row['genes'] = list(np.sort(genes))
            go_row['n_genes'] = len(genes)

            all_go_rows.append(go_row)
            if self.verbose:
                print(go_id, float(len(genes) / float(ref)) * 100,
                      self.ontology[go_id].name, p_value, len(genes), ref,
                      score)

        if n == 0:
            print("No significant p-values")

        cols2 = ['GO_id', 'enrichment_score', 'pvalue', 'genes', 'n_genes']
        return pd.DataFrame(all_go_rows, columns=cols2)

    def create_enrichment_array(self, list_of_exp, labels=None):
        """
        Performs enrichment analysis of list of list of genes

        Parameters
        ----------
        list_of_exp: list_of_list
            A list of samples, where each sample is a list that contains
            species names
        labels: list_like
            list of labels for each sample

        Returns
        -------

        """
        self.num_data_sets = len(list_of_exp)

        assert self.num_data_sets != 0, "Must provide at least one data list"

        if isinstance(labels, list):
            assert (len(labels) == self.num_data_sets,
                    "If you provide a list, must be the same length"
                    " as list_of_exp", AttributeError)

        if labels is None:
            labels = range(0, self.num_data_sets)

        all_data = []
        for n, i in enumerate(list_of_exp):
            tmp_array_3 = self.enrichment_analysis_of_single_sample(i)
            tmp_array_3['sample_index'] = labels[n]
            all_data.append(tmp_array_3)

        results = pd.concat(all_data)

        # add aspect, n_ref, depth, and if in a slim to dataset
        for go_id in results['GO_id'].unique():
            index = results['GO_id'] == go_id
            go_name = self.ontology[go_id].name
            aspect_name = goa_tools[go_id].namespace
            n_reference_genes = len(self.annotations.get_all_genes(go_id))
            depth = self.ontology.term_depth(go_id)
            results.loc[index, 'aspect'] = aspect_name
            results.loc[index, 'ref'] = n_reference_genes
            results.loc[index, 'depth'] = depth
            results.loc[index, 'GO_name'] = go_name
            if go_id in self.slim_set:
                results.loc[index, 'slim'] = True
            else:
                results.loc[index, 'slim'] = False
        if results is None:
            print('No data! Error! Returning nothing')
            return
        self.data = results
        return results

    def create_gene_plots_per_go(self, data=None):
        """ Creates a figure for each GO term in data

        Data should be a result of running create_enrichment_array.
        This function creates a plot of all proteins per term if a term is
        significant and the number of the reference set is larger than 5 and
        the total number of species measured is less than 50.



        Parameters
        ----------
        data : pandas.DataFrame, optional
            previously ran enrichment analysis

        Returns
        -------
        out_array : dict
            dict where keys are pointers to figure locations
        """

        if data is None and self.data is None:
            print("Array is empty. This could be due to no significant species"
                  "provided in list or because analysis has not been run yet.")
            return

        elif data is None:
            data = self.data.copy()
        elif isinstance(data, str):
            data = pd.read_csv(data)
        # get list of all terms
        list_of_go_terms = data['GO_id'].unique()

        # filter data by significance and number of references
        data = data[data['ref'] >= 5]
        data = data[data['pvalue'] < 0.05]

        # create filtered list of species
        list_of_sig_go = data['GO_id'].unique()

        # here we are going to iterate through all sig GO terms and create
        # a list of plots to create. For the HTML side, we need to point to
        # a location
        figure_locations = {}
        plots_to_create = []
        # create plot of genes over time
        for n, i in enumerate(list_of_go_terms):

            # want to plot all species over time
            index = data['GO_id'] == i
            name = data[index]['GO_name'].unique()
            if len(name) > 0:
                name = name[0]
            # want to only plot significant species
            if i not in list_of_sig_go:
                figure_locations[i] = '<a>{0}</a>'.format(name)
                continue

            gene_set = set()
            genes = data[index]['genes']
            for g in genes:
                for j in g:
                    gene_set.add(j)

            # too many genes isn't helpful on plots, so skip them
            if len(gene_set) > 50:
                figure_locations[i] = '<a>{0}</a>'.format(name)
                continue

            save_name = '{0}/Figures/go_{1}_{2}'.format(self.out_dir, i,
                                                        self.savename)
            save_name = save_name.replace(':', '')
            title = "{0} : {1}".format(str(i), name)

            plots_to_create.append(
                    (list(gene_set), save_name, '.', title, True, True))
            out_point = '<a href="Figures/go_{0}_{1}.pdf">{2} ({0})</a>'
            out_point = out_point.format(i, self.savename, name).replace(':',
                                                                         '')
            figure_locations[i] = out_point

        print("Starting to create plots for each GO term")
        # just keeping this code just in case using pathos is a bad idea
        # ultimately, using matplotlib is slow.
        run_seq = False
        run_par = True
        if run_seq:
            st1 = time.time()
            for i in plots_to_create:
                self.exp_data.plot_list_of_genes(i)
            end1 = time.time()
            print("sequential time = {}".format(end1 - st1))

        if run_par:
            st2 = time.time()
            pool = mp.Pool(8)
            pool.map(self.exp_data.plot_list_of_genes, plots_to_create)
            pool.close()
            pool.join()
            end2 = time.time()
            print("parallel time = {}".format(end2 - st2))
        print("Done creating plots for each GO term")

        return figure_locations

    def write_table_to_html(self, data=None, save_name='index'):
        # create plots of everything

        if data is None and self.data is None:
            print("Array is empty. This could be due to no significant species"
                  "provided in list or because analysis has not been run yet.")
            return

        if data is None:
            data = self.data.copy()

        fig_dict = self.create_gene_plots_per_go()

        for i in fig_dict:
            data.loc[data['GO_id'] == i, 'GO_name'] = fig_dict[i]

        tmp = pd.pivot_table(data,
                             index=['GO_id', 'GO_name', 'depth', 'ref', 'slim',
                                    'aspect'],
                             columns='sample_index')

        html_out = os.path.join(self.out_dir, save_name)
        write_single_table(tmp, html_out, 'MAGINE GO analysis')

    # TODO remove and replace with updated heatplots
    def create_heatmaps(self, labels, savename):

        # depracated
        data = self.data
        names, array = self.sort_by_hierarchy(data)
        # creats plots of top and bottom of hierarchy sorted arrays
        plot_heatmap(array, names, labels, start=0, stop=100,
                     savename='%s_top' % savename, )
        plot_heatmap(array, names, labels, start=-100, stop=None,
                     savename='%s_bottom' % savename)
        tmp_array = array[:, :].copy()

        # Centroid clustering of data
        # Creates dendrogram plots
        clustered_array, clustered_names = plot_heatmap_cluster(tmp_array)

        plot_heatmap(clustered_array, clustered_names, labels, start=0,
                     stop=100, savename='%s_top_clustered' % savename)
        plot_heatmap(clustered_array, clustered_names, labels, start=-100,
                     stop=None,
                     savename='%s_bottom_clustered' % savename)

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
        """
        Order array by each sample and print information
        Parameters
        ----------
        savename : str
            name to save file as
        labels : list like
            list of samples
        number : int
            number of top hits to print
        create_plots : bool
            create a pdf of top hits

        Returns
        -------

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
                tmp.append(np.array(
                    self.data_2[self.data_2['GO_id'] == t][score_names])[0])

            tmp = np.array(tmp)[::-1]
            names = np.array(names)[::-1]

            if create_plots:
                plot_heatmap(tmp, names, labels, start=-1 * number, stop=None,
                             savename='top_hits_entry_%i_%s' % (i, savename),)
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
            self.html_pdfs2 = pd.DataFrame(html_pages,
                                           columns=['Top hits per time'])

    def retrieve_top_ranked(self, index, number=20):
        """
        Sorts array by column while excluding similar related GO terms
        Parameters
        ----------
        index : int
            column number to sort GO array by
        number : int
            number of terms wanted in top hits

        Returns
        -------

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
        # return go_terms

        term_to_add = -1 * number
        while counter < number:
            parents = dict(
                [(term, self.get_parents(term, go_terms)) for term in
                 go_terms])
            top_level_terms = [id for id in parents if not parents[id]]

            terms_to_remove = []
            index_to_add = []
            for term in top_level_terms:
                child = get_children(term, go_terms, parents)
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

    def get_parents(self, term, data):
        """
        Get parents of GO term from provided lsit
        Parameters
        ----------
        term : str
            GO ontology term
        data : list
            list of GO terms in which ones wants to find parents of provided
            term

        Returns
        -------

        """
        parents = self.ontology.extract_super_graph([term])
        parents = [id for id in parents if id in data and id != term]
        lst = [set(self.ontology.extract_super_graph([id])) - set([id]) for id
               in parents]
        c = set.union(*lst) if lst else set()
        parents = [t for t in parents if t not in c]
        return parents

    def sort_by_hierarchy(self, data):
        """
        Returns numpy array of GO terms by GO hierarchy
        Parameters
        ----------
        data : array
            Array of GO terms where first row is GO id, other rows enrichment
            values

        Returns
        -------

        """
        names = data[:, 0].copy()
        array = data[:, 1:].astype(np.float32).copy()
        parents = dict(
            [(term, self.get_parents(term, names)) for term in names])
        top_level_terms = [id for id in parents if not parents[id]]
        term_list = []
        visited = []

        def _collect(go_term, parent):
            """
            Collects all children terms to a hierarchy
            Modified from bio.orange.go

            Parameters
            ----------
            go_term : str
                GO term
            parent : int or None
                Counter to provide GO hierarchy

            Returns
            -------

            """
            term_list.append((go_term, self.global_go[go_term], parent))
            parent = len(term_list) - 1
            for c in get_children(go_term, names, parents):
                if c in visited:
                    continue
                else:
                    visited.append(c)
                    _collect(c, parent)

        for topTerm in top_level_terms:
            _collect(topTerm, None)
        term_list = np.array(term_list)
        list_sorted = list(term_list[:, 0])
        tmp_array = []
        for i in list_sorted:
            ind = np.where(names == i)
            tmp_array.append(ind[0][0])
        array = array[tmp_array]
        names = names[tmp_array]
        return names, array


def get_children(term, data, parents):
    """

    Parameters
    ----------
    term : str
        GO term
    data : array

    parents :
        parents of terms

    Returns
    -------
    out : list
        list of children of GO term provided

    """
    return [id for id in data if term in parents[id]]


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
