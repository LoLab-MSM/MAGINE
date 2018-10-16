"""
GO analysis function using orange bioinformatics
"""

import os

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection
from statsmodels.stats.proportion import binom_test

from magine.data.storage import network_data_dir
from magine.enrichment.deprecated.databases.gene_ontology import \
    download_and_process_go
from magine.html_templates.html_tools import write_filter_table
from magine.plotting.species_plotting import plot_genes_by_ont

try:
    import cPickle as pickle
except:  # python3 doesnt have cPickle
    import pickle

pd.set_option('display.max_colwidth', -1)

evidence_codes = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC', 'IEA']


def pivot_table_for_export(data, save_name=None):
    """
    creates a pivot table of combined data that is in MAGINE format

    Parameters
    ----------
    data : pd.DataFrame
        output from magine analysis
    save_name : str
        save name for excel merged format

    Returns
    -------

    """
    index = ['GO_id', 'GO_name', 'depth', 'ref', 'aspect']
    tmp = pd.pivot_table(data, index=index, columns='sample_index',
                         aggfunc='first'
                         )
    tmp = tmp[['pvalue', 'enrichment_score', 'genes', 'n_genes', ]]
    if save_name:
        tmp.to_excel('{}.xlsx'.format(save_name),
                     merge_cells=True)
    return tmp


def write_table_to_html_with_figures(data, exp_data, save_name='index',
                                     out_dir=None, run_parallel=True):
    # create plots of everything
    if isinstance(data, str):
        data = pd.read_csv(data)

    fig_dict, to_remove = plot_genes_by_ont(
        data, save_name, out_dir, exp_data, run_parallel=run_parallel
    )

    for i in fig_dict:
        data.loc[data['GO_id'] == i, 'GO_name'] = fig_dict[i]

    data = data[~data['GO_id'].isin(to_remove)]

    tmp = pivot_table_for_export(data)

    html_out = save_name
    if out_dir is not None:
        html_out = os.path.join(out_dir, html_out)
    print("Saving to : {}".format(html_out))

    html_out = save_name + '_filter'
    if out_dir is not None:
        html_out = os.path.join(out_dir, html_out)

    write_filter_table(tmp, html_out)


class GoAnalysis(object):
    """
    Go analysis class.
    Uses the orangecontrib.bio package
    """

    def __init__(self, species='hsa', output_directory='tmp',
                 reference=None, verbose=False, experimental_data=None,
                 save_name='tmp', slim_name='goslim_pir'):

        self.exp_data = experimental_data
        self.data = None
        self.magine_go = MagineGO(species)
        self.savename = save_name
        self.out_dir = output_directory
        self.reference = reference
        self.verbose = verbose

        # options can be goslim_pir goslim_generic goslim_chembl
        self.slim_name = slim_name

        if reference is not None:
            self.gene_annotations = \
                set(reference).intersection(self.magine_go.gene_to_go.keys())
        else:
            self.gene_annotations = set(self.magine_go.gene_to_go.keys())
        self.number_of_total_reference_genes = len(self.gene_annotations)

        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)
        if not os.path.exists(os.path.join(self.out_dir, 'Figures')):
            os.mkdir(os.path.join(self.out_dir, 'Figures'))
        self.top_hits = []
        self.num_data_sets = 0
        self.created_go_pds = set()

    def _calculate_enrichment_single_sample(self, gene_list):
        """
        Performs enrichment analysis of list of genes

        Parameters
        ----------
        gene_list : array like
            List of genes to perform analysis
        """
        # checks first to see if all genes are annotated in GO
        genes_present = set()
        genes_missing = set()
        for i in gene_list:
            if type(i) != str:
                continue
            elif i in self.gene_annotations:
                genes_present.add(i)
            else:
                genes_missing.add(i)
                split_name = i.split(',')
                if len(split_name) > 1:
                    for n in split_name:
                        if n in self.gene_annotations:
                            genes_present.add(n)
        n_genes = len(genes_present)
        if len(genes_missing) > 0:
            print("Genes not in GO = {}".format(genes_missing))
        print("Number of genes given = {0}."
              " Number of genes in GO = {1}".format(len(gene_list), n_genes))

        # res = self.__calculate_enrichment(genes_present,
        # aspect=['F', 'P', 'C'])
        res = self.magine_go.calculate_enrichment(genes_present,
                                                  reference=self.reference)
        n = 0
        column_names = dict(GO_id=None, enrichment_score=None, pvalue=None,
                            genes=None, n_genes=None, aspect=None,
                            ref=None, depth=None, GO_name=None, slim=None)

        all_go_rows = []
        for go_id, (genes, p_value, ref) in res.items():
            n += 1
            go_row = column_names.copy()

            expected_value = n_genes * float(
                ref) / self.number_of_total_reference_genes
            enrichment = float(len(genes)) / expected_value

            go_row['GO_id'] = go_id
            go_row['enrichment_score'] = enrichment
            go_row['pvalue'] = p_value
            go_row['genes'] = list(np.sort(genes))
            go_row['n_genes'] = len(genes)
            go_row['aspect'] = self.magine_go.go_aspect[go_id]
            go_row['ref'] = len(self.magine_go.go_to_gene[go_id])
            go_row['depth'] = self.magine_go.go_depth[go_id]
            go_row['GO_name'] = self.magine_go.go_to_name[go_id]

            all_go_rows.append(go_row)
            if self.verbose:  # pragma: no cover
                print(go_id, float(len(genes) / float(ref)) * 100,
                      self.magine_go.go_to_name[go_id],
                      p_value, len(genes), ref,
                      enrichment)

        if n == 0:
            print("No significant p-values")

        return pd.DataFrame(all_go_rows, columns=column_names.keys())

    def calculate_enrichment(self, list_of_exp, labels=None):
        """
        Performs enrichment analysis of list of list of genes

        Parameters
        ----------
        list_of_exp: list_of_list or list
            Can be a single list of genes or a list of lists, where each sample
             is a list that contains species names.
        labels: list_like
            list of labels for each sample

        Returns
        -------
        pandas.DataFrame

        """
        self.num_data_sets = len(list_of_exp)

        if not all(isinstance(el, list) for el in list_of_exp):
            print("Running single GO enrichment since single list provided.")
            return self._calculate_enrichment_single_sample(list_of_exp)

        if self.num_data_sets == 0:
            print("Must provide at least one data list")
            print("Returning nothing")
            return

        if isinstance(labels, list):
            assert len(labels) == self.num_data_sets, \
                "If you provide a list, must be the same length" \
                " as list_of_exp"

        if labels is None:
            labels = range(0, self.num_data_sets)

        all_data = []
        for n, i in enumerate(list_of_exp):
            tmp = self._calculate_enrichment_single_sample(i)
            tmp['sample_index'] = labels[n]
            all_data.append(tmp)

        results = pd.concat(all_data)

        if results is None:  # pragma: no cover
            print('No data! Error! Returning nothing')
            return
        self.data = results
        return results

    def write_table_to_html(self, save_name='index', out_dir=None,
                            run_parallel=False):
        """
        Creates a table of all the plots of all genes for each GO term.
        
        Uses last calculated enrichment array.
        
        Parameters
        ----------
        save_name : str
            name of html output file
        out_dir : str, optional
            output path for all plots
        run_parallel : bool
            Create plots in parallel

        Returns
        -------

        """
        if out_dir is None:
            out_dir = self.out_dir
        print("BaseData", self.data)
        print("exp_data", self.exp_data)
        write_table_to_html_with_figures(
                self.data, self.exp_data, save_name, out_dir=out_dir,
                run_parallel=run_parallel
        )


class MagineGO(object):
    def __init__(self, species='hsa'):
        dirname = network_data_dir
        gene_to_go_name = os.path.join(dirname,
                                       '{}_gene_to_go.p'.format(species))
        go_to_gene_name = os.path.join(dirname,
                                       '{}_goids_to_genes.p'.format(
                                           species))
        go_to_go_name = os.path.join(dirname,
                                     '{}_goids_to_goname.p'.format(
                                         species))
        go_depth = os.path.join(dirname, '{}_godepth.p'.format(species))
        go_aspect = os.path.join(dirname, '{}_go_aspect.p'.format(species))
        for i in [gene_to_go_name, go_to_gene_name, go_to_go_name,
                  go_depth,
                  go_aspect]:
            if not os.path.exists(i):
                download_and_process_go(species=species)

        self.gene_to_go = pickle.load(open(gene_to_go_name, 'rb'))
        self.go_to_gene = pickle.load(open(go_to_gene_name, 'rb'))
        self.go_to_name = pickle.load(open(go_to_go_name, 'rb'))
        self.go_depth = pickle.load(open(go_depth, 'rb'))
        self.go_aspect = pickle.load(open(go_aspect, 'rb'))

    def calculate_enrichment(self, genes, reference=None,
                             evidence_codes=None,
                             aspect=None, use_fdr=True):
        """

        Parameters
        ----------
        genes : list
            list of genes
        reference : list
            reference list of species to calculate enrichment
        evidence_codes : list
            GO evidence codes
        use_fdr : bool
            Correct for multiple hypothesis testing

        Returns
        -------

        """

        # TODO check for alias for genes
        genes = set(genes)
        # TODO add aspects
        term_reference = self.go_to_gene.keys()
        aspect_dict = {
            'P': 'biological_process',
            'C': 'cellular_component',
            'F': 'molecular_function'
        }
        if aspect is None:
            term_reference = self.go_to_gene
            gene_reference = self.gene_to_go
        else:
            term_reference = dict()
            gene_reference = dict()

        if aspect is not None:
            for i in aspect:
                if i not in ['P', 'C', 'F']:
                    print("Error: Aspects are only 'P', 'C', and 'F' \n")
                    quit()
            for i in ['P', 'C', 'F']:
                if i in aspect:
                    term_reference = None

        # TODO add reference
        if reference:
            # TODO check for reference alias
            reference = set(reference)
            reference.intersection_update(set(self.gene_to_go.keys()))
        else:
            reference = set(self.gene_to_go.keys())

        # TODO add evidence_codes

        terms = set()
        for i in genes:
            if i in self.gene_to_go:
                for t in self.gene_to_go[i]:
                    terms.add(t)

        n_genes = len(genes)
        n_ref = float(len(reference))
        res = {}
        for term in terms:

            all_annotated_genes = set(self.go_to_gene[term])
            mapped_genes = genes.intersection(all_annotated_genes)
            n_mapped_genes = len(mapped_genes)

            if n_ref > len(all_annotated_genes):
                mapped_reference_genes = \
                    reference.intersection(all_annotated_genes)
            else:
                mapped_reference_genes = \
                    all_annotated_genes.intersection(reference)

            n_mapped_ref = len(mapped_reference_genes)

            prob = float(n_mapped_ref) / n_ref

            p_value = binom_test(n_mapped_genes, n_genes, prob, 'larger')

            res[term] = ([i for i in mapped_genes], p_value, n_mapped_ref)
        if use_fdr:
            res = sorted(res.items(), key=lambda x: x[1][1])
            fdr = fdrcorrection([p for _, (_, p, _) in res],
                                is_sorted=True)
            values = fdr[1]
            res = dict([(index, (genes, p, ref))
                        for (index, (genes, _, ref)), p in
                        zip(res, values)])
        return res

# All code below is old and will be removed. Keeping for just a bit longer.
'''  
# from magine.ontology.go_from_goatools import go as goa_tools
# from orangecontrib.bio import go
# from orangecontrib.bio.go import evidenceDict
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
from statsmodels.stats.proportion import binom_test
from collections import defaultdict
    def __calculate_enrichment(self, genes, reference=None, evi_codes=None,
                               aspect=None, use_fdr=True, ):
        rev_genes_dict = self.annotations.get_gene_names_translator(genes)
        genes = set(rev_genes_dict.keys())
        if reference:
            ref_genes_dict = self.annotations.get_gene_names_translator(
                    reference)
            reference = set(ref_genes_dict.keys())
        else:
            reference = self.annotations.gene_names

        if aspect is None:
            aspects_set = set(["P", "C", "F"])
        elif isinstance(aspect, str):
            aspects_set = set([aspect])
        else:
            aspects_set = aspect

        evi_codes = set(evi_codes or evidenceDict.keys())

        annotations_dict = defaultdict(set)
        for gene in genes:
            for ann in self.annotations.gene_annotations[gene]:
                if ann.Evidence_Code in evi_codes and ann.Aspect in aspects_set:
                    annotations_dict[ann.GO_ID].add(ann)
        ref_annotations = set([ann for gene in reference for ann in
                               self.annotations.gene_annotations[gene] if
                               ann.Evidence_Code in evi_codes and ann.Aspect in aspects_set])

        terms = annotations_dict.keys()
        filtered_terms = [term for term in terms if
                          term in self.annotations.ontology]

        if len(terms) != len(filtered_terms):
            term_diff = set(terms) - set(filtered_terms)
            warnings.warn(
                    "%s terms in the annotations were not found in the "
                    "ontology." % ",".join(map(repr, term_diff)),
                    UserWarning)

        terms = self.annotations.ontology.extract_super_graph(filtered_terms)
        n_genes = len(genes)
        n_ref = float(len(reference))

        res = {}
        for i, term in enumerate(terms):
            # all_annotated_genes = self.annotations.get_all_genes(term)
            # 46 % of time is all_annotations
            all_annotations = self.annotations.get_all_annotations(
                term).intersection(ref_annotations)
            # 39.1 % of time is all_annotated_genes
            all_annotated_genes = set(
                    [ann.geneName for ann in all_annotations])

            mapped_genes = genes.intersection(all_annotated_genes)
            n_mapped_genes = len(mapped_genes)
            if len(reference) > len(all_annotated_genes):
                mapped_reference_genes = reference.intersection(
                        all_annotated_genes)
            else:
                mapped_reference_genes = all_annotated_genes.intersection(
                        reference)
            n_mapped_ref = len(mapped_reference_genes)
            prob = float(n_mapped_ref) / n_ref

            p_value = binom_test(n_mapped_genes, n_genes, prob, 'larger')

            res[term] = ([rev_genes_dict[g] for g in mapped_genes],
                         p_value, n_mapped_ref)
        if use_fdr:
            res = sorted(res.items(), key=lambda x: x[1][1])
            fdr = fdrcorrection0([p for _, (_, p, _) in res],
                                 is_sorted=True)
            values = fdr[1]
            res = dict([(index, (genes, p, ref))
                        for (index, (genes, _, ref)), p in
                        zip(res, values)])

        return res

    # TODO remove and replace with updated heatplots
    def __create_heatmaps(self, labels, savename):
        from magine.plotting.heatmaps import plot_heatmap, plot_heatmap_cluster
        # depracated
        data = self.data
        names, array = self.__sort_by_hierarchy(data)
        # creats plots of top and bottom of hierarchy sorted arrays
        plot_heatmap(array, names, labels, self.magine_go.go_to_name, start=0,
                     stop=100, savename='%s_top' % savename, )
        plot_heatmap(array, names, labels, self.magine_go.go_to_name,
                     start=-100, stop=None, savename='%s_bottom' % savename)
        tmp_array = array[:, :].copy()

        # Centroid clustering of data
        # Creates dendrogram plots
        clustered_array, clustered_names = plot_heatmap_cluster(
                tmp_array, names, savename="{}_dendrogram".format(savename),
                out_dir=self.out_dir)

        plot_heatmap(clustered_array, clustered_names, labels,
                     self.magine_go.go_to_name, start=0, stop=100,
                     savename='%s_top_clustered' % savename)
        plot_heatmap(clustered_array, clustered_names, labels,
                     self.magine_go.go_to_name, start=-100, stop=None,
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
        self.__print_ranked_over_time(savename=savename, labels=labels)

    def __print_ranked_over_time(self, savename=None, labels=None, number=10,
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
        from magine.plotting.heatmaps import plot_heatmap
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
            terms = self.__retrieve_top_ranked(i, number)
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

    def __retrieve_top_ranked(self, index, number=20):
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

        names, tmp = _sort_data_by_index(self.data, index)

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
                    [(term, self.__get_parents(term, go_terms)) for term in
                     go_terms])
            top_level_terms = [id for id in parents if not parents[id]]

            terms_to_remove = []
            index_to_add = []
            for term in top_level_terms:
                child = _get_children(term, go_terms, parents)
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

    def __get_parents(self, term, data):
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

    def __sort_by_hierarchy(self, data):
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
                [(term, self.__get_parents(term, names)) for term in names])
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
            for c in _get_children(go_term, names, parents):
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


def _get_children(term, data, parents):
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


def _sort_data_by_index(data, index=0):
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
'''
