"""
GO analysis function using orange bioinformatics
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
from orangecontrib.bio import go

pd.set_option('display.max_colwidth', -1)

evidence_codes = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC', 'IEA']


# noinspection PyUnresolvedReferences
class GoAnalysis:
    """
    Go analysis class.
    Uses the orangecontrib.bio package
    """

    def __init__(self, species='hsa', slim=False, output_directory='tmp',
                 reference=None, metric='enrichment',
                 verbose=False, experimental_data=None, save_png=False,
                 save_name='tmp'):

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
            # options can be goslim_pir goslim_generic goslim_chembl
            self.ontology.set_slims_subset(slim)
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

    def enrichment_analysis_of_single_sample(self, gene_list, aspect='F',
                                             num=0, slims_only=False):
        """
        Performs enrichment analysis of list of genes

        Parameters
        ----------
        gene_list : array like
            List of genes to perform analysis
        aspect : str or list like
            Can be str of list containing F, P, C
                F : Molecular Function
                P : Biological Process
                C : Cellular Component

        num : int
            sample number used for pandas dataframe
        slims_only : bool
            Use slim terms
        """
        # checks first to see if all genes are annotated in GO
        genes_present = set()
        for i in gene_list:
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
        print("Number of genes given = {0}."
              " Number of genes in GO = {1}".format(len(gene_list),
                                                    number_of_genes))
        res = self.annotations.get_enriched_terms(genes_present,
                                                  slims_only=self.slim,
                                                  aspect=aspect,
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
            # if p_value > 0.05:
            #     continue
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
            tmp_entry.append(p_value)
            tmp_entry.append(list(np.sort(genes)))
            tmp_entry.append(len(genes))
            tmp_entry.append(self.ontology.term_depth(go_id))
            tmp_entry.append(ref)
            sorted_list_2.append(tmp_entry)
            # End of pandas upgrade space
            if self.verbose:
                print(go_id, float(len(genes) / float(ref)) * 100,
                      self.ontology[go_id].name, p_value, len(genes), ref,
                      score)

        if n == 0:
            print("No significant p-values")
        cols = ['GO_name', 'GO_id', 'score_{0}'.format(num),
                'pvalue_{0}'.format(num), 'genes_{0}'.format(num),
                'num_{0}'.format(num), 'depth', 'ref']
        cols2 = ['GO_name', 'GO_id', 'enrichment_score', 'pvalue', 'genes',
                 'n_genes', 'depth', 'ref']
        gene_list = pd.DataFrame(sorted_list_2, columns=cols)
        gene_list2 = pd.DataFrame(sorted_list_2, columns=cols2)
        return sorted_list[:n, :], gene_list, gene_list2

    def create_enrichment_array(self, list_of_exp, aspect, use_slims=False):
        """
        Performs enrichment analysis of list of list of genes

        Parameters
        ----------
        list_of_exp: array like
            List of genes

        aspect : str or list like
            Can be str of list containing F, P, C
                F : Molecular Function
                P : Biological Process
                C : Cellular Component

        use_slims : bool
            Use slim terms from GO

        Returns
        -------

        """

        results = []
        self.num_data_sets = len(list_of_exp)
        assert self.num_data_sets != 0, "Must provide at least one data list"
        all_data = []
        for n, i in enumerate(list_of_exp):
            tmp_array, tmp_array_2, tmp_array_3 = self.enrichment_analysis_of_single_sample(
                i, aspect, n, slims_only=use_slims)
            tmp_array_3['sample_index'] = n
            results.append(tmp_array)
            all_data.append(tmp_array_3)
            if n == 0:
                data_set_2 = tmp_array_2
            else:
                data_set_2 = pd.merge(data_set_2, tmp_array_2,
                                      on=['GO_id', 'GO_name'], how='outer')

        self.data_2 = data_set_2.fillna(0)
        self.data_3 = pd.concat(all_data).fillna(0)

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
            print('Warning! No significant species!!!')
            return None

        return data_set

    def plot_heatmap(self, enrichment_array, names_col, labels, start=0,
                     stop=None, savename='tmp'):
        """
        General function to create a heatmap of enrichment array data

        Parameters
        ----------
        enrichment_array : array like
            numpy array of enrichment values
        names_col : array like
            array of names
        labels : list
            list of column names
        start : int
            starting index for size of array to be shown
        stop : int or None
            stopping index for size of array to be shown
        savename : str
            name of figure to be saved as

        Returns
        -------

        """

        if stop is None:
            matrix = enrichment_array[start:, :]
        else:
            matrix = enrichment_array[start:stop, :]

        size_of_data = np.shape(enrichment_array)[1]
        length_matrix = len(matrix)
        fig = plt.figure(figsize=(14, 20))
        ax1 = fig.add_subplot(111)

        im = ax1.imshow(matrix, aspect=.25, interpolation='nearest',
                        extent=(0, size_of_data, 0, length_matrix + 1),
                        origin='lower')

        names_2 = []
        for i in names_col[start:stop]:
            names_2.append(self.global_go[i])

        x_ticks = np.linspace(.5, size_of_data - .5, size_of_data)
        y_ticks = np.linspace(.5, length_matrix + .5, length_matrix)

        plt.yticks(y_ticks, names_2, fontsize=16)
        if labels:
            plt.xticks(x_ticks, labels, fontsize=16, rotation='90')
        else:
            print("Provide labels")
        plt.colorbar(im, fraction=0.046, pad=0.04)
        plt.savefig(os.path.join(self.out_dir, 'Figures', '%s.pdf' % savename),
                    dpi=150, bbox_inches='tight')
        if self.save_png:
            plt.savefig(
                os.path.join(self.out_dir, 'Figures', '%s.png' % savename),
                dpi=150, bbox_inches='tight')
        plt.close()

    def analyze_data(self, data, aspect='F', savename='tmp', labels=None,
                     analyze=True, slim=None):
        """

        Parameters
        ----------
        data : list of array like
            list of samples to perform analysis
            samples are list like of gene names
        aspect : str or list like
            Can be str of list containing F, P, C
                F : Molecular Function
                P : Biological Process
                C : Cellular Component
        analyze : bool
            Perform heatplot analysis
        savename : str
            name to save plots to if analyze is True
        labels : list
            list of sample labels

        """
        use_slims = False
        if slim is not None:
            self.ontology.set_slims_subset(slim)
            use_slims = True

        data = self.create_enrichment_array(data, aspect, use_slims)
        if data is None:
            print('No data! Error! Returning nothing')
            return
        self.data_2.to_csv('{0}.csv'.format(savename))
        self.data_3.to_csv('{0}_all_data.csv'.format(savename))
        self._savename = savename
        self.data = data

        self.names, self.array = self.sort_by_hierarchy(data)

        self.go_terms = data[:, 0]

        if not analyze:
            return
        savename += '_%s' % aspect
        self.savename = savename
        # creats plots of top and bottom of hierarchy sorted arrays
        self.plot_heatmap(self.array, self.names, labels, start=0, stop=100,
                          savename='%s_top' % savename, )
        self.plot_heatmap(self.array, self.names, labels, start=-100,
                          stop=None,
                          savename='%s_bottom' % savename)
        tmp_array = self.array[:, :].copy()
        # Centroid clustering of data
        # Creates dendrogram plots
        fig = plt.figure()
        axdendro = fig.add_axes([0.09, 0.1, 0.2, 0.8])
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

        self.plot_heatmap(tmp_array, names_sorted, labels, start=0, stop=100,
                          savename='%s_top_clustered' % savename)
        self.plot_heatmap(tmp_array, names_sorted, labels, start=-100,
                          stop=None, savename='%s_bottom_clustered' % savename)

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
                self.plot_heatmap(tmp, names, labels, start=-1 * number,
                                  stop=None
                                  , savename='top_hits_entry_%i_%s' % (
                    i, savename),
                                  )
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

    def export_to_html(self, labels, html_name='tmp'):
        """

        Parameters
        ----------
        labels : list
            list of sample labels for html table
        html_name : str
            name of html page to save as

        Returns
        -------

        """

        real_names = []
        if self.array is None:
            print("Array is empty. This could be due to no signficant species"
                  "provided in list or because analysis has not been run yet.")
            return
        html_array = self.array.copy()

        # replace GO numbers with GO name
        # create plot of genes over time
        # processes = []
        print("Starting to create plots for each GO term")
        for n, i in enumerate(self.names):

            gene_set = set()
            for k in range(self.num_data_sets):
                genes = list(self.data_2[self.data_2['GO_id'] == i][
                                 'genes_{0}'.format(k)])

                for g in genes:
                    if type(g) == list:
                        for e in g:
                            gene_set.add(e)
            save_name = '{0}/Figures/go_{1}_{2}.pdf'.format(self.out_dir, i,
                                                            self.savename).replace(
                ':', '')
            title = "{0} : {1}".format(str(i), self.global_go[i])

            if len(gene_set) > 50:
                real_names.append('<a>{0}</a>'.format(self.global_go[i]))
            else:
                if save_name not in self.created_go_pds:
                    self.created_go_pds.add(save_name)
                    # processes.append((list(gene_set), save_name, '.', title, False, True))

                    # self.exp_data.plot_list_of_genes(list(gene_set), save_name,
                    #                                  title=title)
                real_names.append(
                    '<a href="Figures/go_{0}_{1}.pdf">{2} ({0})</a>'.format(i,
                        self.savename, self.global_go[i]).replace(':', ''))
        print("Done creating plots for each GO term")
        # pool = mp.Pool(4)
        # pool.map(self.exp_data.plot_list_of_genes, processes)
        # pool.close()
        # pool.join()
        # turn it into a pandas dataframe
        d = pd.DataFrame(data=html_array, index=real_names, columns=labels)

        header = """
                  <html>
                  <head>
                  <link rel="stylesheet" href="//cdn.datatables.net/1.10.13/css/jquery.dataTables.min.css"/>
                  </head>
                  <body>
                  <div style="float: left">
                  """
        footer = """
    <script src="https://code.jquery.com/jquery-3.1.1.min.js" crossorigin="anonymous"></script>
    <script src="//cdn.datatables.net/1.10.13/js/jquery.dataTables.min.js"></script>
    <script>
      $(document).ready(function(){
            $('table').DataTable({
                searchDelay:250
            });
      });
    </script>
    </body>
    </html>
    """
        # write out the html file
        with open(os.path.join(self.out_dir, '%s.html' % html_name), 'w') as f:
            f.write(header)
            # f.write(self.html_pdfs.to_html(classes='sortable', escape=False,
            #                                justify='left'))
            # f.write('\n</div>\n')
            # f.write(self.html_pdfs2.to_html(classes='sortable', escape=False,
            #                                 justify='left'))
            f.write(d.to_html(classes='sortable', float_format='{0:.4}'.format,
                              escape=False, justify='left'))
            f.write(footer)

    def export_to_html_2(self, labels, html_name='tmp'):
        """

        Parameters
        ----------
        labels : list
            list of sample labels for html table
        html_name : str
            name of html page to save as

        Returns
        -------

        """

        real_names = []
        if self.array is None:
            print("Array is empty. This could be due to no signficant species"
                  "provided in list or because analysis has not been run yet.")
            return

        data = self.data_3
        tmp = pd.pivot_table(data, index=['GO_id', 'GO_name'],
                             columns='sample_index')
        list_of_sig_go_all = data['GO_id'].unique()
        data = data[data['ref'] >= 5]
        data = data[data['pvalue'] < 0.05]
        list_of_sig_go = data['GO_id'].unique()

        processes = []
        # create plot of genes over time
        for n, i in enumerate(list_of_sig_go_all):
            # want to only plot significant species
            if i not in list_of_sig_go:
                real_names.append('<a>{0}</a>'.format(self.global_go[i]))
                continue

            # want to plot all species over time
            gene_set = set()
            for k in range(self.num_data_sets):
                genes = list(self.data_2[self.data_2['GO_id'] == i][
                                 'genes_{0}'.format(k)])
                # can return a list of lists...
                for g in genes:
                    if type(g) == list:
                        for e in g:
                            gene_set.add(e)

            # too many genes isn't helpful on plots, so skip them
            if len(gene_set) > 50:
                real_names.append('<a>{0}</a>'.format(self.global_go[i]))
                continue

            save_name = '{0}/Figures/go_{1}_{2}'.format(self.out_dir, i,
                self.savename).replace(':', '')
            title = "{0} : {1}".format(str(i), self.global_go[i])

            if save_name not in self.created_go_pds:
                self.created_go_pds.add(save_name)
                processes.append(
                    (list(gene_set), save_name, '.', title, True, True))
            real_names.append(
                '<a href="Figures/go_{0}_{1}.pdf">{2} ({0})</a>'.format(i,
                    self.savename, self.global_go[i]).replace(':', ''))
        # 95.6470000744 for parallel, 4 processors
        # 200.871999979 for sequential, 4 processors
        print(len(processes))
        print("Starting to create plots for each GO term")
        # st1 = time.time()
        # for i in processes:
        #     self.exp_data.plot_list_of_genes(i)
        # end1 = time.time()
        # print("sequential time = {}".format(end1 - st1))
        st2 = time.time()
        # pool = mp.Pool(8)
        # pool.map(self.exp_data.plot_list_of_genes, processes)
        # pool.close()
        # pool.join()
        end2 = time.time()
        print("parallel time = {}".format(end2 - st2))
        print("Done creating plots for each GO term")
        # turn it into a pandas dataframe
        # d = pd.DataFrame(data=html_array, index=real_names, columns=labels)

        print(np.shape(tmp))
        print(np.shape(real_names))
        tmp['GO_id'] = real_names

        def color_negative_red(val):
            """
            Takes a scalar and returns a string with
            the css property `'color: red'` for negative
            strings, black otherwise.
            """
            color = 'red' if val < 0 else 'black'
            return 'color: %s' % color

        # html = (tmp.style.applymap(color_negative_red).render())
        # print(html)
        html_out = os.path.join(self.out_dir, '%s.html' % html_name)
        write_single_table(tmp, html_out, 'MAGINE GO analysis')


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
