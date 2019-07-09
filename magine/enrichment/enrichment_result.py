import itertools
import operator
from collections import OrderedDict
from itertools import combinations

import numpy as np
import pandas as pd

from magine.data.base import BaseData
from magine.plotting.heatmaps import cluster_distance_mat

# Will be OK in Python 2
try:
    basestring
# Allows isinstance(foo, basestring) to work in Python 3
except:
    basestring = str

sig = 'significant'


def load_enrichment_csv(file_name, **args):
    """ Load data into EnrichmentResult data class

    Parameters
    ----------
    file_name : str

    Returns
    -------
    EnrichmentResult

    """
    d = pd.read_csv(file_name, **args)
    return EnrichmentResult(d)


class EnrichmentResult(BaseData):

    def __init__(self, *args, **kwargs):
        super(EnrichmentResult, self).__init__(*args, **kwargs)
        self._index = 'term_name'
        self._identifier = 'term_name'
        self._value_name = 'combined_score'
        self._sample_id_name = 'sample_id'

    @property
    def _constructor(self):
        return EnrichmentResult

    def filter_rows(self, column, options, inplace=False):
        """
        Filters a pandas dataframe provides a column and filter selection.

        Parameters
        ----------
        column : str
        options : str, list
            Can be a single entry or a list
        inplace : bool
            Filter inplace
        Returns
        -------
        pd.DataFrame
        """

        new_data = self.copy()
        valid_opts = sorted(new_data[column].unique())
        if isinstance(options, str):
            if options not in valid_opts:
                print('{} not in {}'.format(options, valid_opts))
            else:
                new_data = new_data[new_data[column] == options]
        elif isinstance(options, list):
            for i in options:
                if i not in valid_opts:
                    print('{} not in {}'.format(i, valid_opts))
            new_data = new_data[new_data[column].isin(options)]

        if inplace:
            self._update_inplace(new_data)
        else:
            return new_data

    def filter_multi(self, p_value=None, combined_score=None, db=None,
                     sample_id=None, category=None, rank=None, inplace=False):
        """
        Filters an enrichment array.

        This is an aggregate function that allows ones to filter an entire
        dataframe with a single function call.

        Parameters
        ----------
        p_value : float
            filters all values less than or equal
        combined_score : float
            filters all values greater than or equal
        db : str, list
        sample_id : str, list
        category : str, list
        rank : int
        inplace : bool
            Filter inplace

        Returns
        -------
        new_data : EnrichmentResult
        """
        new_data = self.copy()
        if p_value is not None:
            if not isinstance(p_value, (int, float)):
                raise AssertionError("p_value must be a float or int")
            new_data = new_data[new_data['adj_p_value'] <= p_value]
        if combined_score is not None:
            if not isinstance(combined_score, (int, float)):
                raise AssertionError("combined_score must be a float or int")
            new_data = new_data[new_data['combined_score'] >= combined_score]
        if isinstance(rank, (int, float)):
            new_data = new_data[new_data['rank'] <= rank]
        if db is not None:
            new_data.filter_rows('db', db, inplace=True)
        if sample_id is not None:
            new_data.filter_rows('sample_id', sample_id, inplace=True)
        if category is not None:
            new_data.filter_rows('category', category, inplace=True)
        if inplace:
            self._update_inplace(new_data)
        else:
            return new_data

    def term_to_genes(self, term):
        """ Get set of genes of provides term(s)

        Parameters
        ----------
        term : str, list

        Returns
        -------
        set

        """
        if isinstance(term, basestring):
            genes = self[self['term_name'].isin([term])]['genes']
        else:
            genes = self[self['term_name'].isin(term)]['genes']
        return set(itertools.chain.from_iterable(genes.str.split(',').values))

    def term_to_genes_dict(self, term_list=None):
        """

        Parameters
        ----------
        term_list : list

        Returns
        -------
        OrderedDict

        """
        if term_list is None:
            term_list = set(self['term_name'].values)
        elif isinstance(term_list, basestring):
            term_list = [term_list]
        gene_to_term = {}
        for term in term_list:
            gene_list = self.term_to_genes(term)
            for g in gene_list:
                if g not in gene_to_term:
                    gene_to_term[g] = set()
                gene_to_term[g].add(term)
        term_to_gene = {}
        for i, j in gene_to_term.items():
            name = ','.join(sorted(j))
            if name not in term_to_gene:
                term_to_gene[name] = set()
            term_to_gene[name].add(i)
        return OrderedDict(
            sorted(term_to_gene.items(), key=operator.itemgetter(0))
        )

    def all_genes_from_df(self):
        """ Returns all genes from gene columns in a set

        Returns
        -------
        set
        """
        return set(
            itertools.chain.from_iterable(self['genes'].str.split(',').values)
        )

    def filter_based_on_words(self, words, inplace=False):
        """ Filter term_name based on key terms

        Parameters
        ----------
        words : list, str
            List of words to use to keep rows in dataframe
        inplace : bool
            Filter the dataframe in place or return filtered copy

        Returns
        -------
        pandas.DataFrame

        """
        if isinstance(words, str):
            words = [words]
        df = self.copy()
        df = df[df['term_name'].str.lower().str.contains('|'.join(words))]
        if inplace:
            self._update_inplace(df)
        else:
            return df

    def find_similar_terms(self, term, level='sample', remove_subset=True):
        """ Calculates similarity of all other terms to given term

        Parameters
        ----------
        term : str
        level : str
            Sample or dataframe level, flattens all terms to one set of genes
        remove_subset : bool
            If any term is a subset of the other term, a score of 1 will be
            used instead of jaccard index.


        Returns
        -------
        pd.DataFrame
        """

        rest_of_df = self[~(self['term_name'] == term)].copy()
        first_genes = self.term_to_genes(term)

        if level == 'dataframe':
            vals = [[i, self.term_to_genes(i)] for i in
                    rest_of_df['term_name'].unique()]
        else:
            vals = [[i, j.split(',')] for i, j in
                    rest_of_df[['term_name', 'genes']].values]

        dist_m = [[i, jaccard_index(first_genes, j, remove_subset)]
                  for i, j in vals]

        df = pd.DataFrame(dist_m, columns=['term_name', 'similarity_score'])
        df.sort_values('similarity_score', inplace=True, ascending=False)
        return df

    def show_terms_below(self, term, level='dataframe', threshold=.7,
                         remove_subset=True):
        """
        Find terms that were removed by remove_redundant

        Parameters
        ----------
        term : str
        level : str
        threshold : float
        remove_subset : bool

        Returns
        -------
        EnrichmentResult
        """
        temp_df = self.copy()
        # calculate similarity of term to all terms
        sim_terms = temp_df.find_similar_terms(term,
                                               remove_subset=remove_subset,
                                               level=level)
        # gather terms that are highly similar
        sim_terms = sim_terms.loc[sim_terms.similarity_score >= threshold]
        high_similar_terms = set(sim_terms.term_name.values)

        # this shows all terms that remained after filtering
        # Will be used to compare against
        term_kept = temp_df.remove_redundant(
            threshold=threshold,
            level=level,
            inplace=False,
            sort_by='combined_score'
        )

        term_kept = set(term_kept.term_name.values)
        terms_removed = high_similar_terms.difference(term_kept)
        # Adding the original term to show similarity
        terms_removed.add(term)
        return temp_df.loc[temp_df.term_name.isin(terms_removed)]

    def remove_redundant(self, threshold=0.75, verbose=False, level='sample',
                         sort_by='combined_score', inplace=False):
        """
        Calculate similarity between all term sets and removes redundant terms.

        Parameters
        ----------
        threshold : float, default 0.75
        verbose : bool, default False
            Print similarity scores and removed terms.
        level : {'sample', 'dataframe'}, default 'sample'
            Level to filter dataframe. 'sample' will pivot the dataframe and
            filter each group of 'sample_id' individually. 'dataframe' will
            merge all genes that share the same 'term_name'.
        sort_by : {'combined_score', 'rank', 'adj_p_value', 'n_genes'},
            default 'combined_score'
            Keyword to sort the dataframe. The scoring starts at the top term and
            compares to all the lower terms. Options are
        inplace : bool
            Filter the dataframe in place or return filtered copy

        Returns
        -------
        pandas.DataFrame

        """

        if sort_by in ('rank', 'adj_p_value'):
            ascending = True
        else:
            ascending = False

        self.sort_values(sort_by, inplace=True, ascending=ascending)
        data_copy = self.copy()
        if 'sample_id' not in data_copy.columns or level == 'dataframe':
            to_keep = data_copy.unique_terms(threshold, verbose, level=level)
        else:
            to_keep = set()
            for i in sorted(data_copy['sample_id'].unique()):
                tmp = data_copy[data_copy['sample_id'] == i]
                to_keep.update(
                    tmp.unique_terms(threshold, verbose, level=level)
                )

        data_copy = data_copy[(data_copy['term_name'].isin(to_keep))]
        print("Number of rows went from {} to {}"
              "".format(len(self.term_name.unique()),
                        len(data_copy.term_name.unique()))
              )
        if inplace:
            self._update_inplace(data_copy)
        else:
            return data_copy

    def unique_terms(self, threshold=0.75, verbose=False, level='dataframe'):
        """

        Parameters
        ----------
        threshold : float
        verbose : bool
        level : str, {'dataframe', 'each'}

        Returns
        -------

        """
        if level == 'dataframe':
            names = self['term_name'].unique()
            scores = self._get_distance_all()
        else:
            names = self['term_name'].values
            scores = self._get_distance_each()
        to_remove, to_keep = set(), set()

        n_dim = len(names)
        ind = 0
        for i, term_1 in enumerate(names):
            if term_1 not in to_remove:
                to_keep.add(term_1)
            else:
                for j in range(n_dim):
                    if i >= j:
                        continue
                    ind += 1
                continue
            if verbose:
                print("Finding matches for {}".format(term_1))
            for j, term_2 in enumerate(names):
                if i >= j:
                    continue
                score = scores[ind]
                ind += 1
                if score > threshold:
                    to_remove.add(term_2)
                if verbose:
                    print("\tScore for {} is {:.3f}".format(term_2, score))
                    if verbose:
                        print("\t\tRemoving {}".format(term_2))

        return to_keep

    def dist_matrix(self, figsize=(8, 8), level='dataframe'):
        """ Create a distance matrix of all term similarity

        Parameters
        ----------
        figsize : tuple
            Size of figure
        level : str, {'dataframe', 'each'}
            How to treats term_name to genes. Dataframe compresses all genes
            from all sample_ids into same term. 'each' treats each term_name
            individually.

        Returns
        -------
        matplotlib.Figure

        """
        mat, names = self.calc_dist(level)
        return cluster_distance_mat(mat, names, figsize)

    def calc_dist(self, level='datafame'):
        if level == 'each':
            names = self['term_name'].values
            scores = self._get_distance_each()
        else:
            names = self['term_name'].unique()
            scores = self._get_distance_all()
        n_dim = len(names)

        mat = np.ones((n_dim, n_dim), dtype=float)
        ind = 0
        for i in range(n_dim):
            for j in range(n_dim):
                if i >= j:
                    continue
                mat[i, j] = scores[ind]
                mat[j, i] = scores[ind]
                ind += 1
        return mat, names

    def _get_distance_each(self):
        vals = [set(i.split(',')) for i in self['genes'].values]
        return list(map(_score, combinations(vals, 2)))

    def _get_distance_all(self):
        vals = [self.term_to_genes(i) for i in self['term_name'].unique()]
        return list(map(_score, combinations(vals, 2)))


def _score(vals):
    return jaccard_index(vals[0], vals[1])


def jaccard_index(set1, set2, remove_subset=True):
    """
    Computes the similarity between two sets.
        https://en.wikipedia.org/wiki/Jaccard_index

    Parameters
    ----------
    set1 : set
    set2 : set
    remove_subset : bool
        If a set is a subset of the other, return 1.


    Returns
    -------
    index : float


    References
    ----------
    .. [1] `Wikipedia entry for the Jaccard index
           <https://en.wikipedia.org/wiki/Jaccard_index>`_
    """
    union = len(set1.union(set2))
    max_size = max(len(set1), len(set2))
    if union == max_size and remove_subset:
        return 1.
    return float(len(set1.intersection(set2))) / float(union)
