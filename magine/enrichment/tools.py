import itertools

import numpy as np
import pandas as pd

from magine.data import Data
from magine.plotting.heatmaps import cluster_distance_mat


def jaccard_index(first_set, second_set):
    """
    Computes the similarity between two sets.
        https://en.wikipedia.org/wiki/Jaccard_index

    Parameters
    ----------
    first_set : set
    second_set : set



    Returns
    -------
    index : float


    References
    ----------
    .. [1] `Wikipedia entry for the Jaccard index
           <https://en.wikipedia.org/wiki/Jaccard_index>`_
    """
    set1 = set(first_set)
    set2 = set(second_set)

    return float(len(set1.intersection(set2))) / len(set1.union(set2))


def load_enrichment_csv(file):
    dataframe = pd.read_csv(file)
    return EnrichmentData(dataframe)


class EnrichmentData(Data):

    def __init__(self, *args, **kwargs):
        super(EnrichmentData, self).__init__(*args, **kwargs)
        self._index = 'term_name'

    @property
    def _constructor(self):
        return EnrichmentData

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
        valid_opts = list(new_data[column].unique())
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

        """
        new_data = self.copy()
        if p_value is not None:
            assert isinstance(p_value, (int, float))
            new_data = new_data[new_data['adj_p_value'] <= p_value]
        if combined_score is not None:
            assert isinstance(combined_score, (int, float))
            new_data = new_data[new_data['combined_score'] >= combined_score]
        if isinstance(rank, (int, float)):
            new_data = new_data[new_data['rank'] <= rank]
        if db is not None:
            new_data = filter_rows(new_data, 'db', db)
        if sample_id is not None:
            new_data = filter_rows(new_data, 'sample_id', sample_id)
        if category is not None:
            new_data = filter_rows(new_data, 'category', category)
        if inplace:
            self._update_inplace(new_data)
        else:
            return new_data

    def term_to_genes(self, term):
        genes = self[self['term_name'] == term]['genes']
        return set(itertools.chain.from_iterable(genes.str.split(',').values))

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

    def find_similar_terms(self, term):

        rest_of_df = self[~(self['term_name'] == term)]
        first_genes = self.term_to_genes(term)
        array = rest_of_df[['term_name', 'genes']].values
        dist_m = [None] * len(array)
        for n, index in enumerate(array):
            dist_m[n] = [index[0],
                         jaccard_index(first_genes, set(index[1].split(',')))]

        df = pd.DataFrame(dist_m, columns=['term_name', 'similarity_score'])
        df.sort_values('similarity_score', inplace=True, ascending=False)
        return df

    def all_genes_from_df(self):
        """ Returns all genes from gene columns in a set

        Returns
        -------
        set
        """
        return set(
            itertools.chain.from_iterable(self['genes'].str.split(',').values)
        )

    def remove_redundant(self, threshold=0.75, verbose=False, level='sample',
                         sort_by='combined_score', inplace=False):
        """
        Calculate similarity between all term sets and removes redundant terms.


        Parameters
        ----------
        threshold : float, default 0.75
        verbose : bool, default False
            Print similarity scores and removed terms.
        level : {'sample', 'dataframe', 'overall'}, default 'sample'
            Level to filter dataframe. 'sample' will pivot the dataframe and
            filter each group of 'sample_id' individually. 'dataframe' will merge
            all genes that share the same 'term_name'. 'overall' will consider each
            term_name individually.
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
        data_copy = self.copy()
        if sort_by in ('rank', 'adj_p_value'):
            ascending = True
        else:
            ascending = False

        data_copy.sort_values(sort_by, inplace=True, ascending=ascending)

        if level == 'sample':
            terms_to_remove = set()
            sample_ids = list(sorted(data_copy['sample_id'].unique()))
            for i in sample_ids:
                tmp = data_copy[data_copy['sample_id'] == i]
                redundant_terms = _terms_to_remove(tmp, threshold, verbose)
                terms_to_remove.update(redundant_terms)
        elif level == 'dataframe':
            terms_to_remove = _calculate_similarity(data_copy, threshold,
                                                    verbose)
        else:
            terms_to_remove = _terms_to_remove(data_copy, threshold, verbose)

        data_copy = data_copy[~(data_copy['term_name'].isin(terms_to_remove))]
        print("Number of rows went from {} to {}".format(self.shape[0],
                                                         data_copy.shape[0]))
        if inplace:
            self._update_inplace(data_copy)
        else:
            return data_copy

    def dist_matrix(self, fig_size=(8, 8), level='dataframe'):
        """ Create a distance matrix of all term similarity

        Parameters
        ----------
        fig_size : tuple
            Size of figure
        level : str, {'dataframe', 'each'}
            How to treats term_name to genes. Dataframe compresses all genes
            from all sample_ids into same term. 'each' treats each term_name
            individually.

        Returns
        -------
        matplotlib.Figure

        """
        if level == 'each':
            names = self['term_name'].values
        else:
            names = self['term_name'].unique()
        n_dim = len(names)
        mat = np.empty((n_dim, n_dim), dtype=float)

        if level == 'each':
            array = self['genes'].values
            for i, ac in enumerate(array):
                first_genes = set(ac.split(','))
                for j, bc in enumerate(array):
                    if i > j:
                        continue
                    mat[i, j] = jaccard_index(first_genes, set(bc.split(',')))
                    mat[j, i] = mat[i, j]

        else:
            for i, n1 in enumerate(names):
                first_genes = self.term_to_genes(n1)
                for j, n2 in enumerate(names):
                    if i > j:
                        continue
                    mat[i, j] = jaccard_index(first_genes,
                                              self.term_to_genes(n2))
                    mat[j, i] = mat[i, j]

        return cluster_distance_mat(mat, names, fig_size)


def term_to_genes(df, term):
    genes = df[df['term_name'] == term]['genes']
    return set(itertools.chain.from_iterable(genes.str.split(',').values))


def _terms_to_remove(data, threshold=0.75, verbose=False):
    array = data[['term_name', 'genes']].values
    lookup = dict()
    to_remove = set()
    for j in range(len(data) - 2):
        top_term_name = array[j, 0]
        top = set(array[j, 1].split(','))
        if verbose:
            print("Finding matches for {}".format(array[j, 0]))
        if top_term_name in to_remove:
            continue
        for i in range(j + 1, len(data)):
            term_name = array[i, 0]

            if top_term_name == term_name:
                continue
            if term_name in lookup:
                new_set = lookup[term_name]
            else:
                new_set = set(array[i, 1].split(','))
                lookup[term_name] = new_set

            if len(new_set.difference(top)) == 0:
                to_remove.add(term_name)
                continue
            score = jaccard_index(top, new_set)
            if verbose:
                print("\tScore for {} is {:.3f}".format(array[i, 0], score))
            if score > threshold:
                to_remove.add(term_name)
                if verbose:
                    print("\t\tRemoving {}".format(term_name))
    return to_remove


def _calculate_similarity(df, threshold=0.75, verbose=False):
    names = df['term_name'].unique()
    to_remove = set()
    visited = set()
    for i, term_1 in enumerate(names):
        first_genes = term_to_genes(df, term_1)
        if verbose:
            print("Checking {}".format(term_1))
        for j, term_2 in enumerate(names):
            if (term_1, term_2) in visited or (term_2, term_1) in visited:
                continue
            else:
                visited.add((term_1, term_2))
                visited.add((term_2, term_1))
            if i >= j:
                continue
            second_genes = term_to_genes(df, term_2)
            if len(first_genes.difference(second_genes)) == 0:
                to_remove.add(term_2)
                if verbose:
                    print("\t'{}' is a subset of '{}'".format(term_1, term_2))
                continue
            else:
                score = jaccard_index(first_genes, second_genes)
                if verbose:
                    print("\tScore for {} is {:.3f}".format(term_2, score))
                if score > threshold:
                    to_remove.add(term_2)
                    if verbose:
                        print("\t\tRemoving {}".format(term_2))
    return to_remove


if __name__ == '__main__':
    a = {0, 1, 2}
    b = {0, 1, 3}
    print(b.difference(a))
    jaccard_index(a, b)
