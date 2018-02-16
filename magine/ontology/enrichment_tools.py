import itertools


def jaccard_index(first_set, second_set):
    """
    Computes the similarity between two sets.
        https://en.wikipedia.org/wiki/Jaccard_index

    Parameters
    ----------
    first_set : :obj: `set`
    second_set : :obj: `set`



    Returns
    -------
    index : float


    References
    ----------
    .. [1] `Wikipedia entry for the Jaccard index
           <https://en.wikipedia.org/wiki/Jaccard_index>`_
    """
    first_set = set(first_set)
    second_set = set(second_set)

    return float(len(first_set.intersection(second_set))) / len(
        first_set.union(second_set))


def filter_db(local_df, options, column):
    copy_df = local_df.copy()
    valid_opts = list(local_df[column].unique())
    if isinstance(options, str):
        if options not in valid_opts:
            print('{} not in {}'.format(options, valid_opts))

        copy_df = local_df[local_df[column] == options]
    elif isinstance(options, list):
        for i in options:
            if i not in valid_opts:
                print('{} not in {}'.format(i, valid_opts))
        copy_df = local_df[local_df[column].isin(options)]
    return copy_df


def filter_dataframe(df, p_value=None, combined_score=None, db=None,
                     sample_id=None, category=None):
    copy_df = df.copy()
    if p_value is not None:
        assert isinstance(p_value, float)
        copy_df = copy_df[copy_df['adj_p_value'] <= p_value]
    if combined_score is not None:
        assert isinstance(combined_score, float)
        copy_df = copy_df[copy_df['combined_score'] >= combined_score]

    copy_df = filter_db(copy_df, db, 'db')
    copy_df = filter_db(copy_df, sample_id, 'sample_id')
    copy_df = filter_db(copy_df, category, 'category')
    copy_df.sort_values('combined_score', inplace=True, ascending=False)
    return copy_df


def term_to_genes(df, term):
    genes = df[df['term_name'] == term]['genes']
    return set(itertools.chain.from_iterable(genes.str.split(',').as_matrix()))


def filter_based_on_words(df, words):
    return df[df['term_name'].str.lower().str.contains('|'.join(words))]


def filter_similar_terms(data, threshold=0.75, verbose=False):
    data.sort_values('combined_score', inplace=True, ascending=False)
    data_copy = data.copy()

    to_remove = _terms_to_remove(data_copy, threshold, verbose)

    data_copy = data_copy[~data_copy['term_name'].isin(to_remove)]
    print("Number of rows went from {} to {}".format(data.shape[0],
                                                     data_copy.shape[0]))
    return data_copy


def remove_redundant(data, threshold=0.75, verbose=False):
    data.sort_values('combined_score', inplace=True, ascending=False)

    data_copy = data.copy()
    sample_ids = list(sorted(data['sample_id'].unique()))

    to_remove = set()
    for i in sample_ids:
        tmp = data_copy[data_copy['sample_id'] == i]
        redundant_terms = _terms_to_remove(tmp, threshold, verbose)
        to_remove.update(redundant_terms)

    data_copy = data_copy[~data_copy['term_name'].isin(to_remove)]
    print("Number of rows went from {} to {}".format(data.shape[0],
                                                     data_copy.shape[0]))
    return data_copy


def _terms_to_remove(data, threshold=0.75, verbose=False):

    array = data[['term_name', 'genes']].as_matrix()

    to_remove = set()
    for j in range(len(data) - 2):
        top_term_name = array[j, 0]
        top = set(array[j, 1].split(','))
        if verbose:
            print("Finding matches for {}".format(array[j, 0]))

        for i in range(j + 1, len(data)):
            term_name = array[i, 0]
            if top_term_name == term_name:
                continue
            new_set = set(array[i, 1].split(','))
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


if __name__ == '__main__':
    a = {0, 1, 2}
    b = {0, 1, 3}
    b = {0, 1}
    print(b.difference(a))
    jaccard_index(a, b)
