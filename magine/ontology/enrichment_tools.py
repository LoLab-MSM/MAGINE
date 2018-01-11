import itertools


def jaccard_index(first_set, second_set):
    """

    Parameters
    ----------
    first_set : :obj: `set`
    second_set : :obj: `set`

    Computes the similarity between two sets.
        https://en.wikipedia.org/wiki/Jaccard_index

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

    return float(len(first_set.intersection(second_set))) / len(first_set.union(second_set))


def filter_db(local_df, options, column):
    copy_df = local_df.copy()
    valid_opts = local_df[column].unique()
    if options is not None:
        if isinstance(options, str):
            assert options in valid_opts
            copy_df = local_df[local_df[column] == options]
        elif isinstance(options, list):
            for i in options:
                assert i in valid_opts
            copy_df = local_df[local_df[column].isin(options)]
    return copy_df


def filter_dataframe(df, p_value=0.05, db=None, sample_id=None, category=None):

    copy_df = df.copy()
    copy_df = copy_df[copy_df['adj_p_value'] <= p_value]
    copy_df = copy_df[copy_df['combined_score'] >= 0.0]

    copy_df = filter_db(copy_df, db, 'db')
    copy_df = filter_db(copy_df, sample_id, 'sample_id')
    copy_df = filter_db(copy_df, category, 'category')
    copy_df.sort_values('combined_score', inplace=True, ascending=False)
    return copy_df


def term_to_genes(df, term):
    genes = df[df['term_name'] == term]['genes']
    return set(itertools.chain.from_iterable(genes.str.split(',').as_matrix()))


def find_similar_terms(data, verbose=False):

    data.sort_values('combined_score', inplace=True, ascending=False)
    data_copy = data.copy()

    array = data.as_matrix()

    to_remove = set()
    for j in range(len(data)-2):
        top_term_name = array[j, 0]
        if top_term_name in to_remove:
            continue
        top = set(array[j, 4].split(','))
        if verbose:
            print("Finding matches for {}".format(array[j, 0]))

        for i in range(j+1, len(data)):
            term_name = array[i, 0]
            if top_term_name == term_name:
                continue
            new_set = set(array[i, 4].split(','))
            score = jaccard_index(top, new_set)
            if score > .75:
                to_remove.add(term_name)
                if verbose:
                    print("\t{} with jaccard index of "
                          "{}".format(term_name, score))

    data_copy = data_copy[~data_copy['term_name'].isin(to_remove)]
    print("Number of rows went from {} to {}".format(data.shape[0], data_copy.shape[0]))
    return data_copy


def filter_based_on_words(df, words):
    return df[df['term_name'].str.lower().str.contains('|'.join(words))]


if __name__ == '__main__':
    a = {0 ,1, 2}
    b = {0, 1, 3}
    jaccard_index(a, b)
