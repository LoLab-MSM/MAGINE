import types

import matplotlib.pyplot as plt
import pandas as pd
from wordcloud import STOPWORDS, WordCloud

# Will be OK in Python 2
try:
    basestring
# Allows isinstance(foo, basestring) to work in Python 3
except:
    basestring = str

process_dbs = [
    'GO_Biological_Process_2017',
    'Humancyc_2016',
    'Reactome_2016',
    'KEGG_2016',
    'NCI-Nature_2016',
    'Panther_2016',
    'WikiPathways_2016',
]

basic_words = {'signaling', 'signalling', 'receptor', 'events', 'protein',
               'proteins', 'regulation', 'interactions', 'via', 'signal',
               'mediated', 'pathway', 'activity', 'complex', 'positive',
               'mrna', 'cellular', 'viral', 'host', 'processing', 'activation',
               'rrna', 'network', 'rna', 'cancer', 'disease', 'cascade',
               'transcript', 'influenza', 'beta', 'pathways', 'gene', 'hiv',
               'downstream', 'activated', 'target', 'dn', 'up'
               }

basic_words.update(STOPWORDS)


def word_cloud_from_array(enrichment_array, sample_ids, category=None,
                          database_list=None, save_name=None, p_value=0.05):
    """
    Creates a word cloud for each sample_id


    Parameters
    ----------
    enrichment_array : magine.enrichment.tools.EnrichmentResults
    sample_ids : list
    category : str, list, optional
    database_list : str, list, optional
    save_name : str

    Returns
    -------

    """
    samples = []
    all_samples = []
    for i in sample_ids:
        sample = enrichment_array.filter_multi(
            p_value=p_value, db=database_list, sample_id=i, category=category
        )

        all_samples.append(sample)
        if len(sample) != 0:
            if save_name is not None:
                output = "{}_{}_wordcloud".format(save_name, i)
            else:
                output = None
            hits_1 = create_wordcloud(sample, save_name=output)
            df1 = pd.DataFrame(
                list(hits_1.word_dict.items()), columns=['words', 'counts']
            )
            df1.sort_values('counts', ascending=False, inplace=True)
            df1['sample'] = i
            samples.append(df1.copy())

    df = pd.concat(samples)
    samples = df['sample'].unique()
    df = pd.pivot_table(df, index='words', columns=['sample']).fillna(0)
    df.columns = df.columns.droplevel()
    df['sum'] = df[samples].sum(axis=1)
    df.sort_values('sum', ascending=False, inplace=True)
    # print("\nSorted by sum of all")
    # print(df.sort_values('sum', ascending=False).head(25))
    return all_samples, df


def create_wordcloud(df, save_name=None):
    """
    Creates a word cloud based on enrichment array

    Must have column 'term_name'.
    It takes this column, flattens it into a large string.
    Then we pass it to the python package wordcloud, which generates the
    wordcloud.

    It returns the figure with a save method and a dictionary of counts.



    Parameters
    ----------
    df : pd.DataFrame
    save_name : str

    Returns
    -------

    """
    data = df.apply(_cleanup_term_name, axis=1)
    text = ' '.join(data)
    # Generate a word cloud image
    wc = WordCloud(margin=0, background_color=None, mode='RGBA',
                   # min_count=1,
                   width=800, height=600, collocations=True,
                   stopwords=basic_words)
    wordcloud = wc.generate(text)
    word_dict = wc.process_text(text)

    def plot(self, save_name=None, figsize=(8, 5)):
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        ax.imshow(self, interpolation='bilinear')
        plt.xticks([])
        plt.yticks([])
        plt.axis("off")
        if save_name is not None:
            plt.savefig('{}.png'.format(save_name), bbox_inches='tight',
                        dpi=150)
            plt.title(save_name)
        return fig

    wordcloud.plot = types.MethodType(plot, wordcloud)

    if save_name is not None:
        wordcloud.plot(save_name)

    wordcloud.word_dict = word_dict

    df1 = pd.DataFrame(list(word_dict.items()),
                       columns=['words', 'counts'])

    df1.sort_values('counts', ascending=False, inplace=True)
    wordcloud.data = df1

    return wordcloud


def _cleanup_term_name(row):
    if not isinstance(row['term_name'], basestring):
        print(row)
    x = row['term_name'].split('_')[0].lower()
    x = ' ' + x
    x = x.replace(' p53', ' tp53')
    x = x.replace('  ', ' ')
    if x[0] == ' ':
        x = x[1:]
    if x[-1] == ' ':
        x = x[:-1]
    return x
