import itertools

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

from magine.networks.ontology_network import OntologyNetworkGenerator
from magine.networks.visualization.cytoscapejs_tools import viewer as cyjs
from magine.networks.visualization.util_networkx import from_networkx
from magine.plotting.wordcloud_mod import WordCloud

all_dbs = [
    'GO_Biological_Process_2017',
    'GO_Molecular_Function_2017',
    'GO_Cellular_Component_2017',
    'KEGG_2016',
    'NCI-Nature_2016',
    'Panther_2016',
    'WikiPathways_2016',
    'BioCarta_2016',
    'Humancyc_2016',
    'Reactome_2016',
    'KEA_2015',
    'ChEA_2016',
    'DrugMatrix',
    'Drug_Perturbations_from_GEO_2014',
]

process_dbs = [
    'GO_Biological_Process_2017',
    'Humancyc_2016',
    'Reactome_2016',
    'KEGG_2016',
    'NCI-Nature_2016',
    'Panther_2016',
    'WikiPathways_2016',
]

pd.set_option('display.width', 10000)
pd.set_option('max_colwidth', 100)


def cleanup_term_name(row):
    x = row['term_name'].split('_')[0].lower()
    x = ' ' + x
    x = x.replace(' p53', ' tp53')
    x = x.replace('  ', ' ')
    if x[0] == ' ':
        x = x[1:]
    if x[-1] == ' ':
        x = x[:-1]
    return x


def filter_dataframe(df, p_value=0.05, db=None, sample_id=None, category=None):
    def filter(local_df, options, column):
        copy_df = local_df.copy()
        if options is not None:
            if isinstance(options, str):
                assert options in df[column].unique()
                copy_df = local_df[local_df[column] == options]
            elif isinstance(options, list):
                for i in options:
                    assert i in df[column].unique()
                copy_df = local_df[local_df[column].isin(options)]
        return copy_df

    copy_df = df.copy()
    copy_df = copy_df[copy_df['adj_p_value'] <= p_value]
    copy_df = copy_df[copy_df['combined_score'] >= 0.0]
    copy_df = filter(copy_df, db, 'db')
    copy_df = filter(copy_df, sample_id, 'sample_id')
    copy_df = filter(copy_df, category, 'category')
    copy_df.sort_values('combined_score', inplace=True, ascending=False)
    return copy_df


basic_words = {'signaling', 'signalling', 'receptor', 'events', 'protein',
               'proteins', 'regulation', 'interactions', 'via', 'signal',
               'mediated', 'pathway', 'activity', 'complex', 'positive',
               'mrna',
               'cellular', 'viral', 'host', 'processing', 'activation',
               'rrna', 'network', 'rna', 'cancer', 'disease', 'cascade',
               'transcript', 'influenza', 'beta', 'pathways', 'gene', 'hiv',
               'downstream', 'activated', 'target',
               }


# basic_words = set()


def create_wordcloud(df, save_name=None):
    data = df.apply(cleanup_term_name, axis=1)
    text = ' '.join(data)
    # Generate a word cloud image
    wc = WordCloud(margin=0, background_color=None, mode='RGBA', min_count=1,
                   width=800, height=600, collocations=True,
                   stopwords=basic_words)
    wordcloud = wc.generate(text)
    word_dict = wc.process_text(text)
    # for i,j in sorted(word_dict.items(), key=lambda p:p[1], reverse=True):
    #     print("{} : {}".format(i, j))
    # Display the generated image:
    # the matplotlib way:
    if save_name is not None:
        plt.figure()
        plt.imshow(wordcloud, interpolation='bilinear')
        plt.xticks([])
        plt.yticks([])
        plt.axis("off")
        plt.savefig('{}.png'.format(save_name), bbox_inches='tight', dpi=150)
        plt.title(save_name)
    plt.show()
    plt.close()
    return word_dict


def term_to_genes(df, term):
    genes = df[df['term_name'] == term]['genes']
    return set(itertools.chain.from_iterable(genes.str.split(',').as_matrix()))


def create_subnetwork(terms, df, network, save_name=None, draw_png=False,
                      cytoscape_js=False):
    df = df[df['term_name'].isin(terms)].copy()
    df['combined_score'] = np.abs(df['combined_score'])
    df['combined_score'] = np.log2(df['combined_score'])
    df.loc[df['combined_score'] > 150, 'combined_score'] = 150

    term_dict = dict()
    label_dict = dict()
    for i in terms:
        genes = df[df['term_name'] == i]['genes']
        genes = set(
            itertools.chain.from_iterable(genes.str.split(',').as_matrix()))
        term_dict[i] = genes
        label_dict[i] = i
        print(i, len(genes))
    all_genes = set()
    for i, j in term_dict.items():
        all_genes.update(j)
    """
    ns = NetworkSubgraphs(network)
    print("Creating subnetwork")
    g = ns.shortest_paths_between_lists(all_genes,
                                        save_name='{}_subgraph'.format(save_name),
                                        single_path=False,
                                      draw=False)


    print("Plotting subnetwork")
    create_igraph_figure(g, 'test_igraph')
    """

    ong = OntologyNetworkGenerator(molecular_network=network)
    print("Looking for direct edges")
    term_g, molecular_g = ong.create_network_from_list(
        terms, term_dict, label_dict, save_name=save_name, draw=draw_png)

    if cytoscape_js:
        return display_graph(molecular_g)
    else:
        return term_g, molecular_g


def display_graph(graph):
    new_nodes = []
    for i, data in graph.nodes(data=True):
        if 'termName' in data:
            graph.node[i]['parent'] = data['termName']
            new_nodes.append(data['termName'])
    for each in new_nodes:
        graph.add_node(each, )
    g_cyjs = from_networkx(graph)
    cyjs.render(g_cyjs, style='Directed', layout_algorithm='cose-bilkent')


def find_hits_array(enrichment_array, cateogry, sample_ids,
                    database_list=process_dbs):
    samples = []
    all_samples = []
    for i in sample_ids:
        sample_1 = filter_dataframe(enrichment_array, p_value=0.05,
                                    db=database_list,
                                    category=cateogry, sample_id=i)
        all_samples.append(sample_1)
        if len(sample_1) != 0:
            hits_1 = create_wordcloud(sample_1,
                                      save_name="{}_wordcloud".format(i))

            df1 = pd.DataFrame(hits_1.items(), columns=['words', 'counts'])
            df1.sort_values('counts', ascending=False, inplace=True)
            df1['sample'] = i
            samples.append(df1.copy())

    df = pd.concat(samples)
    samples = df['sample'].unique()
    df = pd.pivot_table(df, index='words', columns=['sample']).fillna(0)
    df.columns = df.columns.droplevel()
    df['sum'] = df[samples].sum(axis=1)
    df.sort_values('sum', ascending=False, inplace=True)
    print("\nSorted by sum of all")
    print(df.sort_values('sum', ascending=False).head(25))
    # for sample in samples:
    #     print("\nSorted by {}".format(sample))
    #     print(df.sort_values(sample, ascending=False).head(10))
    return all_samples, df


def filter_based_on_words(df, words):
    return df[df['term_name'].str.lower().str.contains('|'.join(words))]


if __name__ == '__main__':
    enrichment_array = pd.read_csv(
        '../figure3_scripts/all_cisplatin_out.csv.gz', index_col=0)
    up_rna = filter_dataframe(enrichment_array, category='rna_up',
                              db=process_dbs)
    down_rna = filter_dataframe(enrichment_array, category='rna_down',
                                db=process_dbs)
    up_proteomics = filter_dataframe(enrichment_array, db=process_dbs,
                                     category='proteomics_down')

    down_proteomics = filter_dataframe(enrichment_array, db=process_dbs,
                                       category='proteomics_down')

    # create_wordcloud(up_proteomics, 'proteomics_up_all')
    # create_wordcloud(down_proteomics, 'proteomics_down_all')

    # create_wordcloud(up_rna, 'rna_up_all')
    # create_wordcloud(down_rna, 'rna_down_all')

    prot_enriched, words_array = find_hits_array(enrichment_array,
                                                 'proteomics_up',
                                                 ['01hr', '06hr', '24hr',
                                                  '48hr'])
    print(words_array.sort_values('01hr', ascending=False))

    # added terms for 1hr
    keywords_1 = [
        'ner',
        'nucleotide excision',
        'cell cycle',
        'dna'
    ]
    #  added terms from 6hr
    keywords_6 = [
        'tp53',  # 7 hits
        'dna damage',  # 5 hits
        'cell death',  # 3 hits
        'apoptosis',  # 3 hits
    ]

    keywords_24 = [
        'dna damage',  # 12 hits
        'apoptosis',  # 8 hits
        'damage response',  # 7 hits
        'p38',
        'growth factor',
    ]

    keywords_48 = [
        'apoptotic',  # 6 hits
        'apoptosis',  # 5 hits
        'cell proliferation ',  # 4 hits
        'cell cycle',  # 4 hits
        'damage response',  #
    ]
    # filtered = filter_based_on_words(sample, keywords_1)
    # print(filtered.head(25))
    # print("There are {} terms with selected keywords".format(filtered.shape[0]))

    filtered = filter_based_on_words(prot_enriched[0], keywords_1)
    print(filtered.head(25))
    print(
    "There are {} terms with selected keywords".format(filtered.shape[0]))
    filtered = filter_based_on_words(prot_enriched[1], keywords_6)
    print(filtered.head(25))
    print(
    "There are {} terms with selected keywords".format(filtered.shape[0]))
    filtered = filter_based_on_words(prot_enriched[2], keywords_24)
    print(filtered.head(25))
    print(
    "There are {} terms with selected keywords".format(filtered.shape[0]))

    filtered = filter_based_on_words(prot_enriched[3], keywords_48)
    print(filtered.head(25))
    print(
    "There are {} terms with selected keywords".format(filtered.shape[0]))

    # quit()

    selected = [
        'Cell Cycle_Homo sapiens_R-HSA-1640170',  # 1hr
        'DNA Repair_Homo sapiens_R-HSA-73894',  # 1hr
        # 'Nucleotide Excision Repair_Homo sapiens_R-HSA-5696398'

        # 'cellular response to DNA damage stimulus ', # 6hr
        # 'Regulation of TP53 Activity_Homo sapiens_R-HSA-5633007', # 6hr
        # 'Apoptosis_Homo sapiens_R-HSA-109581', # 6hr

        # 'Caspase Cascade in Apoptosis_Homo sapiens_b9d3ef2e-618f-11e5-8ac5-06603eb7f303', # 24hr
        # 'Intrinsic Pathway for Apoptosis_Homo sapiens_R-HSA-109606',  # 24hr
        # 'G2/M DNA damage checkpoint_Homo sapiens_R-HSA-69473',
        # '',

        # 'negative regulation of apoptotic process ', # 48hr
        # 'Apoptotic execution  phase_Homo sapiens_R-HSA-75153',  # 48hr
        # 'Apoptosis Modulation and Signaling_Homo sapiens_WP1772',  # 48hr

    ]
    # network = nx.read_gml('../network_related/network_painted2.gml')
    network = nx.read_gpickle('../network_related/cisplatin_network.p')

    proteomics = enrichment_array[
        enrichment_array['category'].str.contains('proteomics')]
    create_subnetwork(selected, proteomics, network, 'test')
