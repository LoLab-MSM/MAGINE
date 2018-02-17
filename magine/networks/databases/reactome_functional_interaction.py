import io
import os
import sys

import networkx as nx
import pandas as pd

from magine.data.storage import network_data_dir

if sys.version_info[0] == 3:
    from urllib.request import urlopen
else:
    from urllib import urlopen

_p_name = os.path.join(network_data_dir, 'reactome_fi.p')


def load_reactome_fi():
    """
    Load reactome functional interaction network
    Returns
    -------

    """

    if not os.path.exists(_p_name):
        print("Downloading Reactome Functional interaction network!")
        download_reactome_functional_interaction()
        assert os.path.exists(_p_name), "Error downloading reactome FI. "
    tmp_graph = nx.read_gpickle(_p_name)
    print("Reactome network has {} nodes and {} edges".format(
        len(tmp_graph.nodes()), len(tmp_graph.edges()))
    )
    return tmp_graph


_maps = {

    'inhibited': 'inhibit',
    'inhibite': 'inhibit',
    'inhibition': 'inhibit',

    'expression': 'expression',
    'expressed': 'expression',

    'repressed': 'repression',

    'catalyzed': 'catalyze',

    'complex': 'binding',
    'dissociation': 'binding',

    'activated': 'activate',
    'activation': 'activate',

    'phosphorylated': 'phosphorylate',
    'phosphorylation': 'phosphorylate',

    'dephosphorylation': 'dephosphorylate',
    'dephosphorylated': 'dephosphorylate',

    'ubiquitinated': 'ubiquitinate',
    'ubiquitination': 'ubiquitinate',
    'glycosylation': 'glycosylate',
    'methylation': 'methylate',

    'binding/association': 'binding',

}


def standardize_edge_types(row):
    name = row['Annotation']
    name = name.replace(':', ' ')
    name = name.replace(';', ' ')
    name = name.replace(',', ' ')
    name = name.replace('state change', 'stateChange')
    name = set(name.split(' '))
    to_remove = ['by', 'regulates', 'regulated', 'input',
                 'PPrel', 'PCrel', 'GErel', 'ECrel', 'interaction',
                 '',
                 ]
    for i in to_remove:
        if i in name:
            name.remove(i)
    for k, v in _maps.items():
        if k in name:
            name.remove(k)
            name.add(v)
    name = '|'.join(sorted(name))
    return name


def old_download_reactome_functional_interaction():
    """
    Downloads reactome functional interaction network


    Returns
    -------

    """
    # url = 'http://reactomews.oicr.on.ca:8080/caBigR3WebApp2015/FIsInGene_031516_with_annotations.txt.zip'
    url = 'http://reactomews.oicr.on.ca:8080/caBigR3WebApp2016/FIsInGene_022717_with_annotations.txt.zip'

    table = pd.read_csv(io.BytesIO(urlopen(url).read()), compression='zip',
                        delimiter='\t', error_bad_lines=False, encoding='utf-8'
                        )

    table = table[table['Direction'] != '-']
    table = table[~table['Annotation'].str.contains('indirect effect')]
    table = table[~table['Annotation'].str.contains('predicted')]
    table = table[~table['Annotation'].str.contains('compound')]
    x = set(table['Gene1'])
    y = set(table['Gene2'])
    genes = x.union(y)
    from magine.mappings.gene_mapper import GeneMapper
    gm = GeneMapper()
    missing_uniprot = set(i for i in genes if i not in gm.gene_name_to_uniprot)
    table = table[~table['Gene1'].isin(missing_uniprot)]
    table = table[~table['Gene2'].isin(missing_uniprot)]

    table = table.as_matrix()

    g = nx.DiGraph()
    added_genes = set()

    def _add_node(node):
        if node not in added_genes:
            g.add_node(node, databaseSource='ReactomeFI',
                       speciesType='gene')
            added_genes.add(node)

    for r in table:
        gene1, gene2, inter, mod, score = _check_rows2(r)
        if gene1 == gene2:
            continue
        else:
            if inter != '?':
                _add_node(gene1)
                _add_node(gene2)
                if mod != '':
                    inter = inter + '|' + mod

                g.add_edge(gene1, gene2, interactionType=inter,
                           score=score, databaseSource='ReactomeFI')
                if r[3] in _both:
                    g.add_edge(gene2, gene1, interactionType=inter,
                               score=score, databaseSource='ReactomeFI')

    print("Reactome network has {} nodes "
          "and {} edges".format(len(g.nodes()), len(g.edges())))

    # nx.write_gml(g, path)
    nx.write_gpickle(g, _p_name)


def download_reactome_functional_interaction():
    url = 'http://reactomews.oicr.on.ca:8080/caBigR3WebApp2016/FIsInGene_022717_with_annotations.txt.zip'

    table = pd.read_csv(io.BytesIO(urlopen(url).read()), compression='zip',
                        delimiter='\t', error_bad_lines=False, encoding='utf-8'
                        )
    table = table[table['Direction'] != '-']
    table = table[~table['Annotation'].str.contains('indirect effect')]
    table = table[~table['Annotation'].str.contains('predicted')]
    table = table[~table['Annotation'].str.contains('compound')]
    x = set(table['Gene1'])
    y = set(table['Gene2'])
    genes = x.union(y)
    from magine.mappings.gene_mapper import GeneMapper
    gm = GeneMapper()
    missing_uniprot = set(i for i in genes if i not in gm.gene_name_to_uniprot)
    table = table[~table['Gene1'].isin(missing_uniprot)]
    table = table[~table['Gene2'].isin(missing_uniprot)]

    table['source'] = table['Gene1']
    table['target'] = table['Gene2']
    table['databaseSource'] = 'ReactomeFI'
    table['interactionType'] = table.apply(standardize_edge_types, axis=1)

    rev_cols = table['Direction'].isin(_reverse)
    table.loc[rev_cols, ['source', 'target']] = \
        table.loc[rev_cols, ['target', 'source']].values

    protein_graph = nx.from_pandas_dataframe(
        table,
        'source',
        'target',
        edge_attr=['interactionType', 'databaseSource'],
        create_using=nx.DiGraph()
    )

    table = table.as_matrix(['source', 'target'])
    added_genes = set()

    def _add_node(node):
        if node not in added_genes:
            protein_graph.add_node(node, databaseSource='ReactomeFI',
                                   speciesType='gene')
            added_genes.add(node)

    # add names to graph
    for r in table:
        _add_node(r[0])
        _add_node(r[1])
    print("Reactome network has {} nodes and {} edges"
          "".format(len(protein_graph.nodes()), len(protein_graph.edges())))

    nx.write_gpickle(protein_graph, _p_name)

def _new_download_reactome_functional_interaction(verbose=False):
    """
    Downloads reactome functional interaction network


    Returns
    -------

    """
    # url = 'http://reactomews.oicr.on.ca:8080/caBigR3WebApp2015/FIsInGene_031516_with_annotations.txt.zip'
    url = 'http://reactomews.oicr.on.ca:8080/caBigR3WebApp2016/FIsInGene_022717_with_annotations.txt.zip'

    table = pd.read_csv(io.BytesIO(urlopen(url).read()), compression='zip',
                        delimiter='\t', error_bad_lines=False, encoding='utf-8'
                        )

    print(table.dtypes)
    print(table.shape)
    table = table[table['Direction'] != '-']
    print(table.shape)
    table = table[~table['Annotation'].str.contains('indirect effect')]
    print(table.shape)
    table = table[~table['Annotation'].str.contains('predicted')]
    print(table.shape)
    table = table[~table['Annotation'].str.contains('compound')]
    print(table.shape)
    if verbose:
        for i, row in table.groupby(['Annotation', 'Direction']):
            if i[1] == '<-':
                print(i)
        # int_types = set(
        #             itertools.chain.from_iterable(
        #                     table['Annotation'].str.split(';').as_matrix()))
        # int_types = {i[1:] if i[0] == ' ' else i for i in int_types}
        # int_types = {i[7:] if i[:7] == 'PPrel: ' else i for i in int_types}
        # int_types = {i[7:] if i[:7] == 'GErel: ' else i for i in int_types}
        # int_types = {i[7:] if i[:7] == 'PCrel: ' else i for i in int_types}
        # int_types = {i[7:] if i[:7] == 'ECrel: ' else i for i in int_types}
        quit()
        directions = set(table['Direction'].unique())
        # for i in directions:
        #     print(i)

        for i in sorted(int_types):
            print(i)
    quit()

    x = set(table['Gene1'])
    y = set(table['Gene2'])
    genes = x.union(y)
    from magine.mappings.gene_mapper import GeneMapper
    gm = GeneMapper()
    missing_uniprot = set()
    for i in genes:
        if i not in gm.gene_name_to_uniprot:
            missing_uniprot.add(i)

    table = table[~table['Gene1'].isin(missing_uniprot)]
    table = table[~table['Gene2'].isin(missing_uniprot)]

    print(table.shape)
    quit()
    table = table.as_matrix()

    g = nx.DiGraph()
    added_genes = set()

    def _add_node(node):
        if node not in added_genes:
            g.add_node(node, databaseSource='ReactomeFI',
                       speciesType='gene')
            added_genes.add(node)

    for r in table:
        edges = _check_rows(r)
        for i in edges:
            gene1, gene2, inter, mod, score = i
            if gene1 == gene2:
                continue
            else:
                if inter != '?':
                    _add_node(gene1)
                    _add_node(gene2)
                    if mod != '':
                        inter = inter + '|' + mod

                    g.add_edge(gene1, gene2, interactionType=inter,
                               score=score,
                               databaseSource='ReactomeFI')

    print("Reactome network has {} nodes "
          "and {} edges".format(len(g.nodes()), len(g.edges())))

    # nx.write_gml(g, path)
    nx.write_gpickle(g, _p_name)


_reverse = {"<-", "|-"}
_forward = {"->", "->"}
_both = {'<->', '<-|', '|->', '|-|'}


def _return_direction(text, gene1, gene2):
    if text in _reverse:
        return gene2, gene1
    else:
        return gene1, gene2


def _check_rows2(row):
    interaction_type = None
    annotation = row[2]
    g1, g2 = _return_direction(row[3], row[0], row[1])
    score = str(row[4])
    modification = ''
    if 'repress' in annotation:
        interaction_type = 'repression'
    elif 'express' in annotation:
        interaction_type = 'expression'
    elif 'reaction' in annotation:
        interaction_type = 'reaction'
    elif 'complex' in annotation:
        interaction_type = 'complex'
    elif 'catalyze' in annotation:
        interaction_type = 'catalyze'
    elif 'deactiv' in annotation:
        interaction_type = 'deactivate'
    elif 'activ' in annotation:
        interaction_type = 'activate'
    elif 'inhib' in annotation:
        interaction_type = 'inhibit'
    elif 'binding' in annotation:
        interaction_type = 'binding'
    elif 'input' in annotation:
        interaction_type = '?'
    elif 'indirect effect' in annotation:
        interaction_type = 'indirect'
    elif 'compound' in annotation:
        interaction_type = 'compound'
    elif 'state change' in annotation:
        interaction_type = 'state change'
    elif 'interaction' in annotation:
        interaction_type = 'binding'
    elif 'PPrel' in annotation:
        interaction_type = 'binding'
    if 'dephosphoryl' in annotation:
        modification = 'dephosphorylation'
        if interaction_type is None:
            interaction_type = 'dephosphorylation'
    elif 'phosphoryl' in annotation:
        modification = 'phosphorylation'
        if interaction_type is None:
            interaction_type = 'phosphorylation'
    elif 'ubiquiti' in annotation:
        modification = 'ubiquitination'
        if interaction_type is None:
            interaction_type = 'ubiquitination'
    if interaction_type is None:
        interaction_type = '?'

    return g1, g2, interaction_type, modification, score


if __name__ == '__main__':
    download_reactome_functional_interaction()
