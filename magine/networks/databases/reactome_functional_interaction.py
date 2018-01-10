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


def download_reactome_functional_interaction(verbose=False):
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
    missing_uniprot = set()
    for i in genes:
        if i not in gm.gene_name_to_uniprot:
            missing_uniprot.add(i)
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

    print("Reactome network has {} nodes "
          "and {} edges".format(len(g.nodes()), len(g.edges())))

    # nx.write_gml(g, path)
    nx.write_gpickle(g, _p_name)


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


def _check_rows(row):
    interaction_type = None
    annotation = row[2]
    g1 = row[0]
    g2 = row[1]
    direction = row[3]
    score = str(row[4])
    modification = ''
    _forward_activate = {'activate; catalyze',
                         'activate',
                         'activate; catalyze; complex',
                         }
    if annotation in _forward_activate:
        return [[g1, g2, 'activate', modification, score]]

    _forward_reverse_activate = {'activate; activated by; catalyze',
                                 'activate; activated by',
                                 'activate; activated by; catalyze; complex; input',
                                 'activate; activated by; catalyze; complex',
                                 'activate; activated by; catalyzed by',
                                 'activate; activated by; catalyzed by; complex; input',
                                 'activate; activated by; complex',
                                 'activate; activated by; complex; input',
                                 'activate; activated by; complex; input; reaction',
                                 'activate; activated by; reaction',
                                 }
    if annotation in _forward_reverse_activate:
        return [[g1, g2, 'activate', modification, score],
                [g2, g1, 'activate', modification, score]]
    if annotation == 'expression regulated by':
        return [[g2, g1, 'expression', modification, score]]
    if annotation == 'expression regulates':
        return [[g1, g2, 'expression', modification, score]]

    first = []
    second = []

    def find_hits():
        inh_b = True
        exp_b = True
        cat_b = True
        act_b = True
        if 'inhibited by' in annotation:
            inh_b = True
        if 'expression regulated by' in annotation:
            exp_b = True
        if 'expression by' in annotation:
            exp_b = True
        if 'catalyzed by' in annotation:
            cat_b = True
        if 'activated by' in annotation:
            act_b = True
        opts = [inh_b, exp_b, cat_b, act_b]
        if sum(opts) > 1:
            print(row)

    find_hits()
    if 'inhibited by' in annotation:
        second = [g2, g1, 'inhibit', modification, score]
    if 'expression regulated by' in annotation:
        second = [g2, g1, 'expression', modification, score]
    if 'expression by' in annotation:
        second = [g2, g1, 'expression', modification, score]
    if 'catalyzed by' in annotation:
        second = [g2, g1, 'catalyze', modification, score]
    if 'activated by' in annotation:
        second = [g2, g1, 'activate', modification, score]
    if 'dissociation' in annotation:
        second = [g2, g1, 'dissociation', modification, score]

    if 'expression regulates' in annotation:
        first = [g1, g2, 'expression', modification, score]
    if 'inhibit' in annotation:
        first = [g1, g2, 'inhibit', modification, score]
    if 'activate' in annotation:
        first = [g1, g2, 'activate', modification, score]
    if 'catalyze' in annotation:
        first = [g1, g2, 'catalyze', modification, score]
    if 'GErel: expression' in annotation:
        first = [g1, g2, 'expression', modification, score]

    if direction in ['<-', '|-']:
        return [second]
    elif direction in ['->', '-|']:
        return [first]
    else:
        return [first, second]
    g1, g2 = _return_direction(row[3], row[0], row[1])

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

    return [[g1, g2, interaction_type, modification, score]]


if __name__ == '__main__':
    download_reactome_functional_interaction(verbose=True)
