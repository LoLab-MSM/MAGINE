import io
import os
from magine.data.storage import network_data_dir
import networkx as nx
import pandas as pd
import sys


if sys.version_info[0] == 3:
    from urllib.request import urlopen
else:
    from urllib import urlopen


def load_reactome_fi():
    """
    Load reactome functional interaction network
    Returns
    -------

    """
    p_name = os.path.join(network_data_dir, 'reactome_fi.p')

    if not os.path.exists(p_name):
        print("Downloading Reactome Functional interaction network!")
        download_reactome_functional_interaction(p_name)
        assert os.path.exists(p_name), "Error downloading reactome FI. "
    return nx.read_gpickle(p_name)


def download_reactome_functional_interaction(path):
    """
    Downloads reactome functional interaction network


    Returns
    -------

    """
    url = 'http://reactomews.oicr.on.ca:8080/caBigR3WebApp2015/FIsInGene_031516_with_annotations.txt.zip'

    table = pd.read_csv(io.BytesIO(urlopen(url).read()), compression='zip',
                        delimiter='\t', error_bad_lines=False,  encoding='utf-8'
                        )
    table = table.as_matrix()
    g = nx.DiGraph()
    added_genes = set()

    def _add_node(node):
        if node not in added_genes:
            g.add_node(node, databaseSource='ReactomeFI',
                       speciesType='gene')
            added_genes.add(node)

    for r in table:
        gene1, gene2, inter, prediction, mod, score = _check_rows2(r)
        if gene1 == gene2:
            continue
        else:
            _add_node(gene1)
            _add_node(gene2)

            g.add_edge(gene1, gene2, interactionType=inter,
                       ptm=mod, score=score,
                       databaseSource='ReactomeFI')
    ge = set(g.nodes())
    g1 = set(table[:, 0])
    g2 = set(table[:, 1])
    g_all = set()
    g_all.update(g1)
    g_all.update(g2)

    for i in g_all:
        if i not in ge:
            print(i)

    print("Reactome network has {} nodes "
          "and {} edges".format(len(g.nodes()), len(g.edges())))

    # nx.write_gml(g, path)
    nx.write_gpickle(g, path)


_reverse = {"<-", "?-"}
_forward = {"->", "->"}


def _return_direction(text, gene1, gene2):
    if text in _reverse:
        return gene2, gene1
    else:
        return gene1, gene2


def _check_rows2(row):
    # print(row)
    prediction = ''
    interaction_type = None
    annotation = row[2]
    g1, g2 = _return_direction(row[3], row[0], row[1])
    score = row[4]
    modification = ''
    if 'predicted' in annotation:
        prediction = 'True'
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

    return g1, g2, interaction_type, prediction, modification, score
