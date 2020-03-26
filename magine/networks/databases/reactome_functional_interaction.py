import io
import logging
import os
import urllib.request

import networkx as nx
import pandas as pd

from magine.data.storage import network_data_dir
from magine.logging import get_logger

_p_name = os.path.join(network_data_dir, 'reactome_fi.p.gz')
logger = get_logger(__name__, log_level=logging.INFO)


def load_reactome_fi():
    """
    Load reactome functional interaction network

    Returns
    -------
    pandas.DataFrame
    """

    if not os.path.exists(_p_name):
        download_reactome_fi()
        if not os.path.exists(_p_name):
            raise FileNotFoundError("Error downloading reactome FI. ")
    tmp_graph = nx.read_gpickle(_p_name)

    n_n, n_e = len(tmp_graph.nodes()), len(tmp_graph.edges())
    logger.info("Reactome : {} nodes and {} edges".format(n_n, n_e))
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
    name = name.replace(':', ' ').replace(';', ' ').replace(',', ' ')
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


def download_reactome_fi():
    """ Downloads reactome functional interaction network """

    logger.info("Downloading Reactome Functional interaction network")

    url = 'http://cpws.reactome.org/caBigR3WebApp2019/FIsInGene_020720_with_annotations.txt.zip'

    req = urllib.request.Request(url)
    with urllib.request.urlopen(req) as response:
        # the_page = response.read()
        table = pd.read_csv(
            io.BytesIO(response.read()), compression='zip', delimiter='\t',
            error_bad_lines=False, encoding='utf-8'
        )
    table = table[table['Direction'] != '-']
    table = table[~table['Annotation'].str.contains('indirect effect')]
    table = table[~table['Annotation'].str.contains('predicted')]
    table = table[~table['Annotation'].str.contains('compound')]
    genes = set(table['Gene1'])
    genes.update(set(table['Gene2']))
    from magine.mappings.gene_mapper import GeneMapper
    gm = GeneMapper()
    missing_uniprot = set(i for i in genes if i not in gm.gene_name_to_uniprot)
    table = table[~table['Gene1'].isin(missing_uniprot)]
    table = table[~table['Gene2'].isin(missing_uniprot)]

    table['source'] = table['Gene1']
    table['target'] = table['Gene2']
    table['databaseSource'] = 'ReactomeFI'
    rev_cols = table['Direction'].isin(_reverse)

    table.loc[rev_cols, ['source', 'target']] = \
        table.loc[rev_cols, ['target', 'source']].values

    table['interactionType'] = table.apply(standardize_edge_types, axis=1)
    protein_graph = nx.from_pandas_edgelist(
        table,
        'source',
        'target',
        edge_attr=['interactionType', 'databaseSource'],
        create_using=nx.DiGraph()
    )
    species = set(table['source'].unique()).union(
        set(table['target'].unique())
    )

    # add names to graph
    for node in species:
        protein_graph.add_node(node, databaseSource='ReactomeFI',
                               speciesType='gene')

    logger.info("Done downloading Reactome functional interaction network")
    nx.write_gpickle(protein_graph, _p_name)


_reverse = {"<-", "|-"}
_forward = {"->", "->"}
_both = {'<->', '<-|', '|->', '|-|'}


if __name__ == '__main__':
    download_reactome_fi()
