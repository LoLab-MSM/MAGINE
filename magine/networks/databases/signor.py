import os

import networkx as nx
import pandas as pd

from magine.data.storage import network_data_dir
from magine.networks.standards import edge_standards

_p_name = os.path.join(network_data_dir, 'signor.p.gz')


def download():
    col_names = [
        'ENTITYA', 'TYPEA', 'IDA', 'DATABASEA', 'ENTITYB', 'TYPEB', 'IDB',
        'DATABASEB', 'EFFECT', 'MECHANISM', 'RESIDUE', 'SEQUENCE', 'TAX_ID',
        'CELL_DATA', 'TISSUE_DATA', 'MODULATOR_COMPLEX', 'TARGET_COMPLEX',
        'MODIFICATIONA', 'MODASEQ', 'MODIFICATIONB', 'MODBSEQ', 'PMID',
        'DIRECT', 'SENTENCE', 'SIGNOR_ID', 'NA', 'NA', 'NA']

    table = pd.read_csv('https://signor.uniroma2.it/getData.php?organism=9606',
                        names=col_names, delimiter='\t', index_col=None,
                        error_bad_lines=False, encoding='utf-8'
                        )
    table.to_csv('signor.csv.gz', compression='gzip', index=False,
                 encoding='utf-8')

    table = pd.read_csv('signor.csv.gz', compression='gzip', encoding='utf-8',
                        index_col=None)

    # filter out non direct
    table = table[table['DIRECT'] == 't']

    # Filter out non descriptive
    table = table[~table['MECHANISM'].isnull()]

    # Drop SIGNOR edges, these are generally complexes
    table = table[~(table['DATABASEA'] == 'SIGNOR')]
    table = table[~(table['DATABASEB'] == 'SIGNOR')]

    # Not sure what they mean, so will remove. Ideally other DBs have this info
    table = table[~(table['MECHANISM'] == 'post transcriptional regulation')]

    def map_to_activate_inhibit(row):
        effect = ''
        mechanism = row['MECHANISM']
        if 'down-regulates' in row['EFFECT']:
            effect = 'inhibit'
        elif 'up-regulates' in row['EFFECT']:
            effect = 'activate'
        if mechanism in edge_standards:
            mechanism = edge_standards[mechanism]
        elif mechanism == 'transcriptional regulation':
            if effect == 'inhibit':
                mechanism = 'repression'
            elif effect == 'activate':
                mechanism = 'expression'
        if effect == '':
            return mechanism
        else:
            return "|".join([effect, mechanism])

    # relabel edge types
    table['interactionType'] = table.apply(map_to_activate_inhibit, axis=1)
    table['databaseSource'] = 'SIGNOR'
    table['pmid'] = 'PMID'

    table['source'] = table['ENTITYA']
    table['target'] = table['ENTITYB']

    protein_graph = nx.from_pandas_edgelist(
        table,
        'source',
        'target',
        edge_attr=['interactionType', 'databaseSource'],
        create_using=nx.DiGraph()
    )

    table = table[['source', 'target']].values
    added_genes = set()

    def _add_node(node):
        if node not in added_genes:
            protein_graph.add_node(node, databaseSource='SIGNOR',
                                   speciesType='gene')
            added_genes.add(node)

    # add names to graph
    for r in table:
        _add_node(r[0])
        _add_node(r[1])

    nx.write_gpickle(protein_graph, _p_name)


def load_signor():
    """
    Load reactome functional interaction network
    Returns
    -------

    """

    if not os.path.exists(_p_name):
        print("Downloading Reactome Functional interaction network!")
        download()
        assert os.path.exists(_p_name), "Error downloading reactome FI. "
    tmp_graph = nx.read_gpickle(_p_name)
    print("SIGNOR network has {} nodes and {} edges".format(
        len(tmp_graph.nodes()), len(tmp_graph.edges()))
    )
    return tmp_graph


if __name__ == '__main__':
    download()
