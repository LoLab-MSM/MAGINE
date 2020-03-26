import logging
import os

import networkx as nx
import pandas as pd

from magine.data.storage import network_data_dir
from magine.networks.standards import edge_standards

_p_name = os.path.join(network_data_dir, 'signor.p.gz')
from magine.logging import get_logger

logger = get_logger(__name__, log_level=logging.INFO)


def download_signor():
    logger.info("Downloading SIGNOR")
    col_names = [
        'ENTITYA', 'TYPEA', 'IDA', 'DATABASEA', 'ENTITYB', 'TYPEB', 'IDB',
        'DATABASEB', 'EFFECT', 'MECHANISM', 'RESIDUE', 'SEQUENCE', 'TAX_ID',
        'CELL_DATA', 'TISSUE_DATA', 'MODULATOR_COMPLEX', 'TARGET_COMPLEX',
        'MODIFICATIONA', 'MODASEQ', 'MODIFICATIONB', 'MODBSEQ', 'PMID',
        'DIRECT', 'SENTENCE', 'SIGNOR_ID', 'NA1', 'NA2', 'NA3']

    table = pd.read_csv('https://signor.uniroma2.it/getData.php?organism=9606',
                        names=col_names, delimiter='\t', index_col=None,
                        error_bad_lines=False, encoding='utf-8'
                        )
    # filter out non direct
    table = table.loc[table['DIRECT'] == 't']

    # Filter out non descriptive
    table = table.loc[~table['MECHANISM'].isnull()]

    # Drop SIGNOR edges, these are generally complexes
    table = table[~(table['DATABASEA'] == 'SIGNOR')]
    table = table[~(table['DATABASEB'] == 'SIGNOR')]

    # Not sure what they mean, so will remove. Ideally other DBs have this info
    table = table[~(table['MECHANISM'] == 'post transcriptional regulation')]

    col_a = ['ENTITYA', 'TYPEA', 'IDA', 'DATABASEA']
    col_b = ['ENTITYB', 'TYPEB', 'IDB', 'DATABASEB']
    cols = ['name', 'species_type', 'id', 'db']
    species_a = table[col_a].copy()
    species_b = table[col_b].copy()
    species_a.rename(columns={i: j for i, j in zip(col_a, cols)}, inplace=True)
    species_b.rename(columns={i: j for i, j in zip(col_b, cols)}, inplace=True)
    species_a.drop_duplicates(inplace=True)
    species_b.drop_duplicates(inplace=True)
    all_species = pd.concat([species_a, species_b])
    all_species.drop_duplicates(inplace=True)

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
    table['pmid'] = table['PMID']

    table['source'] = table['ENTITYA']
    table['target'] = table['ENTITYB']

    protein_graph = nx.from_pandas_edgelist(
        table,
        'source',
        'target',
        edge_attr=['interactionType', 'databaseSource'],
        create_using=nx.DiGraph()
    )
    # add names to graph
    for row in all_species.values:
        name, species_type, id_name, db = row
        if species_type != 'protein':
            species_type = 'compound'
        if species_type == 'protein':
            species_type = 'gene'

        protein_graph.add_node(name, databaseSource='SIGNOR',
                               speciesType=species_type)

    nx.write_gpickle(protein_graph, _p_name)
    logger.info("Done downloading SIGNOR")


def load_signor(fresh_download=False):
    """
    Load reactome functional interaction network

    Parameters
    ----------
    fresh_download: bool
        Download fresh network
    verbose : bool

    Returns
    -------
    nx.DiGraph
    """
    if not os.path.exists(_p_name) or fresh_download:
        print("Downloading Signor network!")
        download_signor()
        if not os.path.exists(_p_name):
            raise FileNotFoundError("Error downloading reactome FI. ")
    tmp_graph = nx.read_gpickle(_p_name)

    logger.info("SIGNOR : {} nodes and {} edges".format(len(tmp_graph.nodes),
                                                        len(tmp_graph.edges)))
    return tmp_graph


if __name__ == '__main__':
    download_signor()
