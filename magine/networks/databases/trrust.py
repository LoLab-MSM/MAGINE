import os

import networkx as nx
import pandas as pd

from magine.data.storage import network_data_dir
from magine.logging import get_logger

log = get_logger(__name__)
_p_name = os.path.join(network_data_dir, 'trrust.p.gz')

url = 'http://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv'


def download_trrust():
    log.info("Downloading TRRUST")
    df = pd.read_csv(url,
                     names=['source', 'target', 'interactionType', 'pmid'],
                     delimiter='\t', index_col=None,
                     error_bad_lines=False, encoding='utf-8'
                     )
    int_type = 'interactionType'

    # filter out non d
    df = df[~(df[int_type] == 'Unknown')].copy()

    df.loc[df[int_type] == 'Activation', int_type] = 'activate|expression'
    df.loc[df[int_type] == 'Repression', int_type] = 'inhibit|repression'

    # Merge duplicate rows while maintaining interaction type and pmid
    int_joined = df.groupby(['source', 'target'])[int_type].apply(
        '|'.join).reset_index()
    pmid_joined = df.groupby(['source', 'target'])['pmid'].apply(
        '|'.join).reset_index()

    # merge to get one dataframe
    df = pd.merge(int_joined, pmid_joined, on=['source', 'target'])

    df['databaseSource'] = 'TRRUST'

    protein_graph = nx.from_pandas_edgelist(
        df,
        'source',
        'target',
        edge_attr=[int_type, 'pmid'],
        create_using=nx.DiGraph()
    )

    nodes = set(df[['source', 'target']].values.flatten())

    # add names to graph
    for n in nodes:
        protein_graph.add_node(n, databaseSource='TRRUST', speciesType='gene')
    nx.write_gpickle(protein_graph, _p_name)


def load_trrust(fresh_download=False, verbose=False):
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
        download_trrust()
        if not os.path.exists(_p_name):
            raise AssertionError("Error downloading TRRUST. ")
    tmp_graph = nx.read_gpickle(_p_name)
    if verbose:
        log.info("TRRUST network : {} nodes and {} edges".format(
            len(tmp_graph.nodes), len(tmp_graph.edges)))
    return tmp_graph


if __name__ == '__main__':
    load_trrust(fresh_download=True, verbose=True)
