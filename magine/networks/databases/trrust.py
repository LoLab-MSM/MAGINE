import os

import networkx as nx
import pandas as pd

from magine.data.storage import network_data_dir

_p_name = os.path.join(network_data_dir, 'trrust.p.gz')

url = 'http://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv'


def download_trrust():
    table = pd.read_csv(url,
                        names=['source', 'target', 'interactionType', 'pmid'],
                        delimiter='\t', index_col=None,
                        error_bad_lines=False, encoding='utf-8'
                        )
    print(table.head(10))

    # filter out non d
    table = table[~(table['interactionType'] == 'Unknown')].copy()

    table.loc[table[
                  'interactionType'] == 'Activation', 'interactionType'] = 'activate|expression'
    table.loc[table[
                  'interactionType'] == 'Repression', 'interactionType'] = 'inhibit|repression'
    table = table[~(table['interactionType'] == 'Unknown')].copy()

    table['databaseSource'] = 'TRRUST'

    protein_graph = nx.from_pandas_edgelist(
        table,
        'source',
        'target',
        edge_attr=['interactionType', 'pmid'],
        create_using=nx.DiGraph()
    )

    table = table[['source', 'target']].values
    added_genes = set()

    def _add_node(node):
        if node not in added_genes:
            protein_graph.add_node(node, databaseSource='TRRUST',
                                   speciesType='gene')
            added_genes.add(node)

    # add names to graph
    for r in table:
        _add_node(r[0])
        _add_node(r[1])

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
        print("Downloading TRRUST network!")
        download_trrust()
        assert os.path.exists(_p_name), "Error downloading TRRUST. "
    tmp_graph = nx.read_gpickle(_p_name)
    if verbose:
        print("SIGNOR : {} nodes and {} edges".format(len(tmp_graph.nodes),
                                                      len(tmp_graph.edges)))
    return tmp_graph


if __name__ == '__main__':
    load_trrust(fresh_download=True, verbose=True)
