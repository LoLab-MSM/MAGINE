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

_p_name = os.path.join(network_data_dir, 'reactome_fi.p.gz')


def load_reactome_fi():
    """
    Load reactome functional interaction network
    Returns
    -------

    """

    if not os.path.exists(_p_name):
        print("Downloading Reactome Functional interaction network!")
        download_reactome_fi()
        assert os.path.exists(_p_name), "Error downloading reactome FI. "
    tmp_graph = nx.read_gpickle(_p_name)
    n_n, n_e = len(tmp_graph.nodes()), len(tmp_graph.edges())
    print("Reactome network has {} nodes and {} edges".format(n_n, n_e))
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


def download_reactome_fi():
    """
    Downloads reactome functional interaction network

    Returns
    -------

    """
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


_reverse = {"<-", "|-"}
_forward = {"->", "->"}
_both = {'<->', '<-|', '|->', '|-|'}


if __name__ == '__main__':
    download_reactome_fi()
