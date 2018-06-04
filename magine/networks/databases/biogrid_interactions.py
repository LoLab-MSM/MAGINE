import io
import os
import sys

import networkx as nx
import pandas as pd

import magine.networks.utils as utils
from magine.data.storage import network_data_dir
from magine.mappings import ChemicalMapper

if sys.version_info[0] == 3:
    from urllib.request import urlopen
else:
    from urllib import urlopen

p_name = os.path.join(network_data_dir, 'biogrid.p.gz')


class BioGridDownload(object):
    def __init__(self):
        self.url = 'https://thebiogrid.org/downloads/archives/Latest%20Release/BIOGRID-ALL-LATEST.tab2.zip'
        self.url2 = 'https://thebiogrid.org/downloads/archives/Latest%20Release/BIOGRID-CHEMICALS-LATEST.chemtab.zip'
        self._db_name = 'BioGrid'
        self._cm = ChemicalMapper()

    def _create_chemical_network(self):
        df = pd.read_csv(io.BytesIO(urlopen(self.url2).read()),
                         compression='zip',
                         delimiter='\t',
                         error_bad_lines=False,
                         low_memory=False,
                         encoding='utf8',
                         )

        df = df[df['Organism'] == 'Homo sapiens']

        chem_cols = ['Official Symbol',
                     'Action',
                     'Chemical Name',
                     'ATC Codes',
                     'CAS Number',
                     'Pubmed ID',
                     'Chemical Type',
                     'Chemical Source',
                     'Chemical Source ID',
                     'Interaction Type'
                     ]

        df = df[chem_cols]
        df = df[~df['Action'].isin(['unknown'])]

        df.drop_duplicates(
            subset=['Official Symbol', 'Action', 'Chemical Name'],
            inplace=True
        )

        def convert_to_name(row):
            db = row['Chemical Source']
            id_chem_source = row['Chemical Source ID']
            c_name = row['chemName']
            if db == 'DRUGBANK':
                if id_chem_source in self._cm.drugbank_to_hmdb:
                    new_name = self._cm.drugbank_to_hmdb[id_chem_source][0]
                    return new_name
                elif c_name in self._cm.chem_name_to_hmdb:
                    new_name = \
                        self._cm.chem_name_to_hmdb[c_name][0]
                    return new_name
            return c_name

        def convert_to_hmdb_only(row):
            db = row['Chemical Source']
            id_chem_source = row['Chemical Source ID']
            c_name = row['chemName']
            if db == 'DRUGBANK':
                if id_chem_source in self._cm.drugbank_to_hmdb:
                    name = self._cm.drugbank_to_hmdb[id_chem_source][0]
                    return name
                elif c_name in self._cm.chem_name_to_hmdb:
                    name = self._cm.chem_name_to_hmdb[c_name][0]
                    return name
            return None

        df['databaseSource'] = self._db_name
        df['pubmedId'] = df['Pubmed ID'].astype(str)

        # cleanup names
        df.rename(columns={'Chemical Type': 'chemType',
                           'Official Symbol': 'gene',
                           'Action': 'interactionType'},
                  inplace=True)

        # keep the same info as other databases (store as compound)
        df.loc[df['chemType'] == 'small molecule', 'chemType'] = 'compound'

        # converting to ascii so we can export play with networkx
        df['chemName'] = df['Chemical Name'].astype(str)

        # convert names to HMBD, or keep it the same if HMDB doesnt exist
        df['target'] = df.apply(convert_to_name, axis=1)

        # add HMDB attribute if it exists
        df['hmdbID'] = df.apply(convert_to_hmdb_only, axis=1).astype(str)

        # create network
        chem_g = nx.from_pandas_edgelist(
            df,
            'gene',
            'target',
            edge_attr=['interactionType', 'databaseSource', 'pubmedId'],
            create_using=nx.DiGraph()
        )
        # df.to_csv('biogrid.csv')
        cols = ['gene', 'target', 'chemName', 'chemType', 'hmdbID']
        chem_table = df[cols].values

        nodes_added = set()

        def add_node(node, node_type, chem_name=None, hmdb=None):
            attr = dict()
            attr['speciesType'] = node_type
            attr['databaseSource'] = self._db_name
            if chem_name is not None:
                attr['chemName'] = chem_name
            if hmdb is not None:
                attr['hmdbNames'] = hmdb

            chem_g.add_node(node, **attr)
            nodes_added.add(node)

        # add node names/attributes
        for row in chem_table:
            gene = row[0]
            chemical = row[1]
            chemical_name = row[2]
            chem_typed = row[3]
            hmdb_id = row[4]
            if gene not in nodes_added:
                add_node(gene, 'gene')

            if chemical not in nodes_added:
                add_node(chemical, chem_typed, chemical_name, hmdb_id)
        return chem_g

    def parse_network(self):
        """
        Parses tab delimited file to networkx.DiGraph


        Returns
        -------

        """

        table = pd.read_csv(io.BytesIO(urlopen(self.url).read()),
                            compression='zip',
                            delimiter='\t',
                            encoding='utf8',
                            error_bad_lines=False,
                            low_memory=False)
        # only keep human
        # TODO enable other organisms
        table = table[table['Organism Interactor A'].isin(['9606'])]
        table = table[table['Organism Interactor B'].isin(['9606'])]

        # table.to_csv('biogrid.csv')

        protein_cols = ['Official Symbol Interactor A',
                        'Official Symbol Interactor B',
                        'Modification',
                        'Pubmed ID',
                        'Source Database'
                        ]

        table = table[protein_cols]
        # Remove puring binding and no modification
        table = table[~table['Modification'].isin(['-', 'No Modification'])]

        # clean up names
        table['interactionType'] = table['Modification'].str.lower()
        table['pubmedId'] = table['Pubmed ID'].astype(str)

        # cleanup names
        table.rename(columns={'Official Symbol Interactor A': 'source',
                              'Official Symbol Interactor B': 'target',
                              'Source Database': 'databaseSource'},
                     inplace=True)

        # create graph
        protein_graph = nx.from_pandas_edgelist(
            table,
            'source',
            'target',
            edge_attr=['interactionType', 'databaseSource', 'pubmedId'],
            create_using=nx.DiGraph()
        )

        table = table[['source', 'target']].values

        added_genes = set()

        def _add_node(node):
            protein_graph.add_node(node,
                                   databaseSource=self._db_name,
                                   speciesType='gene')
            added_genes.add(node)

        # add names to graph
        for r in table:
            if r[0] not in added_genes:
                _add_node(r[0])
            if r[1] not in added_genes:
                _add_node(r[1])

        final_graph = utils.compose(protein_graph,
                                    self._create_chemical_network())
        nx.write_gpickle(final_graph, p_name)
        return final_graph


def create_biogrid_network():
    if os.path.exists(p_name):
        g = nx.read_gpickle(p_name)
    else:
        g = BioGridDownload().parse_network()
    nn, ne = len(g.nodes()), len(g.edges())
    print("BIOGRID network has {} nodes and {} edges".format(nn, ne))
    return g


if __name__ == '__main__':
    bgn = BioGridDownload()
    bgn.parse_network()
    # create_biogrid_network()
