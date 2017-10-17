import io
import os
import sys

import networkx as nx
import pandas as pd

from magine.data.storage import network_data_dir
from magine.mappings import ChemicalMapper

if sys.version_info[0] == 3:
    from urllib.request import urlopen
else:
    from urllib import urlopen

_cm = ChemicalMapper()


class BioGridDownload(object):
    def __init__(self):
        self.url = 'https://thebiogrid.org/downloads/archives/Latest%20Release/BIOGRID-ALL-LATEST.tab2.zip'
        self.url2 = 'https://thebiogrid.org/downloads/archives/Latest%20Release/BIOGRID-CHEMICALS-LATEST.chemtab.zip'
        self._reverse = {"<-", "?-"}
        self._forward = {"->", "->"}
        self._db_name = 'BioGrid'

    def _create_chemical_network(self):
        chemical_int = pd.read_csv(io.BytesIO(urlopen(self.url2).read()),
                                   compression='zip',
                                   delimiter='\t',
                                   error_bad_lines=False,
                                   low_memory=False,
                                   encoding='utf-8',
                                   )

        chemical_int = chemical_int[chemical_int['Organism'] == 'Homo sapiens']
        chem_cols = ['Official Symbol',
                     'Action',
                     'Chemical Name',
                     'ATC Codes',
                     'CAS Number',
                     'Pubmed ID',
                     'Chemical Type',
                     'Chemical Source',
                     'Chemical Source ID']

        chemical_int = chemical_int[chem_cols]
        chem_table = chemical_int.as_matrix()
        chem_g = nx.DiGraph()
        nodes_added = set()
        for row in chem_table:
            target = row[0]
            source = row[2]
            name = source.encode('ascii', 'replace')
            action = row[1]
            atc_code = row[3]
            cas_num = row[4]
            pubmed = row[5]
            chem_typed = row[6]
            chem_source = row[7]
            chem_source_id = row[8]
            if chem_typed == 'small molecule':
                node_type = 'compound'
            else:
                node_type = 'biologic'

            if target not in nodes_added:
                chem_g.add_node(target, speciesType='gene',
                                databaseSource='BioGrid')
                nodes_added.add(target)

            def _add_source_node(node):
                if node not in nodes_added:
                    chem_g.add_node(node,
                                    speciesType=node_type,
                                    chemName=name,
                                    databaseSource='BioGrid',
                                    atcCode=atc_code,
                                    chemSource=chem_source,
                                    chemSourceId=chem_source_id,
                                    casNum=cas_num
                                    )
                    nodes_added.add(node)

            def _add_edge(node):
                chem_g.add_edge(node, target,
                                interactionType=action,
                                pubmed=pubmed,
                                databaseSource='BioGrid'
                                )

            if chem_source == 'DRUGBANK':
                if chem_source_id in _cm.drugbank_to_hmdb:
                    source = _cm.drugbank_to_hmdb[chem_source_id]
                    for n in source:
                        _add_source_node(n)
                        _add_edge(n)
                else:
                    _add_source_node(source)
            else:
                _add_source_node(source)

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
                            error_bad_lines=False,
                            low_memory=False)

        table = table[table['Organism Interactor A'].isin(['9606'])]
        table = table[table['Organism Interactor B'].isin(['9606'])]

        protein_cols = ['Official Symbol Interactor A',
                        'Official Symbol Interactor B',
                        'Modification',
                        'Pubmed ID',
                        'Score',
                        'Source Database',
                        'Experimental System Type']
        table = table[protein_cols]

        table = table[~table['Modification'].isin(['-', 'No Modification'])]

        table = table.as_matrix()
        g = self._create_chemical_network()
        added_genes = set()

        for r in table:
            gene1 = r[0]
            gene2 = r[1]
            inter = r[2].lower()
            pubchem = r[3]
            score = r[4]
            source_db = r[5]
            exp_system = r[6]

            def _add_node(node):
                if node not in added_genes:
                    g.add_node(node,
                               databaseSource=source_db,
                               speciesType='gene'
                               )
                    added_genes.add(node)

            _add_node(gene1)
            _add_node(gene2)

            g.add_edge(gene1, gene2,
                       interactionType=inter,
                       pubchem=pubchem,
                       expSystem=exp_system,
                       ptm=inter,
                       score=score,
                       databaseSource=source_db)

        nx.write_gpickle(g, os.path.join(network_data_dir, 'biogrid.p'))
        return g


def create_biogrid_network():
    p_name = os.path.join(network_data_dir, 'biogrid.p')
    if os.path.exists(p_name):
        g = nx.read_gpickle(p_name)
    else:

        bgn = BioGridDownload()
        g = bgn.parse_network()
    print("BIOGRID network has {} nodes and {} edges "
          "".format(len(g.nodes()), len(g.edges())))
    return g


if __name__ == '__main__':
    bgn = BioGridDownload()
    bgn.parse_network()
