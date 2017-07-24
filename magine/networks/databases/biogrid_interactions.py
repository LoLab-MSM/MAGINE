import io
import os
import sys
import networkx as nx
import pandas as pd
if sys.version_info[0] == 3:
    from urllib.request import urlopen
else:
    from urllib import urlopen
from magine.mappings.chemical_mapper import ChemicalMapper



directory = os.path.dirname(__file__)


class BioGridDownload(object):


    def __init__(self):
        self.out_dir = directory
        self.url = 'https://thebiogrid.org/downloads/archives/Latest%20Release/BIOGRID-ALL-LATEST.tab2.zip'
        self.url2 = 'https://thebiogrid.org/downloads/archives/Latest%20Release/BIOGRID-CHEMICALS-LATEST.chemtab.zip'
        self._reverse = {"<-", "?-"}
        self._forward = {"->", "->"}
        self._db_name = 'BioGrid'
        self.cm = ChemicalMapper()
    def _create_chemical_network(self):
        chemical_int = pd.read_csv(io.BytesIO(urlopen(self.url2).read()),
                                   compression='zip',
                                   delimiter='\t',
                                   error_bad_lines=False,
                                   low_memory=False)

        chemical_int = chemical_int[chemical_int['Organism'] == 'Homo sapiens']
        chem_cols = ['Official Symbol', 'Action',
                     'Chemical Name', 'ATC Codes', 'CAS Number', 'Pubmed ID',
                     'Chemical Type', 'Chemical Source', 'Chemical Source ID']

        chemical_int = chemical_int[chem_cols]
        chem_table = chemical_int.as_matrix()
        chem_g = nx.DiGraph()
        nodes_added = set()
        for row in chem_table:
            target = row[0]
            source = row[2]
            action = row[1]
            atc_code = row[3]
            cas_num = row[4]
            pubmed = row[5]
            chem_typed = row[6]
            chem_source = row[7]
            chem_source_id = row[8]
            if chem_source == 'DRUGBANK':
                name = target
                source = self.cm.drugbank_to_hmdb[chem_source_id]
            if source not in nodes_added:
                if chem_typed == 'small molecule':
                    node_type = 'compound'
                else:
                    node_type = 'biologic'
                chem_g.add_node(source, speciesType=node_type,
                                databaseSource='BioGrid',
                                atcCode=atc_code,
                                chemSource=chem_source,
                                chemSourceId=chem_source_id,
                                casNum=cas_num
                                )
                nodes_added.add(source)
            if target not in nodes_added:
                chem_g.add_node(target, speciesType='gene',
                                databaseSource='BioGrid')
                nodes_added.add(target)

            chem_g.add_edge(source, target,
                            attr_dict={
                                'interactionType': action,
                                'pubmed':          pubmed,
                            })
        print(len(chem_g.nodes()))
        print(len(chem_g.edges()))
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

        protein_cols = ['Official Symbol Interactor A',
                        'Official Symbol Interactor B',
                        'Pubmed ID', 'Score', 'Modification',
                        'Source Database',
                        'Experimental System Type']
        table = table[protein_cols]
        print(table['Modification'].unique())
        print(table.shape)
        table = table[~table['Modification'].isin(['-', 'No Modification'])]
        print(table.shape)
        print(table.head(10))
        table = table.as_matrix()
        g = self._create_chemical_network()
        added_genes = set()

        for r in table:
            gene1 = r[0]
            gene2 = r[1]
            inter = r[2]
            pubchem = r[3]
            score = r[4]
            mod = r[4]
            source_db = r[4]
            exp_system = r[4]

            if gene1 == gene2:
                continue
            else:
                if gene1 not in added_genes:
                    g.add_node(gene1, sourceDB=source_db,
                               speciesType='gene')
                    added_genes.add(gene1)
                if gene2 not in added_genes:
                    g.add_node(gene2, sourceDB=source_db,
                               speciesType='gene')
                    added_genes.add(gene2)

                g.add_edge(gene1, gene2, interactionType=inter,
                           pubchem=pubchem, expSystem=exp_system,
                           ptm=mod, score=score, sourceDB=source_db)

        ge = set(g.nodes())
        g1 = set(table[:, 0])
        g2 = set(table[:, 1])
        g_all = set()
        g_all.update(g1)
        g_all.update(g2)

        for i in g_all:
            if i not in ge:
                print(i)

        print("BIOGRID network has {} nodes ".format(len(g.nodes())))
        print("BIOGRID network has {} edges ".format(len(g.edges())))
        nx.write_gml(g, os.path.join(directory, 'biogrid.gml'))


if __name__ == '__main__':
    import time

    st = time.time()
    rfi = BioGridDownload()
    rfi.parse_network()
    end_time = time.time()
    print(end_time - st)
    # 33.48 seconds with set
