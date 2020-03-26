import logging
import math
import os
import zipfile

from defusedxml import cElementTree as ElementTree
import pandas as pd
import requests

from magine.data.storage import id_mapping_dir
from magine.logging import get_logger

logger = get_logger(__name__, log_level=logging.INFO)


def load_hgnc():
    hgnc_name = os.path.join(id_mapping_dir, 'hgnc.gz')
    if not os.path.exists(hgnc_name):
        hgnc = download_hgnc()
        _check_path(hgnc_name)
    else:
        hgnc = pd.read_csv(hgnc_name, low_memory=False)

    return hgnc.loc[hgnc.status == 'Approved']


def load_uniprot():
    # gather data from uniprot
    uniprot_path = os.path.join(id_mapping_dir, 'human_uniprot.csv.gz')
    if not os.path.exists(uniprot_path):
        uniprot = download_uniprot()
        _check_path(uniprot_path)
    else:
        uniprot = pd.read_csv(uniprot_path, low_memory=False)
    return uniprot


def load_ncbi():
    ncbi_name = os.path.join(id_mapping_dir, 'ncbi.gz')
    if not os.path.exists(ncbi_name):
        ncbi = download_ncbi()
        _check_path(ncbi_name)
    else:
        ncbi = pd.read_csv(ncbi_name, low_memory=False)
    return ncbi


def _check_path(path):
    if not os.path.exists(path):
        raise AssertionError()


def download_uniprot(species='hsa'):
    """
    `<https://www.uniprot.org/>`_

    Parameters
    ----------
    species : str
        Species name. Currently only human is supported.
        Please let us know if you need other species


    """
    _url_h = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/{}'

    if species == 'hsa':
        url = _url_h.format('HUMAN_9606_idmapping.dat.gz')
    else:
        raise ValueError("Currently only implemented human uniprot. "
                         "Please contact us for additional databases")

    columns = ['uniprot', 'mapping_type', 'mapping']

    df = _download(url, columns, name='uniprot', save=False, names=columns,
                   compression='gzip')

    uniprot = pd.pivot_table(df, columns='mapping_type', index='uniprot',
                             aggfunc='first')

    uniprot.columns = uniprot.columns.droplevel()
    uniprot.reset_index(inplace=True)
    outfile = os.path.join(id_mapping_dir, 'human_uniprot.csv.gz')
    uniprot.to_csv(outfile, compression='gzip',
                   header=True, index=False)

    return uniprot


def download_hgnc():
    """
    Downloads HGNC and stores it as a pandas.DataFrame
    `<http://www.genenames.org/>`_

    """

    url = 'ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_groups/protein-coding_gene.txt'

    columns = ['symbol', 'uniprot_ids', 'ensembl_gene_id', 'name', 'location',
               'entrez_id', 'ucsc_id', 'vega_id', 'alias_name', 'alias_symbol',
               'status', 'gene_family', 'gene_family_id', 'ena', 'iuphar',
               'cd',
               'refseq_accession', 'ccds_id', 'pubmed_id', 'mgd_id', 'rgd_id',
               'lsdb', 'bioparadigms_slc', 'enzyme_id', 'merops', 'horde_id',
               'pseudogene.org', 'cosmic', 'rna_central_ids', 'omim_id',
               'imgt',
               'intermediate_filament_db',
               ]
    return _download(url, columns, 'hgnc')


def download_ncbi():
    """
    Downloads data from NCBI

    `<https://www.ncbi.nlm.nih.gov/>`_

    """
    url = 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz'

    columns = ['GeneID', 'Symbol', 'description']
    return _download(url, columns, 'ncbi', compression='gzip')


def _download(url, columns, name, save=True, names=None, compression=None):
    logger.info("Downloading from {}.".format(name.upper()))
    if names is None:
        df = pd.read_table(url, delimiter='\t', low_memory=False, verbose=True,
                           compression=compression)
    else:
        df = pd.read_table(url, delimiter='\t', low_memory=False, verbose=True,
                           names=names, compression=compression)
    df = df[columns]
    if save:
        outfile = os.path.join(id_mapping_dir, '{}.gz'.format(name))
        df.to_csv(outfile, compression='gzip', header=True, index=False)
    logger.info("Done downloading from {}.".format(name.upper()))
    return df


def download_hmdb():
    """
    Downloads data from HMDB

    `<http://www.hmdb.ca/>`_

    """
    HMDB().download_db(fresh_download=True)


class HMDB(object):
    """
    Downloads and processes HMDB metabolites database

    `<http://www.hmdb.ca/>`_

    """

    def __init__(self):
        self.id_dir = id_mapping_dir
        self.target_file = 'hmdb_metabolites.zip'
        self.out_name = os.path.join(self.id_dir, 'hmdb_dataframe.csv.gz')
        self.hmdb_file = os.path.join(self.id_dir, self.target_file)

    def load_db(self, fresh_download=False):
        if not os.path.exists(self.out_name) or fresh_download:
            self.download_db(fresh_download)
        return pd.read_csv(self.out_name, low_memory=False, encoding='utf-8')

    def download_db(self, fresh_download):
        """ parse HMDB to Pandas.DataFrame

        """
        out_dir = os.path.join(self.id_dir, 'HMDB')
        # create output directory
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        # download and unzip xml file
        if fresh_download or not len(os.listdir(out_dir)):
            self._download_hmdb()
            # unzips hmdb metabolites file
            logger.info("Unzipping metabolites file")
            zip_ref = zipfile.ZipFile(self.hmdb_file, 'r')
            zip_ref.extractall(out_dir)
            zip_ref.close()
            logger.info("Done unzipping metabolites file")

        # creates dataframe
        logger.info("Parsing metabolites information from files")
        df = pd.DataFrame([
            self._create_dict(e) for ev, e in
            iter(ElementTree.iterparse(
                os.path.join(out_dir, os.listdir(out_dir)[0]),
                events=("start", "end"))
            )
            if e.tag == '{http://www.hmdb.ca}metabolite' and ev == 'end'
        ])

        for i in categories:
            if i in df.columns:
                df[i.split('}')[1]] = df[i]
                del df[i]

        df.to_csv(self.out_name, index=False, encoding='utf-8',
                  compression='gzip')
        logger.info("Done processing HMDB")

    def _download_hmdb(self):
        """ Downloads hmdb metabolites xml file """
        logger.info("Downloading metabolites information from HMDB")
        ur = 'http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip'

        r = requests.get(ur, stream=True)
        file_size = int(r.headers['content-length'])
        logger.info("File size is {} bytes".format(file_size))
        file_size_dl = 0
        block_sz = 1024
        # block_sz = 8024
        v = set()
        milestone_markers = range(0, 101, 10)

        with open(self.hmdb_file, 'wb') as f:
            for chunk in r.iter_content(chunk_size=block_sz):
                file_size_dl += len(chunk)
                percent_done = int(math.floor(file_size_dl * 100. / file_size))
                if percent_done in milestone_markers:
                    if percent_done not in v:
                        logger.info("{}%".format(percent_done))
                        v.add(percent_done)
                if chunk:  # filter out keep-alive new chunks
                    f.write(chunk)

        logger.info("Done downloading {}.".format(ur))

    @staticmethod
    def _create_dict(elem):
        template = {}
        for i in categories:
            n = elem.find(i)
            if n is None:
                # print(i)
                continue
            elif i == '{http://www.hmdb.ca}protein_associations':
                output = []
                for pr in n.findall('{http://www.hmdb.ca}protein'):
                    for gn in pr.findall('{http://www.hmdb.ca}gene_name'):
                        if gn.text is not None:
                            output.append(gn.text)
                template[i] = '|'.join(output)

            elif i == '{http://www.hmdb.ca}synonyms':
                template[i] = '|'.join(
                    [p.text for p in n.findall('{http://www.hmdb.ca}synonym')]
                )

            elif i == '{http://www.hmdb.ca}secondary_accessions':
                accession = n.findall('{http://www.hmdb.ca}accession')
                if len(accession) == 0:
                    accession = ''
                else:
                    accession = '|'.join(sorted([a.text for a in accession]))
                template[i] = accession

            else:
                template[i] = n.text
        elem.clear()
        return template


categories = ['kegg_id', 'name', 'accession', 'chebi_id',
              'chemspider_id', 'biocyc_id', 'synonyms',
              'pubchem_compound_id', 'protein_associations',
              'inchikey', 'iupac_name', 'drugbank_id',
              'chemical_formula', 'smiles', 'metlin_id',
              'average_molecular_weight', 'secondary_accessions',
              # 'normal_concentrations',
              # 'molecular_framework'
              #  'pathways',
              ]

categories = ['{http://www.hmdb.ca}' + i for i in categories]

valid_uniprot_cols = ['uniprot', 'Allergome', 'BioCyc', 'BioGrid', 'BioMuta',
                      'CCDS', 'CRC64', 'ChEMBL', 'ChiTaRS', 'CleanEx', 'DIP',
                      'DMDM', 'DNASU', 'DisProt', 'DrugBank', 'EMBL',
                      'EMBL-CDS', 'ESTHER', 'Ensembl', 'Ensembl_PRO',
                      'Ensembl_TRS', 'GI', 'GeneCards','GeneDB', 'GeneID',
                      'GeneReviews', 'GeneTree', 'GeneWiki', 'Gene_Name',
                      'Gene_ORFName', 'Gene_Synonym', 'GenomeRNAi',
                      'GuidetoPHARMACOLOGY', 'H-InvDB', 'HGNC', 'HOGENOM',
                      'HOVERGEN', 'HPA', 'KEGG', 'KO', 'MEROPS', 'MIM', 'MINT',
                      'NCBI_TaxID', 'OMA', 'Orphanet', 'OrthoDB', 'PATRIC',
                      'PDB', 'PeroxiBase', 'PharmGKB', 'REBASE', 'Reactome',
                      'RefSeq', 'RefSeq_NT', 'STRING', 'SwissLipids', 'TCDB',
                      'TreeFam', 'UCSC', 'UniGene', 'UniParc', 'UniPathway',
                      'UniProtKB-ID', 'UniRef100', 'UniRef50', 'UniRef90',
                      'eggNOG', 'neXtProt']

ncbi_valid_categories = ['GeneID', 'Symbol', 'description']

hgnc_valid_categories = ['symbol', 'uniprot_ids', 'ensembl_gene_id', 'name',
                         'location', 'entrez_id', 'ucsc_id', 'vega_id',
                         'alias_name', 'alias_symbol', 'status', 'gene_family',
                         'gene_family_id', 'ena', 'iuphar', 'cd',
                         'refseq_accession', 'ccds_id', 'pubmed_id', 'mgd_id',
                         'rgd_id', 'lsdb', 'bioparadigms_slc', 'enzyme_id',
                         'merops', 'horde_id', 'pseudogene.org', 'cosmic',
                         'rna_central_ids', 'omim_id', 'imgt',
                         'intermediate_filament_db'
                         ]

if __name__ == '__main__':
    # download_hgnc()
    download_uniprot('hsa')
    # hmdb = HMDB()
    # df = hmdb.load_db(fresh_download=True)
    # print(df.head(10))
    # print(df.columns)
