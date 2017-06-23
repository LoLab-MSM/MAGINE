import os
import tempfile
import time
import zipfile
from xml.etree import cElementTree as ElementTree

import numpy as np
import pandas as pd
import requests

dir_name = os.path.dirname(__file__)


def download_uniprot_db(species='hsa'):
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
        ValueError("Currently only implemented human uniprot. "
                   "Please contain us for additional databases")

    headers = ['uniprot', 'mapping_type', 'mapping']

    print("First time initialzing gene_mapper! \n"
          "Downloading from Uniprot. Might take awhile.")

    human = pd.read_table(url, delimiter='\t', names=headers,
                          compression='gzip', verbose=True)

    outfile = os.path.join(dir_name, '..', 'data', 'human_uniprot.csv.gz')
    human.to_csv(outfile, compression='gzip', columns=headers, header=True)
    return human


def download_hgnc():
    """
    Downloads HGNC and stores it as a pandas.DataFrame
    `<http://www.genenames.org/>`_

    """

    url = 'ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_groups/protein-coding_gene.txt'
    print("First time initializing gene_mapper! \n"
          "Downloading from HGNC. Might take awhile.")
    hgnc = pd.read_table(url, delimiter='\t', low_memory=False)
    hgnc = hgnc[hgnc['status'] == 'Approved']
    hgnc = hgnc[['symbol', 'uniprot_ids', 'ensembl_gene_id', 'name', 'location',
                 'entrez_id', 'ucsc_id', 'mirbase', 'vega_id', 'alias_name',
                 'alias_symbol']]
    outfile = os.path.join(dir_name, '..', 'data', 'hgnc.gz')
    hgnc.to_csv(outfile, compression='gzip', header=True)
    hgnc = pd.read_csv(outfile)
    return hgnc


def download_ncbi():
    """
    Downloads data from NCBI

    `<https://www.ncbi.nlm.nih.gov/>`_

    """
    url = 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz'
    ncbi = pd.read_table(url, delimiter='\t', low_memory=False,
                         compression='gzip')
    cols = ['GeneID', 'Symbol', 'description']
    ncbi = ncbi[cols]
    outfile = os.path.join(dir_name, '..', 'data', 'ncbi.gz')
    ncbi.to_csv(outfile, compression='gzip', header=True)


class HMDB(object):
    """
    Downloads and processes HMDB metabolites database

    `<http://www.hmdb.ca/>`_

    """

    def __init__(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.out_dir = dir_name
        self.target_file = 'hmdb_metabolites.zip'
        self.out_name = os.path.join(self.out_dir, '..', 'data',
                                     'hmdb_dataframe.csv')

    def _download_hmdb(self):
        """ Downloads hmdb metabolites xml file

        :return:
        """

        hmdb_db_url = 'http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip'

        out_path = os.path.join(self.tmp_dir, self.target_file)

        r = requests.get(hmdb_db_url, stream=True)
        response = requests.head(hmdb_db_url)
        file_size = int(response.headers['content-length'])
        print("Downloading: %s Bytes: %s" % (self.target_file, file_size))
        file_size_dl = 0
        block_sz = 1024
        # block_sz = 8024
        v = set()
        milestone_markers = range(0, 101, 10)

        with open(out_path, 'wb') as f:
            for chunk in r.iter_content(chunk_size=block_sz):
                file_size_dl += len(chunk)
                percent_download = int(
                        np.floor(file_size_dl * 100. / file_size))

                if percent_download in milestone_markers:
                    if percent_download not in v:
                        print("{}%".format(percent_download))
                        v.add(percent_download)
                if chunk:  # filter out keep-alive new chunks
                    f.write(chunk)

        print("Downloaded {} and stored {}".format(hmdb_db_url, out_path))
        return

    def _unzip_hmdb(self, out_directory):
        """ unzips hmdb metabolites file

        :return:
        """
        hmdb_file = os.path.join(self.tmp_dir, self.target_file)

        if not os.path.exists(hmdb_file):
            print("Downloading metabolites information from HMDB")
            self._download_hmdb()
        print("Unzipping metabolites file")
        zip_ref = zipfile.ZipFile(hmdb_file, 'r')
        zip_ref.extractall(out_directory)
        zip_ref.close()
        print("Done unzipping metabolites file")

    def setup(self):
        """ parse HMDB to Pandas.DataFrame

        """
        # out_dir = b'c:\users\jamesp~1\\appdata\local\\temp\\tmpnvpdle\HMDB'
        out_dir = os.path.join(self.tmp_dir, 'HMDB')
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
            self._unzip_hmdb(out_dir)

        count = 0

        print("Parsing metabolites information from files")
        tmp_all = [None] * 1000000
        # tmp_all = []
        for i in os.listdir(out_dir):
            print(i)
            filename = os.path.join(out_dir, i)

            # get an iterable
            context = ElementTree.iterparse(filename, events=("start", "end"))

            # turn it into an iterator
            context = iter(context)

            # get the root element
            event, root = next(context)
            for event, elem in context:

                if '}' in elem.tag:
                    elem.tag = elem.tag.split('}', 1)[1]
                if elem.tag == 'metabolite' and event == 'end':
                    # """

                    # if count > 10000:
                    #     break
                    # """
                    template = self._create_dict(elem)

                    tmp_all[count] = template

                    count += 1
                    elem.clear()
                    root.clear()

        df = pd.DataFrame(tmp_all[:count], columns=template.keys())
        df.to_csv(self.out_name, index=False, encoding='utf-8',
                  tupleize_cols=True)
        print("Done processing HMDB")

    def load_db(self):
        df = pd.read_csv(self.out_name, low_memory=False)
        return df

    @staticmethod
    def _create_dict(elem):
        template = {}
        for i in categories:
            n = elem.find(i)
            if i == 'protein_associations':
                output = []
                for pr in n.findall('protein'):
                    for gn in pr.findall('gene_name'):
                        gene_name = gn.text
                        if gene_name is not None:
                            output.append(gene_name)
                template[i] = '|'.join(i for i in output)
            elif i == 'synonyms':
                output = []
                for pr in n.findall('synonym'):
                    output.append(pr.text)
                template[i] = '|'.join(i for i in output)
            elif i == 'ontology':
                bf_all = []
                cl_all = []
                for pr in n.findall('biofunctions'):
                    for bf in pr.findall('biofunction'):
                        bf_text = bf.text
                        if bf_text is not None:
                            bf_all.append(bf_text)
                template['biofunction'] = '|'.join(i for i in bf_all)
                for pr in n.findall('cellular_locations'):
                    for cl in pr.findall('cellular_location'):
                        cl_text = cl.text
                        if cl_text is not None:
                            cl_all.append(cl.text)
                template['cellular_locations'] = '|'.join(i for i in cl_all)

            else:
                output = n.text
                if output is not None:
                    output = output
                    # output = output.encode('utf-8')
                template[i] = output
        return template


if __name__ == '__main__':
    hm = HMDB()
    end = time.time()
    st = time.time()
    human = download_uniprot_db()
    directory = os.path.dirname(__file__)
    categories = ['kegg_id', 'name', 'accession', 'chebi_id', 'chemspider_id',
              'biocyc_id', 'synonyms', 'pubchem_compound_id',
              'protein_associations', 'inchikey', 'iupac_name',
              'ontology',  # 'cellular_location', 'biofunction'
              # 'cellular_location', 'biofunction', 'molecular_framework'
              # 'secondary_accessions',

              # 'normal_concentrations','chemical_formula', 'smiles',
              # 'drugbank_id', 'average_molecular_weight',
              # 'pathways', 'metlin_id',
              ]
