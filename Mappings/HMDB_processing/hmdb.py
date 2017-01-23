"""
HMDB database download and processing
"""

import os
import tempfile
import xml.etree.cElementTree as ElementTree
import zipfile

import numpy as np
import pandas as pd
import requests

from Mappings.xml_to_dictionary import XmlDictConfig, XmlListConfig

directory = os.path.dirname(__file__)

categories = ['kegg_id', 'name', 'accession', 'chebi_id', 'chemspider_id',
              'biocyc_id', 'synonyms', 'pubchem_compound_id',
              'protein_associations', 'inchikey', 'iupac_name',
              # 'secondary_accessions',

              # 'normal_concentrations','chemical_formula', 'smiles',
              # 'drugbank_id', 'average_molecular_weight',
              # 'pathways', 'metlin_id',
              ]

class HMDB:
    """ Downloads and processes HMDB metabolites database

    """

    def __init__(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.out_dir = directory
        self.target_file = 'hmdb_metabolites.zip'

    def download_hmdb(self):
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
        block_sz = 8192
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

    def unzip_hmdb(self, out_directory):
        """ unzips hmdb metabolites file

        :return:
        """
        hmdb_file = os.path.join(self.tmp_dir, self.target_file)
        if not os.path.exists(hmdb_file):
            print("Downloading metabolites information from HMDB")
            self.download_hmdb()
        print("Unzipping metabolites file")
        zip_ref = zipfile.ZipFile(hmdb_file, 'r')
        zip_ref.extractall(out_directory)
        zip_ref.close()
        print("Done unzipping metabolites file")

    def parse_hmdb(self):
        """ parse HMDB to Pandas.DataFrame

        """
        out_dir = os.path.join(self.tmp_dir, 'HMDB')
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
            self.unzip_hmdb(out_dir)

        count = 0

        print("Parsing metabolites information from files")
        tmp_all = []
        for i in os.listdir(out_dir):
            if i.startswith('H'):
                # if count > 100:
                #     continue
                count += 1
                filename = os.path.join(out_dir, i)
                tree = ElementTree.parse(filename)
                xmldict = XmlDictConfig(tree.getroot())
                tmp = []
                for cat in categories:
                    info = xmldict[cat]
                    if cat == 'protein_associations':
                        tmp_list = []
                        if type(info) == str:

                            tmp.append([])
                        else:
                            if type(info['protein']) == XmlDictConfig:
                                tmp_list.append(info['protein']['gene_name'])

                            elif type(info['protein']) == XmlListConfig:
                                for ii in info['protein']:
                                    tmp_list.append(ii['gene_name'])
                            tmp.append(tmp_list)
                    elif cat == 'synonyms':
                        tmp_list = []
                        if type(info) == str:
                            tmp.append([])
                            continue
                        else:
                            if type(info['synonym']) == str:
                                tmp_list.append(info['synonym'])
                            else:
                                for each in info['synonym']:
                                    tmp_list.append(each)
                        tmp.append(tmp_list)
                    else:
                        if info is None:
                            tmp.append([])
                        else:
                            tmp.append(info.encode('ascii', 'ignore'))
                tmp_all.append(tmp)
        df = pd.DataFrame(tmp_all, columns=categories)
        df.to_csv(os.path.join(self.out_dir, 'hmdb_dataframe.csv.gz'),
                  compression='gzip', index=False)
        print("Done processing HMDB")

    def load_db(self):
        df = pd.read_csv(os.path.join(self.out_dir, 'hmdb_dataframe.csv.gz'))
        return df
if __name__ == '__main__':
    hm = HMDB()
    hm.download_hmdb()
    hm.parse_hmdb()
    df = hm.load_db()
    print(df.head(10))
