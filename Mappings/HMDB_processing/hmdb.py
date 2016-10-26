import czipfile as zipfile
import os

import pandas as pd
import urllib2

from Mappings.xml_to_dictionary import ElementTree, XmlDictConfig, \
    XmlListConfig

directory = os.path.dirname(__file__)


def download_hmdb():
    hmdb_db_url = 'http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip'
    file_name = hmdb_db_url.split('/')[-1]
    u = urllib2.urlopen(hmdb_db_url)
    out_path = os.path.join(directory, file_name)
    f = open(out_path, 'wb')
    meta = u.info()
    file_size = int(meta.getheaders("Content-Length")[0])
    print("Downloading: %s Bytes: %s" % (file_name, file_size))
    file_size_dl = 0
    block_sz = 8192
    while True:
        buffer = u.read(block_sz)
        if not buffer:
            break
        file_size_dl += len(buffer)
        f.write(buffer)
        status = r"%10d  [%3.2f%%]" % (
        file_size_dl, file_size_dl * 100. / file_size)
        status = status + chr(8) * (len(status) + 1)
        print(status)
    f.close()


def parse_hmdb():
    out_dir = os.path.join(directory, 'HMDB')

    def unzip_hmdb():

        hmdb_file = os.path.join(directory, 'hmdb_metabolites.zip')
        if not os.path.exists(hmdb_file):
            download_hmdb()
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        zip_ref = zipfile.ZipFile(hmdb_file, 'r')
        zip_ref.extractall(out_dir)
        zip_ref.close()

    count = 0
    categories = ['kegg_id', 'name', 'accession', 'chebi_id', 'chemspider_id',
                  'biocyc_id', 'synonyms', 'pubchem_compound_id',
                  'protein_associations',
                  # 'secondary_accessions',
                  # 'iupac_name',
                  # 'normal_concentrations','chemical_formula', 'smiles',
                  # 'drugbank_id', 'average_molecular_weight',
                  # 'pathways', 'metlin_id',
                  # 'inchikey',
                  ]

    tmp_all = []
    for i in os.listdir(out_dir):
        if i.startswith('H'):
            # if count > 100:
            #     continue
            count += 1
            filename = os.path.join(out_dir, i)
            tree = ElementTree.parse(filename)
            root = tree.getroot()
            xmldict = XmlDictConfig(root)
            tmp = []
            for cat in categories:
                info = xmldict[cat]
                if info is None:
                    tmp.append(info)
                elif cat == 'protein_associations':
                    tmp_list = []
                    if type(info) == str:
                        continue
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
                        continue
                    else:
                        if type(info['synonym']) == str:
                            tmp_list.append(info['synonym'])
                        else:
                            for each in info['synonym']:
                                tmp_list.append(each)
                    tmp.append(tmp_list)
                else:
                    tmp.append(info.encode('ascii', 'ignore'))
            tmp_all.append(tmp)
    df = pd.DataFrame(tmp_all, columns=categories)
    df.to_csv(os.path.join(directory, 'hmdb_dataframe.csv.gz'),
              compression='gzip')


parse_hmdb()
