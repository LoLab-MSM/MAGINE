#!python

import magine.mappings.databases.download_libraries as dl
import magine.networks.databases as nd


def download_id_mapping():
    dl.download_hgnc()
    dl.download_ncbi()
    dl.download_uniprot()


def download_network_dbs():
    nd.download_reactome_functional_interaction()
    nd.download_all_of_kegg()
    dl.HMDB()


if __name__ == '__main__':
    import time

    st = time.time()
    download_id_mapping()
    download_network_dbs()
    et = time.time()
    print(et - st)
