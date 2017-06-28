import magine.mappings.databases.download_libraries as dl
from magine.networks.databases import download_all_of_kegg, \
    download_reactome_functional_interaction
from magine.ontology.databases.gene_ontology import create_dicts_through_orange


def download_id_mapping():
    dl.download_hgnc()
    dl.download_ncbi()
    dl.download_uniprot()


def download_network_dbs():
    download_reactome_functional_interaction()
    download_all_of_kegg()
    dl.HMDB()


def download_go():
    create_dicts_through_orange()


if __name__ == '__main__':
    import time

    st = time.time()
    download_id_mapping()
    download_network_dbs()
    download_go()
    et = time.time()
    print(et - st)
