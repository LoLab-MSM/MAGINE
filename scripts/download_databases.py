#!python
import logging

from magine.data.storage import clear_cached_dbs, create_storage_structure
from magine.logging import get_logger
import magine.mappings.databases.download_libraries as dl
import magine.networks.databases as nd

logger = get_logger(__name__, log_level=logging.INFO)


def download_id_mapping():
    dl.download_hgnc()
    dl.download_ncbi()
    dl.download_uniprot()


def download_network_dbs():
    dl.download_hmdb()
    nd.download_reactome_fi()
    nd.download_signor()
    nd.download_biogrid()


def run():
    import time

    clear_cached_dbs()
    create_storage_structure()
    logger.info("Downloading all network and ID mapping databases")
    st = time.time()
    download_id_mapping()
    download_network_dbs()
    et = time.time()
    logger.info("Took {} seconds".format(et - st))


if __name__ == '__main__':
    run()
