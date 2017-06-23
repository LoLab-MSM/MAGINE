"""Interface to downloading ID mapping databases.

.. topic:: Databases supported

    :URL: https://www.uniprot.org/
    :URL: https://www.ncbi.nlm.nih.gov/
    :URL: http://www.hmdb.ca/

"""
from .download_libraries import download_uniprot, download_hgnc, \
    download_ncbi, HMDB


__all__ = ['download_uniprot', 'download_ncbi', 'download_hgnc', 'HMDB']
