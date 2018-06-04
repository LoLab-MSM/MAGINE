"""Interface to downloading ID mapping databases.

.. topic:: Databases supported

    :URL: https://www.uniprot.org/
    :URL: https://www.ncbi.nlm.nih.gov/
    :URL: http://www.hmdb.ca/

"""
from .download_libraries import load_hgnc, load_uniprot, \
    load_ncbi, HMDB

__all__ = ['load_hgnc', 'load_ncbi', 'load_uniprot', 'HMDB']
