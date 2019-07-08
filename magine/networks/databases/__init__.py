"""
Database downloads.

MAGINE downloads network information from
* Reactome Functional Interactions
* HMDB
* BioGrid
* KEGG
* Signor
With the exception of KEGG and HMDB, all databases are downloaded and processed
with pandas. KEGG is downloaded using Bioservices. HMDB is download in xml
format and processed using lxml parser.

"""
from .biogrid_interactions import load_biogrid_network, download_biogrid
from .hmdb import load_hmdb_network
from .kegg_kgml import load_all_of_kegg, load_kegg_mappings, download_kegg
from .reactome_functional_interaction import load_reactome_fi, \
    download_reactome_fi
from .signor import load_signor, download_signor

__all__ = ['load_reactome_fi', 'download_reactome_fi',
           'load_biogrid_network', 'download_biogrid',
           'load_signor', 'download_signor',
           'load_kegg_mappings', 'load_all_of_kegg',
           'load_hmdb_network']

if __name__ == '__main__':
    download_reactome_fi()
    download_signor()
    download_biogrid()
    download_kegg()
