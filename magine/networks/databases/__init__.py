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
from .biogrid_interactions import download_biogrid, load_biogrid_network
from .hmdb import load_hmdb_network
from .kegg_kgml import download_kegg, load_kegg, load_kegg_mappings
from .reactome_functional_interaction import download_reactome_fi, \
    load_reactome_fi
from .signor import download_signor, load_signor

__all__ = ['load_reactome_fi', 'download_reactome_fi',
           'load_biogrid_network', 'download_biogrid',
           'load_signor', 'download_signor',
           'load_kegg_mappings', 'load_kegg',
           'load_hmdb_network']
