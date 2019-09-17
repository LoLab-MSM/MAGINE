from .biogrid_interactions import download_biogrid
from .hmdb import load_hmdb_network
from .kegg_kgml import download_kegg
from .reactome_functional_interaction import download_reactome_fi
from .signor import download_signor

if __name__ == '__main__':
    download_reactome_fi()
    download_signor()
    download_biogrid()
    download_kegg()
    load_hmdb_network(fresh_download=True)
