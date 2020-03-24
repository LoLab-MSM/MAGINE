from magine.networks.databases.biogrid_interactions import download_biogrid
from magine.networks.databases.hmdb import load_hmdb_network
from magine.networks.databases.kegg_kgml import download_kegg
from magine.networks.databases.reactome_functional_interaction import \
    download_reactome_fi
from magine.networks.databases.signor import download_signor

if __name__ == '__main__':
    download_reactome_fi()
    download_signor()
    download_biogrid()
    download_kegg()
    load_hmdb_network(fresh_download=True)
