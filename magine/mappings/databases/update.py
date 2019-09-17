from .download_libraries import HMDB, download_hgnc, download_ncbi, \
    download_uniprot

if __name__ == '__main__':
    download_hgnc()
    download_ncbi()
    download_uniprot()
    HMDB().download_db()
