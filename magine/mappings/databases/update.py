from magine.mappings.databases.download_libraries \
    import download_hmdb, download_hgnc, download_ncbi, download_uniprot

if __name__ == '__main__':
    download_hgnc()
    download_ncbi()
    download_uniprot()
    download_hmdb()
