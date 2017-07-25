from magine.mappings.gene_mapper import GeneMapper

gm = GeneMapper('hsa')


def test_gene_name_to_ensemble():
    assert (gm.gene_name_to_ensembl['BAX'][0] == 'ENSG00000087088')


def test_gene_name_to_uniprot():
    assert (gm.gene_name_to_uniprot['BAX'][0] == 'Q07812')


def test_uniprot_to_gene_name():
    assert (gm.uniprot_to_gene_name['Q07812'][0] == 'BAX')


def test_kegg_to_gene_name():
    print(gm.kegg_to_uniprot['hsa:581'])
    assert (gm.kegg_to_uniprot['hsa:581'][0] == 'Q07812')
    assert (gm.kegg_to_gene_name['hsa:581'][0] == 'BAX')


if __name__ == '__main__':
    test_uniprot_to_gene_name()
    test_gene_name_to_uniprot()
    test_gene_name_to_ensemble()
    test_kegg_to_gene_name()
