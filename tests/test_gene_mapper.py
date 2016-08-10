from Mappings.gene_mapper import GeneMapper

gm = GeneMapper('hsa')
print(gm.gene_name_to_uniprot['BAX'])
print(gm.gene_name_to_alias_name['BAX'])
print(gm.gene_name_to_ensembl['BAX'])

print(gm.uniprot_to_gene_name)


def test_gene_name_to_ensemble():
    assert (gm.gene_name_to_ensembl['BAX'][0] == 'ENSG00000087088')


def test_gene_name_to_uniprot():
    assert (gm.gene_name_to_uniprot['BAX'][0] == 'Q07812')


def uniprot_to_gene_name():
    assert (gm.uniprot_to_gene_name['Q07812'][0] == 'BAX')


if __name__ == '__main__':
    uniprot_to_gene_name()
    test_gene_name_to_uniprot()
    test_gene_name_to_ensemble()
