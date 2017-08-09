import networkx as nx
from magine.mappings.gene_mapper import GeneMapper
from magine.mappings.chemical_mapper import ChemicalMapper
from magine.networks.network_generator import expand_by_hmdb

gm = GeneMapper()
cm = ChemicalMapper()


def test_synonyms():
    hmdb = cm.check_synonym_dict(term='dodecene', format_name='accession')
    assert (hmdb == 'HMDB59874')


def test_protein_network():
    # item = 'HMDB02865'
    item = 'HMDB42489'
    hit_list = ['PNLIP', 'LIPC', 'LIPA', 'PNLIPRP1', 'PNPLA3', 'LIPF', 'LIPG',
                'CEL', 'DGAT1', 'PNLIPRP2', 'CPT1B', 'CPT1A', 'LPL', 'CPT2',
                'MGLL', 'CES1', 'LIPE', 'MTTP', 'APOA1', 'CETP', 'APOE',
                'APOC3', 'APOB', 'APOA4', 'CD36', 'P4HB', 'MOGAT2', 'PNPLA4',
                'SLC27A1', 'DGAT2', 'MOGAT1', 'MOGAT3', 'PLB1', 'APOA5']
    for protein in cm.hmdb_accession_to_protein[item]:
        assert protein in hit_list


def test_expand_by_hmdb():
    g = nx.DiGraph()
    g.add_edge('PNLIP', 'LIPC')
    new_g = expand_by_hmdb(graph=g,
                           metabolite_list=['HMDB42489'],
                           all_measured=['HMDB02865', 'HMDB59874']
                           )
    for i in new_g.nodes():
        print(i)


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
