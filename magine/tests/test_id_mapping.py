import networkx as nx
from nose.tools import ok_

from magine.mappings.chemical_mapper import ChemicalMapper
from magine.mappings.gene_mapper import GeneMapper

cm = ChemicalMapper()
gm = GeneMapper()


class TestChemicalMapper(object):

    def test_synonyms(self):
        hmdb = cm.check_synonym_dict(term='dodecene',
                                     format_name='main_accession')

        ok_((hmdb == ['HMDB0000933', 'HMDB0059874']))

    def test_protein_network(self):
        item = 'HMDB42489'
        hit_list = ['PNLIP', 'LIPC', 'LIPA', 'PNLIPRP1', 'PNPLA3', 'LIPF',
                    'LIPG', 'CEL', 'DGAT1', 'PNLIPRP2', 'CPT1B', 'CPT1A',
                    'LPL', 'CPT2', 'MGLL', 'CES1', 'LIPE', 'MTTP', 'APOA1',
                    'CETP', 'APOE', 'APOC3', 'APOB', 'APOA4', 'CD36', 'P4HB',
                    'MOGAT2', 'PNPLA4', 'SLC27A1', 'DGAT2', 'MOGAT1',
                    'MOGAT3', 'PLB1', 'APOA5']

        for protein in cm.hmdb_to_protein[item]:
            ok_(protein in hit_list)


class TestGeneMapper(object):
    def test_kegg_to_gene_name(self):
        ok_(gm.kegg_to_gene_name['hsa:581'][0] == 'BAX')

    def test_kegg_to_uniprot(self):
        ok_(gm.kegg_to_uniprot['hsa:581'][0] == 'Q07812')

    def test_gene_name_to_ensemble(self):
        ok_(gm.gene_name_to_ensembl['BAX'][0] == 'ENSG00000087088')

    def test_gene_name_to_uniprot(self):
        ok_(gm.gene_name_to_uniprot['BAX'][0] == 'Q07812')

    def test_uniprot_to_gene_name(self):
        ok_(gm.uniprot_to_gene_name['Q07812'][0] == 'BAX')


def test_kegg_to_hmdb():
    """
    tests kegg compound to hmdb
    :return:
    """
    g = nx.DiGraph()
    g.add_edge('hsa:224', 'hsa:219')
    g.add_edge('hsa:219', 'cpd:C00197')
    g.add_edge('cpd:C00197', 'cpd:C00197')
    g.add_edge('cpd:C15972', 'cpd:C00197')
    g.add_edge('cpd:C15972', 'cpd:C00469')

    change_dict, kegg_short, chem_names = cm.convert_kegg_nodes(g)
    nx.set_node_attributes(g, kegg_short, 'keggName')
    nx.set_node_attributes(g, chem_names, 'chemName')
    nx.set_node_attributes(g, change_dict, 'hmdbNames')

    dict2, kegg_short = gm.convert_kegg_nodes(g)
    nx.set_node_attributes(g, kegg_short, 'keggName')
    change_dict.update(dict2)
    g = nx.relabel_nodes(g, change_dict)
    answer = '(2R)-2-Hydroxy-3-(phosphonatooxy)propanoate'
    ok_(g.node['HMDB0060180']['chemName'] == answer)
    ok_(g.node['HMDB0060180']['keggName'] == 'C00197')


def test_kegg_to_uniprot():
    """
    tests kegg gene to gene name
    :return:
    """
    g = nx.DiGraph()
    g.add_edge('hsa:217', 'hsa:219')
    g.add_edge('hsa:217', 'hsa:223')
    g.add_edge('hsa:501', 'hsa:219')
    g.add_edge('hsa:224', 'hsa:219')
    g.add_node('hsa:857')
    dic, kegg_short = gm.convert_kegg_nodes(g, species='hsa')
    nx.set_node_attributes(g, kegg_short, 'keggName')
    g = nx.relabel_nodes(g, dic)

    valid = {
        'ALDH2': '217',
        'ALDH1B1': '219',
        'ALDH9A1': '223',
        'ALDH7A1': '501',
        'ALDH3A2': '224',
        'CAV1': '857'
    }
    for k, v in valid.items():
        ok_(g.node[k]['keggName'] == v)


if __name__ == '__main__':
    test_kegg_to_hmdb()
