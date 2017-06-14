import networkx as nx

from magine.mappings.chemical_mapper import ChemicalMapper
from magine.networks.network_generator import expand_by_hmdb

cm = ChemicalMapper()


# def test_init():
#     cm.load()
#     cm.reload()
#     cm.print_info()


# """

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


# """
def test_expand_by_hmdb():
    g = nx.DiGraph()
    g.add_edge('PNLIP', 'LIPC')
    new_g = expand_by_hmdb(graph=g,
                           metabolite_list=['HMDB42489'],
                           all_measured=['HMDB02865', 'HMDB59874']
                           )
    for i in new_g.nodes():
        print(i)
    # """
    # if __name__ == "__main__":
    #     test_init()
    # test_synonyms()
    # test_protein_network()

# """
