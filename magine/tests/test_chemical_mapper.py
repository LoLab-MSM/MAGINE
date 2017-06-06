
from magine.mappings.chemical_mapper import ChemicalMapper

cm = ChemicalMapper()


def test_init():
    cm.load()
    cm.reload()


"""

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

if __name__ == "__main__":
    test_init()
    # test_synonyms()
    # test_protein_network()

"""
