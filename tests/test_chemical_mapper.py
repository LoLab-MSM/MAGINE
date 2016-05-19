from Mappings.maps import ChemicalMapper

cm = ChemicalMapper()


def test_synonyms():
    hmdb = cm.check_synonym_dict(term='dodecene', format_name='accession')
    assert (hmdb == 'HMDB59874')


test_synonyms()
