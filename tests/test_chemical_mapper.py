import pygraphviz as pyg

from Mappings.chemical_mapper import ChemicalMapper

cm = ChemicalMapper()
cm.print_info()


def test_synonyms():
    # term_2 = '1-(9Z,12Z,15Z-Octadeatrienoyl)-2-(8Z,11Z,14Z-eicosatrienoyl)-3-(4Z,7Z,10Z,13Z,16Z,19Z-docosahexaenoyl)-glycerol'
    hmdb = cm.check_synonym_dict(term='dodecene', format_name='accession')
    # hmdb = cm.check_synonym_dict(term=term_2, format_name='accession')
    print('here', hmdb)
    assert (hmdb == 'HMDB59874')


def test_protein_network():
    g = pyg.AGraph(directed=True)
    # item = 'HMDB02865'
    item = 'HMDB42489'
    name = cm.hmdb_accession_to_chemical_name[item][0]

    for protein in cm.hmdb_accession_to_protein[item]:
        if protein is None:
            continue
        else:
            g.add_edge(name, protein)
    g.draw('test_oxytocin.png', prog='dot')


if __name__ == "__main__":
    test_synonyms()
    test_protein_network()
