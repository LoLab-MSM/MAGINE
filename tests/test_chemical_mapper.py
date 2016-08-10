import pygraphviz as pyg

from Mappings.maps import cm

cm.print_info()

def test_synonyms():
    hmdb = cm.check_synonym_dict(term='dodecene', format_name='accession')
    assert (hmdb == 'HMDB59874')

def test_protein_network():
    g = pyg.AGraph(directed=True)
    item = 'HMDB02865'
    item = 'HMDB42489'
    name = cm.hmdb_accession_to_chemical_name[item][0]
    for protein in cm.hmdb_accession_to_protein[item][0]:
        if protein is None:
            continue
        else:
            g.add_edge(name, protein)
    g.draw('test_oxytocin.png', prog='dot')


test_synonyms()
test_protein_network()
