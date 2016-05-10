import os
import xml.etree.cElementTree as ET
import cPickle as pickle
"""

For this I downloaded all the xml files from hmdb. They are free. We just need to cite if we use this at all.


"""
#TODO Need to automatically download all HMDB to automate this part.
#TODO need to extract out gene-metabolite pairs to add to network information
hmdb_to_kegg = {}
kegg_to_hmdb = {}
kegg_hmdb_to_name = {}
for i in os.listdir('/home/pinojc/Downloads/HMDB'):
    if i.startswith('H'):
        tree = ET.parse('/home/pinojc/Downloads/HMDB/%s' %i)
        root = tree.getroot()
        kegg_id = None
        name = None
        hmdb_accession = None
        output = ''
        for entry in root.findall('kegg_id'):
            if entry.text is not None:
                kegg_id = entry.text
        for entry in root.findall('name'):
            if entry.text is not None:
                name = entry.text
        for entry in root.findall('accession'):
            if entry.text is not None:
                hmdb_accession = entry.text
        if (kegg_id and hmdb_accession) is not None:
            hmdb_to_kegg[hmdb_accession] = kegg_id
            if kegg_id in kegg_to_hmdb:
                kegg_to_hmdb[kegg_id].append(hmdb_accession)
            else:
                kegg_to_hmdb[kegg_id] = [hmdb_accession]
            kegg_hmdb_to_name[kegg_id] = name


pickle.dump(hmdb_to_kegg, open("../Mappings/hmdb_to_kegg_from_hmdb.p", "wb"))
pickle.dump(kegg_to_hmdb, open("../Mappings/kegg_to_hmdb_from_hmdb.p", "wb"))
pickle.dump(kegg_hmdb_to_name, open("../Mappings/kegg_to_chemical_name_from_hmdb.p", "wb"))