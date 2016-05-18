import os
import xml.etree.cElementTree as ET
import cPickle as pickle
import networkx as nx
import pandas as pd
import numpy as np

"""

For this I downloaded all the xml files from hmdb. They are free. We just need to cite if we use this at all.


"""
#TODO Need to automatically download all HMDB to automate this part.
#TODO need to extract out gene-metabolite pairs to add to network information
hmdb_to_kegg = {}
kegg_to_hmdb = {}
kegg_hmdb_to_name = {}
hmdb_and_related_proteins = {}
all_proteins = []


hmdb_all = {}
tmp_all = []
count = 0
for i in os.listdir('/home/pinojc/Downloads/HMDB'):

    if i.startswith('H'):
        # if count > 10:
        #     continue
        # count+=1

        tree = ET.parse('/home/pinojc/Downloads/HMDB/%s' %i)
        root = tree.getroot()

        kegg_id = None
        name = None
        hmdb_accession = None
        output = ''
        protein_list = []
        hmdb_entry = {}
        tmp = []
        tmp = np.zeros(4,dtype='S120')
        for entry in root.findall('protein_associations'):
            for prot in entry.findall('protein'):
                protein_list.append(prot.find('gene_name').text)

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

        if None in protein_list:
            protein_list.remove(None)

        protein = ''
        for i in protein_list:
            protein+='%s,' % str(i)

        name = name.encode('utf-8')
        tmp[0] = hmdb_accession
        tmp[1] = name
        tmp[2] = kegg_id
        tmp[3] = protein
        tmp_all.append(tmp)
        if len(protein_list) == 0:
            continue
        else:
            hmdb_and_related_proteins[hmdb_accession] = protein_list

dictionary_dataframe = pd.DataFrame(tmp_all,columns=['accession','name','kegg_id','proteins',])
dictionary_dataframe.to_pickle('hmdb_dataframe.p')
print(dictionary_dataframe.dtypes)
print(dictionary_dataframe.name)
quit()
pickle.dump(kegg_hmdb_to_name, open("../Mappings/hmdb_related_genes.p", "wb"))
quit()
g = nx.DiGraph()
for i in hmdb_and_related_proteins:
    for j in hmdb_and_related_proteins[i]:
        if j is None:
            continue
        g.add_edge(i,j)
nx.write_gml(g,'hmdb.gml')
quit()
nx.draw(g,with_labels=True)
import matplotlib.pyplot as plt
plt.savefig('hmdb_interactions.png')
#pickle.dump(hmdb_to_kegg, open("../Mappings/hmdb_to_kegg_from_hmdb.p", "wb"))
#pickle.dump(kegg_to_hmdb, open("../Mappings/kegg_to_hmdb_from_hmdb.p", "wb"))
#pickle.dump(kegg_hmdb_to_name, open("../Mappings/kegg_to_chemical_name_from_hmdb.p", "wb"))