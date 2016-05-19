import numpy as np
import pandas as pd

headers = ["UniProtKB-AC", "UniProtKB-ID", "GeneID (EntrezGene)", "RefSeq", "GI", "PDB", "GO", "UniRef100", "UniRef90",
           "UniRef50", "UniParc", "PIR", "NCBI-taxon", "MIM", "UniGene", "PubMed", "EMBL", "EMBL-CDS", "Ensembl",
           "Ensembl_TRS", "Ensembl_PRO", "Additional PubMed"]

wanted_headers = ["UniProtKB-AC", "UniProtKB-ID", "GeneID (EntrezGene)", "PDB", "GO"]


def create_mouse_dataframe():
    # mouse = pd.read_table('MOUSE_10090_idmapping_selected.tab.gz', delimiter='\t', names=headers)
    # mouse = mouse[wanted_headers]
    # mouse.to_csv('../mouse_uniprot.gz', compression='gzip', columns=wanted_headers, header=True)
    mouse = pd.read_csv('../mouse_uniprot.gz')
    return mouse


mouse = create_mouse_dataframe()


def create_human_dataframe():
    human = pd.read_table('HUMAN_9606_idmapping_selected.tab.gz', delimiter='\t', names=headers)
    human = human[wanted_headers]
    human.to_csv('../human_uniprot.gz', compression='gzip', columns=wanted_headers, header=True)
    human = pd.read_csv('../human_uniprot.gz')
    return human


human = create_human_dataframe()

mouse_string = ''
for i in np.array(mouse["UniProtKB-AC"]):
    mouse_string += '%s\n' % i

with open('query_file_mouse_uniprot.txt', 'w') as f:
    f.write(mouse_string)

human_string = ''
for i in np.array(human["UniProtKB-AC"]):
    human_string += '%s\n' % i

with open('query_file_human_uniprot.txt', 'w') as f:
    f.write(human_string)

uniprot_id_mouse_to_kegg = pd.read_table('mouseID_to_kegg.tab', delimiter='\t', names=["UniProtKB-ID", 'kegg'],
                                         skiprows=1)
print(uniprot_id_mouse_to_kegg.head(1))

print("Shape of original data", np.shape(mouse))
print("Shape of converted data", np.shape(uniprot_id_mouse_to_kegg))
result = uniprot_id_mouse_to_kegg.merge(mouse, how='right')
result = result[["UniProtKB-AC", "UniProtKB-ID", 'kegg']]
print("Shape of combined data", np.shape(result))

result.to_csv('mouse.csv', sheet_name='Sheet1')
j = uniprot_id_mouse_to_kegg["UniProtKB-ID"].astype(str)
k = mouse["UniProtKB-ID"].astype(str)

for i in k:
    if i in j:
        print('here', i)
        # else:
        #     print('here2',i)
