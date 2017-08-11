# -*- coding: utf-8 -*-
"""
Mapping between data identifiers
"""
try:
    import cPickle as pickle
except:
    import pickle as pickle
import networkx as nx
from bioservices import HGNC, KEGG, UniChem, UniProt

from magine.mappings.chemical_mapper import ChemicalMapper
from magine.mappings.gene_mapper import GeneMapper

kegg = KEGG()
uniprot = UniProt()
hugo = HGNC()
chem = UniChem()


# TODO create a database to store all of this
# TODO create a kegg to uniprot identifier dictionary


# creation of a manual dictionary because of kegg to uniprot errors.
# These errors are mostly on KEGG side that link it to an unreviewed Uniprot ID
manual_dict = {'hsa:857'      : 'CAV1',
               'hsa:2250'     : 'FGF5',
               'hsa:5337'     : 'PLD1',
               'hsa:4312'     : 'MMP1',
               'hsa:102723407': 'IGHV4OR15-8',
               'hsa:100132074': 'FOXO6',
               'hsa:728635'   : 'DHRS4L1',
               'hsa:10411'    : 'RAPGEF3',
               'hsa:100101267': 'POM121C',
               'hsa:2768'     : 'GNA12',
               'hsa:2044'     : 'EPHA5',
               'hsa:100533467': 'BIVM-ERCC5',
               'hsa:7403'     : 'KDM6A',
               'hsa:1981': 'EIF4G1',
               'hsa:2906': 'GRIN2D',
               'hsa:4088': 'SMAD3',
               'hsa:6776': 'STAT5A',
               'hsa:182': 'JAG1',
               'hsa:3708': 'ITPR1',
               'hsa:1293': 'COL6A3',
               'hsa:93034': 'NT5C1B',
               }

compound_manual = {'cpd:C07909': 'HMDB15015',
                   'cpd:C16844': 'HMDB01039',
                   'cpd:C00076': 'HMDB00464',
                   'cpd:C00154': 'HMDB01338',
                   'cpd:C01561': 'HMDB03550',
                   'cpd:C04043': 'HMDB03791',
                   'cpd:C01165': 'HMDB02104',
                   'cpd:C00025' : 'HMDB00148',
                   'cpd:C00696': 'HMDB01403',

                   }


def create_gene_dictionaries(network, species='hsa'):
    """
    maps network from kegg to gene names
    
    Parameters
    ----------
    network : networkx.DiGraph
    species : str
        species of genes (HSA, MMU, etc)

    Returns
    -------

    """

    gm = GeneMapper(species)
    # first we will check the pre-created dictionaries

    # Create the dictionary to store all conversions to be returned
    kegg_to_gene_name = {}
    # List to store things not in the initial dictionary
    unknown_genes = set()
    still_missing = set()
    nodes = set(network.nodes())
    # check stores dictionaries
    for i in nodes:
        str_gene = str(i)
        network.node[i]['keggName'] = str_gene
        if str_gene.startswith(species):
            if str_gene in manual_dict:
                kegg_to_gene_name[str_gene] = manual_dict[str_gene]
                print(str_gene, manual_dict[str_gene])
                continue

            tmp_name = str_gene.replace(species + ':', '')
            if str_gene in gm.kegg_to_gene_name:
                if len(gm.kegg_to_gene_name[str_gene]) == 1:
                    kegg_to_gene_name[str_gene] = gm.kegg_to_gene_name[str_gene][0]
                else:
                    for g in gm.kegg_to_gene_name[str_gene]:
                        if g in gm.gene_name_to_uniprot:
                            kegg_to_gene_name[str_gene] = g

            elif int(tmp_name) in gm.ncbi_to_symbol:
                new = gm.ncbi_to_symbol[int(tmp_name)][0]
                if isinstance(new, float):
                    unknown_genes.add(i)
                    continue
                kegg_to_gene_name[str_gene] = str(new)
            else:
                unknown_genes.add(i)
    if len(unknown_genes) == 0:
        return kegg_to_gene_name, 1
    # create string to call uniprot for mapping
    search_string = ''
    for n, i in enumerate(unknown_genes):
        search_string += str(i) + '\t'
    search_string = search_string.rstrip('\t')

    # This is where it gets tricky. Checking to see if there is a uniprot
    # mapping for the species, if not, trying from KEGG side. Sometimes kegg
    # links to a different uniprot, or uniprot links to a different kegg.
    uni_dict = dict(uniprot.mapping("KEGG_ID", "ACC", query=search_string))

    for i in unknown_genes:
        if i in uni_dict:
            for n in uni_dict[i]:
                # print(i, uni_dict[i], n)
                x = uniprot.search("accession:{}".format(n),
                                   columns='genes(PREFERRED),reviewed,id',
                                   limit=1)
                header, data = x.rstrip('\n').split('\n')
                name, review, entry = data.split('\t')
                if n != entry:
                    print(i, n, entry, x, "dont match")
                elif review == 'reviewed':
                    kegg_to_gene_name[i] = name

        else:
            still_missing.add(i)
    print("{} mappings not found from kegg to"
          " gene name".format(len(still_missing)))
    print(still_missing)
    return kegg_to_gene_name, 0


def hugo_mapper(network, species='hsa'):
    """
    Converts all KEGG names to HGNC 
    
    Parameters
    ----------
    network : networkx.DiGraph
    species : str

    Returns
    -------
    dict
    """
    prefix = species+':'
    nodes = set(network.nodes())
    hugo_dict = {}
    not_found = set()
    for i in nodes:
        if str(i).startswith(prefix):
            tmp_name = str(i).replace(prefix, '')
            mapping = hugo.search(tmp_name)
            if 'response' in mapping:
                response = mapping['response']
                if 'numFound' in response:
                    if response['numFound'] == 0:
                        not_found.add(i)
                        continue
                    elif response['numFound'] == 1:
                        docs = response['docs'][0]
                        hugo_dict[i] = docs['symbol']
                        continue
                    else:
                        if 'symbol' in response['docs'][0]:
                            hugo_dict[i] = response['docs'][0]['symbol']
            else:
                not_found.add(i)
    if not_found != 0:
        print("{} mappings not found after HGNC mapping".format(len(not_found)))
        print("{} ".format(not_found))
    return hugo_dict


def create_compound_dictionary(network):
    """
    Maps network from kegg to gene names
    
    Parameters
    ----------
    network : networkx.DiGraph

    Returns
    -------
    dict
    
    """
    cm = ChemicalMapper()
    cpd_to_hmdb = {}  # he
    still_unknown = []
    nodes = set(network.nodes())
    for i in nodes:
        if i.startswith('cpd:'):
            network.node[i]['keggName'] = i
            name_stripped = i.lstrip('cpd:')
            if name_stripped in cm.kegg_to_hmdb_accession:
                mapping = cm.kegg_to_hmdb_accession[name_stripped]
                if type(mapping) == list:
                    names = '|'.join(i for i in mapping)
                    cpd_to_hmdb[i] = names
                    network.node[i]['hmdbNames'] = names
                    chem_names = []
                    for name in mapping:
                        if name in cm.hmdb_accession_to_chemical_name:
                            chem_names += [e for e in cm.hmdb_accession_to_chemical_name[name]]
                    chem_names = '|'.join(e for e in chem_names)
                    network.node[i]['chemName'] = chem_names.encode('ascii',
                                                                    'ignore')

                elif type(mapping) == str:
                    cpd_to_hmdb[i] = mapping
                    chem_n = cm.hmdb_accession_to_chemical_name[mapping]
                    print(chem_n)
                    network.node[i]['chemName'] = chem_n
                else:
                    print('Returned something else...', mapping)
            else:
                still_unknown.append(i)
    if len(still_unknown) == 0:
        return cpd_to_hmdb
    kegg_hmdb = chem.get_mapping("kegg_ligand", "hmdb")
    for i in still_unknown:
        if i.lstrip('cpd') in kegg_hmdb:
            cpd_to_hmdb[i] = kegg_hmdb[i.lstrip('cpd:')][0]
        # else:
        #     print("Cannot find a HMDB mapping for %s " % i)
    return cpd_to_hmdb


def convert_all(network, species='hsa'):
    """ 
    Maps gene names to HGNC and kegg compound to HMDB
    
    Parameters
    ----------
    network : networkx.DiGraph
        network to convert mappings
    species : str   
        species of network (hsa, mmu)

    Returns
    -------

    """

    change_dict = dict()
    change_dict.update(compound_manual)
    renamed_network = network.copy()

    print('Started converting kegg compounds to HMDB')
    dict1 = create_compound_dictionary(renamed_network)
    change_dict.update(dict1)

    print('Started converting kegg genes to HGNC')
    dict2, found_all = create_gene_dictionaries(renamed_network,
                                                species=species)
    change_dict.update(dict2)

    if not found_all:
        print('Started to check for miRNAs')
        dict3 = hugo_mapper(renamed_network, species=species)
        change_dict.update(dict3)
    change_dict = _check_dict_for_int(change_dict)

    renamed_network = nx.relabel_nodes(renamed_network, change_dict)

    return renamed_network


def _check_dict_for_int(dic):
    new_dic = dict()
    for key, value in dic.items():
        if isinstance(value, float):
            continue
        else:
            new_dic[key] = value
    return new_dic

"""
def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield (name, ''.join(seq))


def gather_list_of_mouse_approved():
    approved_dict = set()
    unknown = ''
    with open('/home/pinojc/Downloads/uniprot_sprot.fasta') as fp:
        for name, seq in read_fasta(fp):
            _, acc, acc_id = name.split('|')
            gene_id, rest = acc_id.split(' ', 1)
            if gene_id.endswith('_MOUSE'):
                if acc in mouse_uniprot_to_gene_name:
                    approved_dict[acc] = mouse_uniprot_to_gene_name[acc]
                else:
                    unknown += acc + '\n'
                    # print('No gene name? : {}'.format(acc))
    with open('unknown.txt', 'w')as f:
        f.write(unknown)
    return approved_dict
    # approved_dict = gather_list_of_mouse_approved()
    # acc_to_gene = pandas.read_table('acc_to_geneid.tab')
    #
    # acc_to_gene = acc_to_gene.set_index('From')['To'].to_dict()
"""
