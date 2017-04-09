# -*- coding: utf-8 -*-
"""
Mapping between data identifiers
"""
import cPickle as pickle
import os
from sys import modules

import networkx as nx
from bioservices import HGNC, KEGG, UniChem, UniProt

from chemical_mapper import ChemicalMapper

try:
    kegg = modules['kegg']
except KeyError:
    kegg = KEGG()

uniprot = UniProt()
hugo = HGNC()
chem = UniChem()

# TODO create a database to store all of this
# TODO create a kegg to uniprot identifier dictionary
directory = os.path.dirname(__file__)
# mouse genes
# mouse_uniprot_to_gene_name = pickle.load(open(os.path.join(directory, 'mouse_uniprot_to_gene_mapper.p'), 'rb'))
# human genes
# human_uniprot_to_gene_name = pickle.load(open(os.path.join(directory, 'network_based_UPids_to_genes.p'), 'rb'))
# human_uniprot_to_gene_name_2 = pickle.load(open(os.path.join(directory, 'kegg_to_uniprot_gene_name.p'), 'rb'))
# human_uniprot_to_gene_name = pd.read_table(os.path.join(directory, 'humanID_to_gene_name.tab'))


# creation of a manual dictionary because of kegg to uniprot errors.
# These errors are mostly on KEGG sides that link it to an unreviewed Uniprot ID
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
               'hsa:427': '',
               'hsa:2005': '',
               'hsa:7010': '',
               'hsa:3710': '',
               'hsa:5159': '',
               'hsa:1735': '',





               }

compound_manual = {'cpd:C07909': 'HMDB15015',
                   'cpd:C16844': 'HMDB01039'
                   }


def create_gene_dictionaries(network, species='hsa'):
    """
    maps network from kegg to gene names
    :param network:
    :param species:
    """

    # first we will check the pre-created dictionaries
    if species == 'hsa':
        dictionary = pickle.load(open(os.path.join(directory,
                                                   'human_kegg_mapper.p'),
                                      'rb'))
    elif species == 'mmu':
        dictionary = pickle.load(open(os.path.join(directory,
                                                   'mouse_kegg_mapper.p'),
                                      'rb'))
    # Create the dictionary to store all conversions to be returned
    kegg_to_gene_name = {}
    # List to store things not in the initial dictionary
    unknown_genes = set()
    still_missing = set()
    nodes = set(network.nodes())
    # check stores dictionaries
    for i in nodes:
        if str(i).startswith(species):
            network.node[i]['keggName'] = i
            if i in dictionary:
                kegg_to_gene_name[i] = dictionary[i]
            elif i in manual_dict:
                kegg_to_gene_name[i] = manual_dict[i]
            else:
                if str(i).startswith(species):
                    unknown_genes.add(i)
    if len(unknown_genes) == 0:
        return kegg_to_gene_name
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
    return kegg_to_gene_name


def hugo_mapper(network, species='hsa'):
    """ Converts all MIR from kegg

    :param species:
    :param network:
    :return:

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
    print("{} mappings not found after HGNC mapping".format(len(not_found)))
    return hugo_dict


def create_compound_dictionary(network):
    """ Maps network from kegg to gene names

    :param network:
    :return: dictionary

    """
    cm = ChemicalMapper()

    cm.reload()
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
                    hmdb_options = mapping[0]
                    if len(mapping) == 1:
                        cpd_to_hmdb[i] = hmdb_options
                    else:
                        cpd_to_hmdb[i] = hmdb_options
                        for j in range(1, len(mapping)):
                            hmdb_options += "__%s" % mapping[j]
                    network.node[i]['hmdbNames'] = hmdb_options.rstrip('__')
                    chem_names = ''
                    for name in mapping:
                        if name in cm.hmdb_accession_to_chemical_name:
                            for each in cm.hmdb_accession_to_chemical_name[name]:
                                # print(each)
                                chem_names += "%s__" % each
                    network.node[i]['chemName'] = chem_names.rstrip('__')
                elif type(mapping) == str:
                    cpd_to_hmdb[i] = mapping
                    network.node[i]['chemName'] = cm.hmdb_accession_to_chemical_name[mapping]
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
        else:
            print("Cannot find a HMDB mapping for %s " % i)
    return cpd_to_hmdb


def convert_all(network, species='hsa'):
    """ Maps gene names to HGNC and kegg compound to HMDB

    :param network:
    :param species:
    :return:
    """
    renamed_network = nx.relabel_nodes(network, compound_manual)
    print('Started converting kegg compounds to HMDB')
    dict1 = create_compound_dictionary(renamed_network)
    renamed_network = nx.relabel_nodes(renamed_network, dict1)
    print('Started converting kegg genes to HGNC')
    dict2 = create_gene_dictionaries(renamed_network, species=species)
    renamed_network = nx.relabel_nodes(renamed_network, dict2)
    print('Started to check for miRNAs')
    dict3 = hugo_mapper(renamed_network, species=species)
    renamed_network = nx.relabel_nodes(renamed_network, dict3)

    return renamed_network


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
