import numpy as np
import pandas as pd

from Mappings.gene_mapper import gm

columns = ['uniprot_source', 'ensebml_source', 'ensebml_source_id', 'uniprot_target', 'ensebml_target',
           'ensebml_target_id', 'type', 'database_source', 'else']

import networkx as nx


def parse_reactome_protein_protein():
    data = pd.read_table('homo_sapiens.interactions.txt', sep='\t', names=columns, low_memory=False)
    print(np.shape(data))
    data = data[['uniprot_source', 'uniprot_target', 'type']]
    data = data.drop_duplicates()
    print(np.shape(data))
    not_in = set()
    g = nx.DiGraph()
    for index, i in data.iterrows():
        source = i['uniprot_source'][8:]  # .lstrip('UniProt:')
        target = i['uniprot_target'][8:]  # .lstrip('UniProt:')
        interaction = i['type']
        if source in gm.uniprot_to_gene_name:
            if target in gm.uniprot_to_gene_name:
                # print(gm.uniprot_to_gene_name[source], gm.uniprot_to_gene_name[target])
                g.add_edge(gm.uniprot_to_gene_name[source][0], gm.uniprot_to_gene_name[target][0],
                           interactionType=interaction)
                g.add_edge(gm.uniprot_to_gene_name[target][0], gm.uniprot_to_gene_name[source][0],
                           interactionType=interaction)
            else:
                print(target)
                not_in.add(target)
        else:
            print(source)
            not_in.add(source)
    nx.write_gml(g, 'reactome_network.gml')
    print('Nodes={}'.format(len(g.nodes())))
    print('Edges={}'.format(len(g.edges())))
    print(len(not_in))
parse_reactome_protein_protein()
quit()

interaction_dictionary = {}
interaction_dictionary['inhibited by, dephosphorylation'] = {'action': 'dephosphorylation',
                                                             'effect': 'inhibition',
                                                             'direct': 'True'}

interaction_dictionary['inhibition, phosphorylation'] = {'action': 'phosphorylation',
                                                         'effect': 'inhibition',
                                                         'direct': 'True'}

interaction_dictionary['inhibition, dephosphorylation'] = {'action': 'dephosphorylation',
                                                           'effect': 'inhibition',
                                                           'direct': 'True'}

interaction_dictionary['inhibited by, phosphorylated by'] = {'action': 'phosphorylation',
                                                             'effect': 'inhibition',
                                                             'direct': 'True'}
interaction_dictionary['inhibition, ubiquitination'] = {'action': 'ubiquitination',
                                                        'effect': 'inhibition',
                                                        'direct': 'True'}
interaction_dictionary['inhibited by'] = {'action': 'notSpecified',
                                          'effect': 'inhibition',
                                          'direct': 'True'}
interaction_dictionary['inhibition'] = {'action': 'notSpecified',
                                        'effect': 'inhibition',
                                        'direct': 'True'}

interaction_dictionary['indirect effect; activated by'] = {'action': 'notSpecified',
                                        'effect': 'activation',
                                        'direct': 'False'}

interaction_dictionary['inhibition, indirect effect'] = {'action': 'notSpecified',
                                                         'effect': 'inhibition',
                                                         'direct': 'False'}

interaction_dictionary['phosphorylation'] = {'action': 'phosphorylation',
                                             'effect': 'notSpecified',
                                             'direct': 'True'}
interaction_dictionary['state change'] = {'action': 'stateChange',
                                          'effect': 'notSpecified',
                                          'direct': 'True'}

interaction_dictionary['activation, phosphorylation'] = {'action': 'phosphorylation',
                                                         'effect': 'activation',
                                                         'direct': 'True'}
interaction_dictionary['activation, phosphorylation, binding/association'] = {'action': 'phosphorylation',
                                                                              'effect': 'activation',
                                                                              'direct': 'True'}
interaction_dictionary['activated by, phosphorylation'] = {'action': 'phosphorylation',
                                                           'effect': 'activation',
                                                           'direct': 'True'}
interaction_dictionary['activated by, dephosphorylation'] = {'action': 'dephosphorylation',
                                                           'effect': 'activation',
                                                           'direct': 'True'}
interaction_dictionary['activation, indirect effect'] = {'action': 'indirect',
                                                         'effect': 'activation',
                                                         'direct': 'False'}
interaction_dictionary['activated by, indirect effect'] = {'action': 'indirect',
                                                           'effect': 'activation',
                                                           'direct': 'False'}
interaction_dictionary['activated, indirect effect'] = {'action': 'indirect',
                                                        'effect': 'activation',
                                                        'direct': 'False'}
interaction_dictionary['activation'] = {'action': 'notSpecified',
                                        'effect': 'activation',
                                        'direct': 'True'}
interaction_dictionary['activated by'] = {'action': 'notSpecified',
                                          'effect': 'activation',
                                          'direct': 'True'}

interaction_dictionary['phosphorylated by'] = {'action': 'phosphorylation',
                                               'effect': 'notSpecified',
                                               'direct': 'True'}
interaction_dictionary['binding/association'] = {'action': 'binding/association',
                                                 'effect': 'notSpecified',
                                                 'direct': 'True'}
interaction_dictionary['indirect effect'] = {'action': 'notSpecified',
                                             'effect': 'notSpecified',
                                             'direct': 'False'}
interaction_dictionary['dissociation'] = {'action': 'dissociation',
                                          'effect': 'notSpecified',
                                          'direct': 'True'}
interaction_dictionary['compound'] = {'action': 'dissociation',
                                          'effect': 'notSpecified',
                                          'direct': 'True'}
interaction_dictionary['compound'] = {'action': 'dissociation',
                                          'effect': 'notSpecified',
                                          'direct': 'True'}
interaction_dictionary['activated by; activated by, phosphorylation; activated, phosphorylation, indirect effect'] = {'action': 'phosphorylation',
                                          'effect': 'activation',
                                          'direct': 'False'}
interaction_dictionary['activated by; activated, indirect effect; binding/association'] = {'action': 'binding/association',
                                          'effect': 'activation',
                                          'direct': 'False'}
interaction_dictionary['activated, indirect effect; activated by'] = {'action': 'binding/association',
                                          'effect': 'activation',
                                          'direct': 'False'}
interaction_dictionary['expression by; activated, indirect effect; expression regulated by'] = {'action': 'expression',
                                          'effect': 'activation',
                                          'direct': 'False'}
edge_dict = {'-': 'not directed',
             '->': 'activate',
             '<-': 'activated',
             '|-': 'inhibit',
             '-|': 'inhibit',
             '<->': 'activate',
             '|-|': 'inhibit',}
def get_interaction_type(row, graph):

    edge = {}
    edge['score'] = str(row['Score'])
    annotation = str(row['Annotation'])
    annotation = annotation.replace('PPrel: ','')
    annotation = annotation.replace('PCrel: ', '')
    annotation = annotation.replace('GErel: ', '')
    annotation = annotation.replace('ECrel: ', '')
    gene_1 = str(row['Gene1'])
    gene_2 = str(row['Gene2'])
    if '/' in gene_1:
        return



    if 'indirect effect' in annotation:
        edge['direct'] = 'False'
        annotation = annotation.replace('indirect effect', '')
        annotation = annotation.lstrip(';')
        annotation = annotation.strip()
    if row['Direction'] == '-':
        return
    graph.add_node(row['Gene1'])
    graph.add_node(row['Gene2'])
    def check_action(annotation, edge):
        if 'activation' in annotation or 'activated by' in annotation or 'activate' in annotation:
            for i in ('activation',  'activated by','activate'):
                annotation = annotation.replace(i,'')
                annotation = annotation.strip()
                annotation = annotation.lstrip(';')
            edge['effect'] = 'activation'

        elif 'inhibition' in annotation or 'inhibited by' or 'inhibit' or 'inhibite'in annotation:
            for i in ('inhibition',  'inhibited by','inhibite','inhibit'):
                annotation = annotation.replace(i,'')
                annotation = annotation.strip()
                annotation = annotation.lstrip(';')
            edge['effect'] = 'inhibition'

        elif 'catalyze' or 'catalyzed by' in annotation:
            for i in ('catalyze', 'catalyzed by'):
                annotation = annotation.replace(i, '')
                annotation = annotation.strip()
                annotation = annotation.lstrip(';')
            edge['effect'] = 'catalyze'

        elif 'predicted' in annotation:
            annotation = annotation.replace('predicted', '')
            edge['effect'] = 'notSpecified'
            edge['predicted'] = 'True'

        elif 'complex' in annotation:
            annotation = annotation.replace('complex', '')
            edge['effect'] = 'notSpecified'

        elif 'input' in annotation:
            annotation = annotation.replace('input', '')
            edge['effect'] = 'notSpecified'
        else:
            print("WARNING : annotation not found", annotation)

        if annotation == '':
            pass
        elif 'expression' in annotation:
            edge['action'] = 'expression'
        elif 'dephosphorylation' in annotation:
            edge['action'] = 'dephosphorylation'
        elif 'expression' in annotation:
            edge['action'] = 'expression'
        elif 'phosphorylation' or 'phosphorylated by' in annotation:
            edge['action'] = 'phosphorylation'
        else:
            print("WARNING : annotation not found",annotation)
        return edge
    edge = check_action(annotation, edge)

    if 'by' in annotation:
        graph.add_edge(gene_2, gene_1, **edge)
    else:
        graph.add_edge(gene_1, gene_2, **edge)



def parse_reactome_with_annotations():
    data = pd.read_table('FIsInGene_031516_with_annotations.txt', sep='\t', low_memory=False)
    print(data.dtypes)
    print(data.head(10))
    type_of_annotation = set()
    print(data['Direction'].unique())
    data = data[(data['Direction'] != '<-|') & (data['Direction'] != '|->')]
    print(data['Direction'].unique())
    new_graph = nx.DiGraph()
    for n,row in data.iterrows():
        get_interaction_type(row,new_graph)

    nx.write_gml(new_graph,'reactome_fi.gml')
    quit()
    for i in new_graph.edges(data=True):
        print(i)


parse_reactome_with_annotations()
