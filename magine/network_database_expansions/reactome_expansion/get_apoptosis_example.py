# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 11:19:59 2016

@author: pinojc
"""

from bioservices_local.reactome_changes import Reactome

reactome = Reactome()

apoptosis_complex = reactome.pathway_complexes(109581)

# complexs = []
# for i in apoptosis_complex:
#     if i['schemaClass'] == 'Complex':
#         # print(i)
#         cont = True
#         if len(complexs) == 0:
#             continue
#         else:
#             print(complexs)
#             complexs = []
#     elif i['schemaClass'] == 'EntityWithAccessionedSequence':
#         species = i['displayName']
#         complexs.append(species)
# species = reactome.pathway_participants(109581)
# for i in species:
#     print(i)

# res = reactome.biopax_exporter(109581)
# print(res)
# x = reactome.get_all_reactions()
# for i in x:
#     print(i)

# x = reactome.bioservices_get_reactants_input_output('R-HSA-141409')
#
# for i in x.keys():
#     print(i, x[i])
#     if i =='input':
#         for each in x[i]:
#             print(each['displayName'])

apoptosis = 109581
necrosis = 5218859
x = reactome.get_pathway_events(necrosis)
import pygraphviz as pyg

g = pyg.AGraph(directed=True)


def shorten_compartment_name(string):
    name_dict = {'plasma membrane'                  : 'PlasmaMem',
                 'cytosol'                          : 'CYTO',
                 'mitochondrial outer membrane'     : 'MOM',
                 'nucleoplasm'                      : 'NucPlasm',
                 'endosome membrane'                : 'EndosomeMem',
                 'mitochondrial intermembrane space': 'MIM',
                 'nuclear envelope'                 : 'NucEnv',
                 'endoplasmic reticulum membrane'   : 'ERM',
                 '['                                : r'\n['}
    for s in name_dict:
        if s in string:
            string = string.replace(s, name_dict[s])
    return string


shapes = {'Protein'          : 'box',
          'Chemical Compound': 'oval',
          'Complex'          : 'note'}


def add_to_graph(sample, list_of_species, graph):
    if type(sample) == int:
        # print('output type = int')
        # outputs.add(each)
        return
    if sample['className'] != 'Protein':
        if sample['className'] != 'Chemical Compound':
            print('Is a {}'.format(sample['className']))
    name = sample['dbId']
    graph.add_node(name,
                   label=shorten_compartment_name(sample['displayName']),
                   shape=shapes[sample['className']])

    list_of_species[name] = shorten_compartment_name(sample['displayName'])


def get_reaction_info(reaction):

    # print("Reaction : {}".format(reaction))
    if 'className' not in reaction:
        # print("No class name")
        return
    if reaction['className'] != 'Reaction':
        print("Not a reaction")
        print(reaction)
        # y = reactome.get_pathway_events(reaction['dbId'])
        # print(y)
        # for i in y:
        #     get_reaction_info(i)
        # print('{} name : {}'.format(reaction['className'], reaction['displayName']))
        return
    print(
        '{} name : {}'.format(reaction['className'], reaction['displayName']))
    # for i in reaction.keys():
    #     print('{} : {}'.format(i, reaction[i]))

    y = reactome.get_reaction_info(reaction['stId'])
    # print('Reaction info : {}'.format(y))
    for i in y.keys():
        print('{} : {}'.format(i, y[i]))

    if ('input' or 'output') not in y:
        print("No input or output")
        return
    catalyst = False
    cat_all = []
    if 'catalystActivity' in y:
        # print('Catalyst reaction')
        # print(y['catalystActivity'])
        # print('number of cat {}'.format(len(y['catalystActivity'])))
        catalyst = y['catalystActivity'][0]['dbId']
        test = reactome.get_event_participating_phys_entities(catalyst)
        # print('TEST', test)
        for n in test:
            if 'displayName' in n:
                catalyst = n['dbId']
                # if n['className'] == 'Set':
                #     print('set of species')
                #     print(n['stId'])
                #     x = reactome.get_event_participants(n['stId'])
                #     print('HERE', x)
                #     quit()
                cat_all.append(catalyst)
                g.add_node(catalyst,
                           label=shorten_compartment_name(n['displayName']),
                           shape='box', fillcolor='grey', style='filled')
    if len(cat_all) > 1:
        print(cat_all)
        quit()
    print("number of inputs : {}".format(len(y['input'])))
    print("number of outputs : {}".format(len(y['output'])))

    inputs = {}
    outputs = {}

    for each in y['input']:
        add_to_graph(each, inputs, g)

    for each in y['output']:
        add_to_graph(each, outputs, g)

    for i in inputs:
        print("Input : {}".format(inputs[i]))
    for i in outputs:
        print("Output : {}".format(outputs[i]))
    # z = reactome.get_event_participants(y['input'][0]['consumedByEvent'][0])
    # print('event participant : {}'.format(z))
    # for i in z:
    #     print(i)
    #     p = i['refEntities']
        #     print(i['peDbId'])
        #     label = shorten_compartment_name(i['displayName'])
        #     g.add_node(i['peDbId'], label=label, shape='box')
        # #     for each in p:
        # #        print(each['dbId'])
        # #        print(each['identifier'])
        # #
    # z = reactome.get_event_participants(y['output'][0]['producedByEvent'][0])
    # print('event participant : {}'.format(z))
    # for i in z:
    #     print(i)
    #     p = i['refEntities']
    #     print(i['peDbId'])
        #     label = shorten_compartment_name(i['displayName'])
        #     g.add_node(i['peDbId'], label=label, shape='box')
        # for each in p:
        #    print('dbid',each['dbId'])
        #    print('id',each['identifier'])

    for each in inputs:
        # if catalyst:
        #     g.add_edge(catalyst, each)
        g.add_edge(each, 'rxn_{}'.format(reaction['dbId']))
    for n in outputs:
        g.add_edge('rxn_{}'.format(reaction['dbId']), n)


get_reaction_info(x[3])
g.draw('test.pdf', prog='dot')
g.write('test.dot')
quit()
samples = [149, 150, 151, 152, 153, 154, 155, 156]
for i in range(len(x)):
    # for i in samples:
    # for i in range(11, 20):
    print('ITERATION NUMBER {}\n'.format(i))
    get_reaction_info(x[i])
g.draw('test.pdf', prog='dot')
print("Number of nodes :{}".format(len(g.nodes())))
print("Number of edges :{}".format(len(g.edges())))
quit()
