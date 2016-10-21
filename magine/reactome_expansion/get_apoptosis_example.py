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

x = reactome.get_pathway(109581)

import pygraphviz as pyg

g = pyg.AGraph(directed=True)


def shorten_compartment_name(string):
    name_dict = {'plasma membrane'                  : 'PM',
                 'cytosol'                          : 'CYTO',
                 'mitochondrial outer membrane'     : 'MOM',
                 'nucleoplasm'                      : 'NP',
                 'endosome membrane'                : 'EndoMem',
                 'mitochondrial intermembrane space': 'MIM'}
    for i in name_dict:
        if i in string:
            string = string.replace(i, name_dict[i])
    return string


def get_reaction_info(reaction):
    # print("Reaction : {}".format(reaction))
    if 'className' not in reaction:
        # print("No class name")
        return
    if reaction['className'] != 'Reaction':
        # print("Not a reaction")
        # print('{} name : {}'.format(reaction['className'], reaction['displayName']))
        return
    print(
    '{} name : {}'.format(reaction['className'], reaction['displayName']))
    for i in reaction.keys():
        print('{} : {}'.format(i, reaction[i]))
    print('\n')
    y = reactome.get_reaction_info(reaction['stId'])
    # print('Reaction info : {}'.format(y))
    for i in y.keys():
        print('{} : {}'.format(i, y[i]))
    inputs = []
    outputs = []
    if ('input' or 'output') not in y:
        print("No input or output")
        return
    catalyst = False
    if 'catalystActivity' in y:
        print('Catalyst reaction')
        print(y['catalystActivity'])
        catalyst = y['catalystActivity'][0]['dbId']
        test = reactome.get_event_participating_phys_entities(catalyst)
        print('TEST', test)
        for n in test:
            if 'displayName' in n:
                catalyst = n['displayName']
                catalyst = shorten_compartment_name(catalyst)

    for each in y['input']:
        if type(each) == int:
            print('input type = int')
            continue
        # for k in each.keys():
        #     print('Input {} : {}'.format(k, each[k]))
        # print('input keys {}'.format(each.keys()))
        print('Inputs : {}'.format(each['displayName']))
        name = each['displayName']
        name = shorten_compartment_name(name)
        inputs.append(name)
    for each in y['output']:
        if type(each) == int:
            print('output type = int')
            continue
        # for k in each.keys():
        #     print('Ouput {} : {}'.format(k, each[k]))
        # for k in each.keys():
        #     print('{} : {}'.format(k, each[k]))
        print('Output : {}'.format(each['displayName']))
        name = each['displayName']
        name = shorten_compartment_name(name)
        outputs.append(name)
    # for each in y['compartment']:
    #     print('compartment : {}'.format(each['displayName']))

    # z = reactome.get_event_participants(y['input'][0]['consumedByEvent'][0])
    # print('event participant : {}'.format(z))
    # for i in z:
    #     print(i)
    #     print(i['peDbId'])
    #     p = i['refEntities']
    #     for each in p:
    #        print(each['dbId'])
    #        print(each['identifier'])
    #
    # z = reactome.get_event_participants(y['output'][0]['producedByEvent'][0])
    # print('event participant : {}'.format(z))
    # for i in z:
    #     print(i)
    #     p = i['refEntities']
    #     print(i['peDbId'])
    #     for each in p:
    #        print(each['dbId'])
    #        print(each['identifier'])

    for each in inputs:
        if catalyst:
            g.add_edge(catalyst, each)
        for n in outputs:
            g.add_edge(each, n)


# get_reaction_info(x[155])
# g.draw('test.pdf', prog='dot')
# quit()
# samples = [154, 155, 156]
for i in range(len(x)):
    # for i in samples:
    # for i in range(10):
    print(i)
    get_reaction_info(x[i])
g.draw('test.pdf', prog='dot')
print("Number of nodes :{}".format(len(g.nodes())))
print("Number of edges :{}".format(len(g.edges())))
quit()
