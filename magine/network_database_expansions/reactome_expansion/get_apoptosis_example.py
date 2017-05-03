# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 11:19:59 2016

@author: pinojc
"""

import networkx as nx
import pandas as pd
import pygraphviz as pyg

from magine.bioservices_local.reactome_changes import Reactome

pd.set_option('expand_frame_repr', False)
pd.set_option('display.width', 10000)
pd.set_option('display.max_rows', 5000)
pd.set_option('display.max_columns', 5000)
pd.options.display.max_colwidth = 5000

reactome = Reactome()

verbose = True
"""

reaction
    name
    input
        list of species, compartment
    output
        list of species, compartment


    dbId -> Protein, Chemical Compound, Complex, Set, OtherEntity, Genes and Transcripts, DNA Sequence

    Complex
        Protein, Chemical Compound
"""

name_dict = {
    'plasma membrane':                   'PlasmaMem',
    'cytosol':                           'CYTO',
    'mitochondrial outer membrane':      'MOM',
    'nucleoplasm':                       'NucPlasm',
    'endosome membrane':                 'EndosomeMem',
    'mitochondrial intermembrane space': 'MIM',
    'nuclear envelope':                  'NucEnv',
    'endoplasmic reticulum membrane':    'ERM',
    '[':                                 r'\n['
}


shapes = {
    'Protein':               'box',
    'Chemical Compound':     'oval',
    'Complex':               'note',
    'Set':                   'rectangle',
    'OtherEntity':           'note',
    'Genes and Transcripts': 'note',
    'DNA Sequence':          'note',
    'RNA Sequence':          'note'
}


def shorten_name(string):
    for s in name_dict:
        if s in string:
            string = string.replace(s, name_dict[s])
    return string


def add_to_graph(sample, list_of_species, graph):
    if type(sample) == int:
        # print('output type = int')
        # outputs.add(each)
        return
    if sample['className'] not in shapes:
        if verbose:
            print('Is a {}'.format(sample['className']))
        return
    name = sample['dbId']
    display = shorten_name(sample['displayName'])
    graph.add_node(name,
                   label=display,
                   displayName=display,
                   dbId=name,
                   shape=shapes[sample['className']])

    list_of_species[name] = display


def extract_list_from_key(keyword, y):
    input_list = []
    for i in y[keyword]:
        if isinstance(i, int):
            input_list.append(i)
        else:
            if 'dbId' in i.keys():
                input_list.append(i['dbId'])
            else:
                if verbose:
                    print("No dbID and not an int")
    return input_list


def get_entity_info(species):
    x = reactome.get_entity_info(species)

    info_needed = ['referenceType', 'className', 'startCoordinate',
                   'endCoordinate', 'referenceEntity', 'compartment',
                   'referenceEntity', 'displayName', 'hasModifiedResidue']

    dont_need = ['inDisease', 'name', 'dbId', 'stId', 'inferredTo',
                 'speciesName', 'species', 'schemaClass', 'isChimeric',
                 'literatureReference']

    need_to_investigate = ['hasComponent', 'literatureReference',
                           'hasMember', 'summation', 'hasCandidate']
    entity_info = dict()
    for n in x:
        if n == 'databaseObject':
            if 'className' in x[n]:
                c_name = x[n]['className']
                if c_name == 'Reaction':
                    continue

            if 'referenceEntity' in x[n]:
                info = x[n]['referenceEntity']
                if 'dbId' in info:
                    entity_info['parent_dbid'] = info['dbId']
                if 'databaseName' in info:
                    entity_info['parent_db_name'] = info['databaseName']
                if 'identifier' in info:
                    entity_info['parent_identifier'] = info['identifier']
            if 'hasModifiedResidue' in x[n]:
                info = x[n]['hasModifiedResidue']
                mod_residues = []
                mods = []
                for r in info:
                    if 'coordinate' in r:
                        mod_residues.append(r['coordinate'])
                    if 'psiMod' in r:
                        ms = r['psiMod']
                        if isinstance(ms, int):
                            mods.append(ms)
                        elif isinstance(ms, dict):
                            if 'dbId' in ms:
                                mods.append(ms['dbId'])
                            else:
                                for m_type in ms:
                                    print(m_type, ms[m_type])
            for m in x[n]:
                if m in info_needed:
                    entity_info[m] = x[n][m]
                elif m not in dont_need:
                    print(species, m, x[n][m])
    return entity_info


def count_entites(list_of_object):
    count_dict = {}
    name_dict = {}
    for i in list_of_object:
        if i in count_dict:
            count_dict[i] += 1
        else:
            count_dict[i] = 1
    for i in count_dict:
        x = reactome.get_entity_attribute(i, 'displayName')
        name_dict[i] = x
        # get_entity_info(i)
    return name_dict, count_dict


def create_graph(inputs, outputs, class_names, catalyst, rxn_name, graph,
                 prev_events):
    for each in inputs:
        name = each
        display = inputs[each]
        graph.add_node(name,
                       label=display,
                       displayName=display,
                       dbId=name,
                       shape=shapes[class_names[each]])
        if catalyst:
            if each == catalyst:
                pass
            else:
                graph.add_edge(catalyst, rxn_name)
        graph.add_edge(name, rxn_name)
        for i in prev_events:
            graph.add_edge(i, name)

    for each in outputs:
        name = each
        display = outputs[each]
        graph.add_node(name,
                       label=display,
                       displayName=display,
                       dbId=name,
                       shape=shapes[class_names[each]])
        graph.add_edge(rxn_name, each)


def create_graph_from_df(inputs, outputs, class_names, catalyst, rxn_name,
                         graph, prev_events):
    for each in inputs:
        name = each
        display = inputs[each]
        graph.add_node(name,
                       label=display,
                       displayName=display,
                       dbId=name,
                       shape=shapes[class_names[each]])
        if catalyst:
            if each == catalyst:
                pass
            else:
                graph.add_edge(catalyst, rxn_name)
        graph.add_edge(name, rxn_name)
        for i in prev_events:
            graph.add_edge(i, name)

    for each in outputs:
        name = each
        display = outputs[each]
        graph.add_node(name,
                       label=display,
                       displayName=display,
                       dbId=name,
                       shape=shapes[class_names[each]])
        graph.add_edge(rxn_name, each)


def _extract_info_from_reactants_or_products(r_dict, in_out):
    dict_of_info = dict()
    return_info = dict()
    for each in r_dict[in_out]:
        if isinstance(each, int):
            if each in return_info:
                return_info[each]['counter'] += 1
                continue
            else:
                if verbose:
                    print(each, "Not in inputs")
                continue
        dict_of_info['id'] = each['dbId']
        dict_of_info['counter'] = 1
        dict_of_info['species_type'] = each['className']
        dict_of_info['displayname'] = each['displayName']
        return_info[each['dbId']] = dict_of_info
    return return_info


def get_reaction_info(reaction, graph):

    if 'className' not in reaction:
        return None, None

    # ensure class is a reaction
    if reaction['className'] != 'Reaction':
        if verbose:
            print("Not a reaction")
            print(reaction)
        return None, None

    # get all reaction info
    y = reactome.get_reaction_info(reaction['dbId'])

    # check to see if compartment exists
    if 'compartment' not in y:
        if verbose:
            print("No compartment")
        return None, None

    # check to make sure there is input and output
    if ('input' or 'output') not in y:
        if verbose:
            print("No input or output")
        return None, None

    # get reaction id and name
    rxn_name = reaction['dbId']
    rxn_display_name = reaction['displayName']

    if verbose:
        print("Reaction {} : {} ".format(rxn_name, rxn_display_name))
        for i in y:
            print("\t{} : {}".format(i, y[i]))

    # get compartment info
    comp = y['compartment']
    compartments = []

    for c in comp:
        if 'name' in c:
            compartments.append(c['name'])
    prev_events = []
    if 'precedingEvent' in y:
        prec_event = y['precedingEvent']
        if verbose:
            print("\tPreceeding event = {}".format(prec_event))
        for i in prec_event:
            if 'dbId' in i:
                prev_event_id = i['dbId']
                prev_events.append(prev_event_id)

    catalyst = False
    cat_all = []
    if 'catalystActivity' in y:
        catalyst = y['catalystActivity'][0]['dbId']
        if verbose:
            print(catalyst, 1)
        test = reactome.get_event_participating_phys_entities(catalyst)
        # print('TEST', test)
        for n in test:
            if verbose:
                for j in n:
                    print('\t\t {} : {}'.format(j, n[j]))
        for n in test:
            if 'displayName' in n:
                catalyst = n['dbId']
                if verbose:
                    print(catalyst, 2)
                cat_all.append(catalyst)
                graph.add_node(catalyst,
                               label=shorten_name(n['displayName']),
                               shape='box', fillcolor='grey', style='filled')

    inputs = _extract_info_from_reactants_or_products(y, 'input')
    outputs = _extract_info_from_reactants_or_products(y, 'output')

    if verbose:
        print('Inputs : {}'.format(inputs))
        print('Outputs : {}'.format(outputs))

    reaction_info = dict()
    reaction_info['inputs'] = inputs
    reaction_info['outputs'] = outputs
    reaction_info['name'] = rxn_display_name
    reaction_info['id'] = rxn_name
    reaction_info['compartment'] = compartments
    reaction_info['catalyst'] = catalyst
    reaction_info['prev_events'] = prev_events
    # return reaction_info, rxn_name
    # quit()

    class_names = {}
    inputs = {}
    outputs = {}
    for each in y['input']:
        if isinstance(each, int):
            if each in inputs:
                continue
        id = each['dbId']
        display = shorten_name(each['displayName'])
        inputs[id] = display
        class_names[id] = each['className']

    for each in y['output']:
        if isinstance(each, int):
            if verbose:
                print(type(each))
            if each in outputs:
                continue
        else:
            id = each['dbId']
            display = shorten_name(each['displayName'])
            outputs[id] = display
            class_names[id] = each['className']

    if len(cat_all) > 1:
        if verbose:
            print(cat_all)
        quit()

    graph.add_node(rxn_name, displayName=rxn_display_name,
                   label=rxn_display_name)
    create_graph(inputs, outputs, class_names, catalyst, rxn_name, graph,
                 prev_events)

    return reaction_info, rxn_name


def extract_pathways(pathway_events):
    pathways = []
    for i in pathway_events:
        if verbose:
            print(i)
        if i['className'] == 'Pathway':
            pathways.append(i['dbId'])
    return pathways


if __name__ == "__main__":

    save_name = 'apoptosis'
    save_name = 'bad_activation'
    save_name = 'apoptosis'
    save_name = 'test'

    apoptosis = 109581
    simple_bax_trans = 114294
    necrosis = 5218859
    bad_activation = 111447

    g = pyg.AGraph(directed=True)
    g.graph_attr['rankdir'] = 'LR'
    g.graph_attr['splines'] = 'true'

    # apoptosis_complex = reactome.pathway_complexes(109581)
    x = reactome.get_pathway_events(195721)

    pathways = extract_pathways(x)
    print(pathways)
    # quit()

    reactions = []
    all_info = []
    all_df = []


    def extract_list_of_events(events):
        for i in events[:10]:
            info, rxn = get_reaction_info(i)

            # df = pd.DataFrame(info)
            if info is None:
                continue
            all_info.append(info)
            reactions.append(rxn)


    extract_list_of_events(x)

    pathways = [111457]
    for i in pathways:
        x = reactome.get_pathway_events(i)
        extract_list_of_events(x)

    df = pd.DataFrame(all_info)

    df.to_csv('reactome_df_{}.csv'.format(save_name))
    # """
    # df = pd.read_csv('reactome_df.csv')

    translocations = df[df['compartment'].map(len) > 1]
    print("Translocations found")
    print(translocations.shape)

    gml_g = nx.nx_agraph.from_agraph(g)
    nx.write_gml(gml_g, '{}.gml'.format(save_name), )
    print("Created GML file!")
    node_species = []
    for i in g.nodes():
        if i not in reactions:
            node_species.append(i)

    # g.add_subgraph(reactions, name='cluster1', label='Reactions')
    # g.add_subgraph(species, name='cluster2', label='Species')

    g.write('{}.dot'.format(save_name), )
    # g.draw('{}.png'.format(save_name), prog='dot')
