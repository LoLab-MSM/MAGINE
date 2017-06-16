import textwrap

import networkx as nx
import pandas as pd
import pygraphviz as pyg

import \
    magine.network_database_expansions.reactome.reactome_api as reactome_api


def test_reactome():
    save_name = 'test_apoptosis'

    # apoptosis = 109581
    apoptosis = 111457
    apoptosis = 114294
    #
    g = pyg.AGraph(directed=True)
    # g.graph_attr['rankdir'] = 'LR'
    g.graph_attr['splines'] = 'true'

    x = reactome_api._reactome.get_pathway_events(apoptosis)
    reactions = []
    all_info = []

    def extract_list_of_events(events):
        for i in events:
            info, rxn = reactome_api.get_reaction_info(i, g)
            if info is None:
                continue
            all_info.append(info)
            reactions.append(rxn)

    x = reactome_api._reactome.get_pathway_events(apoptosis)
    print("{} has {} events".format(save_name, len(x)))
    print("Extracting events")
    extract_list_of_events(x)

    df = pd.DataFrame(all_info)

    df.to_csv('reactome_df_{}.csv'.format(save_name))
    # """
    # df = pd.read_csv('reactome_df.csv')
    if 'compartment' in df.columns:
        translocations = df[df['compartment'].map(len) > 1]
        if len(translocations) > 0:
            print("Translocations found")
            print(translocations)
            print(translocations.shape)
    for i in g.nodes():
        ent_dict = reactome_api.get_entity_info(i)
        node = g.get_node(i)
        if 'displayName' not in ent_dict:
            lab = node.attr['label']
            print(lab, textwrap.wrap(str(lab), 40))
            node.attr['label'] = '\n'.join(textwrap.wrap(str(lab), 20))
            continue

        print("Node and attributes = {} : {} ".format(i, node.attr))

        # if 'label' not in node.attr:
        #     node.attr['label'] = reactome_api.shorten_name(
        #             ent_dict['displayName'])
        # else:
        node.attr['label'] = reactome_api.shorten_name(ent_dict['displayName'])

    gml_g = nx.nx_agraph.from_agraph(g)
    nx.write_gml(gml_g, '{}.gml'.format(save_name))
    print("Created GML file!")

    g.write('{}.dot'.format(save_name), )
    g.draw('{}.png'.format(save_name), prog='dot')


if __name__ == '__main__':
    test_reactome()
