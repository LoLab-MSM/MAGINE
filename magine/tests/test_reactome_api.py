import networkx as nx
import pandas as pd
import pygraphviz as pyg

import \
    magine.network_database_expansions.reactome_expansion.get_apoptosis_example as reactome_api


def test_reactome():
    save_name = 'test_apoptosis'

    apoptosis = 109581

    g = pyg.AGraph(directed=True)
    g.graph_attr['rankdir'] = 'LR'
    g.graph_attr['splines'] = 'true'

    x = reactome_api.reactome.get_pathway_events(apoptosis)

    pathways = reactome_api.extract_pathways(x)
    pathways = [111457]
    print("Pathways found: {}".format(pathways))

    reactions = []
    all_info = []

    def extract_list_of_events(events):
        for i in events[:10]:
            info, rxn = reactome_api.get_reaction_info(i, g)
            if info is None:
                continue
            all_info.append(info)
            reactions.append(rxn)

    #
    for i in pathways:
        x = reactome_api.reactome.get_pathway_events(i)
        extract_list_of_events(x)

    df = pd.DataFrame(all_info)

    df.to_csv('reactome_df_{}.csv'.format(save_name))
    # """
    # df = pd.read_csv('reactome_df.csv')

    translocations = df[df['compartment'].map(len) > 1]
    if len(translocations) > 0:
        print("Translocations found")
        print(translocations)
        print(translocations.shape)

    gml_g = nx.nx_agraph.from_agraph(g)
    nx.write_gml(gml_g, '{}.gml'.format(save_name))
    print("Created GML file!")
    species = []
    for i in g.nodes():
        if i not in reactions:
            species.append(i)

    g.write('{}.dot'.format(save_name), )
    g.draw('{}.png'.format(save_name), prog='dot')


if __name__ == '__main__':
    test_reactome()
