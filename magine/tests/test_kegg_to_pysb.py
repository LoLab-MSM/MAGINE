import os

from magine.networks.databases.kgml_to_networkx_parser import kgml_to_graph
from magine.networks.pysb_conversion.to_pysb import translate

rel_dir = os.path.dirname(__file__)

def test_kegg_to_pysb():
    # nodes, edges = KGML2Graph("hsa04210.xml")
    path = os.path.join(rel_dir, 'KEGG')
    g, x, y = kgml_to_graph("hsa04071.xml", output_dir=path)
    nodes = g.nodes(data=True)
    edges = g.edges(data=True)
    gene_monomers = []
    compound_monomers = []

    output = """
import pysb.macros as macros
from pysb import *
def catalyze(enz, sub, product, klist):
    return macros.catalyze(enz, 'bf', sub, 'bf', product, klist)
Model()
"""
    added = set()
    monomers = ''
    initials = ''
    counter = 0
    rules = ''
    parameters = ""
    for i, info in nodes:
        if i in added:
            continue
        added.add(i)
        i = i.replace(':', '')
        if info['speciesType'] == 'gene':
            gene_monomers.append(i)

            monomers += 'Monomer("{}", ["gene", "state", "bf"], ' \
                        '{{"state": ["I", "A"], ' \
                        '"gene": ["on", "off", "protein"]}})\n'.format(i)
            initials += 'Initial({0}(bf=None, gene="on", state="A"), ' \
                        '{0}_0)\n'.format(i)
            initials += 'Initial({0}(bf=None, gene="protein", state="A"), ' \
                        '{0}_0 )\n'.format(i)
            initials += 'Initial({0}(bf=None, gene="on", state="I"), ' \
                        '{0}_0 )\n'.format(i)
            initials += 'Initial({0}(bf=None, gene="off", state="I"), ' \
                        '{0}_0 )\n'.format(i)
            # print 'Monomer("%s",["gene","bs"],{"gene":["on","off"]})'% i[0]

        elif info['speciesType'] == 'compound':
            compound_monomers.append(i)
            # print 'Monomer("%s",["bf"])'% i[0].strip("+")
            monomers += 'Monomer("{0}", ["bf"])\n'.format(i)
            initials += 'Initial({0}(bf=None), {0}_0 )\n'.format(i)
        parameters += 'Parameter("{0}_0", 1)\n'.format(i)

    for info in edges:
        s_1 = info[0]
        s_2 = info[1]
        s_1 = s_1.replace(':', '')
        s_2 = s_2.replace(':', '')
        int_type = info[2]['interactionType']

        rules, parameters = translate(s_1, s_2, int_type, counter, parameters,
                                      rules, gene_monomers, compound_monomers)
        counter += 1
    output += monomers
    output += parameters
    output += initials
    output += rules
    with open('pysb_model.py', 'w') as File:
        File.write(output)
