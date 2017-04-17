# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 20:36:51 2015

@author: pinojc
"""
from bioservices import *

from magine.networks.kgml_to_networkx_parser import kgml_to_graph

kegg = KEGG(verbose=True)
kegg.TIMEOUT = 100

arrowType = {'activation': 'onormal', 'indirect effect': 'odiamondodiamond', \
             'expression': 'normal', 'inhibition': 'tee', \
             'binding/association': 'curve', 'phosphorylation': 'dot', \
             'missing interaction': 'odiamond', 'compound': 'dotodot', \
             'dissociation': 'diamond', 'ubiquitination': 'oldiamond', \
             'state change': 'teetee', 'dephosphorylation': 'onormal', \
             'repression': 'obox'}

colorType = {'activation': 'black', 'indirect effect': 'black', \
             'expression': 'black', 'inhibition': 'red', \
             'binding/association': 'green', 'phosphorylation': 'blue', \
             'missing interaction': 'pink', 'compound': 'pink', \
             'dissociation': 'pink', 'ubiquitination': 'red', \
             'state change': 'red', 'dephosphorylation': 'red', \
             'repression': 'red'}

shapeType = {'gene': 'oval', 'compound': 'box'}

"""
ECrel	enzyme-enzyme relation, indicating two enzymes catalyzing successive reaction steps
PPrel	protein-protein interaction, such as binding and modification
GErel	gene expression interaction, indicating relation of transcription factor and target gene product
PCrel	protein-compound interaction
maplink link to another map]

"""


def translate(X, Y, option, param_counter, parameters, rules, gene_monomers,
              compound_monomers):
    X = str(X)
    Y = str(Y)
    split = option.split("_")
    print(split)
    if split[0] == 'inhibition':
        rate1 = "kf%s" % str(param_counter)
        rate2 = "kr%s" % str(param_counter)
        rate3 = "kc%s" % str(param_counter)
        if X in gene_monomers:
            enzyme = '%s(gene="protein",state="A")' % X
            sub = '%s(gene="protein",state="A")' % Y
            prod = '%s(gene="protein",state="I")' % Y
        else:
            enzyme = '%s(bf=None)' % X
            sub = '%s(gene="protein",state="A")' % Y
            prod = '%s(gene="protein",state="I")' % Y

        if Y in compound_monomers:
            enzyme = '%s(bf=None)' % Y
            sub = '%s(gene="protein",state="A")' % X
            prod = '%s(gene="protein",state="I")' % X

        tmp1 = 'catalyze(%s,%s,%s,[%s,%s,%s] )\n' % (
            enzyme, sub, prod, rate1, rate2, rate3)
        rules += tmp1
        parameters += 'Parameter("%s",1)\n' % rate1
        parameters += 'Parameter("%s",1)\n' % rate2
        parameters += 'Parameter("%s",1)\n' % rate3

    elif split[0] == 'expression':

        rate1 = "kf%s" % str(param_counter)
        rate2 = "kr%s" % str(param_counter)
        rate3 = "kc%s" % str(param_counter)
        rate4 = "k_synth%s" % str(param_counter)
        rule_name = '"%s_expression_%s"' % (Y, str(param_counter))
        tmp1 = 'catalyze(%s(gene="protein"),%s(gene="off"),%s(gene="on",state="A"),[%s,%s,%s] )\n' % (
            X, Y, Y, rate1, rate2, rate3)
        tmp2 = 'Rule(%s,%s(bf=None,gene="on",state="A") !> %s(bf=None,gene="on",state="A") + %s(bf=None,gene="protein",state="A"),%s)\n' % (
            rule_name, Y, Y, Y, rate2)
        rules += tmp1
        rules += tmp2
        parameters += 'Parameter("%s", 1)\n' % rate1
        parameters += 'Parameter("%s", 1)\n' % rate2
        parameters += 'Parameter("%s", 1)\n' % rate3
        parameters += 'Parameter("%s", 1)\n' % rate4

    elif split[0] == 'activation':

        rate1 = "kf%s" % str(param_counter)
        rate2 = "kr%s" % str(param_counter)
        rate3 = "kc%s" % str(param_counter)

        if X in gene_monomers:
            enzyme = '%s(gene="protein", state="A")' % X
            sub = '%s(gene="protein", state="I")' % Y
            prod = '%s(gene="protein", state="A")' % Y
        else:

            enzyme = '%s(bf=None)' % X
            sub = '%s(gene="protein", state="I")' % Y
            prod = '%s(gene="protein", state="A")' % Y
        if Y in compound_monomers:
            enzyme = '%s(bf=None)' % Y
            sub = '%s(gene="protein", state="I")' % X
            prod = '%s(gene="protein", state="A")' % X

        tmp1 = 'catalyze(%s, %s, %s, [%s, %s, %s] )\n' % (
            enzyme, sub, prod, rate1, rate2, rate3)

        rules += tmp1
        parameters += 'Parameter("%s", 1)\n' % rate1
        parameters += 'Parameter("%s", 1)\n' % rate2
        parameters += 'Parameter("%s", 1)\n' % rate3

    elif split[0] == 'binding/association':

        rate1 = "kf%s" % str(param_counter)
        rate2 = "kr%s" % str(param_counter)
        rule_name = '"%s_binds_%s_%s"' % (X, Y, str(param_counter))

        tmp1 = 'Rule(%s, %s(bf=None, gene="protein") +%s(bf=None,gene="protein")!> %s(bf=1,gene="protein") ** %s(bf=1,gene="protein"),%s,%s)\n' % (
            rule_name, X, Y, X, Y, rate1, rate2)

        tmp1 = tmp1.replace("**", "%")
        rules += tmp1
        parameters += 'Parameter("%s",1)\n' % rate1
        parameters += 'Parameter("%s",1)\n' % rate2

    else:
        print(split[0], 'not found')
    return rules, parameters

    # catalyze(enzyme, e_site, substrate, s_site, product, klist):


# synthesize(A(x=None, y='e'), 1e-4)

if __name__ == '__main__':
    # nodes, edges = KGML2Graph("hsa04210.xml")
    g, x, y = kgml_to_graph("hsa04210.xml", output_dir='.')
    print(g, x, y)
    print('made it here')
    nodes = g.nodes(data=True)
    print(nodes)
    edges = g.edges(data=True)
    pysb_file = ''
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
