import networkx as nx
import pandas as pd

url = 'http://reactomews.oicr.on.ca:8080/caBigR3WebApp2015/FIsInGene_031516_with_annotations.txt.zip'
count = 0
directions = set()
genes = set()
annotation = set()


def return_direction(text, gene1, gene2):
    if text in ("<-", "?-"):
        return gene2, gene1
    if text in ("->", "->"):
        return gene1, gene2
    else:
        return gene1, gene2


table = pd.read_table('FIsInGene_031516_with_annotations.txt')


def check_row(row):
    prediction = False
    interaction_type = None
    annotation = row['Annotation']
    modification = ''
    if 'predicted' in annotation:
        prediction = True
    if 'repress' in annotation:
        interaction_type = 'repression'
    elif 'express' in annotation:
        interaction_type = 'expression'
    elif 'reaction' in annotation:
        interaction_type = 'reaction'
    elif 'complex' in annotation:
        interaction_type = 'complex'
    elif 'catalyze' in annotation:
        interaction_type = 'catalyze'
    elif 'deactiv' in annotation:
        interaction_type = 'deactivate'
    elif 'activ' in annotation:
        interaction_type = 'activate'
    elif 'inhib' in annotation:
        interaction_type = 'inhibit'
    elif 'binding' in annotation:
        interaction_type = 'binding'
    elif 'input' in annotation:
        interaction_type = '?'
    elif 'indirect effect' in annotation:
        interaction_type = 'indirect'
    elif 'compound' in annotation:
        interaction_type = 'compound'
    elif 'state change' in annotation:
        interaction_type = 'state change'
    elif 'interaction' in annotation:
        interaction_type = 'binding'
    elif 'PPrel' in annotation:
        interaction_type = 'binding'
    if 'dephosphoryl' in annotation:
        modification = 'dephosphorylation'
        if interaction_type is None:
            interaction_type = 'dephosphorylation'
    elif 'phosphoryl' in annotation:
        modification = 'phosphorylation'
        if interaction_type is None:
            interaction_type = 'phosphorylation'
    elif 'ubiquiti' in annotation:
        modification = 'ubiquitination'
        if interaction_type is None:
            interaction_type = 'ubiquitination'
    if interaction_type is None:
        print(annotation)
        interaction_type = '?'
        if 'ECrel' in annotation:
            print(row)
    g1, g2 = return_direction(row['Direction'], row['Gene1'], row['Gene2'])
    return g1, g2, interaction_type, prediction, modification


g = nx.DiGraph()
added_genes = set()
for i, row in table.iterrows():
    score = row['Score']
    gene1, gene2, interaction, prediction, modification = check_row(row)
    if prediction:
        if gene1 not in added_genes:
            g.add_node(gene1, sourceDB='ReactomeFI', speciesType='gene')
            added_genes.add(gene1)
        if gene2 not in added_genes:
            g.add_node(gene2, sourceDB='ReactomeFI', speciesType='gene')
            added_genes.add(gene2)
        g.add_edge(gene1, gene2, interactionType=interaction, ptm=modification,
                   score=score, predection='true', sourceDB='ReactomeFI')
        # print(gene1, gene2, interaction, prediction, modification, score)
    elif gene1 == gene2:
        continue
    else:
        if gene1 not in added_genes:
            g.add_node(gene1, sourceDB='ReactomeFI', speciesType='gene')
            added_genes.add(gene1)
        if gene2 not in added_genes:
            g.add_node(gene2, sourceDB='ReactomeFI', speciesType='gene')
            added_genes.add(gene2)

        g.add_edge(gene1, gene2, interactionType=interaction, ptm=modification,
                   score=score, sourceDB='ReactomeFI')

        # print(gene1, gene2, interaction, prediction)
ge = set(g.nodes())
g1 = set(table['Gene1'].unique())
g2 = set(table['Gene2'].unique())
g_all = set(list(g1) + list(g2))
for i in g_all:
    if i not in ge:
        print(i)

print(len(g.nodes()))
print(len(g.edges()))
nx.write_gml(g, 'reactome_fi.gml')
quit()
with open('FIsInGene_031516_with_annotations.txt', 'r') as f:
    f.readline()
    alist = [line.rstrip() for line in f]
    for line in alist:
        if count > 100:
            continue
        Gene1, Gene2, Annotation, Direction, Score = line.split('\t')
        print(Gene1, Gene2, Annotation, Direction, Score)
        directions.add(Direction)
        annotation.add(Annotation)
        genes.add(Gene1)
        genes.add(Gene2)
        count += 1
for each in annotation:
    print(each)
for i in [directions, genes, annotation]:
    print(len(i))
