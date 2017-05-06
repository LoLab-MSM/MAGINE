import os

from statsmodels.sandbox.stats.multicomp import fdrcorrection0
from statsmodels.stats.proportion import binom_test

try:
    import cPickle as pickle
except:
    import pickle

dirname = os.path.dirname(__file__)


class MagineGO(object):
    def __init__(self, species='hsa'):
        gene_to_go_name = os.path.join(dirname,
                                       '{}_gene_to_go.p'.format(species))
        go_to_gene_name = os.path.join(dirname,
                                       '{}_goids_to_genes.p'.format(species))
        go_to_go_name = os.path.join(dirname,
                                     '{}_goids_to_goname.p'.format(species))
        go_depth = os.path.join(dirname, '{}_godepth.p'.format(species))
        go_slim = os.path.join(dirname, '{}_go_slims.p'.format(species))
        go_aspect = os.path.join(dirname, '{}_go_aspect.p'.format(species))
        for i in [gene_to_go_name, go_to_gene_name, go_to_go_name, go_depth,
                  go_aspect]:
            if not os.path.exists(i):
                create_dicts_through_orange(species=species)

        self.gene_to_go = pickle.load(open(gene_to_go_name, 'rb'))
        self.go_to_gene = pickle.load(open(go_to_gene_name, 'rb'))
        self.go_to_name = pickle.load(open(go_to_go_name, 'rb'))
        self.go_depth = pickle.load(open(go_depth, 'rb'))
        self.go_slim = pickle.load(open(go_slim, 'rb'))
        self.go_aspect = pickle.load(open(go_aspect, 'rb'))

    def calculate_enrichment(self, genes, reference=None, evidence_codes=None,
                             aspect=None,
                             use_fdr=True, ):
        """
    
        Parameters
        ----------
        genes : list
            list of genes
        reference : list
            reference list of species to calculate enrichment
        evidence_codes : list
            GO evidence codes
        use_fdr : bool
            Correct for multiple hypothesis testing
    
        Returns
        -------
    
        """

        # TODO check for alias for genes
        genes = set(genes)
        # TODO add aspects
        term_reference = self.go_to_gene.keys()
        aspect_dict = {
            'P': 'biological_process',
            'C': 'cellular_component',
            'F': 'molecular_function'
        }
        if aspect is None:
            term_reference = self.go_to_gene
            gene_reference = self.gene_to_go
        else:
            term_reference = dict()
            gene_reference = dict()

        if aspect is not None:
            for i in aspect:
                if i not in ['P', 'C', 'F']:
                    print("Error: Aspects are only 'P', 'C', and 'F' \n")
                    quit()
            for i in ['P', 'C', 'F']:
                if i in aspect:
                    term_reference = None

        # TODO add reference
        if reference:
            # TODO check for reference alias
            reference = set(reference)
            reference.intersection_update(set(self.gene_to_go.keys()))
        else:
            reference = set(self.gene_to_go.keys())

        # TODO add evidence_codes

        terms = set()
        for i in genes:
            if i in self.gene_to_go:
                for t in self.gene_to_go[i]:
                    terms.add(t)

        n_genes = len(genes)
        n_ref = float(len(reference))
        res = {}
        for term in terms:

            all_annotated_genes = set(self.go_to_gene[term])
            mapped_genes = genes.intersection(all_annotated_genes)
            n_mapped_genes = len(mapped_genes)

            if n_ref > len(all_annotated_genes):
                mapped_reference_genes = \
                    reference.intersection(all_annotated_genes)
            else:
                mapped_reference_genes = \
                    all_annotated_genes.intersection(reference)

            n_mapped_ref = len(mapped_reference_genes)

            prob = float(n_mapped_ref) / n_ref

            p_value = binom_test(n_mapped_genes, n_genes, prob, 'larger')

            res[term] = ([i for i in mapped_genes], p_value, n_mapped_ref)
        if use_fdr:
            res = sorted(res.items(), key=lambda x: x[1][1])
            fdr = fdrcorrection0([p for _, (_, p, _) in res], is_sorted=True)
            values = fdr[1]
            res = dict([(index, (genes, p, ref))
                        for (index, (genes, _, ref)), p in zip(res, values)])
        return res


def create_dicts_through_orange(species='hsa', rev=None, rev_ass=None):
    print("Creating GO files")
    from orangecontrib.bio import go
    ont = go.Ontology(rev=rev)
    ann = go.Annotations(species, ontology=ont, rev=rev_ass)
    go_to_gene = dict()
    gene_to_go = dict()
    goid_to_name = dict()
    go_aspect = dict()
    go_depth = dict()
    go_slims = ont.named_slims_subset('goslim_pir')

    go_to_gene_name = os.path.join(dirname,
                                   '{}_goids_to_genes.p'.format(species))
    go_to_go_name = os.path.join(dirname,
                                 '{}_goids_to_goname.p'.format(species))
    gene_to_go_name = os.path.join(dirname, '{}_gene_to_go.p'.format(species))
    go_depth_name = os.path.join(dirname, '{}_godepth.p'.format(species))
    go_slim_name = os.path.join(dirname, '{}_go_slims.p'.format(species))
    go_aspect_name = os.path.join(dirname, '{}_go_aspect.p'.format(species))
    for i in ont.terms:
        go_to_gene[i] = ann.get_all_genes(i)
        goid_to_name[i] = ont[i].name
        go_depth[i] = ont.term_depth(i)
        go_aspect[i] = ont[i].namespace

    pickle.dump(go_to_gene, open(go_to_gene_name, 'wb'))
    pickle.dump(goid_to_name, open(go_to_go_name, 'wb'))
    pickle.dump(go_depth, open(go_depth_name, 'wb'))
    pickle.dump(go_slims, open(go_slim_name, 'wb'))
    pickle.dump(go_aspect, open(go_aspect_name, 'wb'))

    for i in go_to_gene:
        term = i
        genes = go_to_gene[i]
        for g in genes:
            if g in gene_to_go:
                gene_to_go[g].add(term)
            else:
                gene_to_go[g] = set()
                gene_to_go[g].add(term)
    pickle.dump(gene_to_go, open(gene_to_go_name, 'wb'))
    print("Done creating GO files")


if __name__ == '__main__':
    create_dicts_through_orange()
    go = MagineGO('hsa')
    terms = go.calculate_enrichment(['BAX'],
                                    reference=['BAX', "LL", 'BCL2', 'CASP1'])
    # 'GO:0001569', (['BAX'], 0.78977569118414181, 1))
    for t in terms:
        print(t, terms[t])
