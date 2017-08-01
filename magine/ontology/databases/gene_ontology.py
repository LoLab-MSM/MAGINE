import gzip
import os
from collections import defaultdict

import requests
import wget
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
from statsmodels.stats.proportion import binom_test

from magine.data.storage import id_mapping_dir
from magine.data.storage import network_data_dir

try:
    import cPickle as pickle
except:
    import pickle


class MagineGO(object):
    def __init__(self, species='hsa'):
        dirname = network_data_dir
        gene_to_go_name = os.path.join(dirname,
                                       '{}_gene_to_go.p'.format(species))
        go_to_gene_name = os.path.join(dirname,
                                       '{}_goids_to_genes.p'.format(species))
        go_to_go_name = os.path.join(dirname,
                                     '{}_goids_to_goname.p'.format(species))
        go_depth = os.path.join(dirname, '{}_godepth.p'.format(species))
        go_aspect = os.path.join(dirname, '{}_go_aspect.p'.format(species))
        for i in [gene_to_go_name, go_to_gene_name, go_to_go_name, go_depth,
                  go_aspect]:
            if not os.path.exists(i):
                download_and_process_go(species=species)

        self.gene_to_go = pickle.load(open(gene_to_go_name, 'rb'))
        self.go_to_gene = pickle.load(open(go_to_gene_name, 'rb'))
        self.go_to_name = pickle.load(open(go_to_go_name, 'rb'))
        self.go_depth = pickle.load(open(go_depth, 'rb'))
        self.go_aspect = pickle.load(open(go_aspect, 'rb'))

    def calculate_enrichment(self, genes, reference=None, evidence_codes=None,
                             aspect=None, use_fdr=True):
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


def download_and_process_go(species='hsa'):
    print("Creating GO files")
    from goatools import obo_parser
    obo_file = os.path.join(id_mapping_dir, 'go.obo')
    if not os.path.exists(obo_file):
        download_current_go()
    go = obo_parser.GODag(obo_file)
    gene_to_go, go_to_gene, goid_to_name = create_annotations()

    go_aspect = dict()
    go_depth = dict()

    dirname = network_data_dir

    go_to_gene_name = os.path.join(dirname,
                                   '{}_goids_to_genes.p'.format(species))
    go_to_go_name = os.path.join(dirname,
                                 '{}_goids_to_goname.p'.format(species))
    gene_to_go_name = os.path.join(dirname, '{}_gene_to_go.p'.format(species))
    go_depth_name = os.path.join(dirname, '{}_godepth.p'.format(species))
    go_aspect_name = os.path.join(dirname, '{}_go_aspect.p'.format(species))
    for i in go_to_gene.keys():
        go_depth[i] = go[i].depth
        go_aspect[i] = go[i].namespace

    pickle.dump(go_to_gene, open(go_to_gene_name, 'wb'))
    pickle.dump(goid_to_name, open(go_to_go_name, 'wb'))
    pickle.dump(go_depth, open(go_depth_name, 'wb'))
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


def download_ncbi_gene_file(out_dir, force_dnld=False):
    """Download a file from NCBI Gene's ftp server."""
    out_name = os.path.join(out_dir, "gene2go")
    if not os.path.exists(out_name) or force_dnld:
        tmp_out = os.path.join(out_dir, 'tmp.gz')

        if os.path.exists(tmp_out):
            os.remove(tmp_out)

        if os.path.exists(out_name):
            os.remove(out_name)
        fin_ftp = "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"

        wget.download(fin_ftp, tmp_out)

        with gzip.open(tmp_out, 'rb') as zstrm:
            print("\n  READ:  {F}\n".format(F=tmp_out))
            with open(out_name, 'wb') as ostrm:
                ostrm.write(zstrm.read())
                print("  WROTE: {F}\n".format(F=out_name))


def read_ncbi_gene2go(fin_gene2go, taxids=None, **kws):
    """Read NCBI's gene2go. Return gene2go data for user-specified taxids."""
    from magine.mappings.gene_mapper import GeneMapper
    gm = GeneMapper()

    id2gos = defaultdict(set)
    go2term = dict()
    go2genes = defaultdict(set)

    evs = kws.get('evidence_set', None)
    if taxids is None:  # Default taxid is Human
        taxids = {9606}
    with open(fin_gene2go) as ifstrm:
        for line in ifstrm:
            if line[0] != '#':  # Line contains data. Not a comment
                line = line.rstrip()  # chomp
                flds = line.split('\t')
                if len(flds) >= 5:
                    taxid_curr, geneid, go_id, evidence, qualifier, go_term = flds[
                                                                              :6]
                    taxid_curr = int(taxid_curr)
                    if taxid_curr in taxids and qualifier != 'NOT' and evidence != 'ND':

                        if evs is None or evidence in evs:
                            geneid = int(geneid)
                            symbol = gm.ncbi_to_symbol[geneid]
                            if len(symbol) == 1:
                                symbol = symbol[0]
                            else:
                                print(symbol)
                            go2genes[go_id].add(symbol)
                            id2gos[symbol].add(go_id)
                            go2term[go_id] = go_term

    return id2gos, go2genes, go2term


def download_current_go(redownload=False):
    go_url = 'http://purl.obolibrary.org/obo/go.obo'
    target_file = 'go.obo'
    out_path = os.path.join(id_mapping_dir, target_file)
    if not os.path.exists(out_path) or redownload:
        print("GO exists, default behavior is to re-download.")

        r = requests.get(go_url, stream=True)
        response = requests.head(go_url)
        print(response.headers)
        # file_size = int(response.headers['content-length'])
        print("Downloading ontology file")
        file_size_dl = 0
        block_sz = 1024
        # block_sz = 8024

        with open(out_path, 'wb') as f:
            for chunk in r.iter_content(chunk_size=block_sz):
                file_size_dl += len(chunk)

                if chunk:  # filter out keep-alive new chunks
                    f.write(chunk)

        print("Downloaded {} and stored {}".format(go_url, out_path))


def download_annotations(redownload=False):
    annotations_file = 'http://geneontology.org/gene-associations/goa_human.gaf.gz'
    annotations_target_file = 'goa_human.gaf.gz'
    out_path = os.path.join(id_mapping_dir, annotations_target_file)
    real_path = os.path.join(id_mapping_dir, 'annotations.gaf')
    if not os.path.exists(out_path) or redownload:
        r = requests.get(annotations_file, stream=True)
        response = requests.head(annotations_file)
        print(response.headers)
        # file_size = int(response.headers['content-length'])
        print("Downloading ontology file")
        file_size_dl = 0
        block_sz = 1024
        # block_sz = 8024

        with open(out_path, 'wb') as f:
            for chunk in r.iter_content(chunk_size=block_sz):
                file_size_dl += len(chunk)

                if chunk:  # filter out keep-alive new chunks
                    f.write(chunk)

        with gzip.open(out_path, 'rb') as f:
            file_content = f.read()
        with open(real_path, 'wb') as f:
            f.write(file_content)


def create_annotations():
    out_path = os.path.join(id_mapping_dir, 'gene2go')
    download_ncbi_gene_file(id_mapping_dir, force_dnld=True)
    id2gos, go2genes, go2term = read_ncbi_gene2go(out_path, [9606])

    return id2gos, go2genes, go2term


if __name__ == '__main__':
    # create_annotations()
    download_and_process_go()
    go = MagineGO('hsa')
    terms = go.calculate_enrichment(['BAX'],
                                    reference=['BAX', "LL", 'BCL2', 'CASP1'])
    #    'GO:0001569', (['BAX'], 0.78977569118414181, 1))
    for t in terms:
        print(t, terms[t])
