import os
from collections import defaultdict

import pandas as pd
import requests
from magine.enrichment.ontology_analysis import MagineGO

from magine.data.storage import id_mapping_dir
from magine.data.storage import network_data_dir

try:
    import cPickle as pickle
except:
    import pickle


def download_and_process_go(species='hsa'):
    print("Creating GO files")
    from goatools import obo_parser
    obo_file = os.path.join(id_mapping_dir, 'go.obo')
    if not os.path.exists(obo_file):
        download_current_go()
    go = obo_parser.GODag(obo_file)
    gene_to_go, go_to_gene, goid_to_name = download_ncbi_gene_file()

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


def download_ncbi_gene_file(tax_ids=None):
    """ Downloads gene2go associations files from ncbi

    Parameters
    ----------
    tax_ids : list
        list of tax ids to consider, default [9606], human
    force_dnld

    Returns
    -------

    """

    from magine.mappings.gene_mapper import GeneMapper
    gm = GeneMapper()

    gene2go_url = "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"
    print("Downloading gene2go assocations file")
    r = pd.read_table(gene2go_url, compression='gzip')

    id2gos = defaultdict(set)
    go2term = dict()
    go2genes = defaultdict(set)

    if tax_ids is None:  # Default taxid is Human
        tax_ids = {9606}
    for line in r.values:
        t_id, gene_id, go_id, evidence, qual, go_term = line[:6]

        t_id = int(t_id)
        if t_id in tax_ids and qual != 'NOT' and evidence != 'ND':
            gene_id = int(gene_id)
            symbol = gm.ncbi_to_symbol[gene_id]
            if len(symbol) == 1:
                symbol = symbol[0]
            else:
                print(symbol)
            go2genes[go_id].add(symbol)
            id2gos[symbol].add(go_id)
            go2term[go_id] = go_term

    return id2gos, go2genes, go2term


def download_current_go(redownload=True):
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


if __name__ == '__main__':
    # create_annotations()
    download_and_process_go()
    go = MagineGO('hsa')
    terms = go.calculate_enrichment(['BAX'],
                                    reference=['BAX', "LL", 'BCL2', 'CASP1'])
    #    'GO:0001569', (['BAX'], 0.78977569118414181, 1))
    for t in terms:
        print(t, terms[t])
