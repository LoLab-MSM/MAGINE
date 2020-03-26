import logging
import os

import networkx as nx

from magine.data.storage import network_data_dir
from magine.logging import get_logger

logger = get_logger(__name__, log_level=logging.INFO)


def load_hmdb_network(fresh_download=False):
    """ Create HMDB network containing all metabolite-protein interactions

    Parameters
    ----------
    fresh_download : bool
        Download fresh copy from HMDB
    verbose : bool

    Returns
    -------
    nx.DiGraph
    """
    out_name = os.path.join(network_data_dir, 'hmdb_graph.p.gz')

    if not fresh_download and os.path.exists(out_name):
        tmp_graph = nx.read_gpickle(out_name)
    else:
        from magine.mappings.chemical_mapper import ChemicalMapper

        cm = ChemicalMapper()

        tmp_graph = nx.DiGraph()

        def _add_node(node, node_type):
            attrs = {'databaseSource': 'HMDB', 'speciesType': node_type}
            if node_type == 'compound':
                if node in cm.hmdb_to_chem_name:
                    attrs['chemName'] = sorted(cm.hmdb_to_chem_name[node])[0]
            tmp_graph.add_node(node, **attrs)

        for source, genes in cm.hmdb_main_to_protein.items():
            if source == '':
                continue
            _add_node(source, 'compound')
            for target in genes:
                if target == '':
                    continue
                _add_node(target, 'gene')
                tmp_graph.add_edge(source, target, interactionType='chemical',
                                   databaseSource='HMDB')
        nx.write_gpickle(tmp_graph, out_name)

    logger.info("HMDB : {} nodes and {} edges".format(
        len(tmp_graph.nodes), len(tmp_graph.edges))
    )

    return tmp_graph


if __name__ == '__main__':
    load_hmdb_network(True, True)
