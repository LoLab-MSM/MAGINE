import itertools
import os
from random import randint

import networkx as nx

from magine.networks.network_tools import export_to_dot, networkx_to_igraph


class OntologyNetworkGenerator(object):
    """ Generates ontology network from molecular networks.

    Nodes are the ontology terms
    Edges between nodes are determined based on molecular network graphs.




    """

    def __init__(self, molecular_network=None):

        self.mol_network = molecular_network
        self._nodes = None
        self._edges = None
        self.molecular_network = None

    @property
    def edges(self):
        if self._edges is None:
            self._edges = set(self.mol_network.edges())
        return self._edges

    @property
    def nodes(self):
        if self._nodes is None:
            self._nodes = set(self.mol_network.nodes())
        return self._nodes

    def _count_neighbors(self, term_a, term_b):
        """
        Calculate the number of direct edges between species of two terms


        Parameters
        ----------
        term_a : list_like
            list of species
        term_b : list_like
            list of species

        Returns
        -------
        int, int, list_like
            number of edges from A to B
            number of edges from B to A
            genes responsible for edges
        """
        term_a = set(term_a)
        term_b = set(term_b)
        genes_in_both_go = set()
        # calculate edges between A and B
        a_to_b, gene_in = self._determine_edges(term_1=term_a, term_2=term_b)
        genes_in_both_go.update(gene_in)

        # calculate edges between B and A
        b_to_a, gene_in = self._determine_edges(term_1=term_b, term_2=term_a)
        genes_in_both_go.update(gene_in)

        return a_to_b, b_to_a, genes_in_both_go

    def _determine_edges(self, term_1, term_2):
        """
        calculate the number of neighbors that connect between two terms

        Parameters
        ----------
        term_1
        term_2

        Returns
        -------

        """
        counter = 0
        genes_in_go = set()
        term_1_good = {i for i in term_1 if i in self.nodes}
        term_2_good = {i for i in term_2 if i in self.nodes}

        for i in term_1_good:
            for j in term_2_good:
                if i != j:
                    if (i, j) in self.edges:
                        genes_in_go.add(i)
                        genes_in_go.add(j)
                        counter += 1
        return counter, genes_in_go

    def create_network_from_list(self, list_of_ontology_terms,
                                 ont_to_species_dict, ont_to_label_dict,
                                 save_name=None,
                                 draw=False, threshold=0, out_dir=None,
                                 merge_edges=False):
        """
        Creates a GO level network from list of GO terms

        Parameters
        ----------

        list_of_ontology_terms : list_list
            list of GO terms
        ont_to_species_dict : dict
            dictionary, where keys are ontology terms and values are species
            for that term
        ont_to_label_dict: dict
            dictionary where keys are ontology terms and values are labels
        save_name : str
            name to save network
        draw : bool
            create a go_graph of network
        threshold : int
            integer threshold of number of neighbors between two GO terms
            default = 0
        out_dir : str
            output directory
        merge_edges : bool
            merge the edges between GO nodes

        Returns
        -------
        networkx.DiGraph

        """
        if out_dir is not None:
            out_path = os.path.join(out_dir, 'Network_files')
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
                os.mkdir(out_path)
        else:
            out_path = '.'

        # make sure a network exists
        if self.mol_network is None:
            print("Must provide a network! Returning None")
            return None
        go_graph = nx.DiGraph()
        mol_net = nx.DiGraph()
        all_genes = set()
        list_of_go_terms = set(list_of_ontology_terms)
        sp_to_term = dict()
        sp_to_label = dict()

        for i in itertools.combinations(list_of_go_terms, 2):
            term1 = i[0]
            term2 = i[1]
            term_1 = set(ont_to_species_dict[term1])
            term_2 = set(ont_to_species_dict[term2])
            label_1 = ont_to_label_dict[term1]
            label_2 = ont_to_label_dict[term2]

            a_to_b, b_to_a, genes_in_edges = self._count_neighbors(term_1,
                                                                   term_2)

            # add to graph if at least one edge found between terms
            if a_to_b or b_to_a:
                x = self.mol_network.subgraph(genes_in_edges)
                mol_net.add_edges_from(x.edges(data=True))
                for gene in genes_in_edges:
                    if gene in sp_to_term:
                        if gene in term_1:
                            sp_to_term[gene].add(term1)
                            sp_to_label[gene].add(label_1)
                        if gene in term_2:
                            sp_to_term[gene].add(term2)
                            sp_to_label[gene].add(label_2)
                    else:
                        if gene in term_1:
                            sp_to_term[gene] = {term1}
                            sp_to_label[gene] = {label_1}
                        if gene in term_2:
                            sp_to_term[gene] = {term2}
                            sp_to_label[gene] = {label_2}

                all_genes.update(genes_in_edges)
            else:
                print('No edges between {} and {}'.format(label_1, label_2))
            go_graph.add_node(label_1, term=term1, label=label_1)
            go_graph.add_node(label_2, term=term2, label=label_2)

            if merge_edges:
                if a_to_b > threshold or b_to_a > threshold:
                    if a_to_b > b_to_a:
                        go_graph.add_edge(label_2, label_1,
                                          label=str(a_to_b + b_to_a),
                                          weight=a_to_b + b_to_a, dir='both',
                                          weightAtoB=b_to_a, weightBtoA=a_to_b)
                    else:
                        go_graph.add_edge(label_1, label_2,
                                          label=str(a_to_b + b_to_a),
                                          weight=a_to_b + b_to_a, dir='both',
                                          weightAtoB=a_to_b, weightBtoA=b_to_a)
            else:
                if a_to_b > threshold:
                    go_graph.add_edge(label_1, label_2, label=str(a_to_b),
                                      weight=a_to_b)

                if b_to_a > threshold:
                    go_graph.add_edge(label_2, label_1, label=str(b_to_a),
                                      weight=b_to_a)
        # print(sp_to_term)
        for i in sp_to_term:
            # print(i)
            labels = sp_to_label[i]
            terms = sp_to_term[i]
            assert len(labels) == len(terms), \
                'len(labels) should equal len(terms)'
            mol_net.node[i]['termName'] = ','.join(i for i in sorted(labels))
            mol_net.node[i]['terms'] = ','.join(i for i in sorted(terms))

        self.molecular_network = mol_net
        if save_name is not None:
            if out_dir is not None:
                save_name = os.path.join(out_path, save_name)

            nx.write_gml(mol_net, '{}_subgraph.gml'.format(save_name))
            nx.write_gml(go_graph, '{0}.gml'.format(save_name))

            if draw:
                export_to_dot(go_graph, save_name)
                # export_to_dot(mol_net, save_name + '_subgraph')

                plot(mol_net, save_name + '_subgraph_igraph')

        return go_graph, mol_net


def plot(mol_net, save_name):
    try:
        import cairo
    except ImportError:
        print("Please install pycairo to use igraph plotting")
        return

    try:
        import igraph
        from igraph.drawing.text import TextDrawer
        from igraph.drawing.colors import color_name_to_rgba
    except ImportError:
        print("No igraph, cannot use plotting function")
        return

    g = networkx_to_igraph(mol_net)
    print(g)
    cl = igraph.VertexClustering(g).FromAttribute(g, attribute='termName')
    membership = cl.membership

    if membership is not None:
        gcopy = g.copy()
        edges = []
        edges_colors = []
        for edge in g.es():
            if membership[edge.tuple[0]] != membership[edge.tuple[1]]:
                # edges.append(edge)
                # edges_colors.append("gray")
                edges_colors.append("black")
            else:
                edges_colors.append("black")
        gcopy.delete_edges(edges)
        layout = gcopy.layout("kk")
        # layout = gcopy.layout("tree")
        # layout = gcopy.layout("drl")
        # layout = gcopy.layout("fr")
        g.es["color"] = edges_colors
    visual_style = dict()
    visual_style["vertex_label_dist"] = 0
    visual_style["vertex_shape"] = "circle"
    visual_style["edge_color"] = g.es["color"]
    # visual_style["bbox"] = (4000, 2500)
    visual_style["vertex_size"] = 50
    visual_style["layout"] = layout
    visual_style["bbox"] = (2000, 2000)
    visual_style["margin"] = 100
    # visual_style["node_label"] = g.vs["label"]
    for vertex in g.vs():
        vertex["label"] = vertex['name']

    colors = []
    for i in range(0, max(membership) + 1):
        # colors.append(k_colors[randint(0, len(k_colors))])
        colors.append('#%06X' % randint(0, 0xFFFFFF))
    _groups = dict()
    for vertex in g.vs():
        _groups[vertex.index] = colors[membership[vertex.index]]
        vertex["color"] = colors[membership[vertex.index]]
    _mark_groups = dict()
    # print(_groups)
    for i in _groups:
        co = _groups[i]
        if co in _mark_groups:
            _mark_groups[co].append(i)
        else:
            _mark_groups[co] = [i]
    mark_groups = dict()
    for i in _mark_groups:
        mark_groups[tuple(_mark_groups[i])] = i
    # print(_mark_groups)
    # print(mark_groups)

    visual_style["vertex_color"] = g.vs["color"]

    # igraph.plot(cl, "{}.pdf".format(save_name), mark_groups=True, **visual_style)
    plot = igraph.Plot("{}2.png".format(save_name), bbox=(3000, 2200),
                       background="white")
    plot.add(cl, mark_groups=mark_groups, **visual_style)
    # plot = igraph.plot(cl, mark_groups=True, **visual_style)
    plot.redraw()
    # Grab the surface, construct a drawing context and a TextDrawer
    ctx = cairo.Context(plot.surface)
    ctx.set_font_size(36)
    # drawer = TextDrawer(ctx, "Test title", halign=TextDrawer.CENTER)
    # drawer.draw_at(0, 40, width=600)
    labels = dict()
    for vertex in g.vs():
        labels[colors[membership[vertex.index]]] = vertex['termName'].replace(
            ',', '\n')

    for n, i in enumerate(colors):
        text_drawer = TextDrawer(ctx, text=labels[i], halign=TextDrawer.LEFT)
        text_drawer.draw_at(x=2100, y=100 + (n * 200), width=900, wrap=True)
        ctx.set_source_rgba(*color_name_to_rgba(i))
        ctx.arc(2000, 100 + (n * 200), 100 / 2, 0, 2 * 3.14)
        ctx.fill()
        ctx.set_source_rgba(0, 0, 0, 1)

    # Make the plot draw itself on the Cairo surface
    # Save the plot
    plot.save()
