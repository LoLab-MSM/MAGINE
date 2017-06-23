from bioservices.services import REST
import pandas as pd


class Reactome(REST):
    """
    Reactome api modified from bioservices

    """

    _url = "http://reactomews.oicr.on.ca:8080/ReactomeRESTfulAPI/RESTfulWS"

    def __init__(self, verbose=True, cache=False):
        super(Reactome, self).__init__("Reactome(URL)", url=Reactome._url,
                                       verbose="ERROR", cache=cache)
        self.debugLevel = verbose
        self.test = 2

        # for buffering
        self._list_pathways = []
        self._content_url = 'http://reactome.org/ContentService/data/'

    def get_reaction_info(self, reaction):
        """Fetch information from the reaction HTML page

        .. note:: draft version
        """
        q = 'query/{}'.format(reaction)
        res = self.http_get(self._content_url + q, frmt='json')
        return res

    def get_entity_info(self, entity):
        """Fetch information from the reaction HTML page

        .. note:: draft version
        """
        q = 'query/{}/extended'.format(entity)
        res = self.http_get(self._content_url+q, frmt='json')
        return res

    def get_entity_attribute(self, entity, attribute):
        """Fetch information from the reaction HTML page

        .. note:: draft version
        """
        q = 'query/{}/{}'.format(entity, attribute)
        res = self.http_get(self._content_url+q, frmt='string')

        return res

    def get_event_participants(self, event):
        q = 'event/{}/participants'.format(event)
        res = self.http_get(self._content_url+q, frmt='json')
        return res

    def get_event_participating_phys_entities(self, event):
        q = 'event/{}/participatingPhysicalEntities'.format(event)
        res = self.http_get(self._content_url+q, frmt='json')
        return res

    def get_pathway(self, pathway_id):
        """Fetch information from the reaction HTML page

        .. note:: draft version
        """
        q = 'pathway/{}/Complex'.format(pathway_id)
        res = self.http_get(self._content_url +q, frmt='json')

        return res

    def get_pathway_events(self, pathway_id):
        """Fetch information from the reaction HTML page

        .. note:: draft version
        """
        q = 'pathway/{}/containedEvents'.format(pathway_id)
        res = self.http_get(self._content_url + q, frmt='json')
        return res


'''
    def _download_list_pathways(self):
        url = "http://www._reactome.org/download/current/ReactomePathways.txt"
        if len(self._list_pathways) == 0:
            res = self.session.get(url)
            if res.status_code == 200:
                res = res.text  # content does not work in python 3.3
                res = res.strip()
                self._list_pathways = [x.split("\t") for x in res.split("\n")]
            else:
                self.logging.error("could not fetch the pathways")
        return self._list_pathways

    def get_list_pathways(self):
        """Return list of pathways from _reactome website

        :return: list of lists. Each sub-lis contains 3 items: _reactome pathway
            identifier, description and species

        """
        res = self._download_list_pathways()
        return res

    def get_species(self):
        """Return list of species from all pathways"""
        res = set([x[2] for x in self.get_list_pathways()])
        return res

    def pathway_participants(self, identifier):
        """Get list of pathway participants for a pathway using

        :param int identifier: Pathway database identifier
        :return: list of fully encoded PhysicalEntity objects in the pathway
            (in JSON)

        """
        res = self.http_get("pathwayParticipants/{0}".format(identifier),
                            frmt='json')
        return res

    def pathway_complexes(self, identifier):
        """Get complexes belonging to a pathway

        :param int identifier: Pathway database identifier
        :return: list of all PhysicalEntity objects that participate in the
            Pathway.(in JSON)


        """
        res = self.http_get("pathwayComplexes/{0}".format(identifier),
                            frmt="json")
        return res

    def query_by_id(self, classname, identifier):
        """Get Reactome Database for a specific object.


        :param str classname: e.g. Pathway
        :param int identifier: database identifier or stable identifier if available

        It returns a full object, including full class information about
        all the attributes of the returned object. For example, if the object has
        one PublicationSource attribute, it will return a full PublicationSource
        object within the returned object.
        """
        url = "queryById/{0}/{1}".format(classname, identifier)
        res = self.http_get(url, frmt='json')
        return res

    def query_by_ids(self, classname, identifiers):
        """

        :param str classname: e.g. Pathway
        :param list identifiers: list of strings or int

        .. warning:: not sure the wrapping is correct
        """

        identifiers = self.devtools.list2string(identifiers)
        url = "queryByIds/{0}".format(classname)
        res = self.http_post(url, frmt="json", data=identifiers)
        # headers={'Content-Type': "application/json"})
        return res

    def query_hit_pathways(self, query):
        """Get pathways that contain one or more genes passed in the query list.

        In the Reactome data model, pathways are organized in a
        hierarchical structure. The returned pathways in this method are pathways
        having detailed manually drawn pathway diagrams. Currently only human
        pathways will be returned from this method.


        """
        identifiers = self.devtools.list2string(query)
        res = self.http_post("queryHitPathways", frmt='json', data=identifiers)
        return res

    def query_pathway_for_entities(self, identifiers):
        """Get pathway objects by specifying an array of PhysicalEntity database identifiers.


        The returned Pathways should
        contain the passed EventEntity objects. All passed EventEntity database
        identifiers should be in the database.

        """
        identifiers = self.devtools.list2string(identifiers, space=False)
        url = "pathwayForEntities"
        res = self.http_post(url, frmt='json', data={'ID': identifiers})
        return res

    def species_list(self):
        """Get the list of species used Reactome"""
        url = "speciesList"
        res = self.http_get(url, frmt='json')
        return res

    def SBML_exporter(self, identifier):
        """Get the SBML XML text of a pathway identifier

        :param int identifier: Pathway database identifier
        :return: SBML object in XML format as a string


        """
        url = "sbmlExporter/{0}".format(identifier)
        res = self.http_get(url, frmt='xml')
        return res

    def get_all_reactions(self):
        """Return list of reactions from the Pathway"""
        res = self.get_list_pathways()
        return [x[0] for x in res]

    def bioservices_get_reactants_from_reaction_identifier(self, reaction):
        """Fetch information from the reaction HTML page

        .. note:: draft version
        """
        res = self.http_get(
                'http://www._reactome.org/content/detail/%s' % reaction)
        res = res.content

        try:
            reactants = [x for x in res.split("\n") if '<title>' in x]
            reactants = reactants[0].split("|")[1].strip().strip('</title>')
        except  Exception as err:
            print('Could not interpret title')
            return res

        if reactants.count(':') == 1:
            reactants = reactants.split(":")
        else:
            pass

        return reactants
'''

pd.set_option('expand_frame_repr', False)
pd.set_option('display.width', 10000)
pd.set_option('display.max_rows', 5000)
pd.set_option('display.max_columns', 5000)
pd.options.display.max_colwidth = 5000

_reactome = Reactome()

verbose = True


name_dict = {
    'plasma membrane':                   'PlasmaMem',
    'cytosol':                           'CYTO',
    'mitochondrial outer membrane':      'MOM',
    'nucleoplasm':                       'NucPlasm',
    'endosome membrane':                 'EndosomeMem',
    'mitochondrial intermembrane space': 'MIM',
    'nuclear envelope':                  'NucEnv',
    'endoplasmic reticulum membrane':    'ERM',
    '[':                                 r'\n['
}


shapes = {
    'Protein':               'box',
    'Chemical Compound':     'oval',
    'Complex':               'note',
    'Set':                   'rectangle',
    'OtherEntity':           'note',
    'Genes and Transcripts': 'note',
    'DNA Sequence':          'note',
    'RNA Sequence':          'note'
}


def shorten_name(string):
    for s in name_dict:
        if s in string:
            string = string.replace(s, name_dict[s])
    return string


def add_to_graph(sample, list_of_species, graph):
    if type(sample) == int:
        # print('output type = int')
        # outputs.add(each)
        return
    if sample['className'] not in shapes:
        if verbose:
            print('Is a {}'.format(sample['className']))
        return
    name = sample['dbId']
    display = shorten_name(sample['displayName'])
    graph.add_node(name,
                   label=display,
                   displayName=display,
                   dbId=name,
                   shape=shapes[sample['className']])

    list_of_species[name] = display


def extract_list_from_key(keyword, y):
    input_list = []
    for i in y[keyword]:
        if isinstance(i, int):
            input_list.append(i)
        else:
            if 'dbId' in i.keys():
                input_list.append(i['dbId'])
            else:
                if verbose:
                    print("No dbID and not an int")
    return input_list


def get_entity_info(species):
    x = _reactome.get_entity_info(species)

    info_needed = ['referenceType', 'className', 'startCoordinate',
                   'endCoordinate', 'referenceEntity', 'compartment',
                   'referenceEntity', 'displayName', 'hasModifiedResidue']

    dont_need = ['inDisease', 'name', 'dbId', 'stId', 'inferredTo',
                 'speciesName', 'species', 'schemaClass', 'isChimeric',
                 'literatureReference']

    need_to_investigate = ['hasComponent', 'literatureReference',
                           'hasMember', 'summation', 'hasCandidate']
    entity_info = dict()
    for n in x:
        if n == 'databaseObject':
            if 'className' in x[n]:
                c_name = x[n]['className']
                if c_name == 'Reaction':
                    continue

            if 'referenceEntity' in x[n]:
                info = x[n]['referenceEntity']
                if 'dbId' in info:
                    entity_info['parent_dbid'] = info['dbId']
                if 'databaseName' in info:
                    entity_info['parent_db_name'] = info['databaseName']
                if 'identifier' in info:
                    entity_info['parent_identifier'] = info['identifier']
            if 'hasModifiedResidue' in x[n]:
                info = x[n]['hasModifiedResidue']
                mod_residues = []
                mods = []
                for r in info:
                    if 'coordinate' in r:
                        mod_residues.append(r['coordinate'])
                    if 'psiMod' in r:
                        ms = r['psiMod']
                        if isinstance(ms, int):
                            mods.append(ms)
                        elif isinstance(ms, dict):
                            if 'dbId' in ms:
                                mods.append(ms['dbId'])
                            else:
                                for m_type in ms:
                                    print(m_type, ms[m_type])
            for m in x[n]:
                if m in info_needed:
                    entity_info[m] = x[n][m]
                elif m not in dont_need:
                    print(species, m, x[n][m])
    return entity_info


def count_entites(list_of_object):
    count_dict = {}
    name_dict = {}
    for i in list_of_object:
        if i in count_dict:
            count_dict[i] += 1
        else:
            count_dict[i] = 1
    for i in count_dict:
        x = _reactome.get_entity_attribute(i, 'displayName')
        name_dict[i] = x
        # get_entity_info(i)
    return name_dict, count_dict


def create_graph(inputs, outputs, class_names, catalyst, rxn_name, graph,
                 prev_events):
    for each in inputs:
        name = each
        display = inputs[each]
        graph.add_node(name,
                       label=display,
                       displayName=display,
                       dbId=name,
                       shape=shapes[class_names[each]])
        if catalyst:
            if each == catalyst:
                pass
            else:
                graph.add_edge(catalyst, rxn_name)
        graph.add_edge(name, rxn_name)
        for i in prev_events:
            graph.add_edge(i, name)

    for each in outputs:
        name = each
        display = outputs[each]
        graph.add_node(name,
                       label=display,
                       displayName=display,
                       dbId=name,
                       shape=shapes[class_names[each]])
        graph.add_edge(rxn_name, each)


def create_graph_from_df(inputs, outputs, class_names, catalyst, rxn_name,
                         graph, prev_events):
    for each in inputs:
        name = each
        display = inputs[each]
        graph.add_node(name,
                       label=display,
                       displayName=display,
                       dbId=name,
                       shape=shapes[class_names[each]])
        if catalyst:
            if each == catalyst:
                pass
            else:
                graph.add_edge(catalyst, rxn_name)
        graph.add_edge(name, rxn_name)
        for i in prev_events:
            graph.add_edge(i, name)

    for each in outputs:
        name = each
        display = outputs[each]
        graph.add_node(name,
                       label=display,
                       displayName=display,
                       dbId=name,
                       shape=shapes[class_names[each]])
        graph.add_edge(rxn_name, each)


def _extract_info_from_reactants_or_products(r_dict, in_out):
    dict_of_info = dict()
    return_info = dict()
    for each in r_dict[in_out]:
        if isinstance(each, int):
            if each in return_info:
                return_info[each]['counter'] += 1
                continue
            else:
                if verbose:
                    print(each, "Not in inputs")
                continue
        dict_of_info['id'] = each['dbId']
        dict_of_info['counter'] = 1
        dict_of_info['species_type'] = each['className']
        dict_of_info['displayname'] = each['displayName']
        return_info[each['dbId']] = dict_of_info
    return return_info


def get_reaction_info(reaction, graph):

    if 'className' not in reaction:
        return None, None

    # ensure class is a reaction
    if reaction['className'] != 'Reaction':
        if verbose:
            print("Not a reaction")
            print(reaction)
        return None, None

    # get all reaction info
    y = _reactome.get_reaction_info(reaction['dbId'])

    # check to see if compartment exists
    if 'compartment' not in y:
        if verbose:
            print("No compartment")
        return None, None

    # check to make sure there is input and output
    if ('input' or 'output') not in y:
        if verbose:
            print("No input or output")
        return None, None

    # get reaction id and name
    rxn_name = reaction['dbId']
    rxn_display_name = reaction['displayName']

    if verbose:
        print("Reaction {} : {} ".format(rxn_name, rxn_display_name))
        for i in y:
            print("\t{} : {}".format(i, y[i]))

    # get compartment info
    comp = y['compartment']
    compartments = []

    for c in comp:
        if 'name' in c:
            compartments.append(c['name'])
    prev_events = []
    if 'precedingEvent' in y:
        prec_event = y['precedingEvent']
        if verbose:
            print("\tPreceeding event = {}".format(prec_event))
        for i in prec_event:
            if 'dbId' in i:
                prev_event_id = i['dbId']
                prev_events.append(prev_event_id)

    catalyst = False
    cat_all = []
    if 'catalystActivity' in y:
        catalyst = y['catalystActivity'][0]['dbId']
        if verbose:
            print(catalyst, 1)
        test = _reactome.get_event_participating_phys_entities(catalyst)
        # print('TEST', test)
        for n in test:
            if verbose:
                for j in n:
                    print('\t\t {} : {}'.format(j, n[j]))
        for n in test:
            if 'displayName' in n:
                catalyst = n['dbId']
                if verbose:
                    print(catalyst, 2)
                cat_all.append(catalyst)
                graph.add_node(catalyst,
                               label=shorten_name(n['displayName']),
                               shape='box', fillcolor='grey', style='filled')

    inputs = _extract_info_from_reactants_or_products(y, 'input')
    outputs = _extract_info_from_reactants_or_products(y, 'output')

    if verbose:
        print('Inputs : {}'.format(inputs))
        print('Outputs : {}'.format(outputs))

    reaction_info = dict()
    reaction_info['inputs'] = inputs
    reaction_info['outputs'] = outputs
    reaction_info['name'] = rxn_display_name
    reaction_info['id'] = rxn_name
    reaction_info['compartment'] = compartments
    reaction_info['catalyst'] = catalyst
    reaction_info['prev_events'] = prev_events
    # return reaction_info, rxn_name
    # quit()

    class_names = {}
    inputs = {}
    outputs = {}
    for each in y['input']:
        if isinstance(each, int):
            if each in inputs:
                continue
        id = each['dbId']
        display = shorten_name(each['displayName'])
        inputs[id] = display
        class_names[id] = each['className']

    for each in y['output']:
        if isinstance(each, int):
            if verbose:
                print(type(each))
            if each in outputs:
                continue
        else:
            id = each['dbId']
            display = shorten_name(each['displayName'])
            outputs[id] = display
            class_names[id] = each['className']

    if len(cat_all) > 1:
        if verbose:
            print(cat_all)
        quit()

    graph.add_node(rxn_name, displayName=rxn_display_name,
                   label=rxn_display_name)
    create_graph(inputs, outputs, class_names, catalyst, rxn_name, graph,
                 prev_events)

    return reaction_info, rxn_name


def extract_pathways(pathway_events):
    pathways = []
    for i in pathway_events:
        if verbose:
            print(i)
        if i['className'] == 'Pathway':
            pathways.append(i['dbId'])
    return pathways
