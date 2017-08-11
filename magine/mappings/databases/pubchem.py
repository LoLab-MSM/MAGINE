'''
from bioservices.services import REST


__all__ = ["PubChem"]


class PubChem(REST):
    """Interface to the `PubChem <todo>`_ service

    """
    _url = "http://pubchem.ncbi.nlm.nih.gov/rest/pug"

    def __init__(self, verbose=False, cache=False):
        """**Constructor**

        :param verbose: set to False to prevent informative messages
        """
        print("PubChem is not finalised yet. This is currently only a draft version")
        super(PubChem, self).__init__(name="PubChem", url=PubChem._url,
                verbose=verbose, cache=cache)

    def get_compound_by_smiles(self, identifier, frmt='json'):

        res = self.http_get("compound/smiles/" + identifier
                + '/cids/%s' % frmt, frmt=frmt,
                headers=self.get_headers(content=frmt))
        return res

    def get_compound_by_name(self, name, frmt='json'):
        res = self.http_get("compound/name/{}/{}".format(name,frmt),
                            frmt=frmt,
                            headers=self.get_headers(content=frmt))

        return res

    def get_drug_classification(self, name, frmt='json'):
        res = self.http_get("compound/name/{}/classification/{}".format(name, frmt),
                            frmt=frmt,
                            headers=self.get_headers(content=frmt))
        return res


def test_drug_classification():
    pc = PubChem()
    x = pc.get_drug_classification('cis-Diaminedichloroplatinum')
    print(x.keys())
    hie = x['Hierarchies']['Hierarchy']
    for i in hie:
        print(i['SourceName'])
        print(i.keys())
        if i['SourceName'] == 'MeSH':
            for n in i['Node']:
                print(n)
                print(n['Information']['Name'])
                if 'Description' in n['Information']:
                    print(n['Information']['Description'])
                if n['Information']['Name'] == 'Molecular Mechanisms of Pharmacological Action':
                    print(n)

        # print(i.keys())
        # for k in i.keys():
        #     print(k)




def test_drug_chem_properties():
    pc = PubChem()
    x = pc.get_compound_by_name('Tazobactam')
    dtype_dict = dict()
    dtype_dict[1] = 'sval'
    dtype_dict[7] = 'fval'
    dtype_dict[5] = 'ival'
    dtype_dict[16] = 'binary'

    for i in x['PC_Compounds']:
        print(i)
        print(i.keys())
        # [u'count', u'stereo', u'bonds', u'atoms', u'charge', u'coords', u'props', u'id']
        keys = ['name']
        for j in i['props']:
            # print(j)
            print(j['value'][dtype_dict[j['urn']['datatype']]])
            print()
            continue

            if 'name' in j['urn']:
                print(j['urn']['name'])

            if 'label' in j['urn']:
                print(j['urn']['label'])
                print(j['urn'], j['value'])


if __name__ == '__main__':

    # test_drug_chem_properties()
    test_drug_classification()
    quit()
    pc = PubChem()
    # x = pc.get_compound_by_smiles('CC1(C(N2C(S1(=O)=O)CC2=O)C(=O)O)CN3C=CN=N3')
    # print(x)

    # x = pc.get_compound_by_name('Tazobactam')
'''
