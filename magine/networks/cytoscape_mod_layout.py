import requests

IP = 'localhost'
PORT = 1234
VERSION = 'v1'
BASE_URL = 'http://' + IP + ':' + str(PORT) + '/' + VERSION + '/'
HEADERS = {'Content-Type': 'application/json'}

SUID_LIST = 'suid'
url = 'http://' + IP + ':' + str(PORT) + '/' + VERSION + '/'


class LayoutClient(object):
    def __init__(self):
        self.__url = url + 'apply/layouts'

    def get_all(self):
        return requests.get(self.__url).json()

    def get_option(self, name):
        url = self.__url + '/' + name + '/columntypes'
        return requests.get(url).json()

    def get_parameters(self, name):
        url = self.__url + '/' + name + '/parameters'
        return requests.get(url).json()

    def update(self, name, parameters='color'):
        url = self.__url + '/' + name + '/' + str(parameters)
        requests.put(url)

    def apply(self, name='force-directed', network=None, params=None):
        if network is None:
            raise ValueError('Target network is required')

        url = self.__url + '/' + name + '/' + str(network.get_id())

        requests.get(url, params)

    def bundle_edge(self, network=None):
        if network is None:
            raise ValueError('Target network is required')

        url = self.__url + '/edgebundling/' + str(network.get_id())
        requests.get(url)
