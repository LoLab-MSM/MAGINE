from bioservices import PathwayCommons
pc = PathwayCommons()
pc.TIMEOUT =100


def find_binding_partners(gene):
    x = pc.graph(source=gene,
                 target='BAX',
                 kind="neighborhood",
                 frmt="EXTENDED_BINARY_SIF")
    print(x)
    return x
if __name__ == '__main__':
    find_binding_partners('BID')
