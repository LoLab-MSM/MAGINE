from bioservices import PathwayCommons

pc = PathwayCommons()
pc.TIMEOUT =100
x = pc.graph(source="BAX",
             kind="neighborhood",
             frmt="EXTENDED_BINARY_SIF")
print(x)