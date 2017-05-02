from magine.plotting.venn_diagram_maker import create_venn3, create_venn2

x = ['A', 'B', 'C', 'D']
y = ['C', 'D', 'E', 'F']
z = ['D', 'E', 'F', 'N', 'Z']


def test_venn_2():
    create_venn2(x, y, 'X', 'Y', 'test_1')


def test_venn_3():
    create_venn3(x, y,z, 'X', 'Y','z', 'test_1')


