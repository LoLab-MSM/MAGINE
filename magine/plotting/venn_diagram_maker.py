import matplotlib

matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn2


def create_venn3(list1, list2, list3, label1, label2, label3, savename):
    set1 = set(list1)
    set2 = set(list2)
    set3 = set(list3)
    print(len(set3))
    v = venn3([set1, set2, set3], ('%s(%s)' % (label1, str(len(set1))),
                                   '%s(%s)' % (label2, str(len(set2))),
                                   '%s(%s)' % (label3, str(len(set3)))))
    plt.tight_layout()
    plt.savefig("%s.png" % savename, bbox_inches='tight')
    plt.savefig("%s.pdf" % savename, bbox_inches='tight')
    plt.close()


def create_venn2(list1, list2, label1, label2, savename):
    set1 = set(list1)
    set2 = set(list2)
    v = venn2([set1, set2], ('%s(%s)' % (label1, str(len(set1))), '%s(%s)' % (label2, str(len(set2))),))
    plt.tight_layout()
    plt.savefig("%s.png" % savename, bbox_inches='tight')
    plt.savefig("%s.pdf" % savename, bbox_inches='tight')
    plt.close()
