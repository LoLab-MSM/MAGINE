import matplotlib.pyplot as plt
from matplotlib_venn import venn2 as _venn2, venn3 as _venn3


def create_venn3(list1, list2, list3, label1, label2, label3, save_name=None,
                 image_format='png', title=None, ax=None,
                 colors=('g', 'r', 'b')):
    """ Creates a venn digram containing for 3 lists
    
    Parameters
    ----------
    list1 : list_like
    list2 : list_like
    list3 : list_like
    label1 : str
    label2 : str
    label3 : str
    save_name : str
    image_format : str
        default png
    title: str
    ax : matplotlib.axes

    Returns
    -------

    """
    set1 = set(list1)
    set2 = set(list2)
    set3 = set(list3)
    if len(colors) != 3:
        raise Exception("Must provide a list of three colors")
    v = _venn3([set1, set2, set3], ('%s(%s)' % (label1, str(len(set1))),
                                    '%s(%s)' % (label2, str(len(set2))),
                                    '%s(%s)' % (label3, str(len(set3)))),
               set_colors=colors,
               ax=ax)
    plt.tight_layout()
    if title is not None:
        plt.title(title)
    if save_name is not None:
        save_name = "{}.{}".format(save_name, image_format)
        plt.savefig(save_name, bbox_inches='tight')
    return v


def create_venn2(list1, list2, label1, label2, save_name=None, title=None,
                 image_format='png', ax=None):
    """Creates a venn digram containing for 2 lists
    
    Parameters
    ----------
    list1 : list_like
    list2 : list_like
    label1 : str
    label2 : str
    save_name : str
    title : str
    image_format : str, optional
        default png
    ax : matplotlib.axes
    
    Returns
    -------
    
    """
    set1 = set(list1)
    set2 = set(list2)
    v = _venn2([set1, set2], ('%s(%s)' % (label1, str(len(set1))),
                              '%s(%s)' % (label2, str(len(set2))),),
               ax=ax)
    plt.tight_layout()
    if title is not None:
        plt.title(title)
    if save_name is not None:
        save_name = "{}.{}".format(save_name, image_format)
        plt.savefig(save_name, bbox_inches='tight')

    return v
