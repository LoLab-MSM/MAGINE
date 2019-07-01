import os

import jinja2

env = jinja2.Environment(
    loader=jinja2.FileSystemLoader(
        searchpath=os.path.join(os.path.dirname(__file__), 'templates')
    ),
    autoescape=False
)

workflow_template = env.get_template('workflow_template.html')
plotly_template = env.get_template('plotly_template.html')
single_template = env.get_template('single_table_view.html')
filter_template = env.get_template('filter_table.html')

range_ = {'ref', 'depth', 'enrichment_score', 'rank', 'fold_change',
          'p_value', 'adj_p_value', 'combined_score', 'n_genes',
          'z_score', 'pvalue', 'treated_control_fold_change',
          'p_value_group_1_and_group_2'}

chosen_ = {
    'GO_name', 'slim', 'aspect', 'time', 'data_type', 'identifier', 'label',
    'term_name', 'term_id', 'category', 'db',  # enrichr
    'significant_flag', 'project_name', 'compound', 'GO_id',
    'sample_id', 'sample_index'
}
multi_choose = {'genes', 'source', 'significant'}

auto_complete_ = {'protein', 'gene', 'identifier', 'label', }


def _add_filter(column_num, f_type):
    _default = {'column_number': column_num}
    if f_type == 'range':
        _default.update({'filter_type': 'range_number'})
    elif f_type == 'select':
        _default.update({
            'select_type_options': {'width': "150px"},
            'filter_type': 'multi_select',
            'select_type': 'select2'
        })
    elif f_type == 'chosen':
        _default.update({'select_type': "multi_select",
                         'select_type_options': {'width': "150px"},
                         'text_data_delimiter': "chosen"},
                        )
    elif f_type == 'auto_complete':
        _default.update({'filter_type': "auto_complete",
                         'text_data_delimiter': ","})
    elif f_type == 'multichoose':
        _default.update({
            'filter_type': "multi_select",
            'select_type': "select2",
            'select_type_options': {'width': "150px"},
            'text_data_delimiter': ","},
        )
    return _default


def create_yadf_filters(table):
    _format_dict = []
    for n, i in enumerate(table.columns):
        if isinstance(i, (tuple, list)):
            i = i[0]
        if i in range_:
            _format_dict.append(_add_filter(n, 'range'))
        elif i in chosen_:
            _format_dict.append(_add_filter(n, 'select'))
        elif i in auto_complete_:
            _format_dict.append(_add_filter(n, 'complete'))
        elif i in multi_choose:
            _format_dict.append(_add_filter(n, 'multichoose'))
        else:
            print("'{}' not found in yadf".format(i))
    return _format_dict


def write_filter_table(table, save_name):
    """

    Parameters
    ----------
    table : pandas.DataFrame
    save_name : str

    Returns
    -------

    """

    # formats output to less precision and ints rather than floats
    tmp_table = _format_simple_table(table)

    data = tmp_table.to_dict('split')
    data['filters'] = create_yadf_filters(table)
    html_out = filter_template.render({"data": data})
    with open('{}.html'.format(save_name), 'w') as f:
        f.write(html_out)


def _format_simple_table(data):
    """
    formats precession of data for outputs

    Parameters
    ----------
    data : pandas.DataFrame

    Returns
    -------
    pandas.DataFrame, dict
    """
    tmp_table = data.copy()
    pvalue_float_type = ['pvalue', 'p_value', 'p_value_group_1_and_group_2',
                         'adj_p_value']

    float_type = ['z_score', 'combined_score', 'enrichment_score',
                  'treated_control_fold_change']

    int_type = ['n_genes', 'rank']
    for i in data.columns:
        if isinstance(i, (tuple, list)):
            t_i = i[0]
        else:
            t_i = i
        if t_i in float_type:
            tmp_table[i] = tmp_table[i].fillna(0)
            tmp_table[i] = tmp_table[i].astype(float)
            tmp_table[i] = tmp_table[i].round(2)
            tmp_table[i] = tmp_table[i].apply('{:.4g}'.format)
        elif t_i in pvalue_float_type:
            tmp_table[i] = tmp_table[i].fillna(1)
            tmp_table[i] = tmp_table[i].apply('{:.2g}'.format)
        elif t_i in int_type:
            tmp_table[i] = tmp_table[i].fillna(0)
            tmp_table[i] = tmp_table[i].astype(int)
            tmp_table[i] = tmp_table[i].apply('{:,d}'.format)
    return tmp_table


def format_ploty(text, save_name):
    """
    
    Parameters
    ----------
    text : str
        html code to embed in file
    save_name : str
        html output filename

    Returns
    -------

    """

    html_out = plotly_template.render({"plotly_code": text})
    with open('{}'.format(save_name), 'w') as f:
        f.write(html_out)
