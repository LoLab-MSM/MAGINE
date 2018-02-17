import os

import jinja2

env = jinja2.Environment(
    loader=jinja2.FileSystemLoader(
        searchpath=os.path.join(os.path.dirname(__file__), 'templates')
    )
)

workflow_template = env.get_template('workflow_template.html')
plotly_template = env.get_template('plotly_template.html')
single_template = env.get_template('single_table_view.html')
filter_template = env.get_template('filter_table.html')

range_ = {'GO_id', 'ref', 'depth', 'enrichment_score', 'rank',
          'p_value', 'adj_p_value', 'combined_score', 'n_genes',
          'z_score', 'pvalue', 'treated_control_fold_change',
          'p_value_group_1_and_group_2'}

chosen_ = {'GO_name', 'slim', 'aspect', 'term_name', 'term_id', 'genes',
           'significant_flag', 'db', 'sample_id', 'sample_index'}

auto_complete_ = {'protein', 'gene'}


def _add_filter(column_num, f_type):
    _default = {'column_number': column_num}
    if f_type == 'range':
        _default.update({'filter_type': 'range_number'})
    if f_type == 'select':
        _default.update({'column_data_type': 'html',
                         'filter_type': "multi_select",
                         'select_type': "chosen"})
    if f_type == 'chosen':
        _default.update({'column_data_type': 'html',
                         'filter_type': "','",
                         'select_type': "multi_select",
                         'select_type_options': '{width:"150px"}',
                         'text_data_delimiter': "select2"},
                        )
    if f_type == 'auto_complete':
        _default.update({'filter_type': "auto_complete",
                         'text_data_delimiter': "','"})

    return _default


def write_single_table(table, title, save_name=None):
    """

    Parameters
    ----------
    table : pandas.DataFrame
    title : str
    save_name : str

    Returns
    -------

    """

    # formats output to less precision and ints rather than floats
    tmp_table = _format_simple_table(table)
    html_out = single_template.render({"title": title,
                                       "table_name": tmp_table})
    if save_name is None:
        return html_out
    with open('{}.html'.format(save_name), 'w') as f:
        f.write(html_out)


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


def create_yadf_filters(table):
    _format_dict = []
    for n, i in enumerate(table.columns):
        if len(i) == 2:
            i, j = i
        if i in range_:
            _format_dict.append(_add_filter(n, 'range'))
        if i in chosen_:
            _format_dict.append(_add_filter(n, 'select'))
        if i in auto_complete_:
            _format_dict.append(_add_filter(n, 'complete'))
    return _format_dict


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
        if i in float_type:
            tmp_table[i] = tmp_table[i].fillna(0)
            tmp_table[i] = tmp_table[i].astype(float)
            tmp_table[i] = tmp_table[i].round(2)
            tmp_table[i] = tmp_table[i].apply('{:.4g}'.format)
        elif i in pvalue_float_type:
            tmp_table[i] = tmp_table[i].fillna(1)
            tmp_table[i] = tmp_table[i].apply('{:.2g}'.format)
        elif i in int_type:
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
