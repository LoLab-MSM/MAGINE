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
# filter_template = env.get_template('filter_table_base.html')
filter_template = env.get_template('filter_table.html')
enrich_template = env.get_template('enrichment_template.html')

range_number = 'column_number:{},' \
               'filter_type: "range_number"'

auto_complete = 'column_number:{},' \
                'filter_type: "auto_complete",' \
                'text_data_delimiter: ","'

chosen = 'column_number:{}, ' \
         'filter_type: "multi_select",' \
         'select_type: "select2",' \
         'select_type_options: {{width: \'150px\'}}, ' \
         'text_data_delimiter: ","'

html_selector = 'column_number:{},' \
                'column_data_type: "html",' \
                'filter_type: "multi_select",' \
                'select_type: "chosen"'

dict_of_templates = dict(GO_id=range_number,
                         GO_name=chosen,
                         slim=chosen,
                         aspect=chosen,
                         ref=range_number,
                         depth=range_number,
                         enrichment_score=range_number,
                         term_name=chosen,
                         term_id=chosen,
                         rank=range_number,
                         p_value=range_number,
                         adj_p_value=range_number,
                         combined_score=range_number,
                         genes=chosen,
                         n_genes=range_number,
                         z_score=range_number,
                         significant_flag=chosen,
                         data_type=chosen,
                         pvalue=range_number,
                         treated_control_fold_change=range_number,
                         p_value_group_1_and_group_2=range_number,
                         protein=auto_complete,
                         gene=chosen,
                         time=chosen,
                         compound=auto_complete,
                         compound_id=auto_complete,
                         db=chosen,
                         )


def write_single_table(table, save_name, title):

    # formats output to less precision and ints rather than floats
    tmp_table, format_dict = _format_data_table(table)
    html_table = tmp_table.to_html(escape=False,
                                   # na_rep='-',
                                   formatters=format_dict,
                                   )
    template_vars = {"title":      title,
                     "table_name": html_table
                     }

    html_out = single_template.render(template_vars)
    with open('{}.html'.format(save_name), 'w') as f:
        f.write(html_out)


def process_filter_table(table, title):
    """

    Parameters
    ----------
    table : pandas.DataFrame
    save_name : str
    title : str

    Returns
    -------

    """
    """{column_number: 0},
    {column_number: 1, filter_type: "range_number_slider"},
    {column_number: 2, filter_type: "date"},
    {
        column_number:       3,
        filter_type:         "auto_complete",
        text_data_delimiter: ","
        },
    {
        column_number:        4,
        column_data_type:     "html",
        html_data_type:       "text",
        filter_default_label: "Select tag"
        }"""

    out_string = ''
    leave = ['GO_id', 'genes']
    n = 0
    for n, i in enumerate(table.index.names):
        if i in leave:
            continue
        if i not in dict_of_templates:
            print(i)
            continue
        new_string = dict_of_templates[i].format(n)
        out_string += '{' + new_string + '},\n'

    for m, i in enumerate(table.columns):
        if isinstance(i, str):
            new_string = dict_of_templates[i].format(n + m + 1)
            out_string += '{' + new_string + '},\n'
            continue
        elif i[0] in leave:
            continue
        elif i[0] not in dict_of_templates:
            print(i[0])
            continue
        new_string = dict_of_templates[i[0]].format(n + m + 1)
        out_string += '{' + new_string + '},\n'

    # formats output to less precision and ints rather than floats
    tmp_table, format_dict = _format_data_table(table)

    html_table = tmp_table.to_html(escape=False,
                                   # na_rep='-',
                                   formatters=format_dict
                                   )
    template_vars = {"title": title,
                     "table_name": html_table,
                     "filter_table": out_string}
    return template_vars


def write_filter_table(table, save_name, title):
    """

    Parameters
    ----------
    table : pandas.DataFrame
    save_name : str
    title : str

    Returns
    -------

    """
    """{column_number: 0},
    {column_number: 1, filter_type: "range_number_slider"},
    {column_number: 2, filter_type: "date"},
    {
        column_number:       3,
        filter_type:         "auto_complete",
        text_data_delimiter: ","
        },
    {
        column_number:        4,
        column_data_type:     "html",
        html_data_type:       "text",
        filter_default_label: "Select tag"
        }"""

    out_string = ''
    leave = ['GO_id', 'genes']
    n=0
    for n, i in enumerate(table.index.names):
        if i in leave:
            continue
        if i not in dict_of_templates:
            print(i)
            continue
        new_string = dict_of_templates[i].format(n)
        out_string += '{' + new_string + '},\n'

    for m, i in enumerate(table.columns):
        if isinstance(i, str):
            new_string = dict_of_templates[i].format(n + m + 1)
            out_string += '{' + new_string + '},\n'
            continue
        elif i[0] in leave:
            continue
        elif i[0] not in dict_of_templates:
            print(i[0])
            continue
        new_string = dict_of_templates[i[0]].format(n + m + 1)
        out_string += '{' + new_string + '},\n'

    # formats output to less precision and ints rather than floats
    tmp_table, format_dict = _format_data_table(table)
    
    html_table = tmp_table.to_html(escape=False,
                                   # na_rep='-',
                                   formatters=format_dict
                                   )
    template_vars = {"title":        title,
                     "table_name":   html_table,
                     "filter_table": out_string}

    html_out = filter_template.render(template_vars)
    with open('{}.html'.format(save_name), 'w') as f:
        f.write(html_out)


def _format_data_table(data):
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
    format_dict = {}
    pvalue_float_type = ['pvalue', 'p_value_group_1_and_group_2',
                         'adj_p_value']

    float_type = ['z_score', 'combined_score', 'enrichment_score',
                  'treated_control_fold_change']

    int_type = ['n_genes', 'rank']

    for i in tmp_table.columns:
        if i[0] in float_type:
            # format_dict[i] = '{:.2f}'.format
            format_dict[i] = '{:.4g}'.format
            tmp_table[i] = tmp_table[i].fillna(0)
            tmp_table[i] = tmp_table[i].astype(float)
            tmp_table[i] = tmp_table[i].round(2)
        elif i[0] in pvalue_float_type:
            format_dict[i] = '{:.2g}'.format
            tmp_table[i] = tmp_table[i].fillna(1)
            # tmp_table[i] = tmp_table[i].round(4)
        elif i[0] in int_type:
            format_dict[i] = '{:,d}'.format
            tmp_table[i] = tmp_table[i].fillna(0)
            tmp_table[i] = tmp_table[i].astype(int)

    return tmp_table, format_dict


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

    template_vars = {"plotly_code": text}
    html_out = plotly_template.render(template_vars)
    with open('{}'.format(save_name), 'w') as f:
        f.write(html_out)
