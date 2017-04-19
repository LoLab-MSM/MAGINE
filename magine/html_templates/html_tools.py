import os

import jinja2
import pandas as pd
from magine.plotting.species_plotting import create_gene_plots_per_go
from magine.data.formatter import pivot_table_for_export

env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(
                    searchpath=os.path.join(
                            os.path.dirname(__file__),
                            'templates'))
    )

range_number = """column_number:{},
filter_type: "range_number" """

auto_complete = """column_number:{},
    filter_type: "auto_complete",
    text_data_delimiter: "," """

chosen = """column_number:{},
    filter_type: "chosen"
     """  # text_data_delimiter: ","
"""
multi_select
range_number
filter_type:'select'
select_type: 'chosen',
"""
dict_of_templates = dict()
dict_of_templates['GO_id'] = range_number
dict_of_templates['GO_name'] = auto_complete
dict_of_templates['slim'] = auto_complete
dict_of_templates['aspect'] = auto_complete
dict_of_templates['ref'] = range_number
dict_of_templates['depth'] = range_number
dict_of_templates['enrichment_score'] = range_number
dict_of_templates['pvalue'] = range_number
dict_of_templates['n_genes'] = range_number


def write_single_table(table, save_name, title):

    template = env.get_template('single_table_view.html')

    # formats output to less precision and ints rather than floats
    tmp_table, format_dict = format_data_table(table)
    html_table = tmp_table.to_html(escape=False,
                                   # na_rep='-',
                                   formatters=format_dict,
                                   )
    template_vars = {"title":      title,
                     "table_name": html_table
                     }

    html_out = template.render(template_vars)
    with open('{}.html'.format(save_name), 'w') as f:
        f.write(html_out)


def write_table_to_html_with_figures(data, exp_data, save_name='index',
                                     out_dir='Figures'):
    # create plots of everything
    if isinstance(data, str):
        data = pd.read_csv(data)

    fig_dict, to_remove = create_gene_plots_per_go(data, save_name,
                                                   out_dir, exp_data)
    for i in fig_dict:
        data.loc[data['GO_id'] == i, 'GO_name'] = fig_dict[i]

    data = data[~data['GO_id'].isin(to_remove)]

    tmp = pivot_table_for_export(data)

    html_out = os.path.join(out_dir, save_name)
    print("Saving to : {}".format(html_out))

    write_single_table(tmp, html_out, 'MAGINE GO analysis')

    html_out = os.path.join(out_dir, save_name + '_filter')
    write_filter_table(tmp, html_out, 'MAGINE GO analysis')


def write_filter_table(table, save_name, title):
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
    for n, i in enumerate(table.index.names):
        if i in leave:
            continue
        new_string = dict_of_templates[i].format(n)
        out_string += '{' + new_string + '},\n'

    for m, i in enumerate(table.columns):
        if i[0] in leave:
            continue
        new_string = dict_of_templates[i[0]].format(n + m + 1)
        out_string += '{' + new_string + '},\n'
    # print(out_string)
    # formats output to less precision and ints rather than floats
    tmp_table, format_dict = format_data_table(table)
    html_table = tmp_table.to_html(escape=False,
                                   # na_rep='-',
                                   formatters=format_dict
                                   )
    template_vars = {
        "title":        title,
        "table_name":   html_table,
        "filter_table": out_string
        }

    template = env.get_template('filter_table.html')

    html_out = template.render(template_vars)
    with open('{}.html'.format(save_name), 'w') as f:
        f.write(html_out)


def format_data_table(data):
    tmp_table = data.copy()
    format_dict = {}
    for i in tmp_table.columns:
        if i[0] == 'enrichment_score':
            # format_dict[i] = '{:.2f}'.format
            tmp_table[i] = tmp_table[i].fillna(0)
            tmp_table[i] = tmp_table[i].round(2)
        elif i[0] == 'pvalue':
            format_dict[i] = '{:.2g}'.format
            tmp_table[i] = tmp_table[i].fillna(1)
            # tmp_table[i] = tmp_table[i].round(4)
        elif i[0] == 'n_genes':
            format_dict[i] = '{:,d}'.format
            tmp_table[i] = tmp_table[i].fillna(0)
            tmp_table[i] = tmp_table[i].astype(int)
    return tmp_table, format_dict


def format_ploty(text, save_name):

    template = env.get_template('plotly_template.html')

    template_vars = {
        "plotly_code":        text,
    }
    html_out = template.render(template_vars)
    with open('{}'.format(save_name), 'w') as f:
        f.write(html_out)
