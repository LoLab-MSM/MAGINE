import os
import pandas as pd

import jinja2


def write_single_table(table, save_name, title):
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(
        searchpath=os.path.dirname(__file__))
    )
    template = env.get_template('single_table_view.html')
    tmp_table = table.copy()
    tmp_table = tmp_table.fillna(0)

    format_dict = {}
    for i in tmp_table.columns:
        if i[0] == 'enrichment_score':
            format_dict[i] = '{:.2f}'.format
        elif i[0] == 'pvalue':
            format_dict[i] = '{:.2g}'.format
        elif i[0] == 'n_genes':
            format_dict[i] = '{:,d}'.format
            tmp_table[i] = tmp_table[i].astype(int)
    html_table = tmp_table.to_html(escape=False,
                                   na_rep='-',
                                   formatters=format_dict,
                                   justify='left',
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
    from magine.plotting.species_plotting import create_gene_plots_per_go

    fig_dict, to_remove = create_gene_plots_per_go(data, save_name,
                                                   out_dir, exp_data)
    for i in fig_dict:
        data.loc[data['GO_id'] == i, 'GO_name'] = fig_dict[i]

    data = data[~data['GO_id'].isin(to_remove)]

    tmp = pd.pivot_table(data,
                         index=['GO_id', 'GO_name', 'depth', 'ref', 'slim',
                                'aspect'],
                         columns='sample_index')

    html_out = os.path.join(out_dir, save_name)
    print("Saving to : {}".format(html_out))
    write_single_table(tmp, html_out, 'MAGINE GO analysis')


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

    tem = """
    column_number:{},
         filter_type:"{}" """
    tem2 = """
    column_number:{},
         filter_type:"{}",
          text_data_delimiter: ','"""
    out_string = ''
    for n, i in enumerate(table.index.names):
        new_string = tem2.format(n, 'auto_complete')
        out_string += '{' + new_string + '},\n'

    for m, i in enumerate(table.columns):
        new_string = tem.format(n + m + 1, 'range_number_slider')
        out_string += '{' + new_string + '},\n'
    print(out_string)
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(
            searchpath=os.path.dirname(__file__))
    )
    template = env.get_template('filter_table.html')
    table = table.fillna(0)

    template_vars = {
        "title":        title,
        "table_name":   table.to_html(escape=False),
        "filter_table": out_string
        }

    html_out = template.render(template_vars)
    with open('{}.html'.format(save_name), 'w') as f:
        f.write(html_out)