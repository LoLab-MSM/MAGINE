import os

import jinja2


def write_single_table(table, save_name, title):
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(
        searchpath=os.path.dirname(__file__))
    )
    template = env.get_template('single_table_view.html')
    template_vars = {"title": title,
                     "table_name": table.to_html(escape=False)}
    html_out = template.render(template_vars)
    with open('{}.html'.format(save_name), 'w') as f:
        f.write(html_out)
