import os
import tempfile

import magine.html_templates.html_tools as html_tools
from magine.enrichment.enrichment_result import load_enrichment_csv


def test_filter():
    df = load_enrichment_csv(os.path.join(
        os.path.dirname(__file__), 'Data', 'enrichr_test_enrichr.csv')
    )
    html_tools.write_filter_table(df, 'tmp')
    html_tools.write_filter_table(df, tempfile.mkdtemp())
