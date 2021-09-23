Tutorial
========

.. code:: ipython3

    from IPython.display import display, Image
    %matplotlib inline
    import matplotlib.pyplot as plt


.. code:: ipython3

    import pandas as pd
    import seaborn as sns
    import numpy as np

ExperimentalData
----------------

Since MAGINE is built for multi-sample, multi-omics data, it is no
surprise that the data is the most important aspect. Here we should how
to use the ``ExperimentalData`` class. We designed MAGINE data input to
be as flexible as possible, requiring a standard format for 8 columns.
Users are required to format the input files to share the same column
names. Additional columns can still be used for additional tags on the
data.

+-----------------------+-----------------------+-----------------------+
| Column                | Data type             | Description           |
+=======================+=======================+=======================+
| identifier            | string                | HNGC or HMDB ids      |
+-----------------------+-----------------------+-----------------------+
| label                 | string                | any label that you    |
|                       |                       | would like to have as |
|                       |                       | an alterative to      |
|                       |                       | identifier. Useful    |
|                       |                       | for PTMS              |
|                       |                       | (BAX_S(ph)292),       |
|                       |                       | aliases, or chemical  |
|                       |                       | names (Deoxyuridine)  |
+-----------------------+-----------------------+-----------------------+
| species_type          | string                | Options (gene,        |
|                       |                       | protein, or           |
|                       |                       | metabolite) Right now |
|                       |                       | we are trying to work |
|                       |                       | out a more            |
|                       |                       | comprehensive tag to  |
|                       |                       | use. For now, we      |
|                       |                       | suggest using         |
|                       |                       | ``protein`` for all   |
|                       |                       | gene related products |
|                       |                       | and use the source    |
|                       |                       | column to further     |
|                       |                       | distinguish.          |
+-----------------------+-----------------------+-----------------------+
| significant           | bool                  | used to label if the  |
|                       |                       | measurment was        |
|                       |                       | signficant compared   |
|                       |                       | to control.           |
+-----------------------+-----------------------+-----------------------+
| fold_change           | float                 | Assumes non-log space |
+-----------------------+-----------------------+-----------------------+
| p_value               | float                 | HNGC or HMDB ids      |
+-----------------------+-----------------------+-----------------------+
| identifier            | string                | Significance          |
|                       |                       | statistic for fold    |
|                       |                       | change calculation    |
|                       |                       | (can be FDR, BH, etc  |
|                       |                       | corrected.            |
+-----------------------+-----------------------+-----------------------+
| source                | string                | Experimental platform |
|                       |                       | (ie SILAC, PH_SILAC,  |
|                       |                       | LABEL_FREE, RNASEQ)   |
+-----------------------+-----------------------+-----------------------+
| sample_id             | string                | Used to identify      |
|                       |                       | sample (time points,  |
|                       |                       | drug dose, etc)       |
+-----------------------+-----------------------+-----------------------+

For this tutorial, we are going to use our time series multi-omic
response of A549 cells to cisplatin.

The description of the experiments and dataset can be found in `Norris,
Jeremy L., et
al. <https://pubs.acs.org/doi/abs/10.1021/acs.jproteome.6b01004>`__

This file is in the format as show above in the table.

.. code:: ipython3

    # load the experimental data
    from magine.data.experimental_data import load_data
    
    exp_data = load_data('Data/norris_et_al_2017_cisplatin_data.csv.gz', low_memory=False)

About the data
~~~~~~~~~~~~~~

This dataset consists of 6 experimental platforms across 4 time points.
The ExperimentData class is designed to explore these data in a seamless
way. The core of the class is a pandas.DataFrame, however we built
additional functions that we frequently used. We store meta data, such
as the ``sample_ids``, ``exp_methods``.

.. code:: ipython3

    exp_data.sample_ids




.. parsed-literal::

    ['01hr', '06hr', '24hr', '48hr']



.. code:: ipython3

    exp_data.exp_methods




.. parsed-literal::

    ['rna_seq', 'ph_silac', 'label_free', 'silac', 'C18', 'HILIC']



There are many functions built around the ExperimentalData class. This
notebook is meant to demostrate an examples workflow of how they can be
used. Please refer to the documentation for exhuastive information about
the functions.

Summary of data
~~~~~~~~~~~~~~~

Gather number of measured species per time point and platform.

.. code:: ipython3

    exp_data.create_summary_table()




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th>sample_id</th>
          <th>01hr</th>
          <th>06hr</th>
          <th>24hr</th>
          <th>48hr</th>
          <th>Total Unique Across</th>
        </tr>
        <tr>
          <th>source</th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>C18</th>
          <td>522</td>
          <td>227</td>
          <td>653</td>
          <td>685</td>
          <td>1402</td>
        </tr>
        <tr>
          <th>HILIC</th>
          <td>471</td>
          <td>605</td>
          <td>930</td>
          <td>613</td>
          <td>1504</td>
        </tr>
        <tr>
          <th>label_free</th>
          <td>2766</td>
          <td>2742</td>
          <td>2551</td>
          <td>2261</td>
          <td>3447</td>
        </tr>
        <tr>
          <th>ph_silac</th>
          <td>2608</td>
          <td>3298</td>
          <td>3384</td>
          <td>3236</td>
          <td>5113</td>
        </tr>
        <tr>
          <th>rna_seq</th>
          <td>18741</td>
          <td>19104</td>
          <td>19992</td>
          <td>-</td>
          <td>20642</td>
        </tr>
        <tr>
          <th>silac</th>
          <td>2923</td>
          <td>3357</td>
          <td>3072</td>
          <td>3265</td>
          <td>4086</td>
        </tr>
      </tbody>
    </table>
    </div>



Count number of significantly changed (signficant_flag=True).

.. code:: ipython3

    exp_data.create_summary_table(sig=True)




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th>sample_id</th>
          <th>01hr</th>
          <th>06hr</th>
          <th>24hr</th>
          <th>48hr</th>
          <th>Total Unique Across</th>
        </tr>
        <tr>
          <th>source</th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>C18</th>
          <td>522</td>
          <td>227</td>
          <td>653</td>
          <td>685</td>
          <td>1402</td>
        </tr>
        <tr>
          <th>HILIC</th>
          <td>471</td>
          <td>605</td>
          <td>930</td>
          <td>613</td>
          <td>1504</td>
        </tr>
        <tr>
          <th>label_free</th>
          <td>196</td>
          <td>46</td>
          <td>271</td>
          <td>874</td>
          <td>1085</td>
        </tr>
        <tr>
          <th>ph_silac</th>
          <td>514</td>
          <td>888</td>
          <td>1227</td>
          <td>851</td>
          <td>2278</td>
        </tr>
        <tr>
          <th>rna_seq</th>
          <td>73</td>
          <td>1999</td>
          <td>12215</td>
          <td>-</td>
          <td>12340</td>
        </tr>
        <tr>
          <th>silac</th>
          <td>38</td>
          <td>52</td>
          <td>228</td>
          <td>266</td>
          <td>485</td>
        </tr>
      </tbody>
    </table>
    </div>



MAGINE uses the ``identifier`` column as the default index. This keeps
things simple when using the output for other tools (passing to
molecular networks). You can also pass an index argument to calculate
other values. Here, we use the ``label`` column, which contains PTMs of
our protein species. See how ``rna_seq`` values do not change, but there
is an increase number of ``ph_silac`` unique species.

.. code:: ipython3

    exp_data.create_summary_table(sig=True, index='label')




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th>sample_id</th>
          <th>01hr</th>
          <th>06hr</th>
          <th>24hr</th>
          <th>48hr</th>
          <th>Total Unique Across</th>
        </tr>
        <tr>
          <th>source</th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>C18</th>
          <td>528</td>
          <td>227</td>
          <td>657</td>
          <td>689</td>
          <td>1412</td>
        </tr>
        <tr>
          <th>HILIC</th>
          <td>479</td>
          <td>611</td>
          <td>941</td>
          <td>621</td>
          <td>1521</td>
        </tr>
        <tr>
          <th>label_free</th>
          <td>201</td>
          <td>46</td>
          <td>281</td>
          <td>911</td>
          <td>1149</td>
        </tr>
        <tr>
          <th>ph_silac</th>
          <td>594</td>
          <td>1370</td>
          <td>2414</td>
          <td>1368</td>
          <td>4757</td>
        </tr>
        <tr>
          <th>rna_seq</th>
          <td>73</td>
          <td>1999</td>
          <td>12215</td>
          <td>-</td>
          <td>12340</td>
        </tr>
        <tr>
          <th>silac</th>
          <td>38</td>
          <td>52</td>
          <td>228</td>
          <td>266</td>
          <td>485</td>
        </tr>
      </tbody>
    </table>
    </div>



Filter by category (experimental method)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We can access the input data using the ``.species`` property. This
returns a modified pandas.Datatable.

.. code:: ipython3

    exp_data.species.head(5)




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>identifier</th>
          <th>label</th>
          <th>species_type</th>
          <th>fold_change</th>
          <th>p_value</th>
          <th>significant</th>
          <th>sample_id</th>
          <th>source</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>HOXD1</td>
          <td>HOXD1_rnaseq</td>
          <td>protein</td>
          <td>-520.256762</td>
          <td>0.00102</td>
          <td>True</td>
          <td>06hr</td>
          <td>rna_seq</td>
        </tr>
        <tr>
          <th>1</th>
          <td>MIR7704</td>
          <td>MIR7704_rnaseq</td>
          <td>protein</td>
          <td>-520.256762</td>
          <td>0.00102</td>
          <td>True</td>
          <td>06hr</td>
          <td>rna_seq</td>
        </tr>
        <tr>
          <th>2</th>
          <td>AC078814.1</td>
          <td>AC078814.1_rnaseq</td>
          <td>protein</td>
          <td>-76.022260</td>
          <td>0.00102</td>
          <td>True</td>
          <td>06hr</td>
          <td>rna_seq</td>
        </tr>
        <tr>
          <th>3</th>
          <td>PPM1H</td>
          <td>PPM1H_rnaseq</td>
          <td>protein</td>
          <td>-76.022260</td>
          <td>0.00102</td>
          <td>True</td>
          <td>06hr</td>
          <td>rna_seq</td>
        </tr>
        <tr>
          <th>4</th>
          <td>PLCH1</td>
          <td>PLCH1_rnaseq</td>
          <td>protein</td>
          <td>-17.888990</td>
          <td>0.00102</td>
          <td>True</td>
          <td>06hr</td>
          <td>rna_seq</td>
        </tr>
      </tbody>
    </table>
    </div>



We added attributes to the classs to quickly separate the data based on
various input columns. We use the ``species_type`` and ``source`` column
name to split data into ``compounds``, ``genes`` (includes
``species_type``\ ==\ ``gene``), ``rna`` (includes
``species_type``\ ==\ ``protein``, ``source`` == ``rna``), or
``protein`` (``species_type``\ ==\ ``gene``, ``source`` != ``rna``).
They can be accessed with the “.prefix”, such as

.. code:: ipython3

    exp_data.genes.head(5)




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>identifier</th>
          <th>label</th>
          <th>species_type</th>
          <th>fold_change</th>
          <th>p_value</th>
          <th>significant</th>
          <th>sample_id</th>
          <th>source</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>HOXD1</td>
          <td>HOXD1_rnaseq</td>
          <td>protein</td>
          <td>-520.256762</td>
          <td>0.00102</td>
          <td>True</td>
          <td>06hr</td>
          <td>rna_seq</td>
        </tr>
        <tr>
          <th>1</th>
          <td>MIR7704</td>
          <td>MIR7704_rnaseq</td>
          <td>protein</td>
          <td>-520.256762</td>
          <td>0.00102</td>
          <td>True</td>
          <td>06hr</td>
          <td>rna_seq</td>
        </tr>
        <tr>
          <th>2</th>
          <td>AC078814.1</td>
          <td>AC078814.1_rnaseq</td>
          <td>protein</td>
          <td>-76.022260</td>
          <td>0.00102</td>
          <td>True</td>
          <td>06hr</td>
          <td>rna_seq</td>
        </tr>
        <tr>
          <th>3</th>
          <td>PPM1H</td>
          <td>PPM1H_rnaseq</td>
          <td>protein</td>
          <td>-76.022260</td>
          <td>0.00102</td>
          <td>True</td>
          <td>06hr</td>
          <td>rna_seq</td>
        </tr>
        <tr>
          <th>4</th>
          <td>PLCH1</td>
          <td>PLCH1_rnaseq</td>
          <td>protein</td>
          <td>-17.888990</td>
          <td>0.00102</td>
          <td>True</td>
          <td>06hr</td>
          <td>rna_seq</td>
        </tr>
      </tbody>
    </table>
    </div>



.. code:: ipython3

    exp_data.compounds.head(5)




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>identifier</th>
          <th>label</th>
          <th>species_type</th>
          <th>fold_change</th>
          <th>p_value</th>
          <th>significant</th>
          <th>sample_id</th>
          <th>source</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>128152</th>
          <td>HMDB0036114</td>
          <td>(-)-3-Thujone</td>
          <td>metabolites</td>
          <td>1.6</td>
          <td>2.100000e-02</td>
          <td>True</td>
          <td>06hr</td>
          <td>C18</td>
        </tr>
        <tr>
          <th>128153</th>
          <td>HMDB0001320</td>
          <td>(13E)-11a-Hydroxy-9,15-dioxoprost-13-enoic acid</td>
          <td>metabolites</td>
          <td>88.8</td>
          <td>5.800000e-12</td>
          <td>True</td>
          <td>24hr</td>
          <td>C18</td>
        </tr>
        <tr>
          <th>128154</th>
          <td>HMDB0012113</td>
          <td>(22Alpha)-hydroxy-campest-4-en-3-one</td>
          <td>metabolites</td>
          <td>100.0</td>
          <td>9.500000e-04</td>
          <td>True</td>
          <td>48hr</td>
          <td>HILIC</td>
        </tr>
        <tr>
          <th>128155</th>
          <td>HMDB0010361</td>
          <td>(23S)-23,25-dihdroxy-24-oxovitamine D3 23-(bet...</td>
          <td>metabolites</td>
          <td>-100.0</td>
          <td>1.000000e-12</td>
          <td>True</td>
          <td>48hr</td>
          <td>C18</td>
        </tr>
        <tr>
          <th>128156</th>
          <td>HMDB0011644</td>
          <td>(24R)-Cholest-5-ene-3-beta,7-alpha,24-triol</td>
          <td>metabolites</td>
          <td>1.6</td>
          <td>7.400000e-05</td>
          <td>True</td>
          <td>01hr</td>
          <td>C18</td>
        </tr>
      </tbody>
    </table>
    </div>



Similarily, we can also filter the data by ``source`` using the
``.name``, where name is anything in the ``source`` column. We can get a
list of these by printing ``exp_data.exp_methods``.

.. code:: ipython3

    # prints all the available exp_methods
    exp_data.exp_methods




.. parsed-literal::

    ['rna_seq', 'ph_silac', 'label_free', 'silac', 'C18', 'HILIC']



.. code:: ipython3

    # filters to only the 'label_free' 
    exp_data.label_free.shape




.. parsed-literal::

    (13085, 8)



.. code:: ipython3

    exp_data.label_free.head(5)




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>identifier</th>
          <th>label</th>
          <th>species_type</th>
          <th>fold_change</th>
          <th>p_value</th>
          <th>significant</th>
          <th>sample_id</th>
          <th>source</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>102446</th>
          <td>LIMS1</td>
          <td>LIMS1_lf</td>
          <td>protein</td>
          <td>12.42</td>
          <td>0.00003</td>
          <td>True</td>
          <td>01hr</td>
          <td>label_free</td>
        </tr>
        <tr>
          <th>102447</th>
          <td>SMARCE1</td>
          <td>SMARCE1_lf</td>
          <td>protein</td>
          <td>-2.49</td>
          <td>0.00030</td>
          <td>True</td>
          <td>01hr</td>
          <td>label_free</td>
        </tr>
        <tr>
          <th>102448</th>
          <td>HEXA</td>
          <td>HEXA_lf</td>
          <td>protein</td>
          <td>6.42</td>
          <td>0.00060</td>
          <td>True</td>
          <td>01hr</td>
          <td>label_free</td>
        </tr>
        <tr>
          <th>102449</th>
          <td>SRSF1</td>
          <td>SRSF1_lf</td>
          <td>protein</td>
          <td>-3.21</td>
          <td>0.00060</td>
          <td>True</td>
          <td>01hr</td>
          <td>label_free</td>
        </tr>
        <tr>
          <th>102450</th>
          <td>SF3B1</td>
          <td>SF3B1_lf</td>
          <td>protein</td>
          <td>-1.57</td>
          <td>0.00130</td>
          <td>True</td>
          <td>01hr</td>
          <td>label_free</td>
        </tr>
      </tbody>
    </table>
    </div>



.. code:: ipython3

    exp_data.HILIC.head(5)




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>identifier</th>
          <th>label</th>
          <th>species_type</th>
          <th>fold_change</th>
          <th>p_value</th>
          <th>significant</th>
          <th>sample_id</th>
          <th>source</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>128154</th>
          <td>HMDB0012113</td>
          <td>(22Alpha)-hydroxy-campest-4-en-3-one</td>
          <td>metabolites</td>
          <td>100.0</td>
          <td>0.000950</td>
          <td>True</td>
          <td>48hr</td>
          <td>HILIC</td>
        </tr>
        <tr>
          <th>128157</th>
          <td>HMDB0011644</td>
          <td>(24R)-Cholest-5-ene-3-beta,7-alpha,24-triol</td>
          <td>metabolites</td>
          <td>1.7</td>
          <td>0.000072</td>
          <td>True</td>
          <td>24hr</td>
          <td>HILIC</td>
        </tr>
        <tr>
          <th>128162</th>
          <td>HMDB0012114</td>
          <td>(3S)-3,6-Diaminohexanoate</td>
          <td>metabolites</td>
          <td>-1.9</td>
          <td>0.000030</td>
          <td>True</td>
          <td>06hr</td>
          <td>HILIC</td>
        </tr>
        <tr>
          <th>128164</th>
          <td>HMDB0012114</td>
          <td>(3S)-3,6-Diaminohexanoate</td>
          <td>metabolites</td>
          <td>-3.0</td>
          <td>0.002000</td>
          <td>True</td>
          <td>24hr</td>
          <td>HILIC</td>
        </tr>
        <tr>
          <th>128166</th>
          <td>HMDB0012115</td>
          <td>(3S,5S)-3,5-Diaminohexanoate</td>
          <td>metabolites</td>
          <td>-1.9</td>
          <td>0.000030</td>
          <td>True</td>
          <td>06hr</td>
          <td>HILIC</td>
        </tr>
      </tbody>
    </table>
    </div>



Significant filter
~~~~~~~~~~~~~~~~~~

We can use the ``significant`` column to filter that data to only
contain those species by applying ``.sig`` .

.. code:: ipython3

    exp_data.species.shape




.. parsed-literal::

    (132932, 8)



.. code:: ipython3

    exp_data.species.sig.shape




.. parsed-literal::

    (27288, 8)



.. code:: ipython3

    exp_data.label_free.sig.head(5)




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>identifier</th>
          <th>label</th>
          <th>species_type</th>
          <th>fold_change</th>
          <th>p_value</th>
          <th>significant</th>
          <th>sample_id</th>
          <th>source</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>102446</th>
          <td>LIMS1</td>
          <td>LIMS1_lf</td>
          <td>protein</td>
          <td>12.42</td>
          <td>0.00003</td>
          <td>True</td>
          <td>01hr</td>
          <td>label_free</td>
        </tr>
        <tr>
          <th>102447</th>
          <td>SMARCE1</td>
          <td>SMARCE1_lf</td>
          <td>protein</td>
          <td>-2.49</td>
          <td>0.00030</td>
          <td>True</td>
          <td>01hr</td>
          <td>label_free</td>
        </tr>
        <tr>
          <th>102448</th>
          <td>HEXA</td>
          <td>HEXA_lf</td>
          <td>protein</td>
          <td>6.42</td>
          <td>0.00060</td>
          <td>True</td>
          <td>01hr</td>
          <td>label_free</td>
        </tr>
        <tr>
          <th>102449</th>
          <td>SRSF1</td>
          <td>SRSF1_lf</td>
          <td>protein</td>
          <td>-3.21</td>
          <td>0.00060</td>
          <td>True</td>
          <td>01hr</td>
          <td>label_free</td>
        </tr>
        <tr>
          <th>102450</th>
          <td>SF3B1</td>
          <td>SF3B1_lf</td>
          <td>protein</td>
          <td>-1.57</td>
          <td>0.00130</td>
          <td>True</td>
          <td>01hr</td>
          <td>label_free</td>
        </tr>
      </tbody>
    </table>
    </div>



Filter data to up or down regulated species.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For enrichment analysis, we will want to access up-regulated and
down-regulated species using ``.up`` and ``.down``.

.. code:: ipython3

    exp_data.label_free.up.head(5)




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>identifier</th>
          <th>label</th>
          <th>species_type</th>
          <th>fold_change</th>
          <th>p_value</th>
          <th>significant</th>
          <th>sample_id</th>
          <th>source</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>102446</th>
          <td>LIMS1</td>
          <td>LIMS1_lf</td>
          <td>protein</td>
          <td>12.42</td>
          <td>0.00003</td>
          <td>True</td>
          <td>01hr</td>
          <td>label_free</td>
        </tr>
        <tr>
          <th>102448</th>
          <td>HEXA</td>
          <td>HEXA_lf</td>
          <td>protein</td>
          <td>6.42</td>
          <td>0.00060</td>
          <td>True</td>
          <td>01hr</td>
          <td>label_free</td>
        </tr>
        <tr>
          <th>102451</th>
          <td>USP15</td>
          <td>USP15_N-term A(ace)2_lf</td>
          <td>protein</td>
          <td>18.78</td>
          <td>0.00270</td>
          <td>True</td>
          <td>01hr</td>
          <td>label_free</td>
        </tr>
        <tr>
          <th>102460</th>
          <td>SBDS</td>
          <td>SBDS_lf</td>
          <td>protein</td>
          <td>2.79</td>
          <td>0.00560</td>
          <td>True</td>
          <td>01hr</td>
          <td>label_free</td>
        </tr>
        <tr>
          <th>102488</th>
          <td>CLIC4</td>
          <td>CLIC4_lf</td>
          <td>protein</td>
          <td>2.03</td>
          <td>0.01880</td>
          <td>True</td>
          <td>01hr</td>
          <td>label_free</td>
        </tr>
      </tbody>
    </table>
    </div>



.. code:: ipython3

    exp_data.label_free.down.head(5)




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>identifier</th>
          <th>label</th>
          <th>species_type</th>
          <th>fold_change</th>
          <th>p_value</th>
          <th>significant</th>
          <th>sample_id</th>
          <th>source</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>102447</th>
          <td>SMARCE1</td>
          <td>SMARCE1_lf</td>
          <td>protein</td>
          <td>-2.49</td>
          <td>0.0003</td>
          <td>True</td>
          <td>01hr</td>
          <td>label_free</td>
        </tr>
        <tr>
          <th>102449</th>
          <td>SRSF1</td>
          <td>SRSF1_lf</td>
          <td>protein</td>
          <td>-3.21</td>
          <td>0.0006</td>
          <td>True</td>
          <td>01hr</td>
          <td>label_free</td>
        </tr>
        <tr>
          <th>102450</th>
          <td>SF3B1</td>
          <td>SF3B1_lf</td>
          <td>protein</td>
          <td>-1.57</td>
          <td>0.0013</td>
          <td>True</td>
          <td>01hr</td>
          <td>label_free</td>
        </tr>
        <tr>
          <th>102452</th>
          <td>CKAP4</td>
          <td>CKAP4_lf</td>
          <td>protein</td>
          <td>-3.26</td>
          <td>0.0030</td>
          <td>True</td>
          <td>01hr</td>
          <td>label_free</td>
        </tr>
        <tr>
          <th>102453</th>
          <td>DDX17</td>
          <td>DDX17_lf</td>
          <td>protein</td>
          <td>-3.08</td>
          <td>0.0034</td>
          <td>True</td>
          <td>01hr</td>
          <td>label_free</td>
        </tr>
      </tbody>
    </table>
    </div>



Extracting by sample (time point)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We also added an index filter to segregate by ``sample_id``.

.. code:: ipython3

    for i in exp_data.sample_ids:
        print(i)
        display(exp_data[i].head(5))


.. parsed-literal::

    01hr
    


.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>identifier</th>
          <th>label</th>
          <th>species_type</th>
          <th>fold_change</th>
          <th>p_value</th>
          <th>significant</th>
          <th>sample_id</th>
          <th>source</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>19160</th>
          <td>GRIK4</td>
          <td>GRIK4_rnaseq</td>
          <td>protein</td>
          <td>77.555651</td>
          <td>0.019824</td>
          <td>True</td>
          <td>01hr</td>
          <td>rna_seq</td>
        </tr>
        <tr>
          <th>19161</th>
          <td>GRIK4_3p_UTR</td>
          <td>GRIK4_3p_UTR_rnaseq</td>
          <td>protein</td>
          <td>77.555651</td>
          <td>0.019824</td>
          <td>True</td>
          <td>01hr</td>
          <td>rna_seq</td>
        </tr>
        <tr>
          <th>19162</th>
          <td>AP001187.9</td>
          <td>AP001187.9_rnaseq</td>
          <td>protein</td>
          <td>-25.455050</td>
          <td>0.019824</td>
          <td>True</td>
          <td>01hr</td>
          <td>rna_seq</td>
        </tr>
        <tr>
          <th>19163</th>
          <td>MIR192</td>
          <td>MIR192_rnaseq</td>
          <td>protein</td>
          <td>-25.455050</td>
          <td>0.019824</td>
          <td>True</td>
          <td>01hr</td>
          <td>rna_seq</td>
        </tr>
        <tr>
          <th>19164</th>
          <td>MIR194-2</td>
          <td>MIR194-2_rnaseq</td>
          <td>protein</td>
          <td>-25.455050</td>
          <td>0.019824</td>
          <td>True</td>
          <td>01hr</td>
          <td>rna_seq</td>
        </tr>
      </tbody>
    </table>
    </div>


.. parsed-literal::

    06hr
    


.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>identifier</th>
          <th>label</th>
          <th>species_type</th>
          <th>fold_change</th>
          <th>p_value</th>
          <th>significant</th>
          <th>sample_id</th>
          <th>source</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>HOXD1</td>
          <td>HOXD1_rnaseq</td>
          <td>protein</td>
          <td>-520.256762</td>
          <td>0.00102</td>
          <td>True</td>
          <td>06hr</td>
          <td>rna_seq</td>
        </tr>
        <tr>
          <th>1</th>
          <td>MIR7704</td>
          <td>MIR7704_rnaseq</td>
          <td>protein</td>
          <td>-520.256762</td>
          <td>0.00102</td>
          <td>True</td>
          <td>06hr</td>
          <td>rna_seq</td>
        </tr>
        <tr>
          <th>2</th>
          <td>AC078814.1</td>
          <td>AC078814.1_rnaseq</td>
          <td>protein</td>
          <td>-76.022260</td>
          <td>0.00102</td>
          <td>True</td>
          <td>06hr</td>
          <td>rna_seq</td>
        </tr>
        <tr>
          <th>3</th>
          <td>PPM1H</td>
          <td>PPM1H_rnaseq</td>
          <td>protein</td>
          <td>-76.022260</td>
          <td>0.00102</td>
          <td>True</td>
          <td>06hr</td>
          <td>rna_seq</td>
        </tr>
        <tr>
          <th>4</th>
          <td>PLCH1</td>
          <td>PLCH1_rnaseq</td>
          <td>protein</td>
          <td>-17.888990</td>
          <td>0.00102</td>
          <td>True</td>
          <td>06hr</td>
          <td>rna_seq</td>
        </tr>
      </tbody>
    </table>
    </div>


.. parsed-literal::

    24hr
    


.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>identifier</th>
          <th>label</th>
          <th>species_type</th>
          <th>fold_change</th>
          <th>p_value</th>
          <th>significant</th>
          <th>sample_id</th>
          <th>source</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>37960</th>
          <td>LHX3</td>
          <td>LHX3_rnaseq</td>
          <td>protein</td>
          <td>202.225343</td>
          <td>0.005180</td>
          <td>True</td>
          <td>24hr</td>
          <td>rna_seq</td>
        </tr>
        <tr>
          <th>37961</th>
          <td>C17orf67</td>
          <td>C17orf67_rnaseq</td>
          <td>protein</td>
          <td>2.571464</td>
          <td>0.000123</td>
          <td>True</td>
          <td>24hr</td>
          <td>rna_seq</td>
        </tr>
        <tr>
          <th>37962</th>
          <td>ALX1</td>
          <td>ALX1_rnaseq</td>
          <td>protein</td>
          <td>-2.572587</td>
          <td>0.000123</td>
          <td>True</td>
          <td>24hr</td>
          <td>rna_seq</td>
        </tr>
        <tr>
          <th>37963</th>
          <td>MIR7844</td>
          <td>MIR7844_rnaseq</td>
          <td>protein</td>
          <td>2.573033</td>
          <td>0.009349</td>
          <td>True</td>
          <td>24hr</td>
          <td>rna_seq</td>
        </tr>
        <tr>
          <th>37964</th>
          <td>TMCC3</td>
          <td>TMCC3_rnaseq</td>
          <td>protein</td>
          <td>2.573033</td>
          <td>0.009349</td>
          <td>True</td>
          <td>24hr</td>
          <td>rna_seq</td>
        </tr>
      </tbody>
    </table>
    </div>


.. parsed-literal::

    48hr
    


.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>identifier</th>
          <th>label</th>
          <th>species_type</th>
          <th>fold_change</th>
          <th>p_value</th>
          <th>significant</th>
          <th>sample_id</th>
          <th>source</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>58025</th>
          <td>TNS3</td>
          <td>TNS3_1188_1197_phsilac</td>
          <td>protein</td>
          <td>-3.837129</td>
          <td>0.049</td>
          <td>True</td>
          <td>48hr</td>
          <td>ph_silac</td>
        </tr>
        <tr>
          <th>58026</th>
          <td>SIPA1L3</td>
          <td>SIPA1L3_S(ph)158_phsilac</td>
          <td>protein</td>
          <td>-5.119600</td>
          <td>0.049</td>
          <td>True</td>
          <td>48hr</td>
          <td>ph_silac</td>
        </tr>
        <tr>
          <th>58027</th>
          <td>TNS3</td>
          <td>TNS3_Y(ph)780_phsilac</td>
          <td>protein</td>
          <td>-4.986421</td>
          <td>0.049</td>
          <td>True</td>
          <td>48hr</td>
          <td>ph_silac</td>
        </tr>
        <tr>
          <th>58028</th>
          <td>FGD6</td>
          <td>FGD6_S(ph)554_phsilac</td>
          <td>protein</td>
          <td>-3.900705</td>
          <td>0.049</td>
          <td>True</td>
          <td>48hr</td>
          <td>ph_silac</td>
        </tr>
        <tr>
          <th>58029</th>
          <td>GPN1</td>
          <td>GPN1_S(ph)312_phsilac</td>
          <td>protein</td>
          <td>2.901199</td>
          <td>0.049</td>
          <td>True</td>
          <td>48hr</td>
          <td>ph_silac</td>
        </tr>
      </tbody>
    </table>
    </div>


.. code:: ipython3

    exp_data['01hr'].head(5)




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>identifier</th>
          <th>label</th>
          <th>species_type</th>
          <th>fold_change</th>
          <th>p_value</th>
          <th>significant</th>
          <th>sample_id</th>
          <th>source</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>19160</th>
          <td>GRIK4</td>
          <td>GRIK4_rnaseq</td>
          <td>protein</td>
          <td>77.555651</td>
          <td>0.019824</td>
          <td>True</td>
          <td>01hr</td>
          <td>rna_seq</td>
        </tr>
        <tr>
          <th>19161</th>
          <td>GRIK4_3p_UTR</td>
          <td>GRIK4_3p_UTR_rnaseq</td>
          <td>protein</td>
          <td>77.555651</td>
          <td>0.019824</td>
          <td>True</td>
          <td>01hr</td>
          <td>rna_seq</td>
        </tr>
        <tr>
          <th>19162</th>
          <td>AP001187.9</td>
          <td>AP001187.9_rnaseq</td>
          <td>protein</td>
          <td>-25.455050</td>
          <td>0.019824</td>
          <td>True</td>
          <td>01hr</td>
          <td>rna_seq</td>
        </tr>
        <tr>
          <th>19163</th>
          <td>MIR192</td>
          <td>MIR192_rnaseq</td>
          <td>protein</td>
          <td>-25.455050</td>
          <td>0.019824</td>
          <td>True</td>
          <td>01hr</td>
          <td>rna_seq</td>
        </tr>
        <tr>
          <th>19164</th>
          <td>MIR194-2</td>
          <td>MIR194-2_rnaseq</td>
          <td>protein</td>
          <td>-25.455050</td>
          <td>0.019824</td>
          <td>True</td>
          <td>01hr</td>
          <td>rna_seq</td>
        </tr>
      </tbody>
    </table>
    </div>



Pivot table to get table across time
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We also provide a function to quickly pivot the data to for easy export.

.. code:: ipython3

    exp_data.label_free.pivoter(
        convert_to_log=False, 
        index='identifier',
        columns='sample_id',
        values=['fold_change', 'p_value']
    ).head(10)




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead tr th {
            text-align: left;
        }
    
        .dataframe thead tr:last-of-type th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr>
          <th></th>
          <th colspan="4" halign="left">fold_change</th>
          <th colspan="4" halign="left">p_value</th>
        </tr>
        <tr>
          <th>sample_id</th>
          <th>01hr</th>
          <th>06hr</th>
          <th>24hr</th>
          <th>48hr</th>
          <th>01hr</th>
          <th>06hr</th>
          <th>24hr</th>
          <th>48hr</th>
        </tr>
        <tr>
          <th>identifier</th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>A2M</th>
          <td>1.040000</td>
          <td>1.140</td>
          <td>51.93</td>
          <td>11.58</td>
          <td>0.514800</td>
          <td>0.44370</td>
          <td>0.24260</td>
          <td>0.11130</td>
        </tr>
        <tr>
          <th>AACS</th>
          <td>-1.100000</td>
          <td>3.740</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>0.281800</td>
          <td>0.26950</td>
          <td>NaN</td>
          <td>NaN</td>
        </tr>
        <tr>
          <th>AAGAB</th>
          <td>1.000000</td>
          <td>-1.150</td>
          <td>1.46</td>
          <td>-2.03</td>
          <td>0.968100</td>
          <td>0.39240</td>
          <td>0.84450</td>
          <td>0.09760</td>
        </tr>
        <tr>
          <th>AAK1</th>
          <td>1.320000</td>
          <td>1.590</td>
          <td>NaN</td>
          <td>1.72</td>
          <td>0.715800</td>
          <td>0.18110</td>
          <td>NaN</td>
          <td>0.95660</td>
        </tr>
        <tr>
          <th>AAMP</th>
          <td>-1.200000</td>
          <td>-1.460</td>
          <td>1.85</td>
          <td>1.78</td>
          <td>0.836800</td>
          <td>0.55420</td>
          <td>0.13640</td>
          <td>0.32460</td>
        </tr>
        <tr>
          <th>AAR2</th>
          <td>NaN</td>
          <td>-1.690</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>0.96510</td>
          <td>NaN</td>
          <td>NaN</td>
        </tr>
        <tr>
          <th>AARS</th>
          <td>0.326667</td>
          <td>-0.035</td>
          <td>-1.44</td>
          <td>-3.12</td>
          <td>0.299867</td>
          <td>0.62425</td>
          <td>0.46725</td>
          <td>0.00045</td>
        </tr>
        <tr>
          <th>AARS2</th>
          <td>1.170000</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>0.253000</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
        </tr>
        <tr>
          <th>AARSD1</th>
          <td>1.210000</td>
          <td>4.070</td>
          <td>-2.05</td>
          <td>NaN</td>
          <td>0.459700</td>
          <td>0.49160</td>
          <td>0.78440</td>
          <td>NaN</td>
        </tr>
        <tr>
          <th>AASDHPPT</th>
          <td>-0.330000</td>
          <td>1.020</td>
          <td>1.07</td>
          <td>-1.11</td>
          <td>0.709600</td>
          <td>0.81160</td>
          <td>0.45290</td>
          <td>0.00070</td>
        </tr>
      </tbody>
    </table>
    </div>



Note that in the previous two examples, we find that there are NaN
values. This is because of our experiental data. We can easy check what
species are not found in all 4 of our label free experiements.

.. code:: ipython3

    print(len(exp_data.label_free.present_in_all_columns(
        index='identifier',
        columns='sample_id',
    ).id_list))


.. parsed-literal::

    Number in index went from 3447 to 1819
    1819
    

This shows that out of the 3447 unique species measured in label-free
proteomics, only 1819 were measured in all time points. What one can do
with this information is dependent on the analysis. We can filter by
requiring a species to be signficantly changed in at least ``n``
samples.

.. code:: ipython3

    print(exp_data.label_free.require_n_sig(n_sig=2).identifier.unique().shape)
    print(exp_data.label_free.require_n_sig(n_sig=3).identifier.unique().shape)
    print(exp_data.label_free.require_n_sig(n_sig=4).identifier.unique().shape)


.. parsed-literal::

    (247,)
    (53,)
    (2,)
    

It is important to note that this class is basically a hopped up
pandas.DataFrame, so the commands can be chained together.

.. code:: ipython3

    exp_data.label_free.require_n_sig(n_sig=3).pivoter(
        convert_to_log=False, 
        index='identifier',
        columns='sample_id',
        values=['fold_change', 'p_value']
    ).head(10)




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead tr th {
            text-align: left;
        }
    
        .dataframe thead tr:last-of-type th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr>
          <th></th>
          <th colspan="4" halign="left">fold_change</th>
          <th colspan="4" halign="left">p_value</th>
        </tr>
        <tr>
          <th>sample_id</th>
          <th>01hr</th>
          <th>06hr</th>
          <th>24hr</th>
          <th>48hr</th>
          <th>01hr</th>
          <th>06hr</th>
          <th>24hr</th>
          <th>48hr</th>
        </tr>
        <tr>
          <th>identifier</th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>COX5A</th>
          <td>-2.185000</td>
          <td>-1.550</td>
          <td>-1.610</td>
          <td>5.480</td>
          <td>0.121850</td>
          <td>0.07850</td>
          <td>0.11900</td>
          <td>2.800000e-03</td>
        </tr>
        <tr>
          <th>DAZAP1</th>
          <td>-1.860000</td>
          <td>1.240</td>
          <td>2.040</td>
          <td>2.740</td>
          <td>0.364050</td>
          <td>0.56310</td>
          <td>0.08610</td>
          <td>1.100000e-03</td>
        </tr>
        <tr>
          <th>DDX3X</th>
          <td>-1.940000</td>
          <td>1.090</td>
          <td>-1.620</td>
          <td>-3.030</td>
          <td>0.089000</td>
          <td>0.70330</td>
          <td>0.05480</td>
          <td>2.700000e-03</td>
        </tr>
        <tr>
          <th>DDX5</th>
          <td>-7.623333</td>
          <td>1.090</td>
          <td>-2.085</td>
          <td>-6.160</td>
          <td>0.094733</td>
          <td>0.59080</td>
          <td>0.00725</td>
          <td>2.900000e-03</td>
        </tr>
        <tr>
          <th>ERH</th>
          <td>-2.340000</td>
          <td>-1.095</td>
          <td>1.440</td>
          <td>2.020</td>
          <td>0.073300</td>
          <td>0.48190</td>
          <td>0.06335</td>
          <td>1.885000e-02</td>
        </tr>
        <tr>
          <th>FUS</th>
          <td>-2.415000</td>
          <td>1.200</td>
          <td>1.940</td>
          <td>2.620</td>
          <td>0.096850</td>
          <td>0.08740</td>
          <td>0.03320</td>
          <td>6.000000e-04</td>
        </tr>
        <tr>
          <th>GPRC5A</th>
          <td>-2.450000</td>
          <td>-1.160</td>
          <td>1.800</td>
          <td>31.455</td>
          <td>0.099200</td>
          <td>0.09660</td>
          <td>0.05520</td>
          <td>4.000000e-04</td>
        </tr>
        <tr>
          <th>H1FX</th>
          <td>-1.660000</td>
          <td>-1.320</td>
          <td>-1.540</td>
          <td>-6.340</td>
          <td>0.095700</td>
          <td>0.00250</td>
          <td>0.01970</td>
          <td>9.345000e-07</td>
        </tr>
        <tr>
          <th>HNRNPA2B1</th>
          <td>-4.365000</td>
          <td>1.360</td>
          <td>2.600</td>
          <td>4.960</td>
          <td>0.050550</td>
          <td>0.09280</td>
          <td>0.05390</td>
          <td>2.000000e-04</td>
        </tr>
        <tr>
          <th>HNRNPC</th>
          <td>-4.543333</td>
          <td>4.170</td>
          <td>2.270</td>
          <td>4.920</td>
          <td>0.070000</td>
          <td>0.15085</td>
          <td>0.12120</td>
          <td>3.800000e-03</td>
        </tr>
      </tbody>
    </table>
    </div>



Visualization
~~~~~~~~~~~~~

We provide commonly used plotting functions. ``.volcano_plot``
``.volcano_by_sample`` ``.plot_histogram`` ``.plot_species``
``.heatmap``

Volcano plots
^^^^^^^^^^^^^

.. code:: ipython3

    exp_data.label_free.volcano_plot();



.. image:: Tutorial_files/Tutorial_48_0.png


.. code:: ipython3

    exp_data.label_free.volcano_by_sample(sig_column=True);



.. image:: Tutorial_files/Tutorial_49_0.png


.. code:: ipython3

    exp_data.label_free.plot_histogram();



.. image:: Tutorial_files/Tutorial_50_0.png


Plotting subset of species
~~~~~~~~~~~~~~~~~~~~~~~~~~

We provide the a few plotting interfaces to explore that subsets of the
data. Basically, you create a list of species and provide it to the
function. It filters based on these and then returns the results.

Time series using plot’y or matplotlib
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    exp_data.label_free.plot_species(['LMNA', 'VDAC1'], plot_type='plotly')



.. raw:: html

    <script type="text/javascript">
    window.PlotlyConfig = {MathJaxConfig: 'local'};
    if (window.MathJax) {MathJax.Hub.Config({SVG: {font: "STIX-Web"}});}
    if (typeof require !== 'undefined') {
    require.undef("plotly");
    requirejs.config({
        paths: {
            'plotly': ['https://cdn.plot.ly/plotly-2.4.1.min']
        }
    });
    require(['plotly'], function(Plotly) {
        window._Plotly = Plotly;
    });
    }
    </script>
    



.. raw:: html

    <div>                            <div id="aa2a77bc-663f-48bb-9218-a2e95921975c" class="plotly-graph-div" style="height:525px; width:100%;"></div>            <script type="text/javascript">                require(["plotly"], function(Plotly) {                    window.PLOTLYENV=window.PLOTLYENV || {};                                    if (document.getElementById("aa2a77bc-663f-48bb-9218-a2e95921975c")) {                    Plotly.newPlot(                        "aa2a77bc-663f-48bb-9218-a2e95921975c",                        [{"hoveron":"points","legendgroup":"group_LMNA_S(ph)22_lf","line":{"color":"rgba(0.12156862745098039,0.4666666666666667,0.7058823529411765,1.)"},"marker":{"color":"rgba(0.12156862745098039,0.4666666666666667,0.7058823529411765,1.)","size":8,"symbol":"circle"},"mode":"lines+markers","name":"LMNA_S(ph)22_lf","showlegend":true,"type":"scatter","visible":true,"x":[0,1,3],"y":[0.02856915219677092,-0.5849625007211562,1.4489009511451278]},{"hoveron":"points","legendgroup":"group_LMNA_S(ph)22_lf","line":{"color":"rgba(0.12156862745098039,0.4666666666666667,0.7058823529411765,1.)"},"marker":{"color":"rgba(0.12156862745098039,0.4666666666666667,0.7058823529411765,1.)","size":12,"symbol":"x-open-dot"},"mode":"markers","name":"LMNA_S(ph)22_lf","showlegend":false,"type":"scatter","visible":true,"x":[1,3],"y":[-0.5849625007211562,1.4489009511451278]},{"hoveron":"points","legendgroup":"group_LMNA_S(ph)392_lf","line":{"color":"rgba(0.6823529411764706,0.7803921568627451,0.9098039215686274,1.)"},"marker":{"color":"rgba(0.6823529411764706,0.7803921568627451,0.9098039215686274,1.)","size":8,"symbol":"circle"},"mode":"lines+markers","name":"LMNA_S(ph)392_lf","showlegend":true,"type":"scatter","visible":true,"x":[0],"y":[-2.7398481026993275]},{"hoveron":"points","legendgroup":"group_LMNA_S(ph)392_lf","line":{"color":"rgba(0.6823529411764706,0.7803921568627451,0.9098039215686274,1.)"},"marker":{"color":"rgba(0.6823529411764706,0.7803921568627451,0.9098039215686274,1.)","size":12,"symbol":"x-open-dot"},"mode":"markers","name":"LMNA_S(ph)392_lf","showlegend":false,"type":"scatter","visible":true,"x":[],"y":[]},{"hoveron":"points","legendgroup":"group_LMNA_lf","line":{"color":"rgba(1.0,0.4980392156862745,0.054901960784313725,1.)"},"marker":{"color":"rgba(1.0,0.4980392156862745,0.054901960784313725,1.)","size":8,"symbol":"circle"},"mode":"lines+markers","name":"LMNA_lf","showlegend":true,"type":"scatter","visible":true,"x":[0,0,1,2,3],"y":[-1.835924074254375,-1.2141248053528473,-0.36737106564852945,0.9781956296816515,1.7484612330040354]},{"hoveron":"points","legendgroup":"group_LMNA_lf","line":{"color":"rgba(1.0,0.4980392156862745,0.054901960784313725,1.)"},"marker":{"color":"rgba(1.0,0.4980392156862745,0.054901960784313725,1.)","size":12,"symbol":"x-open-dot"},"mode":"markers","name":"LMNA_lf","showlegend":false,"type":"scatter","visible":true,"x":[0,0,2,3],"y":[-1.835924074254375,-1.2141248053528473,0.9781956296816515,1.7484612330040354]},{"hoveron":"points","legendgroup":"group_VDAC1_N-term A(ace)2_lf","line":{"color":"rgba(1.0,0.7333333333333333,0.47058823529411764,1.)"},"marker":{"color":"rgba(1.0,0.7333333333333333,0.47058823529411764,1.)","size":8,"symbol":"circle"},"mode":"lines+markers","name":"VDAC1_N-term A(ace)2_lf","showlegend":true,"type":"scatter","visible":true,"x":[0],"y":[-2.920293300211007]},{"hoveron":"points","legendgroup":"group_VDAC1_N-term A(ace)2_lf","line":{"color":"rgba(1.0,0.7333333333333333,0.47058823529411764,1.)"},"marker":{"color":"rgba(1.0,0.7333333333333333,0.47058823529411764,1.)","size":12,"symbol":"x-open-dot"},"mode":"markers","name":"VDAC1_N-term A(ace)2_lf","showlegend":false,"type":"scatter","visible":true,"x":[],"y":[]},{"hoveron":"points","legendgroup":"group_VDAC1_lf","line":{"color":"rgba(0.17254901960784313,0.6274509803921569,0.17254901960784313,1.)"},"marker":{"color":"rgba(0.17254901960784313,0.6274509803921569,0.17254901960784313,1.)","size":8,"symbol":"circle"},"mode":"lines+markers","name":"VDAC1_lf","showlegend":true,"type":"scatter","visible":true,"x":[0,0,1,2,3],"y":[-1.9927684307689242,-1.480265122054463,0.25096157353321874,2.526068811667588,5.8938475581329195]},{"hoveron":"points","legendgroup":"group_VDAC1_lf","line":{"color":"rgba(0.17254901960784313,0.6274509803921569,0.17254901960784313,1.)"},"marker":{"color":"rgba(0.17254901960784313,0.6274509803921569,0.17254901960784313,1.)","size":12,"symbol":"x-open-dot"},"mode":"markers","name":"VDAC1_lf","showlegend":false,"type":"scatter","visible":true,"x":[0,2,3],"y":[-1.9927684307689242,2.526068811667588,5.8938475581329195]}],                        {"hovermode":"closest","showlegend":true,"template":{"data":{"bar":[{"error_x":{"color":"#2a3f5f"},"error_y":{"color":"#2a3f5f"},"marker":{"line":{"color":"#E5ECF6","width":0.5},"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"bar"}],"barpolar":[{"marker":{"line":{"color":"#E5ECF6","width":0.5},"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"barpolar"}],"carpet":[{"aaxis":{"endlinecolor":"#2a3f5f","gridcolor":"white","linecolor":"white","minorgridcolor":"white","startlinecolor":"#2a3f5f"},"baxis":{"endlinecolor":"#2a3f5f","gridcolor":"white","linecolor":"white","minorgridcolor":"white","startlinecolor":"#2a3f5f"},"type":"carpet"}],"choropleth":[{"colorbar":{"outlinewidth":0,"ticks":""},"type":"choropleth"}],"contour":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"contour"}],"contourcarpet":[{"colorbar":{"outlinewidth":0,"ticks":""},"type":"contourcarpet"}],"heatmap":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"heatmap"}],"heatmapgl":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"heatmapgl"}],"histogram":[{"marker":{"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"histogram"}],"histogram2d":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"histogram2d"}],"histogram2dcontour":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"histogram2dcontour"}],"mesh3d":[{"colorbar":{"outlinewidth":0,"ticks":""},"type":"mesh3d"}],"parcoords":[{"line":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"parcoords"}],"pie":[{"automargin":true,"type":"pie"}],"scatter":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatter"}],"scatter3d":[{"line":{"colorbar":{"outlinewidth":0,"ticks":""}},"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatter3d"}],"scattercarpet":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scattercarpet"}],"scattergeo":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scattergeo"}],"scattergl":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scattergl"}],"scattermapbox":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scattermapbox"}],"scatterpolar":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatterpolar"}],"scatterpolargl":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatterpolargl"}],"scatterternary":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatterternary"}],"surface":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"surface"}],"table":[{"cells":{"fill":{"color":"#EBF0F8"},"line":{"color":"white"}},"header":{"fill":{"color":"#C8D4E3"},"line":{"color":"white"}},"type":"table"}]},"layout":{"annotationdefaults":{"arrowcolor":"#2a3f5f","arrowhead":0,"arrowwidth":1},"autotypenumbers":"strict","coloraxis":{"colorbar":{"outlinewidth":0,"ticks":""}},"colorscale":{"diverging":[[0,"#8e0152"],[0.1,"#c51b7d"],[0.2,"#de77ae"],[0.3,"#f1b6da"],[0.4,"#fde0ef"],[0.5,"#f7f7f7"],[0.6,"#e6f5d0"],[0.7,"#b8e186"],[0.8,"#7fbc41"],[0.9,"#4d9221"],[1,"#276419"]],"sequential":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"sequentialminus":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]},"colorway":["#636efa","#EF553B","#00cc96","#ab63fa","#FFA15A","#19d3f3","#FF6692","#B6E880","#FF97FF","#FECB52"],"font":{"color":"#2a3f5f"},"geo":{"bgcolor":"white","lakecolor":"white","landcolor":"#E5ECF6","showlakes":true,"showland":true,"subunitcolor":"white"},"hoverlabel":{"align":"left"},"hovermode":"closest","mapbox":{"style":"light"},"paper_bgcolor":"white","plot_bgcolor":"#E5ECF6","polar":{"angularaxis":{"gridcolor":"white","linecolor":"white","ticks":""},"bgcolor":"#E5ECF6","radialaxis":{"gridcolor":"white","linecolor":"white","ticks":""}},"scene":{"xaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","gridwidth":2,"linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white"},"yaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","gridwidth":2,"linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white"},"zaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","gridwidth":2,"linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white"}},"shapedefaults":{"line":{"color":"#2a3f5f"}},"ternary":{"aaxis":{"gridcolor":"white","linecolor":"white","ticks":""},"baxis":{"gridcolor":"white","linecolor":"white","ticks":""},"bgcolor":"#E5ECF6","caxis":{"gridcolor":"white","linecolor":"white","ticks":""}},"title":{"x":0.05},"xaxis":{"automargin":true,"gridcolor":"white","linecolor":"white","ticks":"","title":{"standoff":15},"zerolinecolor":"white","zerolinewidth":2},"yaxis":{"automargin":true,"gridcolor":"white","linecolor":"white","ticks":"","title":{"standoff":15},"zerolinecolor":"white","zerolinewidth":2}}},"updatemenus":[{"buttons":[{"args":["visible",[true,true,true,true,true,true,true,true,true,true]],"label":"All","method":"restyle"},{"args":["visible",[true,true,true,true,true,true,false,false,false,false]],"label":"LMNA","method":"restyle"},{"args":["visible",[false,false,false,false,false,false,true,true,true,true]],"label":"VDAC1","method":"restyle"}],"x":-0.05,"y":1,"yanchor":"top"}],"xaxis":{"range":[0,3],"showticklabels":true,"tickmode":"array","ticktext":["01hr","06hr","24hr","48hr"],"tickvals":[0,1,2,3],"title":{"text":"Sample index"}},"yaxis":{"title":{"text":"log2fc"}}},                        {"responsive": true}                    ).then(function(){
    
    var gd = document.getElementById('aa2a77bc-663f-48bb-9218-a2e95921975c');
    var x = new MutationObserver(function (mutations, observer) {{
            var display = window.getComputedStyle(gd).display;
            if (!display || display === 'none') {{
                console.log([gd, 'removed!']);
                Plotly.purge(gd);
                observer.disconnect();
            }}
    }});
    
    // Listen for the removal of the full notebook cells
    var notebookContainer = gd.closest('#notebook-container');
    if (notebookContainer) {{
        x.observe(notebookContainer, {childList: true});
    }}
    
    // Listen for the clearing of the current output cell
    var outputEl = gd.closest('.output');
    if (outputEl) {{
        x.observe(outputEl, {childList: true});
    }}
    
                            })                };                });            </script>        </div>


.. code:: ipython3

    exp_data.label_free.plot_species(['LMNA', 'VDAC1'], plot_type='matplotlib');



.. image:: Tutorial_files/Tutorial_54_0.png


Heatplots
^^^^^^^^^

.. code:: ipython3

    exp_data.label_free.heatmap(
        ['LMNA', 'VDAC1'], 
        figsize=(6,4), 
        linewidths=0.01
    );



.. image:: Tutorial_files/Tutorial_56_0.png


Notice that the above plot doesn’t show any of the modifiers of LMBA (no
\_s(ph)22_lf). This is because the default index to pivot plots is the
``identifier`` column. You can set the ``label`` column for plotting by
passing index=\ ``label`` to the function. Note, if you want to filter
the data using the more generic ‘identifier’ column, you just specify
that with subset_index=‘identifier’

.. code:: ipython3

    exp_data.label_free.heatmap(
        ['LMNA', 'VDAC1'], 
        subset_index='identifier', 
        index='label',
        figsize=(6,4), 
        linewidths=0.01
    );



.. image:: Tutorial_files/Tutorial_58_0.png


Examples
~~~~~~~~

Here are a few examples how all the above commands can be chained
together to create plots with varying degrees of critera.

Query 1:
^^^^^^^^

::

   Heatmap of label-free proteomics that are signficantly change in at least 3 time points.

.. code:: ipython3

    lf_sig = exp_data.label_free.require_n_sig(
        index='label', 
        columns='sample_id', 
        n_sig=3
    ).heatmap(
        convert_to_log=True, 
        cluster_row=True, 
        index='label',
        values='fold_change', 
        columns='sample_id', 
        annotate_sig=True, 
        figsize=(8, 12), 
        div_colors=True,
        num_colors=21, 
        linewidths=0.01
    );



.. image:: Tutorial_files/Tutorial_60_0.png


Query 2:
^^^^^^^^

::

   Changes that happen at all 3 timepoints for RNA-seq.

.. code:: ipython3

    exp_data.rna.require_n_sig(n_sig=3, index='label').plot_species(plot_type='plotly');



.. raw:: html

    <script type="text/javascript">
    window.PlotlyConfig = {MathJaxConfig: 'local'};
    if (window.MathJax) {MathJax.Hub.Config({SVG: {font: "STIX-Web"}});}
    if (typeof require !== 'undefined') {
    require.undef("plotly");
    requirejs.config({
        paths: {
            'plotly': ['https://cdn.plot.ly/plotly-2.4.1.min']
        }
    });
    require(['plotly'], function(Plotly) {
        window._Plotly = Plotly;
    });
    }
    </script>
    



.. raw:: html

    <div>                            <div id="c19a0e71-9da0-491b-96c7-15eb7a6b27d3" class="plotly-graph-div" style="height:525px; width:100%;"></div>            <script type="text/javascript">                require(["plotly"], function(Plotly) {                    window.PLOTLYENV=window.PLOTLYENV || {};                                    if (document.getElementById("c19a0e71-9da0-491b-96c7-15eb7a6b27d3")) {                    Plotly.newPlot(                        "c19a0e71-9da0-491b-96c7-15eb7a6b27d3",                        [{"hoveron":"points","legendgroup":"group_AC009487.6_rnaseq","line":{"color":"rgba(0.12156862745098039,0.4666666666666667,0.7058823529411765,1.)"},"marker":{"color":"rgba(0.12156862745098039,0.4666666666666667,0.7058823529411765,1.)","size":8,"symbol":"circle"},"mode":"lines+markers","name":"AC009487.6_rnaseq","showlegend":true,"type":"scatter","visible":true,"x":[0,1,2],"y":[-0.35562200000175564,-0.5685870000046604,-0.9286870000014209]},{"hoveron":"points","legendgroup":"group_AC009487.6_rnaseq","line":{"color":"rgba(0.12156862745098039,0.4666666666666667,0.7058823529411765,1.)"},"marker":{"color":"rgba(0.12156862745098039,0.4666666666666667,0.7058823529411765,1.)","size":12,"symbol":"x-open-dot"},"mode":"markers","name":"AC009487.6_rnaseq","showlegend":false,"type":"scatter","visible":true,"x":[0,1,2],"y":[-0.35562200000175564,-0.5685870000046604,-0.9286870000014209]},{"hoveron":"points","legendgroup":"group_C11orf68_rnaseq","line":{"color":"rgba(0.6823529411764706,0.7803921568627451,0.9098039215686274,1.)"},"marker":{"color":"rgba(0.6823529411764706,0.7803921568627451,0.9098039215686274,1.)","size":8,"symbol":"circle"},"mode":"lines+markers","name":"C11orf68_rnaseq","showlegend":true,"type":"scatter","visible":true,"x":[0,1,2],"y":[0.3330860000036713,0.792480000000864,2.744710000000026]},{"hoveron":"points","legendgroup":"group_C11orf68_rnaseq","line":{"color":"rgba(0.6823529411764706,0.7803921568627451,0.9098039215686274,1.)"},"marker":{"color":"rgba(0.6823529411764706,0.7803921568627451,0.9098039215686274,1.)","size":12,"symbol":"x-open-dot"},"mode":"markers","name":"C11orf68_rnaseq","showlegend":false,"type":"scatter","visible":true,"x":[0,1,2],"y":[0.3330860000036713,0.792480000000864,2.744710000000026]},{"hoveron":"points","legendgroup":"group_CITED2_rnaseq","line":{"color":"rgba(1.0,0.4980392156862745,0.054901960784313725,1.)"},"marker":{"color":"rgba(1.0,0.4980392156862745,0.054901960784313725,1.)","size":8,"symbol":"circle"},"mode":"lines+markers","name":"CITED2_rnaseq","showlegend":true,"type":"scatter","visible":true,"x":[0,1,2],"y":[0.5648999999955752,1.03988000000217,3.0135899999999594]},{"hoveron":"points","legendgroup":"group_CITED2_rnaseq","line":{"color":"rgba(1.0,0.4980392156862745,0.054901960784313725,1.)"},"marker":{"color":"rgba(1.0,0.4980392156862745,0.054901960784313725,1.)","size":12,"symbol":"x-open-dot"},"mode":"markers","name":"CITED2_rnaseq","showlegend":false,"type":"scatter","visible":true,"x":[0,1,2],"y":[0.5648999999955752,1.03988000000217,3.0135899999999594]},{"hoveron":"points","legendgroup":"group_FOXC1_rnaseq","line":{"color":"rgba(1.0,0.7333333333333333,0.47058823529411764,1.)"},"marker":{"color":"rgba(1.0,0.7333333333333333,0.47058823529411764,1.)","size":8,"symbol":"circle"},"mode":"lines+markers","name":"FOXC1_rnaseq","showlegend":true,"type":"scatter","visible":true,"x":[0,1,2],"y":[0.3772070000010807,0.5609149999964104,0.49152899999906363]},{"hoveron":"points","legendgroup":"group_FOXC1_rnaseq","line":{"color":"rgba(1.0,0.7333333333333333,0.47058823529411764,1.)"},"marker":{"color":"rgba(1.0,0.7333333333333333,0.47058823529411764,1.)","size":12,"symbol":"x-open-dot"},"mode":"markers","name":"FOXC1_rnaseq","showlegend":false,"type":"scatter","visible":true,"x":[0,1,2],"y":[0.3772070000010807,0.5609149999964104,0.49152899999906363]},{"hoveron":"points","legendgroup":"group_PLCH1_rnaseq","line":{"color":"rgba(0.17254901960784313,0.6274509803921569,0.17254901960784313,1.)"},"marker":{"color":"rgba(0.17254901960784313,0.6274509803921569,0.17254901960784313,1.)","size":8,"symbol":"circle"},"mode":"lines+markers","name":"PLCH1_rnaseq","showlegend":true,"type":"scatter","visible":true,"x":[0,1,2],"y":[2.4764699999998228,-4.161000000003726,-6.711799999997806]},{"hoveron":"points","legendgroup":"group_PLCH1_rnaseq","line":{"color":"rgba(0.17254901960784313,0.6274509803921569,0.17254901960784313,1.)"},"marker":{"color":"rgba(0.17254901960784313,0.6274509803921569,0.17254901960784313,1.)","size":12,"symbol":"x-open-dot"},"mode":"markers","name":"PLCH1_rnaseq","showlegend":false,"type":"scatter","visible":true,"x":[0,1,2],"y":[2.4764699999998228,-4.161000000003726,-6.711799999997806]},{"hoveron":"points","legendgroup":"group_RBMX_rnaseq","line":{"color":"rgba(0.596078431372549,0.8745098039215686,0.5411764705882353,1.)"},"marker":{"color":"rgba(0.596078431372549,0.8745098039215686,0.5411764705882353,1.)","size":8,"symbol":"circle"},"mode":"lines+markers","name":"RBMX_rnaseq","showlegend":true,"type":"scatter","visible":true,"x":[0,1,2],"y":[1.412020000002594,-0.44532300000365466,-1.4826299999998749]},{"hoveron":"points","legendgroup":"group_RBMX_rnaseq","line":{"color":"rgba(0.596078431372549,0.8745098039215686,0.5411764705882353,1.)"},"marker":{"color":"rgba(0.596078431372549,0.8745098039215686,0.5411764705882353,1.)","size":12,"symbol":"x-open-dot"},"mode":"markers","name":"RBMX_rnaseq","showlegend":false,"type":"scatter","visible":true,"x":[0,1,2],"y":[1.412020000002594,-0.44532300000365466,-1.4826299999998749]},{"hoveron":"points","legendgroup":"group_RP11-639F1.1_rnaseq","line":{"color":"rgba(0.8392156862745098,0.15294117647058825,0.1568627450980392,1.)"},"marker":{"color":"rgba(0.8392156862745098,0.15294117647058825,0.1568627450980392,1.)","size":8,"symbol":"circle"},"mode":"lines+markers","name":"RP11-639F1.1_rnaseq","showlegend":true,"type":"scatter","visible":true,"x":[0,1,2],"y":[2.4764699999998228,-4.161000000003726,-6.711799999997806]},{"hoveron":"points","legendgroup":"group_RP11-639F1.1_rnaseq","line":{"color":"rgba(0.8392156862745098,0.15294117647058825,0.1568627450980392,1.)"},"marker":{"color":"rgba(0.8392156862745098,0.15294117647058825,0.1568627450980392,1.)","size":12,"symbol":"x-open-dot"},"mode":"markers","name":"RP11-639F1.1_rnaseq","showlegend":false,"type":"scatter","visible":true,"x":[0,1,2],"y":[2.4764699999998228,-4.161000000003726,-6.711799999997806]},{"hoveron":"points","legendgroup":"group_SESN1_rnaseq","line":{"color":"rgba(1.0,0.596078431372549,0.5882352941176471,1.)"},"marker":{"color":"rgba(1.0,0.596078431372549,0.5882352941176471,1.)","size":8,"symbol":"circle"},"mode":"lines+markers","name":"SESN1_rnaseq","showlegend":true,"type":"scatter","visible":true,"x":[0,1,2],"y":[-0.5685350000037165,1.4627100000004496,3.484260000006248]},{"hoveron":"points","legendgroup":"group_SESN1_rnaseq","line":{"color":"rgba(1.0,0.596078431372549,0.5882352941176471,1.)"},"marker":{"color":"rgba(1.0,0.596078431372549,0.5882352941176471,1.)","size":12,"symbol":"x-open-dot"},"mode":"markers","name":"SESN1_rnaseq","showlegend":false,"type":"scatter","visible":true,"x":[0,1,2],"y":[-0.5685350000037165,1.4627100000004496,3.484260000006248]},{"hoveron":"points","legendgroup":"group_SNHG1_rnaseq","line":{"color":"rgba(0.5803921568627451,0.403921568627451,0.7411764705882353,1.)"},"marker":{"color":"rgba(0.5803921568627451,0.403921568627451,0.7411764705882353,1.)","size":8,"symbol":"circle"},"mode":"lines+markers","name":"SNHG1_rnaseq","showlegend":true,"type":"scatter","visible":true,"x":[0,1,2],"y":[-1.784870000001268,2.3110999999993025,3.7938299999958955]},{"hoveron":"points","legendgroup":"group_SNHG1_rnaseq","line":{"color":"rgba(0.5803921568627451,0.403921568627451,0.7411764705882353,1.)"},"marker":{"color":"rgba(0.5803921568627451,0.403921568627451,0.7411764705882353,1.)","size":12,"symbol":"x-open-dot"},"mode":"markers","name":"SNHG1_rnaseq","showlegend":false,"type":"scatter","visible":true,"x":[0,1,2],"y":[-1.784870000001268,2.3110999999993025,3.7938299999958955]},{"hoveron":"points","legendgroup":"group_SNORD22_rnaseq","line":{"color":"rgba(0.7725490196078432,0.6901960784313725,0.8352941176470589,1.)"},"marker":{"color":"rgba(0.7725490196078432,0.6901960784313725,0.8352941176470589,1.)","size":8,"symbol":"circle"},"mode":"lines+markers","name":"SNORD22_rnaseq","showlegend":true,"type":"scatter","visible":true,"x":[0,1,2],"y":[-1.784870000001268,2.3110999999993025,3.7938299999958955]},{"hoveron":"points","legendgroup":"group_SNORD22_rnaseq","line":{"color":"rgba(0.7725490196078432,0.6901960784313725,0.8352941176470589,1.)"},"marker":{"color":"rgba(0.7725490196078432,0.6901960784313725,0.8352941176470589,1.)","size":12,"symbol":"x-open-dot"},"mode":"markers","name":"SNORD22_rnaseq","showlegend":false,"type":"scatter","visible":true,"x":[0,1,2],"y":[-1.784870000001268,2.3110999999993025,3.7938299999958955]},{"hoveron":"points","legendgroup":"group_SNORD25_rnaseq","line":{"color":"rgba(0.5490196078431373,0.33725490196078434,0.29411764705882354,1.)"},"marker":{"color":"rgba(0.5490196078431373,0.33725490196078434,0.29411764705882354,1.)","size":8,"symbol":"circle"},"mode":"lines+markers","name":"SNORD25_rnaseq","showlegend":true,"type":"scatter","visible":true,"x":[0,1,2],"y":[-1.784870000001268,2.3110999999993025,3.7938299999958955]},{"hoveron":"points","legendgroup":"group_SNORD25_rnaseq","line":{"color":"rgba(0.5490196078431373,0.33725490196078434,0.29411764705882354,1.)"},"marker":{"color":"rgba(0.5490196078431373,0.33725490196078434,0.29411764705882354,1.)","size":12,"symbol":"x-open-dot"},"mode":"markers","name":"SNORD25_rnaseq","showlegend":false,"type":"scatter","visible":true,"x":[0,1,2],"y":[-1.784870000001268,2.3110999999993025,3.7938299999958955]},{"hoveron":"points","legendgroup":"group_SNORD26_rnaseq","line":{"color":"rgba(0.7686274509803922,0.611764705882353,0.5803921568627451,1.)"},"marker":{"color":"rgba(0.7686274509803922,0.611764705882353,0.5803921568627451,1.)","size":8,"symbol":"circle"},"mode":"lines+markers","name":"SNORD26_rnaseq","showlegend":true,"type":"scatter","visible":true,"x":[0,1,2],"y":[-1.784870000001268,2.3110999999993025,3.7938299999958955]},{"hoveron":"points","legendgroup":"group_SNORD26_rnaseq","line":{"color":"rgba(0.7686274509803922,0.611764705882353,0.5803921568627451,1.)"},"marker":{"color":"rgba(0.7686274509803922,0.611764705882353,0.5803921568627451,1.)","size":12,"symbol":"x-open-dot"},"mode":"markers","name":"SNORD26_rnaseq","showlegend":false,"type":"scatter","visible":true,"x":[0,1,2],"y":[-1.784870000001268,2.3110999999993025,3.7938299999958955]},{"hoveron":"points","legendgroup":"group_SNORD27_rnaseq","line":{"color":"rgba(0.8901960784313725,0.4666666666666667,0.7607843137254902,1.)"},"marker":{"color":"rgba(0.8901960784313725,0.4666666666666667,0.7607843137254902,1.)","size":8,"symbol":"circle"},"mode":"lines+markers","name":"SNORD27_rnaseq","showlegend":true,"type":"scatter","visible":true,"x":[0,1,2],"y":[-1.784870000001268,2.3110999999993025,3.7938299999958955]},{"hoveron":"points","legendgroup":"group_SNORD27_rnaseq","line":{"color":"rgba(0.8901960784313725,0.4666666666666667,0.7607843137254902,1.)"},"marker":{"color":"rgba(0.8901960784313725,0.4666666666666667,0.7607843137254902,1.)","size":12,"symbol":"x-open-dot"},"mode":"markers","name":"SNORD27_rnaseq","showlegend":false,"type":"scatter","visible":true,"x":[0,1,2],"y":[-1.784870000001268,2.3110999999993025,3.7938299999958955]},{"hoveron":"points","legendgroup":"group_SNORD28_rnaseq","line":{"color":"rgba(0.9686274509803922,0.7137254901960784,0.8235294117647058,1.)"},"marker":{"color":"rgba(0.9686274509803922,0.7137254901960784,0.8235294117647058,1.)","size":8,"symbol":"circle"},"mode":"lines+markers","name":"SNORD28_rnaseq","showlegend":true,"type":"scatter","visible":true,"x":[0,1,2],"y":[-1.784870000001268,2.3110999999993025,3.7938299999958955]},{"hoveron":"points","legendgroup":"group_SNORD28_rnaseq","line":{"color":"rgba(0.9686274509803922,0.7137254901960784,0.8235294117647058,1.)"},"marker":{"color":"rgba(0.9686274509803922,0.7137254901960784,0.8235294117647058,1.)","size":12,"symbol":"x-open-dot"},"mode":"markers","name":"SNORD28_rnaseq","showlegend":false,"type":"scatter","visible":true,"x":[0,1,2],"y":[-1.784870000001268,2.3110999999993025,3.7938299999958955]},{"hoveron":"points","legendgroup":"group_SNORD30_rnaseq","line":{"color":"rgba(0.4980392156862745,0.4980392156862745,0.4980392156862745,1.)"},"marker":{"color":"rgba(0.4980392156862745,0.4980392156862745,0.4980392156862745,1.)","size":8,"symbol":"circle"},"mode":"lines+markers","name":"SNORD30_rnaseq","showlegend":true,"type":"scatter","visible":true,"x":[0,1,2],"y":[-1.784870000001268,2.3110999999993025,3.7938299999958955]},{"hoveron":"points","legendgroup":"group_SNORD30_rnaseq","line":{"color":"rgba(0.4980392156862745,0.4980392156862745,0.4980392156862745,1.)"},"marker":{"color":"rgba(0.4980392156862745,0.4980392156862745,0.4980392156862745,1.)","size":12,"symbol":"x-open-dot"},"mode":"markers","name":"SNORD30_rnaseq","showlegend":false,"type":"scatter","visible":true,"x":[0,1,2],"y":[-1.784870000001268,2.3110999999993025,3.7938299999958955]},{"hoveron":"points","legendgroup":"group_SNORD61_rnaseq","line":{"color":"rgba(0.7803921568627451,0.7803921568627451,0.7803921568627451,1.)"},"marker":{"color":"rgba(0.7803921568627451,0.7803921568627451,0.7803921568627451,1.)","size":8,"symbol":"circle"},"mode":"lines+markers","name":"SNORD61_rnaseq","showlegend":true,"type":"scatter","visible":true,"x":[0,1,2],"y":[1.412020000002594,-0.44532300000365466,-1.4826299999998749]},{"hoveron":"points","legendgroup":"group_SNORD61_rnaseq","line":{"color":"rgba(0.7803921568627451,0.7803921568627451,0.7803921568627451,1.)"},"marker":{"color":"rgba(0.7803921568627451,0.7803921568627451,0.7803921568627451,1.)","size":12,"symbol":"x-open-dot"},"mode":"markers","name":"SNORD61_rnaseq","showlegend":false,"type":"scatter","visible":true,"x":[0,1,2],"y":[1.412020000002594,-0.44532300000365466,-1.4826299999998749]},{"hoveron":"points","legendgroup":"group_TNFSF9_rnaseq","line":{"color":"rgba(0.7372549019607844,0.7411764705882353,0.13333333333333333,1.)"},"marker":{"color":"rgba(0.7372549019607844,0.7411764705882353,0.13333333333333333,1.)","size":8,"symbol":"circle"},"mode":"lines+markers","name":"TNFSF9_rnaseq","showlegend":true,"type":"scatter","visible":true,"x":[0,1,2],"y":[0.46080100000002855,1.1222299999968819,0.9594490000001178]},{"hoveron":"points","legendgroup":"group_TNFSF9_rnaseq","line":{"color":"rgba(0.7372549019607844,0.7411764705882353,0.13333333333333333,1.)"},"marker":{"color":"rgba(0.7372549019607844,0.7411764705882353,0.13333333333333333,1.)","size":12,"symbol":"x-open-dot"},"mode":"markers","name":"TNFSF9_rnaseq","showlegend":false,"type":"scatter","visible":true,"x":[0,1,2],"y":[0.46080100000002855,1.1222299999968819,0.9594490000001178]},{"hoveron":"points","legendgroup":"group_Y_RNA_rnaseq","line":{"color":"rgba(0.8588235294117647,0.8588235294117647,0.5529411764705883,1.)"},"marker":{"color":"rgba(0.8588235294117647,0.8588235294117647,0.5529411764705883,1.)","size":8,"symbol":"circle"},"mode":"lines+markers","name":"Y_RNA_rnaseq","showlegend":true,"type":"scatter","visible":true,"x":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],"y":[-0.03979469999437695,-0.3770459999960318,-0.00028178099850017667,0.0042912199985277455,-0.004495210001085427,-0.012645300006003541,-0.053455599996809176,0.07774510000444408,-0.08100109999516826,0.08674340000254914,0.09390679999416822,-0.10215500000555035,-0.10313400000273575,0.6367640000028714,-0.5576310000026997,1.877099999999678,-0.48350000000236854,-1.8092400000008713,-0.4274079999986531,0.36720399999814607,0.2633629999974941,0.2294609999958497,-0.33422200000247915,0.18013500000324925,0.15127500000598973,0.08085280000656944,-0.07281840000661273,-0.06820389999520148,-0.2149290000056653,0.4875230000051198,0.20009700000625938,-1.2051899999987983,-0.8491409999984292,0.0694757000059894,-0.08869969999858511,-0.17830699999615418,0.2719069999989679,0.8901329999985126,-4.070149999996867,-1.8393799999988054,0.3909170000040683,0.48948200000350894,0.7747470000029545,-1.1018199999997569,-0.5363539999961259,2.248810000001098,1.858429999998833]},{"hoveron":"points","legendgroup":"group_Y_RNA_rnaseq","line":{"color":"rgba(0.8588235294117647,0.8588235294117647,0.5529411764705883,1.)"},"marker":{"color":"rgba(0.8588235294117647,0.8588235294117647,0.5529411764705883,1.)","size":12,"symbol":"x-open-dot"},"mode":"markers","name":"Y_RNA_rnaseq","showlegend":false,"type":"scatter","visible":true,"x":[0,1,1,2,2,2,2,2,2,2,2],"y":[0.6367640000028714,-0.5576310000026997,-1.8092400000008713,0.8901329999985126,-4.070149999996867,-1.8393799999988054,0.3909170000040683,0.48948200000350894,0.7747470000029545,-1.1018199999997569,2.248810000001098]},{"hoveron":"points","legendgroup":"group_ZDBF2_rnaseq","line":{"color":"rgba(0.09019607843137255,0.7450980392156863,0.8117647058823529,1.)"},"marker":{"color":"rgba(0.09019607843137255,0.7450980392156863,0.8117647058823529,1.)","size":8,"symbol":"circle"},"mode":"lines+markers","name":"ZDBF2_rnaseq","showlegend":true,"type":"scatter","visible":true,"x":[0,1,2],"y":[-0.39830199999853977,0.6498229999960774,0.9563879999966954]},{"hoveron":"points","legendgroup":"group_ZDBF2_rnaseq","line":{"color":"rgba(0.09019607843137255,0.7450980392156863,0.8117647058823529,1.)"},"marker":{"color":"rgba(0.09019607843137255,0.7450980392156863,0.8117647058823529,1.)","size":12,"symbol":"x-open-dot"},"mode":"markers","name":"ZDBF2_rnaseq","showlegend":false,"type":"scatter","visible":true,"x":[0,1,2],"y":[-0.39830199999853977,0.6498229999960774,0.9563879999966954]}],                        {"hovermode":"closest","showlegend":true,"template":{"data":{"bar":[{"error_x":{"color":"#2a3f5f"},"error_y":{"color":"#2a3f5f"},"marker":{"line":{"color":"#E5ECF6","width":0.5},"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"bar"}],"barpolar":[{"marker":{"line":{"color":"#E5ECF6","width":0.5},"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"barpolar"}],"carpet":[{"aaxis":{"endlinecolor":"#2a3f5f","gridcolor":"white","linecolor":"white","minorgridcolor":"white","startlinecolor":"#2a3f5f"},"baxis":{"endlinecolor":"#2a3f5f","gridcolor":"white","linecolor":"white","minorgridcolor":"white","startlinecolor":"#2a3f5f"},"type":"carpet"}],"choropleth":[{"colorbar":{"outlinewidth":0,"ticks":""},"type":"choropleth"}],"contour":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"contour"}],"contourcarpet":[{"colorbar":{"outlinewidth":0,"ticks":""},"type":"contourcarpet"}],"heatmap":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"heatmap"}],"heatmapgl":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"heatmapgl"}],"histogram":[{"marker":{"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"histogram"}],"histogram2d":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"histogram2d"}],"histogram2dcontour":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"histogram2dcontour"}],"mesh3d":[{"colorbar":{"outlinewidth":0,"ticks":""},"type":"mesh3d"}],"parcoords":[{"line":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"parcoords"}],"pie":[{"automargin":true,"type":"pie"}],"scatter":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatter"}],"scatter3d":[{"line":{"colorbar":{"outlinewidth":0,"ticks":""}},"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatter3d"}],"scattercarpet":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scattercarpet"}],"scattergeo":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scattergeo"}],"scattergl":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scattergl"}],"scattermapbox":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scattermapbox"}],"scatterpolar":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatterpolar"}],"scatterpolargl":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatterpolargl"}],"scatterternary":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatterternary"}],"surface":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"surface"}],"table":[{"cells":{"fill":{"color":"#EBF0F8"},"line":{"color":"white"}},"header":{"fill":{"color":"#C8D4E3"},"line":{"color":"white"}},"type":"table"}]},"layout":{"annotationdefaults":{"arrowcolor":"#2a3f5f","arrowhead":0,"arrowwidth":1},"autotypenumbers":"strict","coloraxis":{"colorbar":{"outlinewidth":0,"ticks":""}},"colorscale":{"diverging":[[0,"#8e0152"],[0.1,"#c51b7d"],[0.2,"#de77ae"],[0.3,"#f1b6da"],[0.4,"#fde0ef"],[0.5,"#f7f7f7"],[0.6,"#e6f5d0"],[0.7,"#b8e186"],[0.8,"#7fbc41"],[0.9,"#4d9221"],[1,"#276419"]],"sequential":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"sequentialminus":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]},"colorway":["#636efa","#EF553B","#00cc96","#ab63fa","#FFA15A","#19d3f3","#FF6692","#B6E880","#FF97FF","#FECB52"],"font":{"color":"#2a3f5f"},"geo":{"bgcolor":"white","lakecolor":"white","landcolor":"#E5ECF6","showlakes":true,"showland":true,"subunitcolor":"white"},"hoverlabel":{"align":"left"},"hovermode":"closest","mapbox":{"style":"light"},"paper_bgcolor":"white","plot_bgcolor":"#E5ECF6","polar":{"angularaxis":{"gridcolor":"white","linecolor":"white","ticks":""},"bgcolor":"#E5ECF6","radialaxis":{"gridcolor":"white","linecolor":"white","ticks":""}},"scene":{"xaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","gridwidth":2,"linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white"},"yaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","gridwidth":2,"linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white"},"zaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","gridwidth":2,"linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white"}},"shapedefaults":{"line":{"color":"#2a3f5f"}},"ternary":{"aaxis":{"gridcolor":"white","linecolor":"white","ticks":""},"baxis":{"gridcolor":"white","linecolor":"white","ticks":""},"bgcolor":"#E5ECF6","caxis":{"gridcolor":"white","linecolor":"white","ticks":""}},"title":{"x":0.05},"xaxis":{"automargin":true,"gridcolor":"white","linecolor":"white","ticks":"","title":{"standoff":15},"zerolinecolor":"white","zerolinewidth":2},"yaxis":{"automargin":true,"gridcolor":"white","linecolor":"white","ticks":"","title":{"standoff":15},"zerolinecolor":"white","zerolinewidth":2}}},"updatemenus":[{"buttons":[{"args":["visible",[true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true]],"label":"All","method":"restyle"},{"args":["visible",[true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false]],"label":"AC009487.6","method":"restyle"},{"args":["visible",[false,false,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false]],"label":"C11orf68","method":"restyle"},{"args":["visible",[false,false,false,false,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false]],"label":"CITED2","method":"restyle"},{"args":["visible",[false,false,false,false,false,false,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false]],"label":"FOXC1","method":"restyle"},{"args":["visible",[false,false,false,false,false,false,false,false,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false]],"label":"PLCH1","method":"restyle"},{"args":["visible",[false,false,false,false,false,false,false,false,false,false,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false]],"label":"RBMX","method":"restyle"},{"args":["visible",[false,false,false,false,false,false,false,false,false,false,false,false,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false]],"label":"RP11-639F1.1","method":"restyle"},{"args":["visible",[false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false]],"label":"SESN1","method":"restyle"},{"args":["visible",[false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false]],"label":"SNHG1","method":"restyle"},{"args":["visible",[false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false]],"label":"SNORD22","method":"restyle"},{"args":["visible",[false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false]],"label":"SNORD25","method":"restyle"},{"args":["visible",[false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false]],"label":"SNORD26","method":"restyle"},{"args":["visible",[false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,false,false,false,false,false,false,false,false,false,false,false,false]],"label":"SNORD27","method":"restyle"},{"args":["visible",[false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,false,false,false,false,false,false,false,false,false,false]],"label":"SNORD28","method":"restyle"},{"args":["visible",[false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,false,false,false,false,false,false,false,false]],"label":"SNORD30","method":"restyle"},{"args":["visible",[false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,false,false,false,false,false,false]],"label":"SNORD61","method":"restyle"},{"args":["visible",[false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,false,false,false,false]],"label":"TNFSF9","method":"restyle"},{"args":["visible",[false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,false,false]],"label":"Y_RNA","method":"restyle"},{"args":["visible",[false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,true]],"label":"ZDBF2","method":"restyle"}],"x":-0.05,"y":1,"yanchor":"top"}],"xaxis":{"range":[0,2],"showticklabels":true,"tickmode":"array","ticktext":["01hr","06hr","24hr"],"tickvals":[0,1,2],"title":{"text":"Sample index"}},"yaxis":{"title":{"text":"log2fc"}}},                        {"responsive": true}                    ).then(function(){
    
    var gd = document.getElementById('c19a0e71-9da0-491b-96c7-15eb7a6b27d3');
    var x = new MutationObserver(function (mutations, observer) {{
            var display = window.getComputedStyle(gd).display;
            if (!display || display === 'none') {{
                console.log([gd, 'removed!']);
                Plotly.purge(gd);
                observer.disconnect();
            }}
    }});
    
    // Listen for the removal of the full notebook cells
    var notebookContainer = gd.closest('#notebook-container');
    if (notebookContainer) {{
        x.observe(notebookContainer, {childList: true});
    }}
    
    // Listen for the clearing of the current output cell
    var outputEl = gd.closest('.output');
    if (outputEl) {{
        x.observe(outputEl, {childList: true});
    }}
    
                            })                };                });            </script>        </div>


Query 3:
^^^^^^^^

-  Heatmap and time series plot of proteins that are consistently down
   regulated at 3 time points.

.. code:: ipython3

    exp_data.proteins.down.require_n_sig(n_sig=3, index='label').plot_species(plot_type='matplotlib');
    exp_data.proteins.down.require_n_sig(n_sig=3, index='label').heatmap(index='label', cluster_row=True, linewidths=0.01);



.. image:: Tutorial_files/Tutorial_64_0.png



.. image:: Tutorial_files/Tutorial_64_1.png


Query 4:
^^^^^^^^

::

   Clustered heatmap of label-free data

.. code:: ipython3

    exp_data.label_free.heatmap(
        linewidths=0.01,
        index='label',
        cluster_row=True, 
        cluster_col=True, 
        min_sig=3, 
        figsize=(12,18)
    );



.. image:: Tutorial_files/Tutorial_66_0.png


Extending to other plots
~~~~~~~~~~~~~~~~~~~~~~~~

Since our exp_data is built off a pandas.DataFrame, we can use other
packages that take that data format. Seaborn is one such tool that
provides some very nice plots.

.. code:: ipython3

    label_free = exp_data.label_free.copy()
    label_free.log2_normalize_df(column='fold_change', inplace=True)
    
    g = sns.PairGrid(label_free,
                     x_vars=('sample_id'),
                     y_vars=('fold_change', 'p_value'),
                     aspect=3.25, height=3.5)
    g.map(
        sns.violinplot, 
        palette="pastel", 
        order=label_free.sample_ids
    );



.. image:: Tutorial_files/Tutorial_68_0.png


Venn diagram comparisons between measurements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    from magine.plotting.venn_diagram_maker import create_venn2, create_venn3

.. code:: ipython3

    lf = exp_data.label_free.sig

.. code:: ipython3

    lf = exp_data.label_free.sig.id_list
    silac = exp_data.silac.sig.id_list
    phsilac = exp_data.ph_silac.sig.id_list
    hilic = exp_data.HILIC.sig.id_list
    rplc = exp_data.C18.sig.id_list
    
    create_venn2(hilic, rplc, 'HILIC', 'RPLC');



.. image:: Tutorial_files/Tutorial_72_0.png


.. code:: ipython3

    create_venn3(lf, silac, phsilac, 'LF', 'SILAC', 'ph-SILAC');



.. image:: Tutorial_files/Tutorial_73_0.png


Networks
--------

Create data driven network
~~~~~~~~~~~~~~~~~~~~~~~~~~

MAGINE generates networks using prior knowledge obtained from network
databases. It starts with ``seed`` species, which are biologically
interesting. This can be from manual curation, based on filtering of
data, or a variety of other means. We locate the seed species in
multiple databases and finds interconnecting edges among them. The goal
of this process was to obtain all the ``known`` biological regulation
among the species. We currently utilize KEGG, Reactome, HMDB, TTRUST,
and BioGrid for node and edge sources.

.. code:: ipython3

    # some imports needed
    from magine.networks.network_generator import build_network
    import magine.networks.utils as utils
    import networkx as nx
    import os


.. parsed-literal::

    2021-09-03 13:59:25.008 - magine - INFO - Logging started on MAGINE version 0.1.5
    2021-09-03 13:59:25.009 - magine - INFO - Log entry time offset from UTC: -7.00 hours
    [33mWARNING [bioservices:UniChem:119]: [0m [34mThe URL (http://www.ebi.ac.uk/unichem/rest) provided cannot be reached.[0m
    [33mWARNING [bioservices:UniChem:119]: [0m [34mUniChem has added new source. Please update the source_ids attribute in bioservices[0m
    

This is done using the ``build_network`` function. Now we will create
the network. We pass the seed (interesting species) and background list
(all things measured that we want to make sure are not deleted from the
network) to the generator as well as flags turning on all of the network
databases.

.. code:: ipython3

    measured = exp_data.species.id_list
    sig_measured = exp_data.species.sig.id_list
    print(len(measured))
    print(len(sig_measured))


.. parsed-literal::

    23725
    15777
    

.. code:: ipython3

    if not os.path.exists('Data/cisplatin_network.p'):
        network = build_network(
            
            # genes seed species
            seed_species=sig_measured, 
            
            # all data measured, used to allow interconnecting nodes that are not in seeds.
            all_measured_list=measured,  
            
            use_biogrid=True,  # expand with biogrid
            use_hmdb=True,  # expand with hmdb
            use_reactome=True,  # expand with reactome
            use_signor=True,  # expand with signor
            trim_source_sink=True,  # remove all source and sink nodes not measured
            save_name='Data/cisplatin_network'
        )
        # add attibutes to graph nodes (measured, measured at which time points, 
        # significantly changed at which time point)
        utils.add_data_to_graph(network, exp_data)
        print("Saving network")
        # write to GML for cytoscape or other program
        nx.write_gml(
            network,
            os.path.join('Data', 'cisplatin_network_w_attributes.gml')
        )
    
        # write to gpickle for fast loading in python
        nx.write_gpickle(
            network,
            os.path.join('Data', 'cisplatin_based_network.p'),
        )
    else:
        # Load the network, note that it is returned above but for time limits,
        # we will just load the generated one.
        network = nx.read_gpickle('Data/cisplatin_based_network.p')
    
    

.. code:: ipython3

    print(network.number_of_nodes())
    print(network.number_of_edges())


.. parsed-literal::

    13308
    181300
    

As you might iMAGINE, the larger number of input nodes and source
databases, the larger the resulting network. 13308 nodes and 181300
edges are too much to manually explore. Thus, we are going to use the
``Subgraph`` Class to being to query the network. We developed multiple
tools to subgraph and explore the network, as well as multiple
visualizations.

Explore subgraphs of network
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    from magine.networks.subgraphs import Subgraph
    from magine.networks.visualization import draw_igraph, draw_graphviz, draw_mpl, draw_cyjs
    net_sub = Subgraph(network)

.. code:: ipython3

    bax_neighbors = net_sub.neighbors(
        'BAX', # node of interest
        upstream=True, # include upstream nodes
        downstream=False,  # include downstream nodes
        include_only=exp_data.species.sig.id_list # limit nodes to only significant changed species
    )

There are multiple ways to visualize the network. draw_igraph, draw_mpl,
draw_graphviz, draw_cyjs. The draw_cyjs provides interactive networks
within the browser. It shouldn’t be used for extremely large networks.
There is also some issues with exporting these views when converting the
notebooks to other formats (such as web documentation).

.. code:: ipython3

    draw_igraph(bax_neighbors, bbox=[400, 400], node_size=25, inline=True)




.. image:: Tutorial_files/Tutorial_86_0.svg



.. code:: ipython3

    draw_graphviz(bax_neighbors, 'fdp')



.. image:: Tutorial_files/Tutorial_87_0.png
   :width: 800px


.. code:: ipython3

    draw_cyjs(bax_neighbors)



.. raw:: html

    <!DOCTYPE html>
    <html>
    <head>
        <meta charset=utf-8/>
        <style type="text/css">
            body {
                font: 14px helvetica neue, helvetica, arial, sans-serif;
            }
    
            #cycaa67223-bbd2-49e2-a62a-c4a40a5c1e9c{
                height: 700px;
                width: 90%;
                border: 5px solid black;
                box-sizing: border-box;
                position: relative;
                top: 5;
                margin-bottom: -700px;
                background: white;
            }
    
        </style>
    
        <script>
    
    
            requirejs.config({
    
                paths: {
                    'popper': 'https://unpkg.com/popper.js@1.14.1/dist/umd/popper',
                    'tippy': 'https://cdnjs.cloudflare.com/ajax/libs/tippy.js/2.3.0/tippy.min',
                    'cytoscape': 'https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.2.10/cytoscape',
                    'cytoscape-popper': 'https://cdn.rawgit.com/cytoscape/cytoscape.js-popper/3ad50859/cytoscape-popper',
                    'jquery': 'https://cdnjs.cloudflare.com/ajax/libs/jquery/2.2.4/jquery.min',
                    'qtip2': 'https://cdnjs.cloudflare.com/ajax/libs/qtip2/2.2.0/basic/jquery.qtip.min',
                    'dagre': 'https://cdn.rawgit.com/cpettitt/dagre/v0.7.4/dist/dagre.min',
                    'cytoscape-dagre': 'https://cdn.rawgit.com/cytoscape/cytoscape.js-dagre/1.5.0/cytoscape-dagre',
                    'cytoscape-cose-bilkent': 'https://cdn.rawgit.com/cytoscape/cytoscape.js-cose-bilkent/1.6.1/cytoscape-cose-bilkent'
                },
                shim: {
                    'cytoscape-popper': {
                        deps: ['cytoscape', 'popper']
                    },
                    'cytoscape-dagre': {
                        deps: ['cytoscape', 'dagre']
                    }
                },
                map: {
                    '*': {
                        'popper.js': 'popper',
                        'webcola': 'cola'
                    }
                }
    
            });
    
    
            require(['cytoscape', 'cytoscape-popper', 'popper', 'tippy', 'jquery',
                    'cytoscape-cose-bilkent', 'cytoscape-dagre', 'dagre'],
                function (cytoscape, cypopper, popper, tippy, jquery, regCose,
                          cydag, dagre) {
                    console.log('Loading Cytoscape.js Module...');
                    window['popper'] = popper;
                    window['tippy'] = tippy;
                    window['cytoscape'] = cytoscape;
                    cypopper(cytoscape);
                    regCose(cytoscape);
                    cydag(cytoscape, dagre);
    
                    function makeTippy(target, text) {
                        return tippy(target.popperRef(),
                            {
                                html: add_tip(text),
                                trigger: 'manual',
                                arrow: true,
                                placement: 'top',
                                hideOnClick: false,
                                interactive: true,
                                multiple: true,
                                sticky: true
                            }).tooltips[0];
                    }
    
                    function add_tip(text) {
                        let div = document.createElement('div');
                        div.innerHTML = "<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=" + text + "' target='_blank'>Gene Card</a>";
                        return div;
                    }
    
                    let cy = window.cy = cytoscape({
                        container: $('#cycaa67223-bbd2-49e2-a62a-c4a40a5c1e9c'),
                        elements: {
                            nodes: [{"data": {"speciesType": "gene", "keggName": "581", "databaseSource": "KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "BAX", "name": "BAX"}}, {"data": {"speciesType": "gene", "keggName": "207", "databaseSource": "BioGrid|HMDB|KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "AKT1", "name": "AKT1"}}, {"data": {"speciesType": "gene", "keggName": "5566", "databaseSource": "BioGrid|HMDB|KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "PRKACA", "name": "PRKACA"}}, {"data": {"speciesType": "gene", "keggName": "5567", "databaseSource": "BioGrid|HMDB|KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "PRKACB", "name": "PRKACB"}}, {"data": {"speciesType": "gene", "keggName": "5590", "databaseSource": "BioGrid|HMDB|KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "PRKCZ", "name": "PRKCZ"}}, {"data": {"speciesType": "gene", "keggName": "2932", "databaseSource": "BioGrid|HMDB|KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "GSK3B", "name": "GSK3B"}}, {"data": {"speciesType": "gene", "keggName": "23411", "databaseSource": "BioGrid|HMDB|KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "SIRT1", "name": "SIRT1"}}, {"data": {"speciesType": "gene", "keggName": "2033", "databaseSource": "BioGrid|HMDB|KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "EP300", "name": "EP300"}}, {"data": {"speciesType": "gene", "keggName": "7157", "databaseSource": "BioGrid|KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "TP53", "name": "TP53"}}, {"data": {"speciesType": "gene", "keggName": "637", "databaseSource": "BioGrid|KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "BID", "name": "BID"}}, {"data": {"speciesType": "gene", "keggName": "598", "databaseSource": "BioGrid|KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "BCL2L1", "name": "BCL2L1"}}, {"data": {"speciesType": "gene", "keggName": "27113", "databaseSource": "KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "BBC3", "name": "BBC3"}}, {"data": {"speciesType": "gene", "keggName": "29108", "databaseSource": "BioGrid|KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "PYCARD", "name": "PYCARD"}}, {"data": {"speciesType": "gene", "keggName": "4170", "databaseSource": "BioGrid|KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "MCL1", "name": "MCL1"}}, {"data": {"speciesType": "gene", "keggName": "10018", "databaseSource": "BioGrid|KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "BCL2L11", "name": "BCL2L11"}}, {"data": {"speciesType": "gene", "keggName": "2934", "databaseSource": "BioGrid|HMDB|KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "GSN", "name": "GSN"}}, {"data": {"speciesType": "gene", "keggName": "51100", "databaseSource": "KEGG|SIGNOR", "color": "white", "id": "SH3GLB1", "name": "SH3GLB1"}}, {"data": {"speciesType": "gene", "keggName": "599", "databaseSource": "KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "BCL2L2", "name": "BCL2L2"}}, {"data": {"speciesType": "gene", "databaseSource": "BioGrid|ReactomeFI|SIGNOR", "color": "white", "id": "YY1", "name": "YY1"}}, {"data": {"speciesType": "gene", "databaseSource": "ReactomeFI", "color": "white", "id": "BOK", "name": "BOK"}}],
                            edges: [{"data": {"databaseSource": "SIGNOR", "interactionType": "activate|phosphorylate", "width": 10.0, "source": "AKT1", "target": "BAX"}}, {"data": {"databaseSource": "KEGG", "interactionType": "inhibit|phosphorylate", "width": 10.0, "source": "PRKACA", "target": "BAX"}}, {"data": {"databaseSource": "KEGG", "interactionType": "inhibit|phosphorylate", "width": 10.0, "source": "PRKACB", "target": "BAX"}}, {"data": {"databaseSource": "SIGNOR", "interactionType": "inhibit|phosphorylate", "width": 10.0, "source": "PRKCZ", "target": "BAX"}}, {"data": {"databaseSource": "SIGNOR", "interactionType": "activate|phosphorylate", "width": 10.0, "source": "GSK3B", "target": "BAX"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "binding|inhibit", "width": 10.0, "source": "SIRT1", "target": "BAX"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "EP300", "target": "BAX"}}, {"data": {"databaseSource": "KEGG|ReactomeFI|SIGNOR", "interactionType": "activate|binding|expression|indirect", "width": 10.0, "source": "TP53", "target": "BAX"}}, {"data": {"databaseSource": "KEGG|SIGNOR", "interactionType": "activate|binding", "width": 10.0, "source": "BID", "target": "BAX"}}, {"data": {"databaseSource": "KEGG|ReactomeFI|SIGNOR", "interactionType": "binding|inhibit", "width": 10.0, "source": "BCL2L1", "target": "BAX"}}, {"data": {"databaseSource": "SIGNOR", "interactionType": "activate|binding", "width": 10.0, "source": "BBC3", "target": "BAX"}}, {"data": {"databaseSource": "SIGNOR", "interactionType": "activate|relocalization", "width": 10.0, "source": "PYCARD", "target": "BAX"}}, {"data": {"databaseSource": "SIGNOR", "interactionType": "binding|inhibit", "width": 10.0, "source": "MCL1", "target": "BAX"}}, {"data": {"databaseSource": "KEGG|ReactomeFI|SIGNOR", "interactionType": "activate|binding", "width": 10.0, "source": "BCL2L11", "target": "BAX"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "inhibit", "width": 10.0, "source": "GSN", "target": "BAX"}}, {"data": {"databaseSource": "SIGNOR", "interactionType": "activate|binding", "width": 10.0, "source": "SH3GLB1", "target": "BAX"}}, {"data": {"databaseSource": "SIGNOR", "interactionType": "binding|inhibit", "width": 10.0, "source": "BCL2L2", "target": "BAX"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "YY1", "target": "BAX"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "BOK", "target": "BAX"}}]
                        },
                        boxSelectionEnabled: true,
                        wheelSensitivity: .25,
    
                    });
    
                    // apply layout
                    cy.layout({"name": "cose-bilkent", "spacingFactor": 1.5, "animate": false}).run();
    
                    // add style
                    cy.style().fromJson([{"selector": "node", "css": {"text-opacity": 1.0, "background-opacity": 1.0, "font-weight": "normal", "color": "rgb(0,153,234)", "border-width": 3.0, "border-color": "rgb(51,51,51)", "font-size": 9, "height": 50, "width": 50, "label": "data(name)", "text-wrap": "wrap", "text-max-width": 50, "shape": "ellipse", "background-color": "data(color)", "text-halign": "center", "font-family": "SansSerif", "text-valign": "center", "border-opacity": 1.0}}, {"selector": "node[speciesType = 'compound']", "css": {"text-opacity": 1.0, "background-opacity": 1.0, "font-weight": "normal", "color": "rgb(0,153,234)", "border-width": 3.0, "border-color": "rgb(51,51,51)", "font-size": 9, "height": 35, "width": 35, "label": "data(chemName)", "text-wrap": "wrap", "text-max-width": 95, "shape": "ellipse", "background-color": "rgb(255,255,255)", "text-halign": "center", "font-family": "SansSerif", "text-valign": "center", "border-opacity": 1.0}}, {"selector": "$node > node", "css": {"font-size": 20, "shape": "ellipse", "padding-right": "10px", "padding-bottom": "10px", "padding-top": "10px", "text-valign": "top", "text-halign": "center", "background-color": "#bbb", "padding-left": "10px"}}, {"selector": "node:selected", "css": {"background-color": "rgb(255,0,102)"}}, {"selector": "edge", "css": {"opacity": 1.0, "text-opacity": 1.0, "font-size": 12, "font-weight": "normal", "width": "2", "color": "rgb(0,0,0)", "line-color": "rgb(51,51,51)", "source-arrow-color": "rgb(0,0,0)", "target-arrow-color": "rgb(51,51,51)", "curve-style": "bezier", "line-style": "solid", "target-arrow-shape": "triangle", "source-arrow-shape": "none", "font-family": "SansSerif"}}, {"selector": "edge[weight>1]", "css": {"opacity": 1.0, "text-opacity": 1.0, "font-size": 12, "font-weight": "normal", "width": "data(width)", "color": "rgb(0,0,0)", "line-color": "rgb(51,51,51)", "source-arrow-color": "rgb(0,0,0)", "target-arrow-color": "rgb(51,51,51)", "curve-style": "bezier", "line-style": "solid", "target-arrow-shape": "triangle", "source-arrow-shape": "none", "font-family": "SansSerif"}}, {"selector": "edge[interactionType *= 'inhibit'],edge[interactionType *= 'deactivat']", "css": {"text-opacity": 1.0, "font-size": 12, "font-weight": "normal", "font-family": "SansSerif", "curve-style": "bezier", "source-arrow-shape": "none", "target-arrow-shape": "tee", "color": "rgb(0,0,0)", "line-color": "rgb(51,51,51)", "source-arrow-color": "rgb(0,0,0)", "target-arrow-color": "rgb(51,51,51)", "line-style": "solid"}}, {"selector": "edge:selected", "css": {"line-color": "rgb(255,0,0)", "label": "data(interactionType)"}}, {"selector": "edge[weight>1]:selected", "css": {"background-color": "yellow", "line-color": "red", "label": "data(label)", "target-arrow-color": "black", "source-arrow-color": "black", "width": "data(width)"}}]).update();
    
                    // Add tippy for each node
                    cy.nodes().forEach(function (n) {
                        n.data()['tip'] = makeTippy(n, n.data('id'));
                    });
    
                    // hide tippy text on click
                    cy.on('tap', 'node', function (evt) {
                        var ele = evt.target;
                        if (ele.data()['tip']['state']['visible']) {
                            ele.data()['tip'].hide();
                        } else {
                            ele.data()['tip'].show();
                        }
                    });
    
                    // put the png data in an img tag
                    let downloadButton = document.createElement("BUTTON");
                    downloadButton.id = 'dbutton';
                    downloadButton.innerHTML = '<i class="fa fa-download" aria-hidden="true"></i>';
                    downloadButton.addEventListener('click', function () {
                        let element = document.createElement('a');
                        element.setAttribute('href', cy.png({scale: 3}));
                        element.setAttribute('download', 'graph.png');
                        element.style.display = 'none';
                        document.body.appendChild(element);
                        element.click();
                        document.body.removeChild(element);
                    });
                    let p = document.getElementById('cycaa67223-bbd2-49e2-a62a-c4a40a5c1e9c');
                    p.parentElement.append(downloadButton);
    
                });
    
    
        </script>
    </head>
    
    <body>
    <div id="cycaa67223-bbd2-49e2-a62a-c4a40a5c1e9c"></div>
    <!-- When only #uuid div is placed on this page,
    the height of output-box on ipynb will be 0px.
    One line below will prevent that. -->
    <div id="dummy"
         style="width:100px;height:700px"></div>
    </body>
    
    </html>


This network can be exanded by a single or list of nodes passed.

.. code:: ipython3

    expand = net_sub.expand_neighbors(bax_neighbors, nodes='BID', downstream=True)

.. code:: ipython3

    draw_igraph(expand, 
                bbox=[500, 500], 
                node_font_size=25,
                font_size=4,
                node_size=25, 
                inline=True, 
                layout='graphopt')




.. image:: Tutorial_files/Tutorial_91_0.svg



.. code:: ipython3

    draw_graphviz(expand, 'sfdp', width=700)



.. image:: Tutorial_files/Tutorial_92_0.png
   :width: 700px


Finding paths between two proteins
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    BAX_to_AKT = net_sub.paths_between_pair('NOTCH1', 'MYC', bidirectional=True)

.. code:: ipython3

    draw_graphviz(BAX_to_AKT)



.. image:: Tutorial_files/Tutorial_95_0.png
   :width: 800px


.. code:: ipython3

    draw_cyjs(BAX_to_AKT)



.. raw:: html

    <!DOCTYPE html>
    <html>
    <head>
        <meta charset=utf-8/>
        <style type="text/css">
            body {
                font: 14px helvetica neue, helvetica, arial, sans-serif;
            }
    
            #cyabc00e71-4be2-477c-abe9-3987594f8f4d{
                height: 700px;
                width: 90%;
                border: 5px solid black;
                box-sizing: border-box;
                position: relative;
                top: 5;
                margin-bottom: -700px;
                background: white;
            }
    
        </style>
    
        <script>
    
    
            requirejs.config({
    
                paths: {
                    'popper': 'https://unpkg.com/popper.js@1.14.1/dist/umd/popper',
                    'tippy': 'https://cdnjs.cloudflare.com/ajax/libs/tippy.js/2.3.0/tippy.min',
                    'cytoscape': 'https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.2.10/cytoscape',
                    'cytoscape-popper': 'https://cdn.rawgit.com/cytoscape/cytoscape.js-popper/3ad50859/cytoscape-popper',
                    'jquery': 'https://cdnjs.cloudflare.com/ajax/libs/jquery/2.2.4/jquery.min',
                    'qtip2': 'https://cdnjs.cloudflare.com/ajax/libs/qtip2/2.2.0/basic/jquery.qtip.min',
                    'dagre': 'https://cdn.rawgit.com/cpettitt/dagre/v0.7.4/dist/dagre.min',
                    'cytoscape-dagre': 'https://cdn.rawgit.com/cytoscape/cytoscape.js-dagre/1.5.0/cytoscape-dagre',
                    'cytoscape-cose-bilkent': 'https://cdn.rawgit.com/cytoscape/cytoscape.js-cose-bilkent/1.6.1/cytoscape-cose-bilkent'
                },
                shim: {
                    'cytoscape-popper': {
                        deps: ['cytoscape', 'popper']
                    },
                    'cytoscape-dagre': {
                        deps: ['cytoscape', 'dagre']
                    }
                },
                map: {
                    '*': {
                        'popper.js': 'popper',
                        'webcola': 'cola'
                    }
                }
    
            });
    
    
            require(['cytoscape', 'cytoscape-popper', 'popper', 'tippy', 'jquery',
                    'cytoscape-cose-bilkent', 'cytoscape-dagre', 'dagre'],
                function (cytoscape, cypopper, popper, tippy, jquery, regCose,
                          cydag, dagre) {
                    console.log('Loading Cytoscape.js Module...');
                    window['popper'] = popper;
                    window['tippy'] = tippy;
                    window['cytoscape'] = cytoscape;
                    cypopper(cytoscape);
                    regCose(cytoscape);
                    cydag(cytoscape, dagre);
    
                    function makeTippy(target, text) {
                        return tippy(target.popperRef(),
                            {
                                html: add_tip(text),
                                trigger: 'manual',
                                arrow: true,
                                placement: 'top',
                                hideOnClick: false,
                                interactive: true,
                                multiple: true,
                                sticky: true
                            }).tooltips[0];
                    }
    
                    function add_tip(text) {
                        let div = document.createElement('div');
                        div.innerHTML = "<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=" + text + "' target='_blank'>Gene Card</a>";
                        return div;
                    }
    
                    let cy = window.cy = cytoscape({
                        container: $('#cyabc00e71-4be2-477c-abe9-3987594f8f4d'),
                        elements: {
                            nodes: [{"data": {"speciesType": "gene", "keggName": "4609", "databaseSource": "BioGrid|KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "MYC", "name": "MYC"}}, {"data": {"speciesType": "gene", "keggName": "26508", "databaseSource": "KEGG|ReactomeFI", "color": "white", "id": "HEYL", "name": "HEYL"}}, {"data": {"speciesType": "gene", "keggName": "3516", "databaseSource": "KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "RBPJ", "name": "RBPJ"}}, {"data": {"speciesType": "gene", "keggName": "3725", "databaseSource": "BioGrid|HMDB|KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "JUN", "name": "JUN"}}, {"data": {"speciesType": "gene", "keggName": "1499", "databaseSource": "BioGrid|HMDB|KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "CTNNB1", "name": "CTNNB1"}}, {"data": {"speciesType": "gene", "keggName": "8312", "databaseSource": "BioGrid|KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "AXIN1", "name": "AXIN1"}}, {"data": {"speciesType": "gene", "keggName": "4851", "databaseSource": "BioGrid|HMDB|KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "NOTCH1", "name": "NOTCH1"}}, {"data": {"speciesType": "gene", "keggName": "11317", "databaseSource": "KEGG|ReactomeFI", "color": "white", "id": "RBPJL", "name": "RBPJL"}}, {"data": {"speciesType": "gene", "keggName": "23493", "databaseSource": "KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "HEY2", "name": "HEY2"}}, {"data": {"speciesType": "gene", "keggName": "23462", "databaseSource": "BioGrid|KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "HEY1", "name": "HEY1"}}, {"data": {"speciesType": "gene", "keggName": "2033", "databaseSource": "BioGrid|HMDB|KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "EP300", "name": "EP300"}}, {"data": {"speciesType": "gene", "keggName": "3280", "databaseSource": "BioGrid|KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "HES1", "name": "HES1"}}, {"data": {"speciesType": "gene", "keggName": "388585", "databaseSource": "KEGG|ReactomeFI", "color": "white", "id": "HES5", "name": "HES5"}}, {"data": {"speciesType": "gene", "keggName": "4088", "databaseSource": "BioGrid|KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "SMAD3", "name": "SMAD3"}}, {"data": {"speciesType": "gene", "keggName": "10524", "databaseSource": "BioGrid|HMDB|KEGG|ReactomeFI|SIGNOR", "color": "white", "id": "KAT5", "name": "KAT5"}}],
                            edges: [{"data": {"databaseSource": "ReactomeFI", "interactionType": "binding|inhibit", "width": 10.0, "source": "MYC", "target": "KAT5"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "binding|inhibit", "width": 10.0, "source": "MYC", "target": "EP300"}}, {"data": {"databaseSource": "SIGNOR", "interactionType": "binding|inhibit", "width": 10.0, "source": "MYC", "target": "SMAD3"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "HEYL", "target": "MYC"}}, {"data": {"databaseSource": "KEGG", "interactionType": "expression", "width": 10.0, "source": "RBPJ", "target": "MYC"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "JUN", "target": "MYC"}}, {"data": {"databaseSource": "KEGG", "interactionType": "expression", "width": 10.0, "source": "CTNNB1", "target": "MYC"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "AXIN1", "target": "MYC"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "NOTCH1", "target": "HEYL"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "activate|expression", "width": 10.0, "source": "NOTCH1", "target": "HES1"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "NOTCH1", "target": "HEY2"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "NOTCH1", "target": "HES5"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "NOTCH1", "target": "HEY1"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "NOTCH1", "target": "RBPJL"}}, {"data": {"databaseSource": "KEGG|ReactomeFI|SIGNOR", "interactionType": "activate|binding|catalyze|inhibit", "width": 10.0, "source": "NOTCH1", "target": "RBPJ"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "NOTCH1", "target": "CTNNB1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "NOTCH1", "target": "AXIN1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "inhibit", "width": 10.0, "source": "NOTCH1", "target": "JUN"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "RBPJL", "target": "MYC"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "HEY2", "target": "MYC"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "HEY1", "target": "MYC"}}, {"data": {"pubmedId": "22100894", "databaseSource": "BIOGRID", "interactionType": "acetylate", "width": 10.0, "source": "EP300", "target": "NOTCH1"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "HES1", "target": "MYC"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "HES5", "target": "MYC"}}, {"data": {"databaseSource": "SIGNOR", "interactionType": "activate|binding", "width": 10.0, "source": "SMAD3", "target": "NOTCH1"}}, {"data": {"databaseSource": "SIGNOR", "interactionType": "acetylate|inhibit", "width": 10.0, "source": "KAT5", "target": "NOTCH1"}}]
                        },
                        boxSelectionEnabled: true,
                        wheelSensitivity: .25,
    
                    });
    
                    // apply layout
                    cy.layout({"name": "cose-bilkent", "spacingFactor": 1.5, "animate": false}).run();
    
                    // add style
                    cy.style().fromJson([{"selector": "node", "css": {"text-opacity": 1.0, "background-opacity": 1.0, "font-weight": "normal", "color": "rgb(0,153,234)", "border-width": 3.0, "border-color": "rgb(51,51,51)", "font-size": 9, "height": 50, "width": 50, "label": "data(name)", "text-wrap": "wrap", "text-max-width": 50, "shape": "ellipse", "background-color": "data(color)", "text-halign": "center", "font-family": "SansSerif", "text-valign": "center", "border-opacity": 1.0}}, {"selector": "node[speciesType = 'compound']", "css": {"text-opacity": 1.0, "background-opacity": 1.0, "font-weight": "normal", "color": "rgb(0,153,234)", "border-width": 3.0, "border-color": "rgb(51,51,51)", "font-size": 9, "height": 35, "width": 35, "label": "data(chemName)", "text-wrap": "wrap", "text-max-width": 95, "shape": "ellipse", "background-color": "rgb(255,255,255)", "text-halign": "center", "font-family": "SansSerif", "text-valign": "center", "border-opacity": 1.0}}, {"selector": "$node > node", "css": {"font-size": 20, "shape": "ellipse", "padding-right": "10px", "padding-bottom": "10px", "padding-top": "10px", "text-valign": "top", "text-halign": "center", "background-color": "#bbb", "padding-left": "10px"}}, {"selector": "node:selected", "css": {"background-color": "rgb(255,0,102)"}}, {"selector": "edge", "css": {"opacity": 1.0, "text-opacity": 1.0, "font-size": 12, "font-weight": "normal", "width": "2", "color": "rgb(0,0,0)", "line-color": "rgb(51,51,51)", "source-arrow-color": "rgb(0,0,0)", "target-arrow-color": "rgb(51,51,51)", "curve-style": "bezier", "line-style": "solid", "target-arrow-shape": "triangle", "source-arrow-shape": "none", "font-family": "SansSerif"}}, {"selector": "edge[weight>1]", "css": {"opacity": 1.0, "text-opacity": 1.0, "font-size": 12, "font-weight": "normal", "width": "data(width)", "color": "rgb(0,0,0)", "line-color": "rgb(51,51,51)", "source-arrow-color": "rgb(0,0,0)", "target-arrow-color": "rgb(51,51,51)", "curve-style": "bezier", "line-style": "solid", "target-arrow-shape": "triangle", "source-arrow-shape": "none", "font-family": "SansSerif"}}, {"selector": "edge[interactionType *= 'inhibit'],edge[interactionType *= 'deactivat']", "css": {"text-opacity": 1.0, "font-size": 12, "font-weight": "normal", "font-family": "SansSerif", "curve-style": "bezier", "source-arrow-shape": "none", "target-arrow-shape": "tee", "color": "rgb(0,0,0)", "line-color": "rgb(51,51,51)", "source-arrow-color": "rgb(0,0,0)", "target-arrow-color": "rgb(51,51,51)", "line-style": "solid"}}, {"selector": "edge:selected", "css": {"line-color": "rgb(255,0,0)", "label": "data(interactionType)"}}, {"selector": "edge[weight>1]:selected", "css": {"background-color": "yellow", "line-color": "red", "label": "data(label)", "target-arrow-color": "black", "source-arrow-color": "black", "width": "data(width)"}}]).update();
    
                    // Add tippy for each node
                    cy.nodes().forEach(function (n) {
                        n.data()['tip'] = makeTippy(n, n.data('id'));
                    });
    
                    // hide tippy text on click
                    cy.on('tap', 'node', function (evt) {
                        var ele = evt.target;
                        if (ele.data()['tip']['state']['visible']) {
                            ele.data()['tip'].hide();
                        } else {
                            ele.data()['tip'].show();
                        }
                    });
    
                    // put the png data in an img tag
                    let downloadButton = document.createElement("BUTTON");
                    downloadButton.id = 'dbutton';
                    downloadButton.innerHTML = '<i class="fa fa-download" aria-hidden="true"></i>';
                    downloadButton.addEventListener('click', function () {
                        let element = document.createElement('a');
                        element.setAttribute('href', cy.png({scale: 3}));
                        element.setAttribute('download', 'graph.png');
                        element.style.display = 'none';
                        document.body.appendChild(element);
                        element.click();
                        document.body.removeChild(element);
                    });
                    let p = document.getElementById('cyabc00e71-4be2-477c-abe9-3987594f8f4d');
                    p.parentElement.append(downloadButton);
    
                });
    
    
        </script>
    </head>
    
    <body>
    <div id="cyabc00e71-4be2-477c-abe9-3987594f8f4d"></div>
    <!-- When only #uuid div is placed on this page,
    the height of output-box on ipynb will be 0px.
    One line below will prevent that. -->
    <div id="dummy"
         style="width:100px;height:700px"></div>
    </body>
    
    </html>


Running enrichment analysis via EnrichR
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MAGINE automates the upload samples to EnrichR and collated the results
into a user friends format (``EnrichmentResult`` Class).

.. code:: ipython3

    from magine.enrichment.enrichr import Enrichr

.. code:: ipython3

    e = Enrichr()

Consistantly down RNA

.. code:: ipython3

    const_dn_lf = exp_data.label_free.up.require_n_sig(n_sig=2, index='label')
    const_dn_lf.heatmap();



.. image:: Tutorial_files/Tutorial_101_0.png


.. code:: ipython3

    # single sample, single background gene set
    enrich_dn_lf = e.run(const_dn_lf.id_list, gene_set_lib= 'Reactome_2016')
    enrich_dn_lf.head(10)

.. code:: ipython3

    enrich_dn_lf.term_name = enrich_dn_lf.term_name.str.split('_').str.get(0)
    enrich_dn_lf.head(10)




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>term_name</th>
          <th>rank</th>
          <th>p_value</th>
          <th>z_score</th>
          <th>combined_score</th>
          <th>adj_p_value</th>
          <th>genes</th>
          <th>n_genes</th>
          <th>db</th>
          <th>significant</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>mrna splicing - major pathway</td>
          <td>1</td>
          <td>0.000003</td>
          <td>17.523290</td>
          <td>224.549447</td>
          <td>0.000275</td>
          <td>FUS,HNRNPA2B1,RBMX,SNRPF,SRSF1,SRSF7</td>
          <td>6</td>
          <td>Reactome_2016</td>
          <td>True</td>
        </tr>
        <tr>
          <th>1</th>
          <td>mrna splicing</td>
          <td>2</td>
          <td>0.000004</td>
          <td>16.245283</td>
          <td>201.396704</td>
          <td>0.000275</td>
          <td>FUS,HNRNPA2B1,RBMX,SNRPF,SRSF1,SRSF7</td>
          <td>6</td>
          <td>Reactome_2016</td>
          <td>True</td>
        </tr>
        <tr>
          <th>2</th>
          <td>apoptotic cleavage of cellular proteins</td>
          <td>3</td>
          <td>0.000004</td>
          <td>43.874380</td>
          <td>543.247522</td>
          <td>0.000275</td>
          <td>GSN,LMNA,PLEC,VIM</td>
          <td>4</td>
          <td>Reactome_2016</td>
          <td>True</td>
        </tr>
        <tr>
          <th>3</th>
          <td>caspase-mediated cleavage of cytoskeletal prot...</td>
          <td>4</td>
          <td>0.000005</td>
          <td>118.642857</td>
          <td>1442.061152</td>
          <td>0.000275</td>
          <td>GSN,PLEC,VIM</td>
          <td>3</td>
          <td>Reactome_2016</td>
          <td>True</td>
        </tr>
        <tr>
          <th>4</th>
          <td>the citric acid (tca) cycle and respiratory el...</td>
          <td>5</td>
          <td>0.000006</td>
          <td>15.243743</td>
          <td>183.651472</td>
          <td>0.000275</td>
          <td>ATP5D,ATP5F1,ATP5L,SLC16A1,UQCRC2,UQCRFS1</td>
          <td>6</td>
          <td>Reactome_2016</td>
          <td>True</td>
        </tr>
        <tr>
          <th>5</th>
          <td>formation of atp by chemiosmotic coupling</td>
          <td>6</td>
          <td>0.000013</td>
          <td>82.120879</td>
          <td>922.112457</td>
          <td>0.000510</td>
          <td>ATP5D,ATP5F1,ATP5L</td>
          <td>3</td>
          <td>Reactome_2016</td>
          <td>True</td>
        </tr>
        <tr>
          <th>6</th>
          <td>apoptotic execution  phase</td>
          <td>7</td>
          <td>0.000015</td>
          <td>30.783752</td>
          <td>341.143188</td>
          <td>0.000510</td>
          <td>GSN,LMNA,PLEC,VIM</td>
          <td>4</td>
          <td>Reactome_2016</td>
          <td>True</td>
        </tr>
        <tr>
          <th>7</th>
          <td>respiratory electron transport, atp synthesis ...</td>
          <td>8</td>
          <td>0.000017</td>
          <td>17.661147</td>
          <td>193.585769</td>
          <td>0.000510</td>
          <td>ATP5D,ATP5F1,ATP5L,UQCRC2,UQCRFS1</td>
          <td>5</td>
          <td>Reactome_2016</td>
          <td>True</td>
        </tr>
        <tr>
          <th>8</th>
          <td>processing of capped intron-containing pre-mrna</td>
          <td>9</td>
          <td>0.000022</td>
          <td>11.958834</td>
          <td>128.249699</td>
          <td>0.000575</td>
          <td>FUS,HNRNPA2B1,RBMX,SNRPF,SRSF1,SRSF7</td>
          <td>6</td>
          <td>Reactome_2016</td>
          <td>True</td>
        </tr>
        <tr>
          <th>9</th>
          <td>mrna splicing - minor pathway</td>
          <td>10</td>
          <td>0.000486</td>
          <td>21.747813</td>
          <td>165.908608</td>
          <td>0.011427</td>
          <td>SNRPF,SRSF1,SRSF7</td>
          <td>3</td>
          <td>Reactome_2016</td>
          <td>True</td>
        </tr>
      </tbody>
    </table>
    </div>



.. code:: ipython3

    from magine.plotting.wordcloud_tools import create_wordcloud
    wc = create_wordcloud(enrich_dn_lf.sig)
    wc.plot();



.. image:: Tutorial_files/Tutorial_104_0.png


Multi-timepoint enrichment analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We can also run enrichment analysis across multiple samples. The
resulting output is similar to the ExperimentalData class, but designed
around functions for futrther enrichment analysis. This class is the
EnrichementResult class. First let’s run the enrichment for each
ph-silac sample. We need a list of those species by sample and the
sample ids to place them in our results. We can access these with :
``exp_data.ph_silac.sig.up_by_sample``,
``exp_data.ph_silac.sig.sample_ids``. Here we will run against the
Reactome_2016 gene set, but you can also run it across a list of gene
sets.

.. code:: ipython3

    ph_silac_enrichment = e.run_samples(
        exp_data.ph_silac.sig.up_by_sample, 
        exp_data.ph_silac.sig.sample_ids,
        gene_set_lib='Reactome_2016'
    )

.. code:: ipython3

    ph_silac_enrichment.head(10)




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>term_name</th>
          <th>rank</th>
          <th>p_value</th>
          <th>z_score</th>
          <th>combined_score</th>
          <th>adj_p_value</th>
          <th>genes</th>
          <th>n_genes</th>
          <th>db</th>
          <th>significant</th>
          <th>sample_id</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>cell cycle_hsa-1640170</td>
          <td>1</td>
          <td>6.097256e-07</td>
          <td>3.317719</td>
          <td>47.477413</td>
          <td>0.00043</td>
          <td>ACD,AKAP9,BRCA1,CDC16,CDC20,CDC7,CLASP2,DCTN1,...</td>
          <td>26</td>
          <td>Reactome_2016</td>
          <td>True</td>
          <td>01hr</td>
        </tr>
        <tr>
          <th>1</th>
          <td>interleukin-2 signaling_hsa-451927</td>
          <td>2</td>
          <td>1.697041e-06</td>
          <td>4.580979</td>
          <td>60.865751</td>
          <td>0.00055</td>
          <td>AKAP9,BRAF,CNKSR2,CUL3,HAVCR2,INPPL1,IRS2,MAPK...</td>
          <td>16</td>
          <td>Reactome_2016</td>
          <td>True</td>
          <td>01hr</td>
        </tr>
        <tr>
          <th>2</th>
          <td>interleukin-3, 5 and gm-csf signaling_hsa-512988</td>
          <td>3</td>
          <td>2.677495e-06</td>
          <td>4.410658</td>
          <td>56.591512</td>
          <td>0.00055</td>
          <td>AKAP9,BRAF,CNKSR2,CUL3,INPPL1,IRS2,MAPK3,MARK3...</td>
          <td>16</td>
          <td>Reactome_2016</td>
          <td>True</td>
          <td>01hr</td>
        </tr>
        <tr>
          <th>3</th>
          <td>interleukin receptor shc signaling_hsa-912526</td>
          <td>4</td>
          <td>5.589099e-06</td>
          <td>4.392809</td>
          <td>53.129672</td>
          <td>0.00055</td>
          <td>AKAP9,BRAF,CNKSR2,CUL3,INPPL1,IRS2,MAPK3,MARK3...</td>
          <td>15</td>
          <td>Reactome_2016</td>
          <td>True</td>
          <td>01hr</td>
        </tr>
        <tr>
          <th>4</th>
          <td>signalling by ngf_hsa-166520</td>
          <td>5</td>
          <td>5.968762e-06</td>
          <td>3.332650</td>
          <td>40.088346</td>
          <td>0.00055</td>
          <td>AKAP13,AKAP9,ARHGEF16,BRAF,CNKSR2,CUL3,HDAC1,I...</td>
          <td>21</td>
          <td>Reactome_2016</td>
          <td>True</td>
          <td>01hr</td>
        </tr>
        <tr>
          <th>5</th>
          <td>mapk family signaling cascades_hsa-5683057</td>
          <td>6</td>
          <td>7.859464e-06</td>
          <td>4.027363</td>
          <td>47.336790</td>
          <td>0.00055</td>
          <td>AKAP9,BRAF,CNKSR2,CUL3,DNAJB1,IRS2,MAPK3,MARK3...</td>
          <td>16</td>
          <td>Reactome_2016</td>
          <td>True</td>
          <td>01hr</td>
        </tr>
        <tr>
          <th>6</th>
          <td>insulin receptor signalling cascade_hsa-74751</td>
          <td>7</td>
          <td>8.967430e-06</td>
          <td>3.982165</td>
          <td>46.280367</td>
          <td>0.00055</td>
          <td>AKAP9,BRAF,CNKSR2,CUL3,INSR,IRS2,MAPK3,MARK3,P...</td>
          <td>16</td>
          <td>Reactome_2016</td>
          <td>True</td>
          <td>01hr</td>
        </tr>
        <tr>
          <th>7</th>
          <td>signaling by interleukins_hsa-449147</td>
          <td>8</td>
          <td>9.805077e-06</td>
          <td>3.453619</td>
          <td>39.829245</td>
          <td>0.00055</td>
          <td>AKAP9,BRAF,CNKSR2,CUL3,HAVCR2,INPPL1,IRS2,MAP3...</td>
          <td>19</td>
          <td>Reactome_2016</td>
          <td>True</td>
          <td>01hr</td>
        </tr>
        <tr>
          <th>8</th>
          <td>signal attenuation_hsa-74749</td>
          <td>9</td>
          <td>1.022501e-05</td>
          <td>43.755556</td>
          <td>502.780826</td>
          <td>0.00055</td>
          <td>INSR,IRS2,MAPK3,SHC1</td>
          <td>4</td>
          <td>Reactome_2016</td>
          <td>True</td>
          <td>01hr</td>
        </tr>
        <tr>
          <th>9</th>
          <td>signaling by fgfr2_hsa-5654738</td>
          <td>10</td>
          <td>1.150341e-05</td>
          <td>3.551082</td>
          <td>40.385980</td>
          <td>0.00055</td>
          <td>AKAP9,BRAF,CNKSR2,CUL3,HNRNPA1,HNRNPM,INSR,IRS...</td>
          <td>18</td>
          <td>Reactome_2016</td>
          <td>True</td>
          <td>01hr</td>
        </tr>
      </tbody>
    </table>
    </div>



Term names from enrichR follow various formats, depending on the library
and when it was created. We provide some tools to clean them up, but
since the EnrichmentResult is based on pandas.Dataframe, we suggest
users clean/shorten names to their liking.

.. code:: ipython3

    # clean up naming of terms
    ph_silac_enrichment.term_name = ph_silac_enrichment.term_name.str.split('_').str.get(0)

.. code:: ipython3

    ph_silac_enrichment.head(10)




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>term_name</th>
          <th>rank</th>
          <th>p_value</th>
          <th>z_score</th>
          <th>combined_score</th>
          <th>adj_p_value</th>
          <th>genes</th>
          <th>n_genes</th>
          <th>db</th>
          <th>significant</th>
          <th>sample_id</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>cell cycle</td>
          <td>1</td>
          <td>6.097256e-07</td>
          <td>3.317719</td>
          <td>47.477413</td>
          <td>0.00043</td>
          <td>ACD,AKAP9,BRCA1,CDC16,CDC20,CDC7,CLASP2,DCTN1,...</td>
          <td>26</td>
          <td>Reactome_2016</td>
          <td>True</td>
          <td>01hr</td>
        </tr>
        <tr>
          <th>1</th>
          <td>interleukin-2 signaling</td>
          <td>2</td>
          <td>1.697041e-06</td>
          <td>4.580979</td>
          <td>60.865751</td>
          <td>0.00055</td>
          <td>AKAP9,BRAF,CNKSR2,CUL3,HAVCR2,INPPL1,IRS2,MAPK...</td>
          <td>16</td>
          <td>Reactome_2016</td>
          <td>True</td>
          <td>01hr</td>
        </tr>
        <tr>
          <th>2</th>
          <td>interleukin-3, 5 and gm-csf signaling</td>
          <td>3</td>
          <td>2.677495e-06</td>
          <td>4.410658</td>
          <td>56.591512</td>
          <td>0.00055</td>
          <td>AKAP9,BRAF,CNKSR2,CUL3,INPPL1,IRS2,MAPK3,MARK3...</td>
          <td>16</td>
          <td>Reactome_2016</td>
          <td>True</td>
          <td>01hr</td>
        </tr>
        <tr>
          <th>3</th>
          <td>interleukin receptor shc signaling</td>
          <td>4</td>
          <td>5.589099e-06</td>
          <td>4.392809</td>
          <td>53.129672</td>
          <td>0.00055</td>
          <td>AKAP9,BRAF,CNKSR2,CUL3,INPPL1,IRS2,MAPK3,MARK3...</td>
          <td>15</td>
          <td>Reactome_2016</td>
          <td>True</td>
          <td>01hr</td>
        </tr>
        <tr>
          <th>4</th>
          <td>signalling by ngf</td>
          <td>5</td>
          <td>5.968762e-06</td>
          <td>3.332650</td>
          <td>40.088346</td>
          <td>0.00055</td>
          <td>AKAP13,AKAP9,ARHGEF16,BRAF,CNKSR2,CUL3,HDAC1,I...</td>
          <td>21</td>
          <td>Reactome_2016</td>
          <td>True</td>
          <td>01hr</td>
        </tr>
        <tr>
          <th>5</th>
          <td>mapk family signaling cascades</td>
          <td>6</td>
          <td>7.859464e-06</td>
          <td>4.027363</td>
          <td>47.336790</td>
          <td>0.00055</td>
          <td>AKAP9,BRAF,CNKSR2,CUL3,DNAJB1,IRS2,MAPK3,MARK3...</td>
          <td>16</td>
          <td>Reactome_2016</td>
          <td>True</td>
          <td>01hr</td>
        </tr>
        <tr>
          <th>6</th>
          <td>insulin receptor signalling cascade</td>
          <td>7</td>
          <td>8.967430e-06</td>
          <td>3.982165</td>
          <td>46.280367</td>
          <td>0.00055</td>
          <td>AKAP9,BRAF,CNKSR2,CUL3,INSR,IRS2,MAPK3,MARK3,P...</td>
          <td>16</td>
          <td>Reactome_2016</td>
          <td>True</td>
          <td>01hr</td>
        </tr>
        <tr>
          <th>7</th>
          <td>signaling by interleukins</td>
          <td>8</td>
          <td>9.805077e-06</td>
          <td>3.453619</td>
          <td>39.829245</td>
          <td>0.00055</td>
          <td>AKAP9,BRAF,CNKSR2,CUL3,HAVCR2,INPPL1,IRS2,MAP3...</td>
          <td>19</td>
          <td>Reactome_2016</td>
          <td>True</td>
          <td>01hr</td>
        </tr>
        <tr>
          <th>8</th>
          <td>signal attenuation</td>
          <td>9</td>
          <td>1.022501e-05</td>
          <td>43.755556</td>
          <td>502.780826</td>
          <td>0.00055</td>
          <td>INSR,IRS2,MAPK3,SHC1</td>
          <td>4</td>
          <td>Reactome_2016</td>
          <td>True</td>
          <td>01hr</td>
        </tr>
        <tr>
          <th>9</th>
          <td>signaling by fgfr2</td>
          <td>10</td>
          <td>1.150341e-05</td>
          <td>3.551082</td>
          <td>40.385980</td>
          <td>0.00055</td>
          <td>AKAP9,BRAF,CNKSR2,CUL3,HNRNPA1,HNRNPM,INSR,IRS...</td>
          <td>18</td>
          <td>Reactome_2016</td>
          <td>True</td>
          <td>01hr</td>
        </tr>
      </tbody>
    </table>
    </div>



The ``EnrichmentResult`` Class shares the same plotting format as the
``ExperimentalData`` Class, so we can resuse a lot of the commands we
demonstrated earlier.

.. code:: ipython3

    ph_silac_enrichment.require_n_sig(n_sig=3, inplace=True)
    ph_silac_enrichment.heatmap(
        figsize=(4,16),
        linewidths=0.01,
        cluster_by_set=False
    );



.. image:: Tutorial_files/Tutorial_112_0.png


.. code:: ipython3

    print(len(ph_silac_enrichment.sig.term_name.unique()))


.. parsed-literal::

    139
    

Term compression
^^^^^^^^^^^^^^^^

The above plot seems really busy. There are 139 enriched terms. If we
look at the top ranked terms, we see that some fo them have similar
descriptions “…. fgfr signaling”. If we look at the gene list, we can
also see that some of the genes are similar. To see if there are
redundant terms that are enriched, we can calculate their similarity
with the Jaccard Index (intersection over union). |width=50|

.. |width=50| image:: https://wikimedia.org/api/rest_v1/media/math/render/svg/eaef5aa86949f49e7dc6b9c8c3dd8b233332c9e7

.. code:: ipython3

    # calculate the Jaccard Index and returns a ranked dataframe of terms and scores.
    # Higher scores means more similar terms
    d = ph_silac_enrichment.find_similar_terms('cell cycle, mitotic')
    display(d.head(20))



.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>term_name</th>
          <th>similarity_score</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>85</th>
          <td>viral messenger rna synthesis</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>490</th>
          <td>ns1 mediated effects on host pathways</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>447</th>
          <td>export of viral ribonucleoproteins from nucleus</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>61</th>
          <td>viral messenger rna synthesis</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>60</th>
          <td>glucose transport</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>59</th>
          <td>mitotic prophase</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>120</th>
          <td>m phase</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>57</th>
          <td>host interactions with influenza factors</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>56</th>
          <td>ns1 mediated effects on host pathways</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>55</th>
          <td>nuclear pore complex (npc) disassembly</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>424</th>
          <td>nuclear envelope breakdown</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>121</th>
          <td>hexose transport</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>495</th>
          <td>host interactions with influenza factors</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>497</th>
          <td>antiviral mechanism by ifn-stimulated genes</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>48</th>
          <td>nuclear envelope breakdown</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>134</th>
          <td>nuclear import of rev protein</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>498</th>
          <td>isg15 antiviral mechanism</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>468</th>
          <td>interactions of vpr with host cellular proteins</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>499</th>
          <td>glucose transport</td>
          <td>1.0</td>
        </tr>
        <tr>
          <th>116</th>
          <td>transport of ribonucleoproteins into the host ...</td>
          <td>1.0</td>
        </tr>
      </tbody>
    </table>
    </div>


We can do this for all terms and view the results in a distance matrix
(plot used for visualization purpose only).

.. code:: ipython3

    ph_silac_enrichment.dist_matrix(figsize=(12, 12));



.. image:: Tutorial_files/Tutorial_117_0.png


We can remove the redundant ones to compress the array. Here, we sort
the terms by combined_score. For each term, we calculate the Jaccard
index with all other terms. If a term falls below above a user defined
threshold, it will be removed in the resulting array. By doing so, we
minimize the total number of terms, while maximizing the information
content of the resulting array.

.. code:: ipython3

    ph_silac_enrichment_slim = ph_silac_enrichment.remove_redundant(level='dataframe')


.. parsed-literal::

    Number of rows went from 139 to 29
    

.. code:: ipython3

    # notive the reduction in size and overlap of terms
    ph_silac_enrichment_slim.dist_matrix();



.. image:: Tutorial_files/Tutorial_120_0.png


.. code:: ipython3

    ph_silac_enrichment_slim.heatmap(
        min_sig=2, 
        figsize=(4,12),
        linewidths=0.01,
        cluster_by_set=True
    );



.. image:: Tutorial_files/Tutorial_121_0.png


Important to known, we can still recover the terms removed based on the
highest level term kept.

.. code:: ipython3

    ph_silac_enrichment.show_terms_below('g2/m checkpoints').heatmap(
        linewidths=0.01, 
        convert_to_log=False,
        figsize=(3, 8));


.. parsed-literal::

    Number of rows went from 139 to 28
    


.. image:: Tutorial_files/Tutorial_123_1.png


.. code:: ipython3

    sorted(ph_silac_enrichment_slim.term_name.unique())




.. parsed-literal::

    ['activation of the ap-1 family of transcription factors',
     'antiviral mechanism by ifn-stimulated genes',
     'apoptotic cleavage of cellular proteins',
     'cellular response to heat stress',
     'cellular senescence',
     'g2/m checkpoints',
     'golgi cisternae pericentriolar stack reorganization',
     'growth hormone receptor signaling',
     'influenza life cycle',
     'late phase of hiv life cycle',
     'major pathway of rrna processing in the nucleolus',
     'mrna splicing - major pathway',
     'negative feedback regulation of mapk pathway',
     'nonsense-mediated decay (nmd)',
     'nuclear envelope breakdown',
     'nuclear import of rev protein',
     'regulation of mrna stability by proteins that bind au-rich elements',
     'rho gtpases activate paks',
     'rho gtpases activate wasps and waves',
     'rrna modification in the nucleus',
     'signal attenuation',
     'sumoylation of dna damage response and repair proteins',
     'sumoylation of dna replication proteins',
     'sumoylation of rna binding proteins',
     'transcriptional regulation by tp53',
     'transport of mature mrna derived from an intron-containing transcript',
     'transport of mature mrna derived from an intronless transcript',
     'trna processing in the nucleus',
     'vpr-mediated nuclear import of pics']



Visualize species of enrichmentt term
'''''''''''''''''''''''''''''''''''''

We can use the EnrichmentResult to extract out any given set of genes
for a term (or the entire array). For a select term, we can extract out
the species of interest to visualize. This makes chaining together the
data, networks, and enrichment output seamlessly.

.. code:: ipython3

    exp_data.ph_silac.heatmap(
        ph_silac_enrichment_slim.sig.term_to_genes('apoptotic cleavage of cellular proteins'),
        subset_index='identifier',
        index='label',
        cluster_row=False,
        sort_row='mean',
        min_sig=2,
        linewidths=0.01,
        figsize=(4, 8),
    );



.. image:: Tutorial_files/Tutorial_126_0.png


.. code:: ipython3

    exp_data.ph_silac.heatmap(
        ph_silac_enrichment_slim.sig.term_to_genes('g2/m checkpoints'),
        subset_index='identifier',
        index='label',
        cluster_row=True,
        sort_row='index',
        min_sig=2,
        linewidths=0.01,
        figsize=(3, 6),
    );



.. image:: Tutorial_files/Tutorial_127_0.png


.. code:: ipython3

    exp_data.ph_silac.heatmap(
        ph_silac_enrichment.sig.term_to_genes('apoptosis'),
        subset_index='identifier',
        index='label',
        cluster_row=True,
        sort_row='index',
        min_sig=2,
        linewidths=0.01,
        figsize=(2,12),
    );



.. image:: Tutorial_files/Tutorial_128_0.png


We can use the ExperimentalData class to filter the data to create lists of genes for further analysis. We take these lists and run enrichment analysis using Enrichr.
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

.. container::

Since this part is time consuming, it is best to do it outside of a
notebook. The code to do so can be found in “run_enrichment.py”. The
results will be a csv file that we will load next.

Annotated Gene set Network (AGN)
--------------------------------

Lastly, we created a function to generate molecular and coarse grain
networks based on enrichment terms. Users can used the compressed
enrichment result terms to generate large scale representations of their
data, or by selecting key terms of importance. Here, we are going to use
3 terms from the compressed enrichment result class.

.. code:: ipython3

    from magine.networks.annotated_set import create_asn
    # selected terms of interest
    terms=['apoptotic cleavage of cellular proteins', 
           'g2/m checkpoints',
           'transcriptional regulation by tp53']
    
    term_net, mol_net = create_asn(
        ph_silac_enrichment_slim, network,
        terms=terms,
        save_name='all_example',
        use_threshold=True,
        use_cytoscape=False, # If you have cytoscape open, this will create a cytoscape session if True
    )


.. parsed-literal::

    Creating ontology network
    

.. code:: ipython3

    draw_cyjs(term_net)



.. raw:: html

    <!DOCTYPE html>
    <html>
    <head>
        <meta charset=utf-8/>
        <style type="text/css">
            body {
                font: 14px helvetica neue, helvetica, arial, sans-serif;
            }
    
            #cydfc98151-0c52-4470-8edf-fbb20cdde6a5{
                height: 700px;
                width: 90%;
                border: 5px solid black;
                box-sizing: border-box;
                position: relative;
                top: 5;
                margin-bottom: -700px;
                background: white;
            }
    
        </style>
    
        <script>
    
    
            requirejs.config({
    
                paths: {
                    'popper': 'https://unpkg.com/popper.js@1.14.1/dist/umd/popper',
                    'tippy': 'https://cdnjs.cloudflare.com/ajax/libs/tippy.js/2.3.0/tippy.min',
                    'cytoscape': 'https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.2.10/cytoscape',
                    'cytoscape-popper': 'https://cdn.rawgit.com/cytoscape/cytoscape.js-popper/3ad50859/cytoscape-popper',
                    'jquery': 'https://cdnjs.cloudflare.com/ajax/libs/jquery/2.2.4/jquery.min',
                    'qtip2': 'https://cdnjs.cloudflare.com/ajax/libs/qtip2/2.2.0/basic/jquery.qtip.min',
                    'dagre': 'https://cdn.rawgit.com/cpettitt/dagre/v0.7.4/dist/dagre.min',
                    'cytoscape-dagre': 'https://cdn.rawgit.com/cytoscape/cytoscape.js-dagre/1.5.0/cytoscape-dagre',
                    'cytoscape-cose-bilkent': 'https://cdn.rawgit.com/cytoscape/cytoscape.js-cose-bilkent/1.6.1/cytoscape-cose-bilkent'
                },
                shim: {
                    'cytoscape-popper': {
                        deps: ['cytoscape', 'popper']
                    },
                    'cytoscape-dagre': {
                        deps: ['cytoscape', 'dagre']
                    }
                },
                map: {
                    '*': {
                        'popper.js': 'popper',
                        'webcola': 'cola'
                    }
                }
    
            });
    
    
            require(['cytoscape', 'cytoscape-popper', 'popper', 'tippy', 'jquery',
                    'cytoscape-cose-bilkent', 'cytoscape-dagre', 'dagre'],
                function (cytoscape, cypopper, popper, tippy, jquery, regCose,
                          cydag, dagre) {
                    console.log('Loading Cytoscape.js Module...');
                    window['popper'] = popper;
                    window['tippy'] = tippy;
                    window['cytoscape'] = cytoscape;
                    cypopper(cytoscape);
                    regCose(cytoscape);
                    cydag(cytoscape, dagre);
    
                    function makeTippy(target, text) {
                        return tippy(target.popperRef(),
                            {
                                html: add_tip(text),
                                trigger: 'manual',
                                arrow: true,
                                placement: 'top',
                                hideOnClick: false,
                                interactive: true,
                                multiple: true,
                                sticky: true
                            }).tooltips[0];
                    }
    
                    function add_tip(text) {
                        let div = document.createElement('div');
                        div.innerHTML = "<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=" + text + "' target='_blank'>Gene Card</a>";
                        return div;
                    }
    
                    let cy = window.cy = cytoscape({
                        container: $('#cydfc98151-0c52-4470-8edf-fbb20cdde6a5'),
                        elements: {
                            nodes: [{"data": {"term": "g2/m checkpoints", "label": "g2/m checkpoints", "color": "red", "sample48hr": 38.75344051566447, "sample24hr": 44.830558085253855, "sample01hr": 1.2162050956428634, "sample06hr": 20.86038580686977, "id": "g2/m checkpoints", "name": "g2/m checkpoints"}}, {"data": {"term": "apoptotic cleavage of cellular proteins", "label": "apoptotic cleavage of cellular proteins", "color": "red", "sample48hr": 151.57345071680933, "sample24hr": 148.5683161707896, "sample01hr": 48.16585615888427, "sample06hr": 38.12803694041321, "id": "apoptotic cleavage of cellular proteins", "name": "apoptotic cleavage of cellular proteins"}}, {"data": {"term": "transcriptional regulation by tp53", "label": "transcriptional regulation by tp53", "color": "red", "sample48hr": 34.529861832956605, "sample24hr": 9.210549049725085, "sample01hr": 6.22853371200079, "sample06hr": 12.766465748111214, "id": "transcriptional regulation by tp53", "name": "transcriptional regulation by tp53"}}],
                            edges: [{"data": {"label": "18", "weight": 18, "pvalue": 1.3992091791086609e-24, "adjPval": 6.996045895543304e-24, "width": 10.0, "source": "g2/m checkpoints", "target": "apoptotic cleavage of cellular proteins"}}, {"data": {"label": "12", "weight": 12, "pvalue": 2.569429978480065e-10, "adjPval": 6.423574946200162e-10, "width": 7.2, "source": "g2/m checkpoints", "target": "transcriptional regulation by tp53"}}, {"data": {"label": "3", "weight": 3, "pvalue": 0.029301930323811457, "adjPval": 0.029301930323811457, "width": 3.0, "source": "apoptotic cleavage of cellular proteins", "target": "transcriptional regulation by tp53"}}, {"data": {"label": "5", "weight": 5, "pvalue": 0.0005950646310354398, "adjPval": 0.0007438307887942997, "width": 3.9333333333333336, "source": "transcriptional regulation by tp53", "target": "apoptotic cleavage of cellular proteins"}}, {"data": {"label": "6", "weight": 6, "pvalue": 6.38682893491284e-05, "adjPval": 0.00010644714891521401, "width": 4.4, "source": "transcriptional regulation by tp53", "target": "g2/m checkpoints"}}]
                        },
                        boxSelectionEnabled: true,
                        wheelSensitivity: .25,
    
                    });
    
                    // apply layout
                    cy.layout({"name": "cose-bilkent", "spacingFactor": 1.5, "animate": false}).run();
    
                    // add style
                    cy.style().fromJson([{"selector": "node", "css": {"text-opacity": 1.0, "background-opacity": 1.0, "font-weight": "normal", "color": "rgb(0,153,234)", "border-width": 3.0, "border-color": "rgb(51,51,51)", "font-size": 9, "height": 50, "width": 50, "label": "data(name)", "text-wrap": "wrap", "text-max-width": 50, "shape": "ellipse", "background-color": "data(color)", "text-halign": "center", "font-family": "SansSerif", "text-valign": "center", "border-opacity": 1.0}}, {"selector": "node[speciesType = 'compound']", "css": {"text-opacity": 1.0, "background-opacity": 1.0, "font-weight": "normal", "color": "rgb(0,153,234)", "border-width": 3.0, "border-color": "rgb(51,51,51)", "font-size": 9, "height": 35, "width": 35, "label": "data(chemName)", "text-wrap": "wrap", "text-max-width": 95, "shape": "ellipse", "background-color": "rgb(255,255,255)", "text-halign": "center", "font-family": "SansSerif", "text-valign": "center", "border-opacity": 1.0}}, {"selector": "$node > node", "css": {"font-size": 20, "shape": "ellipse", "padding-right": "10px", "padding-bottom": "10px", "padding-top": "10px", "text-valign": "top", "text-halign": "center", "background-color": "#bbb", "padding-left": "10px"}}, {"selector": "node:selected", "css": {"background-color": "rgb(255,0,102)"}}, {"selector": "edge", "css": {"opacity": 1.0, "text-opacity": 1.0, "font-size": 12, "font-weight": "normal", "width": "2", "color": "rgb(0,0,0)", "line-color": "rgb(51,51,51)", "source-arrow-color": "rgb(0,0,0)", "target-arrow-color": "rgb(51,51,51)", "curve-style": "bezier", "line-style": "solid", "target-arrow-shape": "triangle", "source-arrow-shape": "none", "font-family": "SansSerif"}}, {"selector": "edge[weight>1]", "css": {"opacity": 1.0, "text-opacity": 1.0, "font-size": 12, "font-weight": "normal", "width": "data(width)", "color": "rgb(0,0,0)", "line-color": "rgb(51,51,51)", "source-arrow-color": "rgb(0,0,0)", "target-arrow-color": "rgb(51,51,51)", "curve-style": "bezier", "line-style": "solid", "target-arrow-shape": "triangle", "source-arrow-shape": "none", "font-family": "SansSerif"}}, {"selector": "edge[interactionType *= 'inhibit'],edge[interactionType *= 'deactivat']", "css": {"text-opacity": 1.0, "font-size": 12, "font-weight": "normal", "font-family": "SansSerif", "curve-style": "bezier", "source-arrow-shape": "none", "target-arrow-shape": "tee", "color": "rgb(0,0,0)", "line-color": "rgb(51,51,51)", "source-arrow-color": "rgb(0,0,0)", "target-arrow-color": "rgb(51,51,51)", "line-style": "solid"}}, {"selector": "edge:selected", "css": {"line-color": "rgb(255,0,0)", "label": "data(interactionType)"}}, {"selector": "edge[weight>1]:selected", "css": {"background-color": "yellow", "line-color": "red", "label": "data(label)", "target-arrow-color": "black", "source-arrow-color": "black", "width": "data(width)"}}]).update();
    
                    // Add tippy for each node
                    cy.nodes().forEach(function (n) {
                        n.data()['tip'] = makeTippy(n, n.data('id'));
                    });
    
                    // hide tippy text on click
                    cy.on('tap', 'node', function (evt) {
                        var ele = evt.target;
                        if (ele.data()['tip']['state']['visible']) {
                            ele.data()['tip'].hide();
                        } else {
                            ele.data()['tip'].show();
                        }
                    });
    
                    // put the png data in an img tag
                    let downloadButton = document.createElement("BUTTON");
                    downloadButton.id = 'dbutton';
                    downloadButton.innerHTML = '<i class="fa fa-download" aria-hidden="true"></i>';
                    downloadButton.addEventListener('click', function () {
                        let element = document.createElement('a');
                        element.setAttribute('href', cy.png({scale: 3}));
                        element.setAttribute('download', 'graph.png');
                        element.style.display = 'none';
                        document.body.appendChild(element);
                        element.click();
                        document.body.removeChild(element);
                    });
                    let p = document.getElementById('cydfc98151-0c52-4470-8edf-fbb20cdde6a5');
                    p.parentElement.append(downloadButton);
    
                });
    
    
        </script>
    </head>
    
    <body>
    <div id="cydfc98151-0c52-4470-8edf-fbb20cdde6a5"></div>
    <!-- When only #uuid div is placed on this page,
    the height of output-box on ipynb will be 0px.
    One line below will prevent that. -->
    <div id="dummy"
         style="width:100px;height:700px"></div>
    </body>
    
    </html>


.. code:: ipython3

    draw_cyjs(mol_net, add_parent=True)



.. raw:: html

    <!DOCTYPE html>
    <html>
    <head>
        <meta charset=utf-8/>
        <style type="text/css">
            body {
                font: 14px helvetica neue, helvetica, arial, sans-serif;
            }
    
            #cyb63e02fc-a04e-4038-add4-84a305eef936{
                height: 700px;
                width: 90%;
                border: 5px solid black;
                box-sizing: border-box;
                position: relative;
                top: 5;
                margin-bottom: -700px;
                background: white;
            }
    
        </style>
    
        <script>
    
    
            requirejs.config({
    
                paths: {
                    'popper': 'https://unpkg.com/popper.js@1.14.1/dist/umd/popper',
                    'tippy': 'https://cdnjs.cloudflare.com/ajax/libs/tippy.js/2.3.0/tippy.min',
                    'cytoscape': 'https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.2.10/cytoscape',
                    'cytoscape-popper': 'https://cdn.rawgit.com/cytoscape/cytoscape.js-popper/3ad50859/cytoscape-popper',
                    'jquery': 'https://cdnjs.cloudflare.com/ajax/libs/jquery/2.2.4/jquery.min',
                    'qtip2': 'https://cdnjs.cloudflare.com/ajax/libs/qtip2/2.2.0/basic/jquery.qtip.min',
                    'dagre': 'https://cdn.rawgit.com/cpettitt/dagre/v0.7.4/dist/dagre.min',
                    'cytoscape-dagre': 'https://cdn.rawgit.com/cytoscape/cytoscape.js-dagre/1.5.0/cytoscape-dagre',
                    'cytoscape-cose-bilkent': 'https://cdn.rawgit.com/cytoscape/cytoscape.js-cose-bilkent/1.6.1/cytoscape-cose-bilkent'
                },
                shim: {
                    'cytoscape-popper': {
                        deps: ['cytoscape', 'popper']
                    },
                    'cytoscape-dagre': {
                        deps: ['cytoscape', 'dagre']
                    }
                },
                map: {
                    '*': {
                        'popper.js': 'popper',
                        'webcola': 'cola'
                    }
                }
    
            });
    
    
            require(['cytoscape', 'cytoscape-popper', 'popper', 'tippy', 'jquery',
                    'cytoscape-cose-bilkent', 'cytoscape-dagre', 'dagre'],
                function (cytoscape, cypopper, popper, tippy, jquery, regCose,
                          cydag, dagre) {
                    console.log('Loading Cytoscape.js Module...');
                    window['popper'] = popper;
                    window['tippy'] = tippy;
                    window['cytoscape'] = cytoscape;
                    cypopper(cytoscape);
                    regCose(cytoscape);
                    cydag(cytoscape, dagre);
    
                    function makeTippy(target, text) {
                        return tippy(target.popperRef(),
                            {
                                html: add_tip(text),
                                trigger: 'manual',
                                arrow: true,
                                placement: 'top',
                                hideOnClick: false,
                                interactive: true,
                                multiple: true,
                                sticky: true
                            }).tooltips[0];
                    }
    
                    function add_tip(text) {
                        let div = document.createElement('div');
                        div.innerHTML = "<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=" + text + "' target='_blank'>Gene Card</a>";
                        return div;
                    }
    
                    let cy = window.cy = cytoscape({
                        container: $('#cyb63e02fc-a04e-4038-add4-84a305eef936'),
                        elements: {
                            nodes: [{"data": {"termName": "g2/m checkpoints,transcriptional regulation by tp53", "terms": "g2/m checkpoints,transcriptional regulation by tp53", "parent": "g2/m checkpoints,transcriptional regulation by tp53", "color": "white", "id": "RAD9A", "name": "RAD9A"}}, {"data": {"termName": "g2/m checkpoints,transcriptional regulation by tp53", "terms": "g2/m checkpoints,transcriptional regulation by tp53", "parent": "g2/m checkpoints,transcriptional regulation by tp53", "color": "white", "id": "TP53", "name": "TP53"}}, {"data": {"termName": "g2/m checkpoints,transcriptional regulation by tp53", "terms": "g2/m checkpoints,transcriptional regulation by tp53", "parent": "g2/m checkpoints,transcriptional regulation by tp53", "color": "white", "id": "NBN", "name": "NBN"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "parent": "g2/m checkpoints", "color": "white", "id": "MCM2", "name": "MCM2"}}, {"data": {"termName": "g2/m checkpoints,transcriptional regulation by tp53", "terms": "g2/m checkpoints,transcriptional regulation by tp53", "parent": "g2/m checkpoints,transcriptional regulation by tp53", "color": "white", "id": "TOPBP1", "name": "TOPBP1"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "parent": "g2/m checkpoints", "color": "white", "id": "PSMD4", "name": "PSMD4"}}, {"data": {"termName": "apoptotic cleavage of cellular proteins", "terms": "apoptotic cleavage of cellular proteins", "parent": "apoptotic cleavage of cellular proteins", "color": "white", "id": "APC", "name": "APC"}}, {"data": {"termName": "apoptotic cleavage of cellular proteins", "terms": "apoptotic cleavage of cellular proteins", "parent": "apoptotic cleavage of cellular proteins", "color": "white", "id": "CTNNB1", "name": "CTNNB1"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "parent": "g2/m checkpoints", "color": "white", "id": "PSMB3", "name": "PSMB3"}}, {"data": {"termName": "apoptotic cleavage of cellular proteins", "terms": "apoptotic cleavage of cellular proteins", "parent": "apoptotic cleavage of cellular proteins", "color": "white", "id": "PTK2", "name": "PTK2"}}, {"data": {"termName": "apoptotic cleavage of cellular proteins", "terms": "apoptotic cleavage of cellular proteins", "parent": "apoptotic cleavage of cellular proteins", "color": "white", "id": "SPTAN1", "name": "SPTAN1"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "parent": "g2/m checkpoints", "color": "white", "id": "PSMF1", "name": "PSMF1"}}, {"data": {"termName": "apoptotic cleavage of cellular proteins", "terms": "apoptotic cleavage of cellular proteins", "parent": "apoptotic cleavage of cellular proteins", "color": "white", "id": "TJP1", "name": "TJP1"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "parent": "g2/m checkpoints", "color": "white", "id": "PSMB10", "name": "PSMB10"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "parent": "g2/m checkpoints", "color": "white", "id": "PSMA5", "name": "PSMA5"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "parent": "g2/m checkpoints", "color": "white", "id": "CDC7", "name": "CDC7"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "parent": "g2/m checkpoints", "color": "white", "id": "MCM7", "name": "MCM7"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "parent": "g2/m checkpoints", "color": "white", "id": "MCM6", "name": "MCM6"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "parent": "g2/m checkpoints", "color": "white", "id": "PSMA3", "name": "PSMA3"}}, {"data": {"termName": "g2/m checkpoints,transcriptional regulation by tp53", "terms": "g2/m checkpoints,transcriptional regulation by tp53", "parent": "g2/m checkpoints,transcriptional regulation by tp53", "color": "white", "id": "YWHAZ", "name": "YWHAZ"}}, {"data": {"termName": "g2/m checkpoints,transcriptional regulation by tp53", "terms": "g2/m checkpoints,transcriptional regulation by tp53", "parent": "g2/m checkpoints,transcriptional regulation by tp53", "color": "white", "id": "ATR", "name": "ATR"}}, {"data": {"termName": "g2/m checkpoints,transcriptional regulation by tp53", "terms": "g2/m checkpoints,transcriptional regulation by tp53", "parent": "g2/m checkpoints,transcriptional regulation by tp53", "color": "white", "id": "SFN", "name": "SFN"}}, {"data": {"termName": "apoptotic cleavage of cellular proteins", "terms": "apoptotic cleavage of cellular proteins", "parent": "apoptotic cleavage of cellular proteins", "color": "white", "id": "LMNB1", "name": "LMNB1"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "parent": "g2/m checkpoints", "color": "white", "id": "UIMC1", "name": "UIMC1"}}, {"data": {"termName": "g2/m checkpoints,transcriptional regulation by tp53", "terms": "g2/m checkpoints,transcriptional regulation by tp53", "parent": "g2/m checkpoints,transcriptional regulation by tp53", "color": "white", "id": "BRCA1", "name": "BRCA1"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "parent": "g2/m checkpoints", "color": "white", "id": "PSMD1", "name": "PSMD1"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "parent": "g2/m checkpoints", "color": "white", "id": "PSMB7", "name": "PSMB7"}}, {"data": {"termName": "g2/m checkpoints,transcriptional regulation by tp53", "terms": "g2/m checkpoints,transcriptional regulation by tp53", "parent": "g2/m checkpoints,transcriptional regulation by tp53", "color": "white", "id": "RAD50", "name": "RAD50"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "parent": "transcriptional regulation by tp53", "color": "white", "id": "MAPK14", "name": "MAPK14"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "parent": "transcriptional regulation by tp53", "color": "white", "id": "JUN", "name": "JUN"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "parent": "transcriptional regulation by tp53", "color": "white", "id": "ATF2", "name": "ATF2"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "parent": "transcriptional regulation by tp53", "color": "white", "id": "TSC2", "name": "TSC2"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "parent": "transcriptional regulation by tp53", "color": "white", "id": "TP53BP2", "name": "TP53BP2"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "parent": "transcriptional regulation by tp53", "color": "white", "id": "RRM2B", "name": "RRM2B"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "parent": "transcriptional regulation by tp53", "color": "white", "id": "TXN", "name": "TXN"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "parent": "transcriptional regulation by tp53", "color": "white", "id": "TAF1L", "name": "TAF1L"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "parent": "transcriptional regulation by tp53", "color": "white", "id": "GPI", "name": "GPI"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "parent": "transcriptional regulation by tp53", "color": "white", "id": "G6PD", "name": "G6PD"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "parent": "transcriptional regulation by tp53", "color": "white", "id": "RFFL", "name": "RFFL"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "parent": "transcriptional regulation by tp53", "color": "white", "id": "HDAC1", "name": "HDAC1"}}, {"data": {"termName": "g2/m checkpoints,transcriptional regulation by tp53", "terms": "g2/m checkpoints,transcriptional regulation by tp53", "parent": "g2/m checkpoints,transcriptional regulation by tp53", "color": "white", "id": "YWHAE", "name": "YWHAE"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "parent": "transcriptional regulation by tp53", "color": "white", "id": "TAF9", "name": "TAF9"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "parent": "transcriptional regulation by tp53", "color": "white", "id": "FAS", "name": "FAS"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "parent": "transcriptional regulation by tp53", "color": "white", "id": "BNIP3L", "name": "BNIP3L"}}, {"data": {"termName": "apoptotic cleavage of cellular proteins", "terms": "apoptotic cleavage of cellular proteins", "parent": "apoptotic cleavage of cellular proteins", "color": "white", "id": "VIM", "name": "VIM"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "parent": "transcriptional regulation by tp53", "color": "white", "id": "MAPKAP1", "name": "MAPKAP1"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "parent": "transcriptional regulation by tp53", "color": "white", "id": "MTOR", "name": "MTOR"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "parent": "transcriptional regulation by tp53", "color": "white", "id": "LAMTOR1", "name": "LAMTOR1"}}, {"data": {"termName": "g2/m checkpoints,transcriptional regulation by tp53", "terms": "g2/m checkpoints,transcriptional regulation by tp53", "parent": "g2/m checkpoints,transcriptional regulation by tp53", "color": "white", "id": "YWHAQ", "name": "YWHAQ"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "parent": "transcriptional regulation by tp53", "color": "white", "id": "SSRP1", "name": "SSRP1"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "parent": "transcriptional regulation by tp53", "color": "white", "id": "PPP2CA", "name": "PPP2CA"}}, {"data": {"termName": "g2/m checkpoints,transcriptional regulation by tp53", "terms": "g2/m checkpoints,transcriptional regulation by tp53", "parent": "g2/m checkpoints,transcriptional regulation by tp53", "color": "white", "id": "YWHAH", "name": "YWHAH"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "parent": "transcriptional regulation by tp53", "color": "white", "id": "RBL1", "name": "RBL1"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "parent": "transcriptional regulation by tp53", "color": "white", "id": "SUPT16H", "name": "SUPT16H"}}, {"data": {"color": "white", "id": "apoptotic cleavage of cellular proteins", "name": "apoptotic cleavage of cellular proteins"}}, {"data": {"color": "white", "id": "g2/m checkpoints", "name": "g2/m checkpoints"}}, {"data": {"color": "white", "id": "g2/m checkpoints,transcriptional regulation by tp53", "name": "g2/m checkpoints,transcriptional regulation by tp53"}}, {"data": {"color": "white", "id": "transcriptional regulation by tp53", "name": "transcriptional regulation by tp53"}}],
                            edges: [{"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "RAD9A", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate|binding", "width": 10.0, "source": "RAD9A", "target": "NBN"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "RAD9A", "target": "MCM2"}}, {"data": {"databaseSource": "SIGNOR", "interactionType": "activate|binding", "width": 10.0, "source": "RAD9A", "target": "TOPBP1"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "TP53", "target": "SFN"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "TP53", "target": "LMNB1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "TP53", "target": "CTNNB1"}}, {"data": {"databaseSource": "KEGG|ReactomeFI|SIGNOR", "interactionType": "activate|expression|indirect", "width": 10.0, "source": "TP53", "target": "FAS"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "TP53", "target": "TSC2"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "TP53", "target": "RRM2B"}}, {"data": {"databaseSource": "KEGG", "interactionType": "expression", "width": 10.0, "source": "TP53", "target": "BNIP3L"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "NBN", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "TOPBP1", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate|binding", "width": 10.0, "source": "TOPBP1", "target": "NBN"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "TOPBP1", "target": "MCM2"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMD4", "target": "APC"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMD4", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMD4", "target": "CTNNB1"}}, {"data": {"databaseSource": "KEGG|SIGNOR", "interactionType": "binding|inhibit", "width": 10.0, "source": "APC", "target": "CTNNB1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate|binding", "width": 10.0, "source": "CTNNB1", "target": "TJP1"}}, {"data": {"databaseSource": "KEGG", "interactionType": "expression", "width": 10.0, "source": "CTNNB1", "target": "JUN"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate|inhibit", "width": 10.0, "source": "CTNNB1", "target": "HDAC1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMB3", "target": "APC"}}, {"data": {"databaseSource": "BIOGRID|ReactomeFI", "interactionType": "catalyze|cleavage", "pubmedId": "20080206", "width": 10.0, "source": "PSMB3", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMB3", "target": "CTNNB1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "binding|catalyze", "width": 10.0, "source": "PTK2", "target": "SPTAN1"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "PTK2", "target": "MAPK14"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMF1", "target": "APC"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMF1", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMF1", "target": "CTNNB1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMB10", "target": "APC"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMB10", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMB10", "target": "CTNNB1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMA5", "target": "APC"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMA5", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMA5", "target": "CTNNB1"}}, {"data": {"databaseSource": "BIOGRID|ReactomeFI|SIGNOR", "interactionType": "activate|binding|catalyze|phosphorylate", "pubmedId": "10846177", "width": 10.0, "source": "CDC7", "target": "MCM7"}}, {"data": {"databaseSource": "BIOGRID|ReactomeFI", "interactionType": "activate|binding|catalyze|phosphorylate", "pubmedId": "10846177", "width": 10.0, "source": "CDC7", "target": "MCM6"}}, {"data": {"databaseSource": "BIOGRID|ReactomeFI|SIGNOR", "interactionType": "activate|binding|catalyze|phosphorylate", "pubmedId": "10846177", "width": 10.0, "source": "CDC7", "target": "MCM2"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMA3", "target": "APC"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMA3", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMA3", "target": "CTNNB1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "YWHAZ", "target": "RAD9A"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "YWHAZ", "target": "ATR"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "YWHAZ", "target": "TOPBP1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "binding|inhibit", "width": 10.0, "source": "YWHAZ", "target": "TSC2"}}, {"data": {"databaseSource": "BIOGRID|KEGG|ReactomeFI|SIGNOR", "interactionType": "activate|binding|phosphorylate", "pubmedId": "11278964", "width": 10.0, "source": "ATR", "target": "BRCA1"}}, {"data": {"databaseSource": "BIOGRID|KEGG|ReactomeFI|SIGNOR", "interactionType": "activate|catalyze|phosphorylate", "pubmedId": "10608806", "width": 10.0, "source": "ATR", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate|binding", "width": 10.0, "source": "ATR", "target": "RAD9A"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "ATR", "target": "CDC7"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "ATR", "target": "MCM7"}}, {"data": {"databaseSource": "ReactomeFI|SIGNOR", "interactionType": "activate|catalyze|phosphorylate", "width": 10.0, "source": "ATR", "target": "MCM6"}}, {"data": {"databaseSource": "BIOGRID|ReactomeFI|SIGNOR", "interactionType": "activate|catalyze|phosphorylate", "pubmedId": "15210935", "width": 10.0, "source": "ATR", "target": "MCM2"}}, {"data": {"databaseSource": "ReactomeFI|SIGNOR", "interactionType": "activate|binding|phosphorylate", "width": 10.0, "source": "ATR", "target": "NBN"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "binding|inhibit", "width": 10.0, "source": "SFN", "target": "TSC2"}}, {"data": {"databaseSource": "SIGNOR", "interactionType": "activate|binding", "width": 10.0, "source": "UIMC1", "target": "BRCA1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "binding|catalyze", "width": 10.0, "source": "BRCA1", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "BRCA1", "target": "PPP2CA"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "BRCA1", "target": "MAPK14"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "binding|expression", "width": 10.0, "source": "BRCA1", "target": "ATF2"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMD1", "target": "APC"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMD1", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMD1", "target": "CTNNB1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMB7", "target": "APC"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMB7", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMB7", "target": "CTNNB1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "RAD50", "target": "TP53"}}, {"data": {"databaseSource": "SIGNOR", "interactionType": "activate|binding", "width": 10.0, "source": "RAD50", "target": "BRCA1"}}, {"data": {"databaseSource": "BIOGRID|KEGG", "interactionType": "activate|phosphorylate", "pubmedId": "25483191", "width": 10.0, "source": "MAPK14", "target": "JUN"}}, {"data": {"databaseSource": "BIOGRID|KEGG|ReactomeFI|SIGNOR", "interactionType": "activate|catalyze|phosphorylate", "pubmedId": "10581258", "width": 10.0, "source": "MAPK14", "target": "TP53"}}, {"data": {"databaseSource": "BIOGRID|KEGG|ReactomeFI|SIGNOR", "interactionType": "activate|catalyze|phosphorylate", "pubmedId": "10581258", "width": 10.0, "source": "MAPK14", "target": "ATF2"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "MAPK14", "target": "TSC2"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "MAPK14", "target": "YWHAZ"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "JUN", "target": "FAS"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "binding|expression", "width": 10.0, "source": "JUN", "target": "TP53"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "JUN", "target": "VIM"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "JUN", "target": "BRCA1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "JUN", "target": "HDAC1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "JUN", "target": "TXN"}}, {"data": {"databaseSource": "KEGG", "interactionType": "expression", "width": 10.0, "source": "ATF2", "target": "VIM"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "binding|expression|inhibit", "width": 10.0, "source": "ATF2", "target": "JUN"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "TSC2", "target": "MAPKAP1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "TSC2", "target": "MTOR"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "TSC2", "target": "LAMTOR1"}}, {"data": {"databaseSource": "SIGNOR", "interactionType": "activate|binding", "width": 10.0, "source": "TP53BP2", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "RRM2B", "target": "TXN"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "TAF1L", "target": "TP53"}}, {"data": {"databaseSource": "KEGG", "interactionType": "chemical", "width": 10.0, "source": "GPI", "target": "G6PD"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "RFFL", "target": "TP53"}}, {"data": {"pubmedId": "11099047", "databaseSource": "BIOGRID", "interactionType": "deacetylate", "width": 10.0, "source": "HDAC1", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "binding|inhibit", "width": 10.0, "source": "YWHAE", "target": "TSC2"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze|inhibit", "width": 10.0, "source": "TAF9", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "binding|inhibit", "width": 10.0, "source": "YWHAQ", "target": "TSC2"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "SSRP1", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "PPP2CA", "target": "SFN"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate|inhibit", "width": 10.0, "source": "PPP2CA", "target": "YWHAE"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate|inhibit", "width": 10.0, "source": "PPP2CA", "target": "YWHAQ"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate|inhibit", "width": 10.0, "source": "PPP2CA", "target": "YWHAH"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate|inhibit", "width": 10.0, "source": "PPP2CA", "target": "YWHAZ"}}, {"data": {"databaseSource": "SIGNOR", "interactionType": "activate|dephosphorylate", "width": 10.0, "source": "PPP2CA", "target": "RBL1"}}, {"data": {"databaseSource": "SIGNOR", "interactionType": "activate|dephosphorylate", "width": 10.0, "source": "PPP2CA", "target": "CTNNB1"}}, {"data": {"databaseSource": "SIGNOR", "interactionType": "activate|dephosphorylate", "width": 10.0, "source": "PPP2CA", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "binding|inhibit", "width": 10.0, "source": "YWHAH", "target": "TSC2"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "SUPT16H", "target": "TP53"}}]
                        },
                        boxSelectionEnabled: true,
                        wheelSensitivity: .25,
    
                    });
    
                    // apply layout
                    cy.layout({"name": "cose-bilkent", "spacingFactor": 1.5, "animate": false}).run();
    
                    // add style
                    cy.style().fromJson([{"selector": "node", "css": {"text-opacity": 1.0, "background-opacity": 1.0, "font-weight": "normal", "color": "rgb(0,153,234)", "border-width": 3.0, "border-color": "rgb(51,51,51)", "font-size": 9, "height": 50, "width": 50, "label": "data(name)", "text-wrap": "wrap", "text-max-width": 50, "shape": "ellipse", "background-color": "data(color)", "text-halign": "center", "font-family": "SansSerif", "text-valign": "center", "border-opacity": 1.0}}, {"selector": "node[speciesType = 'compound']", "css": {"text-opacity": 1.0, "background-opacity": 1.0, "font-weight": "normal", "color": "rgb(0,153,234)", "border-width": 3.0, "border-color": "rgb(51,51,51)", "font-size": 9, "height": 35, "width": 35, "label": "data(chemName)", "text-wrap": "wrap", "text-max-width": 95, "shape": "ellipse", "background-color": "rgb(255,255,255)", "text-halign": "center", "font-family": "SansSerif", "text-valign": "center", "border-opacity": 1.0}}, {"selector": "$node > node", "css": {"font-size": 20, "shape": "ellipse", "padding-right": "10px", "padding-bottom": "10px", "padding-top": "10px", "text-valign": "top", "text-halign": "center", "background-color": "#bbb", "padding-left": "10px"}}, {"selector": "node:selected", "css": {"background-color": "rgb(255,0,102)"}}, {"selector": "edge", "css": {"opacity": 1.0, "text-opacity": 1.0, "font-size": 12, "font-weight": "normal", "width": "2", "color": "rgb(0,0,0)", "line-color": "rgb(51,51,51)", "source-arrow-color": "rgb(0,0,0)", "target-arrow-color": "rgb(51,51,51)", "curve-style": "bezier", "line-style": "solid", "target-arrow-shape": "triangle", "source-arrow-shape": "none", "font-family": "SansSerif"}}, {"selector": "edge[weight>1]", "css": {"opacity": 1.0, "text-opacity": 1.0, "font-size": 12, "font-weight": "normal", "width": "data(width)", "color": "rgb(0,0,0)", "line-color": "rgb(51,51,51)", "source-arrow-color": "rgb(0,0,0)", "target-arrow-color": "rgb(51,51,51)", "curve-style": "bezier", "line-style": "solid", "target-arrow-shape": "triangle", "source-arrow-shape": "none", "font-family": "SansSerif"}}, {"selector": "edge[interactionType *= 'inhibit'],edge[interactionType *= 'deactivat']", "css": {"text-opacity": 1.0, "font-size": 12, "font-weight": "normal", "font-family": "SansSerif", "curve-style": "bezier", "source-arrow-shape": "none", "target-arrow-shape": "tee", "color": "rgb(0,0,0)", "line-color": "rgb(51,51,51)", "source-arrow-color": "rgb(0,0,0)", "target-arrow-color": "rgb(51,51,51)", "line-style": "solid"}}, {"selector": "edge:selected", "css": {"line-color": "rgb(255,0,0)", "label": "data(interactionType)"}}, {"selector": "edge[weight>1]:selected", "css": {"background-color": "yellow", "line-color": "red", "label": "data(label)", "target-arrow-color": "black", "source-arrow-color": "black", "width": "data(width)"}}]).update();
    
                    // Add tippy for each node
                    cy.nodes().forEach(function (n) {
                        n.data()['tip'] = makeTippy(n, n.data('id'));
                    });
    
                    // hide tippy text on click
                    cy.on('tap', 'node', function (evt) {
                        var ele = evt.target;
                        if (ele.data()['tip']['state']['visible']) {
                            ele.data()['tip'].hide();
                        } else {
                            ele.data()['tip'].show();
                        }
                    });
    
                    // put the png data in an img tag
                    let downloadButton = document.createElement("BUTTON");
                    downloadButton.id = 'dbutton';
                    downloadButton.innerHTML = '<i class="fa fa-download" aria-hidden="true"></i>';
                    downloadButton.addEventListener('click', function () {
                        let element = document.createElement('a');
                        element.setAttribute('href', cy.png({scale: 3}));
                        element.setAttribute('download', 'graph.png');
                        element.style.display = 'none';
                        document.body.appendChild(element);
                        element.click();
                        document.body.removeChild(element);
                    });
                    let p = document.getElementById('cyb63e02fc-a04e-4038-add4-84a305eef936');
                    p.parentElement.append(downloadButton);
    
                });
    
    
        </script>
    </head>
    
    <body>
    <div id="cyb63e02fc-a04e-4038-add4-84a305eef936"></div>
    <!-- When only #uuid div is placed on this page,
    the height of output-box on ipynb will be 0px.
    One line below will prevent that. -->
    <div id="dummy"
         style="width:100px;height:700px"></div>
    </body>
    
    </html>


.. code:: ipython3

    draw_cyjs(mol_net, add_parent=False)



.. raw:: html

    <!DOCTYPE html>
    <html>
    <head>
        <meta charset=utf-8/>
        <style type="text/css">
            body {
                font: 14px helvetica neue, helvetica, arial, sans-serif;
            }
    
            #cyb441c7f6-d1e1-4e7f-96c6-add976715bf6{
                height: 700px;
                width: 90%;
                border: 5px solid black;
                box-sizing: border-box;
                position: relative;
                top: 5;
                margin-bottom: -700px;
                background: white;
            }
    
        </style>
    
        <script>
    
    
            requirejs.config({
    
                paths: {
                    'popper': 'https://unpkg.com/popper.js@1.14.1/dist/umd/popper',
                    'tippy': 'https://cdnjs.cloudflare.com/ajax/libs/tippy.js/2.3.0/tippy.min',
                    'cytoscape': 'https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.2.10/cytoscape',
                    'cytoscape-popper': 'https://cdn.rawgit.com/cytoscape/cytoscape.js-popper/3ad50859/cytoscape-popper',
                    'jquery': 'https://cdnjs.cloudflare.com/ajax/libs/jquery/2.2.4/jquery.min',
                    'qtip2': 'https://cdnjs.cloudflare.com/ajax/libs/qtip2/2.2.0/basic/jquery.qtip.min',
                    'dagre': 'https://cdn.rawgit.com/cpettitt/dagre/v0.7.4/dist/dagre.min',
                    'cytoscape-dagre': 'https://cdn.rawgit.com/cytoscape/cytoscape.js-dagre/1.5.0/cytoscape-dagre',
                    'cytoscape-cose-bilkent': 'https://cdn.rawgit.com/cytoscape/cytoscape.js-cose-bilkent/1.6.1/cytoscape-cose-bilkent'
                },
                shim: {
                    'cytoscape-popper': {
                        deps: ['cytoscape', 'popper']
                    },
                    'cytoscape-dagre': {
                        deps: ['cytoscape', 'dagre']
                    }
                },
                map: {
                    '*': {
                        'popper.js': 'popper',
                        'webcola': 'cola'
                    }
                }
    
            });
    
    
            require(['cytoscape', 'cytoscape-popper', 'popper', 'tippy', 'jquery',
                    'cytoscape-cose-bilkent', 'cytoscape-dagre', 'dagre'],
                function (cytoscape, cypopper, popper, tippy, jquery, regCose,
                          cydag, dagre) {
                    console.log('Loading Cytoscape.js Module...');
                    window['popper'] = popper;
                    window['tippy'] = tippy;
                    window['cytoscape'] = cytoscape;
                    cypopper(cytoscape);
                    regCose(cytoscape);
                    cydag(cytoscape, dagre);
    
                    function makeTippy(target, text) {
                        return tippy(target.popperRef(),
                            {
                                html: add_tip(text),
                                trigger: 'manual',
                                arrow: true,
                                placement: 'top',
                                hideOnClick: false,
                                interactive: true,
                                multiple: true,
                                sticky: true
                            }).tooltips[0];
                    }
    
                    function add_tip(text) {
                        let div = document.createElement('div');
                        div.innerHTML = "<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=" + text + "' target='_blank'>Gene Card</a>";
                        return div;
                    }
    
                    let cy = window.cy = cytoscape({
                        container: $('#cyb441c7f6-d1e1-4e7f-96c6-add976715bf6'),
                        elements: {
                            nodes: [{"data": {"termName": "g2/m checkpoints,transcriptional regulation by tp53", "terms": "g2/m checkpoints,transcriptional regulation by tp53", "color": "white", "id": "RAD9A", "name": "RAD9A"}}, {"data": {"termName": "g2/m checkpoints,transcriptional regulation by tp53", "terms": "g2/m checkpoints,transcriptional regulation by tp53", "color": "white", "id": "TP53", "name": "TP53"}}, {"data": {"termName": "g2/m checkpoints,transcriptional regulation by tp53", "terms": "g2/m checkpoints,transcriptional regulation by tp53", "color": "white", "id": "NBN", "name": "NBN"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "color": "white", "id": "MCM2", "name": "MCM2"}}, {"data": {"termName": "g2/m checkpoints,transcriptional regulation by tp53", "terms": "g2/m checkpoints,transcriptional regulation by tp53", "color": "white", "id": "TOPBP1", "name": "TOPBP1"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "color": "white", "id": "PSMD4", "name": "PSMD4"}}, {"data": {"termName": "apoptotic cleavage of cellular proteins", "terms": "apoptotic cleavage of cellular proteins", "color": "white", "id": "APC", "name": "APC"}}, {"data": {"termName": "apoptotic cleavage of cellular proteins", "terms": "apoptotic cleavage of cellular proteins", "color": "white", "id": "CTNNB1", "name": "CTNNB1"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "color": "white", "id": "PSMB3", "name": "PSMB3"}}, {"data": {"termName": "apoptotic cleavage of cellular proteins", "terms": "apoptotic cleavage of cellular proteins", "color": "white", "id": "PTK2", "name": "PTK2"}}, {"data": {"termName": "apoptotic cleavage of cellular proteins", "terms": "apoptotic cleavage of cellular proteins", "color": "white", "id": "SPTAN1", "name": "SPTAN1"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "color": "white", "id": "PSMF1", "name": "PSMF1"}}, {"data": {"termName": "apoptotic cleavage of cellular proteins", "terms": "apoptotic cleavage of cellular proteins", "color": "white", "id": "TJP1", "name": "TJP1"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "color": "white", "id": "PSMB10", "name": "PSMB10"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "color": "white", "id": "PSMA5", "name": "PSMA5"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "color": "white", "id": "CDC7", "name": "CDC7"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "color": "white", "id": "MCM7", "name": "MCM7"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "color": "white", "id": "MCM6", "name": "MCM6"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "color": "white", "id": "PSMA3", "name": "PSMA3"}}, {"data": {"termName": "g2/m checkpoints,transcriptional regulation by tp53", "terms": "g2/m checkpoints,transcriptional regulation by tp53", "color": "white", "id": "YWHAZ", "name": "YWHAZ"}}, {"data": {"termName": "g2/m checkpoints,transcriptional regulation by tp53", "terms": "g2/m checkpoints,transcriptional regulation by tp53", "color": "white", "id": "ATR", "name": "ATR"}}, {"data": {"termName": "g2/m checkpoints,transcriptional regulation by tp53", "terms": "g2/m checkpoints,transcriptional regulation by tp53", "color": "white", "id": "SFN", "name": "SFN"}}, {"data": {"termName": "apoptotic cleavage of cellular proteins", "terms": "apoptotic cleavage of cellular proteins", "color": "white", "id": "LMNB1", "name": "LMNB1"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "color": "white", "id": "UIMC1", "name": "UIMC1"}}, {"data": {"termName": "g2/m checkpoints,transcriptional regulation by tp53", "terms": "g2/m checkpoints,transcriptional regulation by tp53", "color": "white", "id": "BRCA1", "name": "BRCA1"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "color": "white", "id": "PSMD1", "name": "PSMD1"}}, {"data": {"termName": "g2/m checkpoints", "terms": "g2/m checkpoints", "color": "white", "id": "PSMB7", "name": "PSMB7"}}, {"data": {"termName": "g2/m checkpoints,transcriptional regulation by tp53", "terms": "g2/m checkpoints,transcriptional regulation by tp53", "color": "white", "id": "RAD50", "name": "RAD50"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "color": "white", "id": "MAPK14", "name": "MAPK14"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "color": "white", "id": "JUN", "name": "JUN"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "color": "white", "id": "ATF2", "name": "ATF2"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "color": "white", "id": "TSC2", "name": "TSC2"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "color": "white", "id": "TP53BP2", "name": "TP53BP2"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "color": "white", "id": "RRM2B", "name": "RRM2B"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "color": "white", "id": "TXN", "name": "TXN"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "color": "white", "id": "TAF1L", "name": "TAF1L"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "color": "white", "id": "GPI", "name": "GPI"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "color": "white", "id": "G6PD", "name": "G6PD"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "color": "white", "id": "RFFL", "name": "RFFL"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "color": "white", "id": "HDAC1", "name": "HDAC1"}}, {"data": {"termName": "g2/m checkpoints,transcriptional regulation by tp53", "terms": "g2/m checkpoints,transcriptional regulation by tp53", "color": "white", "id": "YWHAE", "name": "YWHAE"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "color": "white", "id": "TAF9", "name": "TAF9"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "color": "white", "id": "FAS", "name": "FAS"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "color": "white", "id": "BNIP3L", "name": "BNIP3L"}}, {"data": {"termName": "apoptotic cleavage of cellular proteins", "terms": "apoptotic cleavage of cellular proteins", "color": "white", "id": "VIM", "name": "VIM"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "color": "white", "id": "MAPKAP1", "name": "MAPKAP1"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "color": "white", "id": "MTOR", "name": "MTOR"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "color": "white", "id": "LAMTOR1", "name": "LAMTOR1"}}, {"data": {"termName": "g2/m checkpoints,transcriptional regulation by tp53", "terms": "g2/m checkpoints,transcriptional regulation by tp53", "color": "white", "id": "YWHAQ", "name": "YWHAQ"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "color": "white", "id": "SSRP1", "name": "SSRP1"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "color": "white", "id": "PPP2CA", "name": "PPP2CA"}}, {"data": {"termName": "g2/m checkpoints,transcriptional regulation by tp53", "terms": "g2/m checkpoints,transcriptional regulation by tp53", "color": "white", "id": "YWHAH", "name": "YWHAH"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "color": "white", "id": "RBL1", "name": "RBL1"}}, {"data": {"termName": "transcriptional regulation by tp53", "terms": "transcriptional regulation by tp53", "color": "white", "id": "SUPT16H", "name": "SUPT16H"}}],
                            edges: [{"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "RAD9A", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate|binding", "width": 10.0, "source": "RAD9A", "target": "NBN"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "RAD9A", "target": "MCM2"}}, {"data": {"databaseSource": "SIGNOR", "interactionType": "activate|binding", "width": 10.0, "source": "RAD9A", "target": "TOPBP1"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "TP53", "target": "SFN"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "TP53", "target": "LMNB1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "TP53", "target": "CTNNB1"}}, {"data": {"databaseSource": "KEGG|ReactomeFI|SIGNOR", "interactionType": "activate|expression|indirect", "width": 10.0, "source": "TP53", "target": "FAS"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "TP53", "target": "TSC2"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "TP53", "target": "RRM2B"}}, {"data": {"databaseSource": "KEGG", "interactionType": "expression", "width": 10.0, "source": "TP53", "target": "BNIP3L"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "NBN", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "TOPBP1", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate|binding", "width": 10.0, "source": "TOPBP1", "target": "NBN"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "TOPBP1", "target": "MCM2"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMD4", "target": "APC"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMD4", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMD4", "target": "CTNNB1"}}, {"data": {"databaseSource": "KEGG|SIGNOR", "interactionType": "binding|inhibit", "width": 10.0, "source": "APC", "target": "CTNNB1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate|binding", "width": 10.0, "source": "CTNNB1", "target": "TJP1"}}, {"data": {"databaseSource": "KEGG", "interactionType": "expression", "width": 10.0, "source": "CTNNB1", "target": "JUN"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate|inhibit", "width": 10.0, "source": "CTNNB1", "target": "HDAC1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMB3", "target": "APC"}}, {"data": {"databaseSource": "BIOGRID|ReactomeFI", "interactionType": "catalyze|cleavage", "pubmedId": "20080206", "width": 10.0, "source": "PSMB3", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMB3", "target": "CTNNB1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "binding|catalyze", "width": 10.0, "source": "PTK2", "target": "SPTAN1"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "PTK2", "target": "MAPK14"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMF1", "target": "APC"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMF1", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMF1", "target": "CTNNB1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMB10", "target": "APC"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMB10", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMB10", "target": "CTNNB1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMA5", "target": "APC"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMA5", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMA5", "target": "CTNNB1"}}, {"data": {"databaseSource": "BIOGRID|ReactomeFI|SIGNOR", "interactionType": "activate|binding|catalyze|phosphorylate", "pubmedId": "10846177", "width": 10.0, "source": "CDC7", "target": "MCM7"}}, {"data": {"databaseSource": "BIOGRID|ReactomeFI", "interactionType": "activate|binding|catalyze|phosphorylate", "pubmedId": "10846177", "width": 10.0, "source": "CDC7", "target": "MCM6"}}, {"data": {"databaseSource": "BIOGRID|ReactomeFI|SIGNOR", "interactionType": "activate|binding|catalyze|phosphorylate", "pubmedId": "10846177", "width": 10.0, "source": "CDC7", "target": "MCM2"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMA3", "target": "APC"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMA3", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMA3", "target": "CTNNB1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "YWHAZ", "target": "RAD9A"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "YWHAZ", "target": "ATR"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "YWHAZ", "target": "TOPBP1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "binding|inhibit", "width": 10.0, "source": "YWHAZ", "target": "TSC2"}}, {"data": {"databaseSource": "BIOGRID|KEGG|ReactomeFI|SIGNOR", "interactionType": "activate|binding|phosphorylate", "pubmedId": "11278964", "width": 10.0, "source": "ATR", "target": "BRCA1"}}, {"data": {"databaseSource": "BIOGRID|KEGG|ReactomeFI|SIGNOR", "interactionType": "activate|catalyze|phosphorylate", "pubmedId": "10608806", "width": 10.0, "source": "ATR", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate|binding", "width": 10.0, "source": "ATR", "target": "RAD9A"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "ATR", "target": "CDC7"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "ATR", "target": "MCM7"}}, {"data": {"databaseSource": "ReactomeFI|SIGNOR", "interactionType": "activate|catalyze|phosphorylate", "width": 10.0, "source": "ATR", "target": "MCM6"}}, {"data": {"databaseSource": "BIOGRID|ReactomeFI|SIGNOR", "interactionType": "activate|catalyze|phosphorylate", "pubmedId": "15210935", "width": 10.0, "source": "ATR", "target": "MCM2"}}, {"data": {"databaseSource": "ReactomeFI|SIGNOR", "interactionType": "activate|binding|phosphorylate", "width": 10.0, "source": "ATR", "target": "NBN"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "binding|inhibit", "width": 10.0, "source": "SFN", "target": "TSC2"}}, {"data": {"databaseSource": "SIGNOR", "interactionType": "activate|binding", "width": 10.0, "source": "UIMC1", "target": "BRCA1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "binding|catalyze", "width": 10.0, "source": "BRCA1", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "BRCA1", "target": "PPP2CA"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "BRCA1", "target": "MAPK14"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "binding|expression", "width": 10.0, "source": "BRCA1", "target": "ATF2"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMD1", "target": "APC"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMD1", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMD1", "target": "CTNNB1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMB7", "target": "APC"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMB7", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "PSMB7", "target": "CTNNB1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "RAD50", "target": "TP53"}}, {"data": {"databaseSource": "SIGNOR", "interactionType": "activate|binding", "width": 10.0, "source": "RAD50", "target": "BRCA1"}}, {"data": {"databaseSource": "BIOGRID|KEGG", "interactionType": "activate|phosphorylate", "pubmedId": "25483191", "width": 10.0, "source": "MAPK14", "target": "JUN"}}, {"data": {"databaseSource": "BIOGRID|KEGG|ReactomeFI|SIGNOR", "interactionType": "activate|catalyze|phosphorylate", "pubmedId": "10581258", "width": 10.0, "source": "MAPK14", "target": "TP53"}}, {"data": {"databaseSource": "BIOGRID|KEGG|ReactomeFI|SIGNOR", "interactionType": "activate|catalyze|phosphorylate", "pubmedId": "10581258", "width": 10.0, "source": "MAPK14", "target": "ATF2"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "MAPK14", "target": "TSC2"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "MAPK14", "target": "YWHAZ"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "JUN", "target": "FAS"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "binding|expression", "width": 10.0, "source": "JUN", "target": "TP53"}}, {"data": {"databaseSource": "KEGG|ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "JUN", "target": "VIM"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "JUN", "target": "BRCA1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "JUN", "target": "HDAC1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "expression", "width": 10.0, "source": "JUN", "target": "TXN"}}, {"data": {"databaseSource": "KEGG", "interactionType": "expression", "width": 10.0, "source": "ATF2", "target": "VIM"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "binding|expression|inhibit", "width": 10.0, "source": "ATF2", "target": "JUN"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "TSC2", "target": "MAPKAP1"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "TSC2", "target": "MTOR"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "TSC2", "target": "LAMTOR1"}}, {"data": {"databaseSource": "SIGNOR", "interactionType": "activate|binding", "width": 10.0, "source": "TP53BP2", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "RRM2B", "target": "TXN"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "TAF1L", "target": "TP53"}}, {"data": {"databaseSource": "KEGG", "interactionType": "chemical", "width": 10.0, "source": "GPI", "target": "G6PD"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "RFFL", "target": "TP53"}}, {"data": {"pubmedId": "11099047", "databaseSource": "BIOGRID", "interactionType": "deacetylate", "width": 10.0, "source": "HDAC1", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "binding|inhibit", "width": 10.0, "source": "YWHAE", "target": "TSC2"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze|inhibit", "width": 10.0, "source": "TAF9", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "binding|inhibit", "width": 10.0, "source": "YWHAQ", "target": "TSC2"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "SSRP1", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate", "width": 10.0, "source": "PPP2CA", "target": "SFN"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate|inhibit", "width": 10.0, "source": "PPP2CA", "target": "YWHAE"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate|inhibit", "width": 10.0, "source": "PPP2CA", "target": "YWHAQ"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate|inhibit", "width": 10.0, "source": "PPP2CA", "target": "YWHAH"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "activate|inhibit", "width": 10.0, "source": "PPP2CA", "target": "YWHAZ"}}, {"data": {"databaseSource": "SIGNOR", "interactionType": "activate|dephosphorylate", "width": 10.0, "source": "PPP2CA", "target": "RBL1"}}, {"data": {"databaseSource": "SIGNOR", "interactionType": "activate|dephosphorylate", "width": 10.0, "source": "PPP2CA", "target": "CTNNB1"}}, {"data": {"databaseSource": "SIGNOR", "interactionType": "activate|dephosphorylate", "width": 10.0, "source": "PPP2CA", "target": "TP53"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "binding|inhibit", "width": 10.0, "source": "YWHAH", "target": "TSC2"}}, {"data": {"databaseSource": "ReactomeFI", "interactionType": "catalyze", "width": 10.0, "source": "SUPT16H", "target": "TP53"}}]
                        },
                        boxSelectionEnabled: true,
                        wheelSensitivity: .25,
    
                    });
    
                    // apply layout
                    cy.layout({"name": "cose-bilkent", "spacingFactor": 1.5, "animate": false}).run();
    
                    // add style
                    cy.style().fromJson([{"selector": "node", "css": {"text-opacity": 1.0, "background-opacity": 1.0, "font-weight": "normal", "color": "rgb(0,153,234)", "border-width": 3.0, "border-color": "rgb(51,51,51)", "font-size": 9, "height": 50, "width": 50, "label": "data(name)", "text-wrap": "wrap", "text-max-width": 50, "shape": "ellipse", "background-color": "data(color)", "text-halign": "center", "font-family": "SansSerif", "text-valign": "center", "border-opacity": 1.0}}, {"selector": "node[speciesType = 'compound']", "css": {"text-opacity": 1.0, "background-opacity": 1.0, "font-weight": "normal", "color": "rgb(0,153,234)", "border-width": 3.0, "border-color": "rgb(51,51,51)", "font-size": 9, "height": 35, "width": 35, "label": "data(chemName)", "text-wrap": "wrap", "text-max-width": 95, "shape": "ellipse", "background-color": "rgb(255,255,255)", "text-halign": "center", "font-family": "SansSerif", "text-valign": "center", "border-opacity": 1.0}}, {"selector": "$node > node", "css": {"font-size": 20, "shape": "ellipse", "padding-right": "10px", "padding-bottom": "10px", "padding-top": "10px", "text-valign": "top", "text-halign": "center", "background-color": "#bbb", "padding-left": "10px"}}, {"selector": "node:selected", "css": {"background-color": "rgb(255,0,102)"}}, {"selector": "edge", "css": {"opacity": 1.0, "text-opacity": 1.0, "font-size": 12, "font-weight": "normal", "width": "2", "color": "rgb(0,0,0)", "line-color": "rgb(51,51,51)", "source-arrow-color": "rgb(0,0,0)", "target-arrow-color": "rgb(51,51,51)", "curve-style": "bezier", "line-style": "solid", "target-arrow-shape": "triangle", "source-arrow-shape": "none", "font-family": "SansSerif"}}, {"selector": "edge[weight>1]", "css": {"opacity": 1.0, "text-opacity": 1.0, "font-size": 12, "font-weight": "normal", "width": "data(width)", "color": "rgb(0,0,0)", "line-color": "rgb(51,51,51)", "source-arrow-color": "rgb(0,0,0)", "target-arrow-color": "rgb(51,51,51)", "curve-style": "bezier", "line-style": "solid", "target-arrow-shape": "triangle", "source-arrow-shape": "none", "font-family": "SansSerif"}}, {"selector": "edge[interactionType *= 'inhibit'],edge[interactionType *= 'deactivat']", "css": {"text-opacity": 1.0, "font-size": 12, "font-weight": "normal", "font-family": "SansSerif", "curve-style": "bezier", "source-arrow-shape": "none", "target-arrow-shape": "tee", "color": "rgb(0,0,0)", "line-color": "rgb(51,51,51)", "source-arrow-color": "rgb(0,0,0)", "target-arrow-color": "rgb(51,51,51)", "line-style": "solid"}}, {"selector": "edge:selected", "css": {"line-color": "rgb(255,0,0)", "label": "data(interactionType)"}}, {"selector": "edge[weight>1]:selected", "css": {"background-color": "yellow", "line-color": "red", "label": "data(label)", "target-arrow-color": "black", "source-arrow-color": "black", "width": "data(width)"}}]).update();
    
                    // Add tippy for each node
                    cy.nodes().forEach(function (n) {
                        n.data()['tip'] = makeTippy(n, n.data('id'));
                    });
    
                    // hide tippy text on click
                    cy.on('tap', 'node', function (evt) {
                        var ele = evt.target;
                        if (ele.data()['tip']['state']['visible']) {
                            ele.data()['tip'].hide();
                        } else {
                            ele.data()['tip'].show();
                        }
                    });
    
                    // put the png data in an img tag
                    let downloadButton = document.createElement("BUTTON");
                    downloadButton.id = 'dbutton';
                    downloadButton.innerHTML = '<i class="fa fa-download" aria-hidden="true"></i>';
                    downloadButton.addEventListener('click', function () {
                        let element = document.createElement('a');
                        element.setAttribute('href', cy.png({scale: 3}));
                        element.setAttribute('download', 'graph.png');
                        element.style.display = 'none';
                        document.body.appendChild(element);
                        element.click();
                        document.body.removeChild(element);
                    });
                    let p = document.getElementById('cyb441c7f6-d1e1-4e7f-96c6-add976715bf6');
                    p.parentElement.append(downloadButton);
    
                });
    
    
        </script>
    </head>
    
    <body>
    <div id="cyb441c7f6-d1e1-4e7f-96c6-add976715bf6"></div>
    <!-- When only #uuid div is placed on this page,
    the height of output-box on ipynb will be 0px.
    One line below will prevent that. -->
    <div id="dummy"
         style="width:100px;height:700px"></div>
    </body>
    
    </html>


Finally, we can bring it full circle and subset our experimental data to
visualize the nodes in the networks measured values over time.

.. code:: ipython3

    exp_data.ph_silac.heatmap(
        mol_net.nodes,
        subset_index='identifier',
        index='label',
        cluster_row=True,
        sort_row='index',
        min_sig=2,
        linewidths=0.01,
        figsize=(6,12),
    );
        



.. image:: Tutorial_files/Tutorial_136_0.png


Custom workflows
''''''''''''''''

Here we presented some examples of how to use MAGINE. The strength is in
MAGINES ability to explore data, enrichment, and networks within a
single space, allowing back and forth exploration. If you have any
suggestions or would like to contribute workflows or pipelines, please
share on our github issues.
