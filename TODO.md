#TODO

* magine_analysis.py

    * Restructure class to have setters for each data type.
    
    * Make pandas dataframe get lists per sample id, get rid of return_rna, return_proteomics
     
    * Have a sorting method to sort sample_ids, right now we assume they can be floats but if they are strings they should be able to defined the order. 

* Function requests
    * Add -log10 pvalue and log2 fold change to network node attributes.
    
    * For network subgraphs, collapse all unmeasured nodes into the edges.

