"""
Manually constructed mapping to most common description.

"""

edge_standards = {
    'activation'                        : 'activate',
    'activator'                         : 'activate',
    'potentiator'                       : 'activate',
    'oxidation'                         : 'oxidate',
    'palmitoylation'                    : 'palmitoylate',
    'inducer'                           : 'expression',
    'stimulator'                        : 'expression',
    'suppressor'                        : 'repression',
    'transcriptional activation'        : 'expression',
    'transcriptional inhibition'        : 'repression',
    'blocker'                           : 'inhibit',
    'inhibitor'                         : 'inhibit',
    'inhibition'                        : 'inhibit',
    'inhibitor, competitive'            : 'inhibit',

    'proteolytic processing'            : 'cleavage',
    'cleavage'                          : 'cleavage',
    'stabilization'                     : 'stabilization',
    'translation regulation'            : 'translation',

    'hydroxylation'                     : 'hydroxylation',
    'lipidation'                        : 'lipidation',
    'guanine nucleotide exchange factor': 'gtpase-activate',
    'gtpase-activating protein'         : 'gtpase-activate',

    # binding
    'binding'                           : 'binding',
    'binding/association'               : 'binding',
    'binder'                            : 'binding',
    'complex'                           : 'binding',
    'dissociation'                      : 'binding',

    # indirect/missing
    'indirect effect'                   : 'indirect',
    'missing interaction'               : 'indirect',

    'state change'                      : 'stateChange',
    'relocalization'                    : 'relocalization',

    's-nitrosylation'                   : 's-nitrosylation',
    'acetylation'                       : 'acetylate',
    'tyrosination'                      : 'tyrosinate',
    'ubiquitination'                    : 'ubiquitinate',
    'methylation'                       : 'methylate',
    'glycosylation'                     : 'glycosylate',
    'sumoylation'                       : 'sumoylate',
    'ribosylation'                      : 'ribosylate',
    'neddylation'                       : 'neddylate',
    'desumoylation'                     : 'desumoylate',
    'deneddylation'                     : 'deneddylate',
    'demethylation'                     : 'demethylate',
    'deacetylation'                     : 'deacetylate',
    'trimethylation'                    : 'trimethylate',
    'desensitize the target'            : 'inhibit',
    'deubiquitination'                  : 'deubiquitinate',
    'nedd(rub1)ylation'                 : 'neddy(rub1)late',

    'dephosphorylation': 'dephosphorylate',
    'phosphorylation': 'phosphorylate',

    'negative modulator': 'inhibit',
    'inhibitory allosteric modulator': 'allosteric|inhibit',
    'allosteric modulator': 'allosteric|modulate',
    'positive allosteric modulator': 'activate|allosteric',
    'positive modulator': 'activate',


    # chemical related
    'compound': 'chemical',
    'product of': 'chemical',
    'ligand': 'chemical',
    'cofactor': 'chemical',
    'multitarget': 'chemical',
    'small molecule catalysis': 'activate|catalyst|chemical',
    'chemical activation': 'activate|chemical',
    'chemical inhibition': 'chemical|inhibit',
    'partial agonist': 'activate|chemical',
    'inverse agonist': 'activate|chemical',
    'agonist': 'activate|chemical',
    'antagonist': 'inhibit|chemical',
    'partial antagonist': 'inhibit|chemical',
}
