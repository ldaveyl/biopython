# Copyright 2018 by Lucas Davey.  All rights reserved.

"""Bio.SearchIO support for Anarci output formats.

This module adds support for parsing ANARCI outputs. ANARCI is a tool
for numbering amino-acid sequences of antibody and T-cell receptor variable domains.
ANARCI aligns a given sequence to a database of Hidden Markov Models 
that describe the germline sequences of antibody and TCR domain types.

Bio.SearchIO.AnarciIO was tested on ANARCI 1.3

More information on ANARCI is available on its home page at
http://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/ANARCI.php

Supported Formats
=================
  - Numbered - 'anarci-num'   - parsing, indexing
  - Hits     - 'anarci-hits'  - parsing, indexing
  
  
anarci-num
=========
The numbered text output is the standard output of ANARCI. 
This format contains only the most significant hit and, in contrast to anarci-hits,
also the annotation.

    >>> from Bio import SearchIO
    >>> anarci_qresult = next(SearchIO.parse('Anarci/num_anarci_001.txt', 'anarci-num'))
    >>> anarci_qresult
    QueryResult(id='AXG50451.1', 1 hits)

Anarci-num output always has only one Hit, HSP and HSPfragment. 
HSPFragment's query object contains a list of SeqFeatures which contain the antibody
annotation. This was done to make writing to different file formats intuitive.

There are currently two types of SeqFeature objects: scheme and 'ANARCI_hit'.
(scheme being 'kabat', 'chothia', etc.)
SeqFeature objects of type scheme contain the qualifiers: 'residue', 'label'
and anarci_scheme_num. ANARCI_hit contains the entire sequence and some other information.
Every SeqFeature object also contains a FeatureLocation which describes the position before
and after every residue.
    
    >>> # first hit, first hsp, first fragment 
    >>> frag = anarci_qresult[0][0][0]
    >>> # first SeqFeature
    >>> frag.query.features[0]
    SeqFeature(FeatureLocation(ExactPosition(20), ExactPosition(21)), type='kabat')
    >>> # last SeqFeature
    >>> frag.query.features[-1]
    SeqFeature(FeatureLocation(ExactPosition(20), ExactPosition(116)), type='ANARCI_hit')
    
If you are just interested in seeing the annotation, HSPFragment's query object
also provides a letter_annotations attribute.

    >>> frag.query.letter_annotations['kabat'][:10]
    ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10']
    
As of now there is no way to include the hit sequence as this is not included in the ANARCI output.
Therefore HSPFragment doesn't have a hit object and hit_start and hit_end are set to 0 to
prevent a warning.

The numbered format provides the following attributes for each SearchIO object:

+---------------------+-------------------------+---------------------------------+
| Object              | Attribute               | Value                           |
+=====================+=========================+=================================+
| QueryResult         | target                  | 'NCBI protein database'         |
|                     +-------------------------+---------------------------------+
|                     | program                 | 'ANARCI'                        |
|                     +-------------------------+---------------------------------+
|                     | version                 | '1.3'                           |
|                     +-------------------------+---------------------------------+
|                     | id                      | query sequence id               |
|                     +-------------------------+---------------------------------+
|                     | description             | query sequence description      |
|                     +-------------------------+---------------------------------+
|                     | scheme                  | numbering scheme                |
|                     +-------------------------+---------------------------------+
|                     | chain_type              | light (L) or heavy chain (H)    |
|                     +-------------------------+---------------------------------+
|                     | species                 | species of the hit              |
+---------------------+-------------------------+---------------------------------+
| Hit                 | id                      | hit id                          |
|                     +-------------------------+---------------------------------+
|                     | description             | hit description                 |
+---------------------+-------------------------+---------------------------------+
| HSP                 | evalue                  | The e-value of the alignment    |
|                     +-------------------------+---------------------------------+               
|                     | bitscore                | The bit-score of the alignment  |
+---------------------+-------------------------+---------------------------------+
| HSPFragment         | query                   | query sequence                  |
|                     +-------------------------+---------------------------------+ 
|                     | query_start             | index of the first residue      |
|                     +-------------------------+---------------------------------+ 
|                     | query_end               | index of the last residue       |
|                     +-------------------------+---------------------------------+              
|                     | alphabet                | 'ProteinAlphabet()'             |
|                     +-------------------------+---------------------------------+              
|                     | letter_annotations      | antibody annotation             |
+---------------------+-------------------------+---------------------------------+ 

Note that this table only includes the attributres explicitly defined. 
You can still use the the default QueryResult, Hit, HSP and HSPFragment attributes.

If the file contains germline information (--assign_germline flag) 
the following attributes may be available as well:

+---------------------+-------------------------+---------------------------------+
| Object              | Attribute               | Value                           |
+=====================+=========================+=================================+
| QueryResult         | species_germline        | species of the germline         |
|                     +-------------------------+---------------------------------+
|                     | v_gene                  | id of germline v-region         |
|                     +-------------------------+---------------------------------+
|                     | v_identity              | sequence identity of v-region   |
|                     +-------------------------+---------------------------------+
|                     | j_gene                  | id of germline j-region         |
|                     +-------------------------+---------------------------------+
|                     | j_identity              | sequence identity of j-region   |
+---------------------+-------------------------+---------------------------------+
| Hit                 | id                      | hit id                          |
|                     +-------------------------+---------------------------------+
|                     | description             | hit description                 |
|                     +-------------------------+---------------------------------+
|                     | query_description       | query description               |
+---------------------+-------------------------+---------------------------------+
| HSP                 | evalue                  | The e-value of the alignment    |
|                     +-------------------------+---------------------------------+               
|                     | bitscore                | The bit-score of the alignment  |
+---------------------+-------------------------+---------------------------------+
| HSPFragment         | query                   | query sequence                  |
|                     +-------------------------+---------------------------------+ 
|                     | query_start             | index of the first residue      |
|                     +-------------------------+---------------------------------+ 
|                     | query_end               | index of the last residue       |
|                     +-------------------------+---------------------------------+              
|                     | alphabet                | 'ProteinAlphabet()'             |
|                     +-------------------------+---------------------------------+              
|                     | letter_annotations      | antibody annotation             |
+---------------------+-------------------------+---------------------------------+ 

anarci-hits
=========

The hits output is the ouput obtained by setting the -hit_file flag
in ANARCI. This format does not contain any annotation, but it does contain
additional hits and the query sequence

    >>> from Bio import SearchIO
    >>> anarci_qresult = next(SearchIO.parse('Anarci/hits_anarci_001.txt', 'anarci-hits'))
    >>> anarci_qresult
    QueryResult(id='AXG50451.1', 6 hits)
    >>> anarci_qresult.input_sequence[:40]
    'METDTLLLWVLLLWVPGSTGQVKLEESGPGLVNPSQSLSL'
    
Anarci-hits output can have more than one hit, but every hit has only one
HSP and HSPFragment.

The hits format provides the following attributes for each SearchIO object:

+---------------------+-------------------------+---------------------------------+
| Object              | Attribute               | Value                           |
+=====================+=========================+=================================+
| QueryResult         | target                  | 'NCBI protein database'         |
|                     +-------------------------+---------------------------------+
|                     | program                 | 'ANARCI'                        |
|                     +-------------------------+---------------------------------+
|                     | version                 | '1.3'                           |
|                     +-------------------------+---------------------------------+
|                     | id                      | query sequence id               |
|                     +-------------------------+---------------------------------+
|                     | scheme                  | numbering scheme                |
|                     +-------------------------+---------------------------------+
|                     | chain_type              | light (L) or heavy chain (H)    |
|                     +-------------------------+---------------------------------+
|                     | species                 | species of the hit              |
|                     +-------------------------+---------------------------------+
|                     | input_sequence          | input sequence                  |
+---------------------+-------------------------+---------------------------------+
| Hit                 | id                      | hit id                          |
|                     +-------------------------+---------------------------------+
|                     | description             | hit description                 |
|                     +-------------------------+---------------------------------+
|                     | query_description       | query description               |
+---------------------+-------------------------+---------------------------------+
| HSP                 | evalue                  | The e-value of the alignment    |
|                     +-------------------------+---------------------------------+               
|                     | bitscore                | The bit-score of the alignment  |
+---------------------+-------------------------+---------------------------------+
| HSPFragment         | query                   | query sequence                  |
|                     +-------------------------+---------------------------------+ 
|                     | query_start             | index of the first residue      |
|                     +-------------------------+---------------------------------+ 
|                     | query_end               | index of the last residue       |
|                     +-------------------------+---------------------------------+              
|                     | alphabet                | 'ProteinAlphabet()'             |
+---------------------+-------------------------+---------------------------------+ 

Note that this table only includes the attributes explicitly defined. 
You can still use the the default QueryResult, Hit, HSP and HSPFragment attributes.
"""

from .anarci_num import AnarciNumParser, AnarciNumIndexer
from .anarci_hits import AnarciHitsParser, AnarciHitsIndexer


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
