# Copyright 2018 by Lucas Davey.  All rights reserved.

# TODO: specify existence in Bio.Application docstring

"""Command line wrapper for the antibody annotation program Anarci."""

from Bio.Application import _Option, _Switch, AbstractCommandline

class AnarciCommandline(AbstractCommandline):
    """Commandline object for Anarci.
    
    http://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/ANARCI.php
    
    Last checked against version 1.3
    
    Attributes:

     - sequence                 sequence or an input fasta file
     - outfile                  The output file to use. Default is stdout
     - scheme                   Which numbering scheme should be used
     - restrict                 Restrict ANARCI to only recognise certain 
                                types of receptor chains.
     - outfile_hits             Output file for domain hit tables for each sequence
     - use_species              Use a specific species in the germline assignment
     - bit_score_threshold      Change the bit score threshold used to 
                                confirm an alignment should be used
     - help                     show this help message and exit
     - assign_germline          Assign the v and j germlines to the sequence
     - hmmer_path               The path to the directory containing hmmer programs
     - csv                      Write the output in csv format
   
    Example:

    >>> from Bio.Anarci import AnarciCommandline
    >>> cline = AnarciCommandline(sequence='ncbi_query_10.fa', scheme='kabat')
    >>> cline.outfile = 'anarci_outfile.txt'
    >>> cline.ht = 'anarci_outfile_hits.txt'
    >>> cline.bogusparameter = 1995     # Invalid parameter
    Traceback (most recent call last):
        ...
    ValueError: Option name bogusparameter was not found.
    >>> print(cline)
    ANARCI --sequence=ncbi_query_10.fa --outfile=anarci_outfile.txt --scheme=kabat --outfile_hits=anarci_outfile_hits.txt

    You would typically run the command line with cline() or via
    the Python subprocess module, as described in the Biopython tutorial.
    """
    def __init__(self, cmd='ANARCI', **kwargs):
        self.program_name = 'ANARCI'
        self.parameters = [
        _Option(['--sequence', 'sequence', 'i'],
                'A sequence or an input fasta file',
                filename=True),
        _Option(['--outfile', 'outfile', 'o'],
                'The output file to use. Default is stdout',
                filename=True),
        _Option(['--scheme', 'scheme', 's'],
                'Which numbering scheme should be used. i, k, c, m, w \
                and a are shorthand for IMGT, Kabat, Chothia, Martin \
                (Extended Chothia), Wolfguy and Aho respectively. \
                Default IMGT',
                checker_function=lambda x: x in ['kabat','aho','wolfguy','imgt',
                                                'a','c','chothia','i',
                                                'k','m','w','martin']),
        _Option(['--restrict', 'restrict', 'r'],
                'Restrict ANARCI to only recognise certain types of \
                receptor chains.',
                checker_function=lambda x: x in ['ig', 'tr', 'heavy', 'light', 'H', 'K', 'L', 'A', 'B']),
        _Option(['--outfile_hits', 'outfile_hits', 'ht'],
                'Output file for domain hit tables for each sequence. \
                Otherwise not output.',
                filename=True),
        _Option(['--ncpu', 'ncpu', 'p'],
                'number of parallel CPU workers to use for \
                multithreads. Default is to use all that can be found.',
                checker_function=lambda x: isinstance(x, int)),
        _Option(['--use_species', 'use_species'],
                'Use a specific species in the germline assignment.',
                checker_function=lambda x: x in ['alpaca', 'rabbit', 'rhesus',
                                                'pig', 'rat', 'human', 'mouse']),
        _Option(['--bit_score_threshold', 'bit_score_threshold'],
                'Change the bit score threshold used to confirm an \
                alignment should be used.',
                checker_function=lambda x: isinstance(x, int)),
        _Switch(['--help', 'h', 'help'],
                'show this help message and exit'),
        _Switch(['--assign_germline', 'assign_germline'],
                'Assign the v and j germlines to the sequence. The most \
                sequence identical germline is assigned.'),
        _Switch(['--hmmerpath', 'hmmerpath', 'hp'],
                'The path to the directory containing hmmer programs. \
                (including hmmscan)'), 
        _Switch(['--csv', 'csv'], 
                'Write the output in csv format. Outfile must be \
                 specified. A csv file is written for each chain type \
                 <outfile>_<chain_type>.csv. Kappa and lambda are \
                 considered together.')      
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
        
# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
