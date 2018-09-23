# Copyright 2018 by Lucas Davey.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Factory object capable of instantiating SearchIO.parse generator objects."""

import io
import os
import warnings

from .Applications import AnarciCommandline
from Bio import SearchIO, Seq, SeqRecord, SeqIO
from Bio.Alphabet import IUPAC

class AnarciFactory:
    """Factory object to run and parse Anarci.
    
    The concept of this class is to provide an easy way to sequentially run and parse 
    Biopython sequence objects into Anarci output. This gives the added benefit that
    it catches when Anarci crashes (e.g. for non-natural amino acids). AnarciFactory
    also uses 

    >>> from Bio.Anarci import AnarciFactory
    >>> from Bio.SeqRecord import SeqRecord
    >>> from Bio.Seq import Seq
    >>> rec = SeqRecord(Seq('slrlscaasgftfnnctihwvrqapgkgldwvavisydgankydadsvkgrftisrdnsk \
                         ntlylqmnslgsedtavyycasessgsyfdlwgrgtlv'), id='AAD14970.2')
    >>> AF = AnarciFactory(scheme='kabat')
    >>> gen = AF.run(rec)
    >>> print(next(gen))
    
    
    """
    def __init__(self, scheme=None, restrict=None, hmmerpath=None, outfile_hits=None,
                ncpu=None, assign_germline=None, use_species=None, bit_score_threshold=None):

        self.cline = AnarciCommandline()
        
        if scheme:
            self.scheme = scheme
            self.cline.scheme = scheme
        if restrict:
            self.restrict = restrict
            self.cline.restrict = restrict
        if hmmerpath:
            self.hmmerpath = hmmerpath
            self.cline.hmmerpath = hmmerpath
        if outfile_hits:
            self.outfile_hits = outfile_hits
            self.cline.outfile_hits = outfile_hits
        if ncpu:
            self.ncpu = ncpu
            self.cline.ncpu = ncpu
        if assign_germline:
            self.assign_germline = assign_germline
            self.cline.assign_germline = assign_germline
        if use_species:
            self.use_species = use_species
            self.cline.use_species = use_species
        if bit_score_threshold:
            self.bit_score_threshold = bit_score_threshold
            self.cline.bit_score_threshold = bit_score_threshold

    def _validate(self, seq):
        """Check if Sequence object contains only the 20 natural amino acids (PRIVATE)."""           
        # return False if sequence is not valid
        return not all([aa.upper() in IUPAC.IUPACProtein.letters for aa in seq])
 
    def run(self, inp):
        """Yield SearchIO.parse generator object for the Anarci numbered output."""        
        if isinstance(inp, Seq.Seq):
            if self._validate(seq):
                self.cline.sequence = str(inp)
            else:
                warnings.warn('Invalid amino acid present in Seq object: {}'.format(seq.id))
        elif isinstance(inp, SeqRecord.SeqRecord):
            FH = open('temp_file.fa', 'w')
            if self._validate(inp):
                SeqIO.write(inp, 'temp_file.fa', 'fasta')
                self.cline.sequence = 'temp_file.fa'
            else:
                warnings.warn('Invalid amino acid present in SeqRecord object: {}'.format(inp.id))
            FH.close()
        elif isinstance(inp, list):
            FH = open('temp_file.fa', 'a')
            for seq in inp:
                if self._validate(seq):
                    SeqIO.write(seq, FH, 'fasta')
                    self.cline.sequence = 'temp_file.fa'
                else:
                    warnings.warn('Invalid amino acid present in input list of Seq(Record) objects: {}'.format(seq.id))
            FH.close()
        else:       
            raise TypeError('Input sequence type must be Fasta, Seq, SeqRecord or a list of Seq objects.')
        
        print(self.cline)
        # call AnarciCommandline object to get output
        # first element of tuple is actual output, second element is a warning
        anarci_output = io.StringIO(self.cline()[0])
        #print(anarci_output.read())
        # remove temporary file if present
        if os.path.exists('temp_file.fa'):
            os.remove('temp_file.fa')

        return SearchIO.parse(anarci_output, 'anarci-num')

    # def run_hits(self, inp):
        # """Yield SearchIO.parse generator object for the Anarci hits output."""
        # assert self.outfile_hits, 'AnarciFactory object has no outfile_hits file specified'


#TODO: check if AnarciFactory can inherit from AnarciCommandline
#TODO: make __repr__ and __str__, see ^
#TODO: make factory run for hits output format

# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
