# Copyright 2018 by Lucas Davey.  All rights reserved.

"""Bio.SearchIO parser for Anarci hits output formats."""

import re
import warnings

from Bio import BiopythonWarning
from Bio._py3k import _as_bytes, _bytes_to_string
from Bio.Alphabet import generic_protein
from Bio.SearchIO._index import SearchIndexer
from Bio.SearchIO._model import QueryResult, Hit, HSP, HSPFragment
from Bio.SeqFeature import SeqFeature, FeatureLocation

__all__ = ('AnarciHitsParser', 'AnarciHitsIndexer')


class AnarciHitsParser:
    """Parser for the Anarci hits format."""
    
    def __init__(self, handle):
        self.handle = handle
        self.line = self.handle.readline()
        
    def __iter__(self):
        """Iterate over the file handle; yields QueryResult object."""
        # break out if it's an empty file
        if not self.line:
            raise StopIteration

        for qresult in self._parse_qresult():
            yield qresult
            
    def _extract_sequence(self, sequence, start, end):
        """Extracts query sequence from input sequence (PRIVATE)."""
        return sequence[start:end] 

    def _parse_qresult(self):
        """Yield QueryResult objects (PRIVATE)."""
        # state values, used to determine where in the file a line is
        _STATE_HEADER = 0
        _STATE_ANNOTATION = 1
        
        state = _STATE_HEADER
        sequence = ''
        hit_id = ''
        hit_lst = list()
        self.handle.seek(0)

        while True:
            line = self.handle.readline()
            assert line.strip(), 'File is either empty or contains blanklines'
            if not line:
                break
            if line[:2] == '//':
                if hit_id:
                    yield qresult
                else:
                    qresult = QueryResult(hits=[], id=query_id)
                    qresult.program = 'ANARCI'
                    qresult.version = '1.3'
                    qresult.target = 'NCBI protein database'
                    qresult.description = query_description   
                    # warnings.warn('Query could not be annotated: "{}"'.format(query_id), BiopythonWarning)          
                    yield qresult
                # reset for next QueryResult
                hit_lst = list()
                sequence = ''
                hit_id = ''
                state = _STATE_HEADER 

            elif state == _STATE_HEADER:
                # find query id
                if line.startswith('NAME'):
                    line = line.split()
                    query_id = line[1]
                    query_description = ' '.join(line[2:])
                # find sequence
                elif line.startswith('SEQUENCE'):
                    line = line.split()
                    sequence += line[-1]
                # find when header has ended
                elif line.startswith('         id description'):
                    state = _STATE_ANNOTATION
            else:
                # extract statistics
                line = line.split()
                species, chain_type = line[0].split('_')
                evalue = float(line[1])
                score = float(line[2])
                seqstart_index = int(line[4])
                seqend_index = int(line[5])
                hit_id = species + '_' + chain_type
                # make HSPFragment object
                frag = HSPFragment(query=self._extract_sequence(sequence, seqstart_index, seqend_index),
                                   hit_id=hit_id, query_id=query_id, alphabet=generic_protein)
                frag.query_start = int(seqstart_index)
                frag.query_end = int(seqend_index)
                frag.hit_start = 0
                frag.hit_end = 0
                frag.query.description = query_description
                # make HSP object
                hsp = HSP([frag])
                hsp.evalue = evalue
                hsp.bitscore = score
                # make Hit object
                hit = Hit(hsps=[hsp], id=hit_id, query_id=query_id)
                hit.description = ''
                hit.query_description = query_description
                hit_lst.append(hit)
                # make QueryResult object
                qresult = QueryResult(hits=hit_lst)
                qresult.program = 'Anarci'
                qresult.version = '1.3'
                qresult.target = 'NCBI protein database'
                qresult.chain_type = chain_type
                qresult.species = species
                qresult.description = query_description
                setattr(qresult, 'input_sequence', sequence)                

class AnarciHitsIndexer(SearchIndexer):
    """Indexer for the Anarci hits format."""

    _parser = AnarciHitsParser

    def __init__(self, filename):
        SearchIndexer.__init__(self, filename)
        self.line = self._handle.readline()
        if not self.line == _as_bytes('# Hit file for ANARCI\n'):
            raise AssertionError('Wrong format or file')
       
    def __iter__(self):
        """Iterate over the file handle; yields key, start offset, and length."""
        handle = self._handle

        for key, offset, length in self._qresult_index():
            yield _bytes_to_string(key), offset, length

    def _qresult_index(self):
        """Indexer for Anarci hit files (PRIVATE)."""
        handle = self._handle
        # first line was read in __init__ and is skipped
        # because it would otherwise mess up indexing
        start_offset = 0
        
        for line in handle:
            end_offset = handle.tell()
            if not line:
                break
            # find query id
            if line.startswith(_as_bytes('NAME')):
                line = line.split()
                qresult_key = line[1]
            if line[:2] == _as_bytes('//'):
                yield qresult_key, start_offset, end_offset - start_offset
                start_offset = end_offset

    def get_raw(self, offset):
        """Return the raw bytes string of a QueryResult object from the given offset."""
        handle = self._handle
        handle.seek(offset)
        qresult_raw = _as_bytes('')

        while True:
            line = handle.readline()
            if not line:
                break
            qresult_raw += line
            if line[:2] == _as_bytes('//'):
                return qresult_raw

# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
