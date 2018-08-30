# Copyright 2018 by Lucas Davey.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for Anarci numbered output formats."""

import re
import warnings

from Bio import BiopythonWarning
from Bio._py3k import _as_bytes, _bytes_to_string
from Bio.Alphabet import generic_protein
from Bio.SearchIO._index import SearchIndexer
from Bio.SearchIO._model import QueryResult, Hit, HSP, HSPFragment
from Bio.SeqFeature import SeqFeature, FeatureLocation

__all__ = ('AnarciNumParser', 'AnarciNumIndexer')


class AnarciNumParser:
    """Parser for the Anarci numbered format."""

    def __init__(self, handle):
        self.handle = handle
        self.line = self.handle.readline()
        assert not self.line.startswith('# Hit file for ANARCI'), 'Format should be "anarci-hits"'

    def __iter__(self):
        # break out if it's an empty file
        if not self.line:
            raise StopIteration

        for qresult in self._parse_qresult():
            yield qresult

    def _parse_qresult(self):
        """Yield QueryResult objects (PRIVATE)."""
        _RE_STATS = re.compile(r'\|\w+\|\w\|.+\|.+\|\d+\|\d+\|')
        _RE_GERMLINE_STATS = re.compile(r'\|\w+\|.+\|\d\.\d+\|.+\|\d\.\d+\|')
        
        # state values, used to determine where in the file a line is
        _STATE_HEADER = 0
        _STATE_ANNOTATION = 1

        sequence = ''
        annotation = list()
        feature_lst = list()
        state = _STATE_HEADER
        assign_germline = False
        residue_counter = 0

        self.handle.seek(0)
        for line in self.handle:
            assert line.strip(), 'File is either empty or contains blanklines'
            # "//" marks the end of a query result
            if line[:2] == '//':   
                # queryresult needs a sequence or it won't contain any hits
                if sequence:
                    # make HSPFragment object
                    frag = HSPFragment(query=sequence, hit_id=hit_id, query_id=query_id, alphabet=generic_protein)
                    frag.query_start = int(seqstart_index)
                    frag.query_end = residue_counter
                    frag.hit_start = 0
                    frag.hit_end = 0
                    feature = SeqFeature(FeatureLocation(seqstart_index, residue_counter), type='ANARCI_hit')
                    feature.qualifiers['sequence'] = sequence
                    feature.qualifiers['evalue'] = evalue
                    feature.qualifiers['bitscore'] = score
                    feature.qualifiers['chain_type'] = chain_type
                    feature.qualifiers['species'] = species
                    feature_lst.append(feature)
                    frag.query.features = feature_lst
                    frag.query.letter_annotations[scheme] = annotation
                    # make HSP object
                    hsp = HSP([frag])
                    hsp.evalue = evalue
                    hsp.bitscore = score
                    # make Hit object
                    hit = Hit(hsps=[hsp], id=hit_id, query_id=query_id)
                    hit.description = ''
                    # make QueryResult object
                    qresult = QueryResult(hits=[hit])
                    qresult.program = 'ANARCI'
                    qresult.version = '1.3'
                    qresult.target = 'NCBI protein database'
                    qresult.description = query_description
                    setattr(qresult, 'scheme', scheme)
                    setattr(qresult, 'chain_type', chain_type)
                    setattr(qresult, 'species', species)
                    # assign germline attributes if specified                
                    if assign_germline:
                       setattr(qresult, 'species_germline', species_germline)
                       setattr(qresult, 'v_gene', v_gene)
                       setattr(qresult, 'v_identity', v_identity)
                       setattr(qresult, 'j_gene', j_gene)
                       setattr(qresult, 'j_identity', j_identity)
                else: 
                    qresult = QueryResult(hits=[], id=query_id)
                    qresult.program = 'ANARCI'
                    qresult.version = '1.3'
                    qresult.target = 'NCBI protein database'
                    qresult.description = query_description
                    # warnings.warn('Query could not be annotated: "{}"'.format(query_id), BiopythonWarning)            
                yield qresult
                # reset for next query result
                sequence =  ''
                annotation = list()
                feature_lst = list()
                residue_counter = 0
                state = _STATE_HEADER
            else:
                # parse comments
                if line[0] == '#':
                    # find query id, is always on the first line
                    if state == _STATE_HEADER:
                        line = line.split()
                        if line[1:3] == ['Input', 'sequence']:
                            query_id = '<unknown id>' 
                            query_description = ''
                        else:
                            query_id = line[1]
                            query_description = ' '.join(line[2:])
                        state = _STATE_ANNOTATION
                    # find scheme
                    elif line.startswith('# Scheme = '):
                        line = line.split()
                        scheme = line[-1]
                    # find stats
                    elif re.search(_RE_STATS, line):
                        line = line.split('|')
                        species = line[1]
                        chain_type = line[2]
                        evalue = float(line[3])
                        score = float(line[4])
                        seqstart_index = int(line[5])
                        seqend_index = int(line[6])
                        hit_id = species + '_' + chain_type                      
                    # find germline stats  
                    elif re.search(_RE_GERMLINE_STATS, line):
                        assign_germline = True
                        line = line.split('|')
                        species_germline = line[1]
                        v_gene = line[2]
                        v_identity = line[3]
                        j_gene = line[4]
                        j_identity = line[5]
                # moved out of comments: parse sequence and annotation
                else:
                    line = line.split()
                    if line[2] != '-':
                        # if length of "line" list is 4 then there are alternative annotations
                        if len(line) == 4:
                            annot = line[1] + line[2].lower()
                        else:
                            annot = line[1]
                        residue = line[-1].upper()
                        sequence += residue
                        annotation.append(annot)
                        # A residue location is defined as between two integers.
                        # e.g. "R (1, 2)" is the second residue in the sequence.              
                        feature = SeqFeature(FeatureLocation(seqstart_index + residue_counter, 
                                                             seqstart_index + residue_counter + 1), type=scheme)
                        feature.qualifiers['residue'] = residue
                        feature.qualifiers['label'] = residue + annot
                        feature.qualifiers['ANARCI_{}_num'.format(scheme)] = annot 
                        feature_lst.append(feature)
                        residue_counter += 1
                        state = _STATE_ANNOTATION

class AnarciNumIndexer(SearchIndexer):
    """Indexer for the Anarci numbered format."""

    _parser = AnarciNumParser

    def __init__(self, filename):
        SearchIndexer.__init__(self, filename)
        self.line = self._handle.readline()
        assert not self.line == _as_bytes('# Hit file for ANARCI\n'), 'format should be "anarci-hits"'

    def __iter__(self):        
        # break out if it's an empty file
        if not self.line:
            raise StopIteration

        for key, offset, length in self._qresult_index():
            yield _bytes_to_string(key), offset, length

    def _qresult_index(self):
        """Indexer for Anarci numbered files (PRIVATE)."""
        # state values, used to determine where in the file a line is
        _STATE_QUERY_ID = 0
        _STATE_ANNOTATION = 1
        
        handle = self._handle
        handle.seek(0)
        qresult_key = None
        start_offset = 0
        state = _STATE_QUERY_ID

        while True:
            line = handle.readline()
            end_offset = handle.tell()
            if not line:
                break
            if line.startswith(_as_bytes('#')):
                if state == _STATE_QUERY_ID:
                    qresult_key = line.split()[1]
                    state = _STATE_ANNOTATION
            else:
                state = _STATE_QUERY_ID
            if line[:2] == _as_bytes('//'):
                yield qresult_key, start_offset, \
                end_offset - start_offset
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
