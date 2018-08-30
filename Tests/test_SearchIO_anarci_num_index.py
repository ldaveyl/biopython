# Copyright 2018 by Lucas Davey.  All rights reserved.

"""Tests for SearchIO anarci-num indexing."""

import unittest

from search_tests_common import CheckRaw, CheckIndex

class AnarciNumRawCases(CheckRaw):
    """Check Anarci numbering get_raw method.""" 
    
    fmt = 'anarci-num'
    
    def test_num_anarci_positive(self):
        """Test anarci-num raw string retrieval, Anarci 1.3, positive hit"""
        filename = 'Anarci/num_anarci_001.txt'
        raw = """# AXG50451.1 immunoglobulin heavy chain variable region, partial [Cricetulus migratorius]
# ANARCI numbered
# Domain 1 of 1
# Most significant HMM hit
#|species|chain_type|e-value|score|seqstart_index|seqend_index|
#|mouse|H|1.5e-53|171.1|20|135|
# Scheme = kabat
H 1       Q
H 2       V
H 3       K
H 4       L
H 5       E
H 6       E
H 7       S
H 8       G
H 9       P
H 10      G
H 11      L
H 12      V
H 13      N
H 14      P
H 15      S
H 16      Q
H 17      S
H 18      L
H 19      S
H 20      L
H 21      S
H 22      C
H 23      S
H 24      V
H 25      T
H 26      G
H 27      Y
H 28      S
H 29      I
H 30      T
H 31      S
H 32      G
H 33      Y
H 34      G
H 35      W
H 35    A N
H 36      W
H 37      I
H 38      R
H 39      Q
H 40      F
H 41      P
H 42      G
H 43      Q
H 44      K
H 45      V
H 46      E
H 47      W
H 48      M
H 49      G
H 50      F
H 51      I
H 52      Y
H 53      Y
H 54      E
H 55      G
H 56      S
H 57      T
H 58      Y
H 59      Y
H 60      N
H 61      P
H 62      S
H 63      I
H 64      K
H 65      S
H 66      R
H 67      I
H 68      S
H 69      I
H 70      T
H 71      R
H 72      D
H 73      T
H 74      S
H 75      K
H 76      N
H 77      Q
H 78      F
H 79      F
H 80      L
H 81      Q
H 82      V
H 82    A N
H 82    B S
H 82    C V
H 83      T
H 84      T
H 85      E
H 86      D
H 87      T
H 88      A
H 89      T
H 90      Y
H 91      Y
H 92      C
H 93      A
H 94      R
H 95      Q
H 96      T
H 97      G
H 98      Y
H 99      F
H 100     -
H 101     D
H 102     Y
H 103     W
H 104     G
H 105     Q
H 106     G
H 107     T
H 108     M
H 109     V
H 110     T
H 111     V
H 112     S
H 113     S
//
"""
       
        self.check_raw(filename, "AXG50451.1", raw)
        
    def test_num_anarci_negative(self):
        """Test anarci-num raw string retrieval, Anarci 1.3, negative hit"""
        filename = 'Anarci/num_anarci_001.txt'
        raw ="""# NP_001138472.1 integrin alpha-V isoform 3 precursor [Homo sapiens]
//
"""

class AnarciNumIndexCases(CheckIndex):

    fmt = 'anarci-num'
    
    def test_num_anarci(self):
        """Test anarci-num indexing, Anarci 1.3"""
        filename = 'Anarci/num_anarci_001.txt'
        self.check_index(filename, self.fmt)     
        
if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
unittest.main(testRunner=runner)
