# Copyright 2018 by Lucas Davey.  All rights reserved.

"""Tests for SearchIO anarci-num parsing."""

import os
import unittest

from Bio.Alphabet import generic_protein
from Bio.SearchIO import parse

# test case files are in the Anarci directory
TEST_DIR = 'Anarci'
FMT = 'anarci-num'

def get_file(filename):
    """Returns the path of a test file."""
    return os.path.join(TEST_DIR, filename)

class BaseAnarciCases(unittest.TestCase):
    """Check common attributes"""
    def check_common_attrs(self, qresults):
        # check common attributes
        for qresult in qresults:
            # only check if has hits
            if qresult.hits:
                for hit in qresult:
                    self.assertEqual(qresult.id, hit.query_id)
                    for hsp in hit:
                        self.assertEqual(hit.id, hsp.hit_id)
                        self.assertEqual(qresult.id, hsp.query_id)

class AnarciNumCases(BaseAnarciCases):

    def test_num_anarci_001(self):
        """Test parsing anarci-num output (num_anarci_001.txt)"""

        anarci_file = get_file('num_anarci_001.txt')
        qresults = list(parse(anarci_file, FMT))
        self.assertEqual(10, len(qresults))
        self.check_common_attrs(qresults)
        # test first qresult
        qresult = qresults[0]
        self.assertEqual(1, len(qresult.hits))
        self.assertEqual('AXG50451.1', qresult.id)
        self.assertEqual('immunoglobulin heavy chain variable region, partial [Cricetulus migratorius]', qresult.description)
        self.assertEqual('NCBI protein database', qresult.target)
        self.assertEqual('ANARCI', qresult.program)
        self.assertEqual('1.3', qresult.version)
        self.assertEqual('kabat', qresult.scheme)
        self.assertEqual('H', qresult.chain_type)
        self.assertEqual('mouse', qresult.species)
        self.assertEqual(1, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('mouse_H', hit.id)
        self.assertEqual('', hit.description)
        self.assertEqual(1, len(hit)) 
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(116, hsp.aln_span)
        self.assertEqual(1.5e-53, hsp.evalue)
        self.assertEqual(171.1, hsp.bitscore)
        self.assertEqual(20, hsp.query_start)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(116, hsp.query_end)
        self.assertEqual(0, hsp.hit_end)
        self.assertEqual('QVKLEESGPGLVNPSQSLSLSCSVTGYSITSGYGWNWIRQ', str(hsp.query.seq)[:40])
        self.assertEqual('NQFFLQVNSVTTEDTATYYCARQTGYFDYWGQGTMVTVSS', str(hsp.query.seq)[-40:])
        # first qresult, first hit, first hsp, first fragment
        frag = qresult[0].hsps[0][0]
        self.assertEqual(generic_protein, frag.alphabet)
        features = frag.query.features
        self.assertEqual(117, len(features))
        self.assertEqual('kabat', features[0].type)
        self.assertEqual('[20:21]', str(features[0].location))
        self.assertEqual('Q', features[0].qualifiers['residue'])
        self.assertEqual('Q1', features[0].qualifiers['label'])
        self.assertEqual('1', features[0].qualifiers['ANARCI_kabat_num'])
        self.assertEqual('kabat', features[51].type)
        self.assertEqual('[71:72]', str(features[51].location))
        self.assertEqual('I', features[51].qualifiers['residue'])
        self.assertEqual('I51', features[51].qualifiers['label'])
        self.assertEqual('51', features[51].qualifiers['ANARCI_kabat_num'])
        self.assertEqual('kabat', features[88].type)
        self.assertEqual('[108:109]', str(features[88].location))
        self.assertEqual('E', features[88].qualifiers['residue'])
        self.assertEqual('E85', features[88].qualifiers['label'])
        self.assertEqual('85', features[88].qualifiers['ANARCI_kabat_num'])
        self.assertEqual('[20:116]', str(features[-1].location))
        self.assertEqual(116, len(features[-1].qualifiers['sequence']))
        self.assertEqual(1.5e-53, features[-1].qualifiers['evalue'])
        self.assertEqual(171.1, features[-1].qualifiers['bitscore'])
        self.assertEqual('H', features[-1].qualifiers['chain_type'])
        self.assertEqual('mouse', features[-1].qualifiers['species'])
        self.assertEqual('QVKLEESGPGLVNPSQSLSLSCSVTGYSITSGYGWNWIRQ', features[-1].qualifiers['sequence'][:40])
        self.assertEqual('NQFFLQVNSVTTEDTATYYCARQTGYFDYWGQGTMVTVSS', features[-1].qualifiers['sequence'][-40:])
        # test second qresult
        qresult = qresults[1]
        self.assertEqual('NP_001138472.1', qresult.id)
        self.assertEqual('integrin alpha-V isoform 3 precursor [Homo sapiens]', qresult.description)
        self.assertEqual('NCBI protein database', qresult.target)
        self.assertEqual('ANARCI', qresult.program)
        self.assertEqual('1.3', qresult.version)
        self.assertEqual(0, len(qresult))
        # test third qresult
        qresult = qresults[2]
        self.assertEqual('NP_001138471.1', qresult.id)
        self.assertEqual('integrin alpha-V isoform 2 precursor [Homo sapiens]', qresult.description)
        self.assertEqual('NCBI protein database', qresult.target)
        self.assertEqual('ANARCI', qresult.program)
        self.assertEqual('1.3', qresult.version)
        self.assertEqual(0, len(qresult))
        # test fourth qresult
        qresult = qresults[3]
        self.assertEqual('NP_002201.1', qresult.id)
        self.assertEqual('integrin alpha-V isoform 1 preproprotein [Homo sapiens]', qresult.description)
        self.assertEqual('NCBI protein database', qresult.target)
        self.assertEqual('ANARCI', qresult.program)
        self.assertEqual('1.3', qresult.version)
        self.assertEqual(0, len(qresult))
        # test fifth qresult
        qresult = qresults[4]
        self.assertEqual('NP_032428.2', qresult.id)
        self.assertEqual('integrin alpha-V precursor [Mus musculus]', qresult.description)
        self.assertEqual('NCBI protein database', qresult.target)
        self.assertEqual('ANARCI', qresult.program)
        self.assertEqual('1.3', qresult.version)
        self.assertEqual(0, len(qresult))
        # test sixth qresult
        qresult = qresults[5]
        self.assertEqual('sp|P0DPI1.1|BXA1_CLOBH', qresult.id)
        self.assertEqual('RecName: Full=Botulinum neurotoxin type A; Short=B', qresult.description[:50])
        self.assertEqual('NCBI protein database', qresult.target)
        self.assertEqual('ANARCI', qresult.program)
        self.assertEqual('1.3', qresult.version)
        self.assertEqual(0, len(qresult))
        # test seventh qresult
        qresult = qresults[6]
        self.assertEqual('sp|P18640.3|BXC_CBCP', qresult.id)
        self.assertEqual('RecName: Full=Botulinum neurotoxin type C; Short=B', qresult.description[:50])
        self.assertEqual('NCBI protein database', qresult.target)
        self.assertEqual('ANARCI', qresult.program)
        self.assertEqual('1.3', qresult.version)
        self.assertEqual(0, len(qresult))
        # test eigth qresult
        qresult = qresults[7]
        self.assertEqual(1, len(qresult.hits))
        self.assertEqual('sp|P0DOX5.2|IGG1_HUMAN', qresult.id)
        self.assertEqual('RecName: Full=Immunoglobulin gamma-1 heavy chain; AltName: Full=Immunoglobulin gamma-1 heavy chain NIE', qresult.description)
        self.assertEqual('NCBI protein database', qresult.target)
        self.assertEqual('ANARCI', qresult.program)
        self.assertEqual('1.3', qresult.version)
        self.assertEqual('kabat', qresult.scheme)
        self.assertEqual('H', qresult.chain_type)
        self.assertEqual('human', qresult.species)
        self.assertEqual(1, len(qresult))
        # eigth qresult, first hit
        hit = qresult[0]
        self.assertEqual('human_H', hit.id)
        self.assertEqual('', hit.description)
        self.assertEqual(1, len(hit)) 
        # eigth qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(119, hsp.aln_span)
        self.assertEqual(3.9e-58, hsp.evalue)
        self.assertEqual(186.1, hsp.bitscore)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(119, hsp.query_end)
        self.assertEqual(0, hsp.hit_end)
        self.assertEqual('QVQLVQSGGGVVQPGRSLRLSCAASGFTFSRYTIHWVRQA', str(hsp.query.seq)[:40])
        self.assertEqual('YLNMNSLRPEDTAVYYCARIRDTAMFFAHWGQGTLVTVSS', str(hsp.query.seq)[-40:])
        # eigth qresult, first hit, first hsp, first fragment
        frag = qresult[0].hsps[0][0]
        self.assertEqual(generic_protein, frag.alphabet)
        features = frag.query.features
        self.assertEqual(120, len(features))
        self.assertEqual('kabat', features[0].type)
        self.assertEqual('[0:1]', str(features[0].location))
        self.assertEqual('Q', features[0].qualifiers['residue'])
        self.assertEqual('Q1', features[0].qualifiers['label'])
        self.assertEqual('1', features[0].qualifiers['ANARCI_kabat_num'])
        self.assertEqual('kabat', features[51].type)
        self.assertEqual('[51:52]', str(features[51].location))
        self.assertEqual('S', features[51].qualifiers['residue'])
        self.assertEqual('S52', features[51].qualifiers['label'])
        self.assertEqual('52', features[51].qualifiers['ANARCI_kabat_num'])
        self.assertEqual('kabat', features[88].type)
        self.assertEqual('[88:89]', str(features[88].location))
        self.assertEqual('E', features[88].qualifiers['residue'])
        self.assertEqual('E85', features[88].qualifiers['label'])
        self.assertEqual('85', features[88].qualifiers['ANARCI_kabat_num'])
        self.assertEqual('[0:119]', str(features[-1].location))
        self.assertEqual(119, len(features[-1].qualifiers['sequence']))
        self.assertEqual(3.9e-58, features[-1].qualifiers['evalue'])
        self.assertEqual(186.1, features[-1].qualifiers['bitscore'])
        self.assertEqual('H', features[-1].qualifiers['chain_type'])
        self.assertEqual('human', features[-1].qualifiers['species'])
        self.assertEqual('QVQLVQSGGGVVQPGRSLRLSCAASGFTFSRYTIHWVRQA', features[-1].qualifiers['sequence'][:40])
        self.assertEqual('YLNMNSLRPEDTAVYYCARIRDTAMFFAHWGQGTLVTVSS', features[-1].qualifiers['sequence'][-40:])
        # test nineth qresult
        qresult = qresults[8]
        self.assertEqual('sp|P01834.2|IGKC_HUMAN', qresult.id)
        self.assertEqual('RecName: Full=Immunoglobulin kappa constant; AltNa', qresult.description[:50])
        self.assertEqual('NCBI protein database', qresult.target)
        self.assertEqual('ANARCI', qresult.program)
        self.assertEqual('1.3', qresult.version)
        self.assertEqual(0, len(qresult))
        # test tenth qresult
        qresult = qresults[9]
        self.assertEqual(1, len(qresult.hits))
        self.assertEqual('sp|P01615.2|KVD28_HUMAN', qresult.id)
        self.assertEqual('RecName: Full=Immunoglobulin kappa variable 2D-28;', qresult.description[:50])
        self.assertEqual('NCBI protein database', qresult.target)
        self.assertEqual('ANARCI', qresult.program)
        self.assertEqual('1.3', qresult.version)
        self.assertEqual('kabat', qresult.scheme)
        self.assertEqual('K', qresult.chain_type)
        self.assertEqual('pig', qresult.species)
        self.assertEqual(1, len(qresult))
        # tenth qresult, first hit
        hit = qresult[0]
        self.assertEqual('pig_K', hit.id)
        self.assertEqual('', hit.description)
        self.assertEqual(1, len(hit)) 
        # tenth qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(97, hsp.aln_span)
        self.assertEqual(2.2e-52, hsp.evalue)
        self.assertEqual(167.3, hsp.bitscore)
        self.assertEqual(20, hsp.query_start)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(97, hsp.query_end)
        self.assertEqual(0, hsp.hit_end)
        self.assertEqual('DIVMTQSPLSLPVTPGEPASISCRSSQSLLHSNGYNYLDW', str(hsp.query.seq)[:40])
        self.assertEqual('NRASGVPDRFSGSGSGTDFTLKISRVEAEDVGVYYCMQAL', str(hsp.query.seq)[-40:])
        # tenth qresult, first hit, first hsp, first fragment
        frag = qresult[0].hsps[0][0]
        self.assertEqual(generic_protein, frag.alphabet)
        features = frag.query.features
        self.assertEqual(98, len(features))
        self.assertEqual('kabat', features[0].type)
        self.assertEqual('[20:21]', str(features[0].location))
        self.assertEqual('D', features[0].qualifiers['residue'])
        self.assertEqual('D1', features[0].qualifiers['label'])
        self.assertEqual('1', features[0].qualifiers['ANARCI_kabat_num'])
        self.assertEqual('kabat', features[51].type)
        self.assertEqual('[71:72]', str(features[51].location))
        self.assertEqual('L', features[51].qualifiers['residue'])
        self.assertEqual('L47', features[51].qualifiers['label'])
        self.assertEqual('47', features[51].qualifiers['ANARCI_kabat_num'])
        self.assertEqual('kabat', features[88].type)
        self.assertEqual('[108:109]', str(features[88].location))
        self.assertEqual('G', features[88].qualifiers['residue'])
        self.assertEqual('G84', features[88].qualifiers['label'])
        self.assertEqual('84', features[88].qualifiers['ANARCI_kabat_num'])
        self.assertEqual('[20:97]', str(features[-1].location))
        self.assertEqual(2.2e-52, features[-1].qualifiers['evalue'])
        self.assertEqual(167.3, features[-1].qualifiers['bitscore'])
        self.assertEqual('K', features[-1].qualifiers['chain_type'])
        self.assertEqual('pig', features[-1].qualifiers['species'])
        self.assertEqual(97, len(features[-1].qualifiers['sequence']))
        self.assertEqual('DIVMTQSPLSLPVTPGEPASISCRSSQSLLHSNGYNYLDW', features[-1].qualifiers['sequence'][:40])
        self.assertEqual('NRASGVPDRFSGSGSGTDFTLKISRVEAEDVGVYYCMQAL', features[-1].qualifiers['sequence'][-40:])

    def test_anarci_num_002(self):
        """Test parsing anarci-num output with germlines assigned (num_anarci_002.txt)"""
        
        anarci_file = get_file('num_anarci_002.txt')
        qresults = list(parse(anarci_file, FMT))
        self.assertEqual(1, len(qresults))
        self.check_common_attrs(qresults)
        # test first qresult
        qresult = qresults[0]
        self.assertEqual(1, len(qresult.hits))
        self.assertEqual('pdb|6CYF|U', qresult.id)
        self.assertEqual('Chain U, IgG1 antibody, heavy chain fragment', qresult.description)
        self.assertEqual('NCBI protein database', qresult.target)
        self.assertEqual('ANARCI', qresult.program)
        self.assertEqual('1.3', qresult.version)
        self.assertEqual('chothia', qresult.scheme)
        self.assertEqual('H', qresult.chain_type)
        self.assertEqual('human', qresult.species_germline)
        self.assertEqual('IGHV3-23*04', qresult.v_gene)
        self.assertEqual('0.92', qresult.v_identity)
        self.assertEqual('IGHJ6*01', qresult.j_gene)
        self.assertEqual('1.00', qresult.j_identity)
        self.assertEqual(1, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('alpaca_H', hit.id)
        self.assertEqual('', hit.description)
        self.assertEqual(1, len(hit)) 
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(124, hsp.aln_span)
        self.assertEqual(4.2e-63, hsp.evalue)
        self.assertEqual(202.3, hsp.bitscore)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(124, hsp.query_end)
        self.assertEqual(0, hsp.hit_end)
        self.assertEqual('EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMDWVRQA', str(hsp.query.seq)[:40])
        self.assertEqual('SLRAEDTAVYYCAKEEFLPGTHYFYGMDVWGQGTTVTVSS', str(hsp.query.seq)[-40:])
        # first qresult, first hit, first hsp, first fragment
        frag = qresult[0].hsps[0][0]
        self.assertEqual(generic_protein, frag.alphabet)
        features = frag.query.features
        self.assertEqual(125, len(features))
        self.assertEqual('chothia', features[0].type)
        self.assertEqual('[0:1]', str(features[0].location))
        self.assertEqual('E', features[0].qualifiers['residue'])
        self.assertEqual('E1', features[0].qualifiers['label'])
        self.assertEqual('1', features[0].qualifiers['ANARCI_chothia_num'])
        self.assertEqual('chothia', features[51].type)
        self.assertEqual('[51:52]', str(features[51].location))
        self.assertEqual('T', features[51].qualifiers['residue'])
        self.assertEqual('T52', features[51].qualifiers['label'])
        self.assertEqual('52', features[51].qualifiers['ANARCI_chothia_num'])
        self.assertEqual('chothia', features[88].type)
        self.assertEqual('[88:89]', str(features[88].location))
        self.assertEqual('E', features[88].qualifiers['residue'])
        self.assertEqual('E85', features[88].qualifiers['label'])
        self.assertEqual('85', features[88].qualifiers['ANARCI_chothia_num'])
        self.assertEqual('[0:124]', str(features[-1].location))
        self.assertEqual(4.2e-63, features[-1].qualifiers['evalue'])
        self.assertEqual(202.3, features[-1].qualifiers['bitscore'])
        self.assertEqual('H', features[-1].qualifiers['chain_type'])
        self.assertEqual('alpaca', features[-1].qualifiers['species'])
        self.assertEqual(124, len(features[-1].qualifiers['sequence']))
        self.assertEqual('EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMDWVRQA', features[-1].qualifiers['sequence'][:40])
        self.assertEqual('SLRAEDTAVYYCAKEEFLPGTHYFYGMDVWGQGTTVTVSS', features[-1].qualifiers['sequence'][-40:])
        
if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
