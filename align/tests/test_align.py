# -*- coding: utf-8 -*-

import unittest

from align import AlignmentResult, aligner
from align.matrix import DNAFULL, BLOSUM62


class TestGlobal(unittest.TestCase):

    def test_global1(self):
        aln, = aligner('WW', 'WEW', method='global')
        assert aln.seq1 == 'W-W', aln
        assert aln.seq2 == 'WEW', aln

    def test_global2(self):
        aln, = aligner('WW', 'WEW', method='global', gap_open=-100)
        assert aln.seq1 == 'W-W', aln
        assert aln.seq2 == 'WEW', aln

    def test_global3(self):
        aln, = aligner('A', 'A', method='global', gap_open=-7)
        assert aln.seq1 == 'A', aln
        assert aln.seq2 == 'A', aln

    def test_global4(self):
        aln, = aligner('R', 'K', method='global', gap_open=-7)
        assert aln.seq1 == 'R', aln
        assert aln.seq2 == 'K', aln

    def test_global5(self):
        aln, = aligner('R', 'AR', method='global', gap_open=-7)
        assert aln.seq1 == '-R', aln
        assert aln.seq2 == 'AR', aln

    def test_global6(self):
        aln, = aligner('AR', 'R', method='global', gap_open=-7)
        assert aln.seq1 == 'AR', aln
        assert aln.seq2 == '-R', aln

    def test_global7(self):
        aln, = aligner('AR', 'RA', method='global', gap_open=-7)
        assert aln.seq1 == 'AR', aln
        assert aln.seq2 == 'RA', aln

    def test_global8(self):
        aln, = aligner('AR', 'RA', method='global', gap_open=-3)
        assert aln.seq1 == 'AR-', aln
        assert aln.seq2 == '-RA', aln

    def test_global9(self):
        aln, = aligner('RAR', 'RR', method='global', gap_open=-3)
        assert aln.seq1 == 'RAR', aln
        assert aln.seq2 == 'R-R', aln

    def test_global10(self):
        aln, = aligner('RAR', 'RR', method='global', gap_open=-10)
        assert aln.seq1 == 'RAR', aln
        assert aln.seq2 == 'R-R', aln

    def test_global11(self):
        aln, = aligner('RAAR', 'RR', method='global', gap_open=-5)
        assert aln.seq1 == 'RAAR', aln
        assert aln.seq2 == 'R--R', aln

    def test_global12(self):
        aln, = aligner('RLR', 'RER', method='global', gap_open=-9)
        assert aln.seq1 == 'RLR', aln
        assert aln.seq2 == 'RER', aln

    def test_global13(self):
        aln, = aligner('RLR', 'RER', method='global', gap_open=-1)
        assert aln.seq1 == 'RL-R', aln
        assert aln.seq2 == 'R-ER', aln

    def test_global14(self):
        aln, = aligner('RLR', 'REER', method='global', gap_open=-1)
        assert aln.seq1 == 'RL--R', aln
        assert aln.seq2 == 'R-EER', aln

    def test_global15(self):
        aln, = aligner('AGEBAM', 'AGEBAMAM', method='global', gap_open=-6)
        assert aln.seq1 == 'AGEBAM--', aln
        assert aln.seq2 == 'AGEBAMAM', aln

    def test_global16(self):
        aln, = aligner('CPELIRKNCANTH', 'PREKRLICAN', method='global',
                       gap_open=-0.5)
        assert aln.seq1 == 'CP-E--LIRKNCANTH', aln
        assert aln.seq2 == '-PREKRLI---CAN--', aln

    def test_global17(self):
        aln, = aligner('CPEL', 'PREK', method='global', gap_open=-5)
        assert aln.seq1 == 'CP-EL', aln
        assert aln.seq2 == '-PREK', aln

    def test_global18(self):
        aln, = aligner('RLRR', 'RRER', method='global', gap_open=-1)
        assert aln.seq1 == 'RLR-R', aln
        assert aln.seq2 == 'R-RER', aln

    def test_global19(self):
        aln, = aligner('TAAT', 'TAATTC', method='global', matrix=DNAFULL)
        assert aln.seq1 == 'TAAT--', aln
        assert aln.seq2 == 'TAATTC', aln

    def test_global20(self):
        aln, = aligner('WR', 'WRR', method='global')
        assert aln.seq1 == 'WR-', aln
        assert aln.seq2 == 'WRR', aln

    def test_global21(self):
        aln, = aligner('AIP', 'AP', method='global')
        assert aln.seq1 == 'AIP', aln
        assert aln.seq2 == 'A-P', aln

    def test_global22(self):
        aln, = aligner('PAA', 'PA', method='global')
        assert aln.seq1 == 'PAA', aln
        assert aln.seq2 == 'PA-', aln

    def test_global23(self):
        aln, = aligner('TAATTC', 'TAAT', method='global', matrix=DNAFULL,
                       gap_open=-10, gap_extend=-1)
        assert aln.seq1 == 'TAATTC', aln
        assert aln.seq2 == 'TAAT--', aln

    def test_global24(self):
        aln, = aligner('CELECANTH', 'PELICAN', method='global')
        assert aln.seq1 == 'CELECANTH', aln
        assert aln.seq2 == 'PELICAN--', aln

    def test_global25(self):
        aln, = aligner('PELICAN', 'CELECANTH', method='global')
        assert aln.seq1 == 'PELICAN--', aln
        assert aln.seq2 == 'CELECANTH', aln

    def test_global26(self):
        aln, = aligner('AGEBANAN', 'ACEBAN', gap_extend=-1, gap_open=-2,
                       method='global', matrix=BLOSUM62)
        assert aln.seq1 == 'AGEBANAN', aln
        assert aln.seq2 == 'ACEBAN--', aln


class TestGlobalCFE(unittest.TestCase):

    def test_global_cfe1(self):
        aln, = aligner('TAAT', 'TAATTC', method='global_cfe', matrix=DNAFULL)
        assert aln.seq1 == 'TAAT--', aln
        assert aln.seq2 == 'TAATTC', aln

    def test_global_cfe2(self):
        aln, = aligner('TCTAAT', 'TAAT', method='global_cfe', matrix=DNAFULL)
        assert aln.seq1 == 'TCTAAT', aln
        assert aln.seq2 == '--TAAT', aln

    def test_global_cfe3(self):
        aln, = aligner('PAA', 'PA', method='global_cfe')
        assert aln.seq1 == 'PAA', aln
        assert aln.seq2 == 'PA-', aln

    def test_global_cfe4(self):
        alns = set(aligner('AATGAA', 'AATGAATGAA', method='global_cfe',
                           matrix=DNAFULL, n_max_return=None))
        aln1 = AlignmentResult(
            seq1='AATGAA----', seq2='AATGAATGAA', pos1=0, pos2=0, score=30.0)
        aln2 = AlignmentResult(
            seq1='----AATGAA', seq2='AATGAATGAA', pos1=0, pos2=0, score=30.0)
        assert aln1 in alns, alns
        assert aln2 in alns, alns

    def test_global_cfe5(self):
        aln, = aligner(
            'AAAAAAAAAAAAACCTGCGCCCCAAAAAAAAAAAAAAAAAAAA',
            'CCTGCGCACCCC', method="global_cfe")
        assert aln.seq1 == 'AAAAAAAAAAAAACCTGCGC-CCCAAAAAAAAAAAAAAAAAAAA', aln
        assert aln.seq2 == '-------------CCTGCGCACCCC-------------------', aln


class TestLocal(unittest.TestCase):

    def test_local1(self):
        aln, = aligner('TCTAAT', 'TAAT', method='local', matrix=DNAFULL)
        assert aln.seq2 == 'TAAT', aln
        assert aln.seq1 == 'TAAT', aln

    def test_local2(self):
        aln, = aligner('TCTAAT', 'TAATCT', method='local', matrix=DNAFULL)
        assert aln.seq2 == 'TAAT', aln
        assert aln.seq1 == 'TAAT', aln

    def test_local3(self):
        aln, = aligner('A', 'A', method='local')
        assert aln.seq1 == 'A', aln
        assert aln.seq2 == 'A', aln

    def test_local4(self):
        aln, = aligner('RA', 'AR', method='local')
        assert aln.seq1 == 'R', aln
        assert aln.seq2 == 'R', aln

    def test_local5(self):
        aln, = aligner('RRR', 'RR', method='local')
        assert aln.seq1 == 'RR', aln
        assert aln.seq2 == 'RR', aln

    def test_local6(self):
        aln, = aligner('PYNCHAN', 'YNCH', method='local')
        assert aln.seq1 == 'YNCH', aln
        assert aln.seq2 == 'YNCH', aln

    def test_local7(self):
        aln, = aligner('AIP', 'AP', method='local')
        assert aln.seq1 == 'P', aln
        assert aln.seq2 == 'P', aln

    def test_local8(self):
        aln, = aligner('PAA', 'PA', method='local')
        assert aln.seq1 == 'PA', aln
        assert aln.seq2 == 'PA', aln


class TestGlocal(unittest.TestCase):

    def test_glocal1(self):
        aln, = aligner('AAATAATAAA', 'TAAT', method='glocal', matrix=DNAFULL)
        assert aln.seq1 == 'TAAT', aln
        assert aln.seq2 == 'TAAT', aln

    def test_glocal2(self):
        aln, = aligner('AAATAATAAA', 'TATAT', method='glocal', gap_open=-1,
                       matrix=DNAFULL)
        assert aln.seq1 == 'TA-AT', aln
        assert aln.seq2 == 'TATAT', aln

    def test_glocal3(self):
        aln, = aligner('TATATAAA', 'CCTATAT', method='glocal', gap_open=-1,
                       matrix=DNAFULL)
        assert aln.seq1 == '--TATAT', aln
        assert aln.seq2 == 'CCTATAT', aln

    def test_glocal4(self):
        aln, = aligner('CCTATAT', 'TATATAAA', method='glocal', gap_open=-1,
                       matrix=DNAFULL)
        assert aln.seq1 == 'CCTATAT', aln
        assert aln.seq2 == '--TATAT', aln


if __name__ == '__main__':
    unittest.main()
