# -*- coding: utf-8 -*-

import unittest

from align import AlignmentResult, aligner
from align.matrix import DNAFULL, BLOSUM62


class TestGlobal(unittest.TestCase):

    method = 'global'

    def test_global1(self):
        alns = aligner('WW', 'WEW', method=self.method, matrix=BLOSUM62,
                       max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'W-W', aln
        assert aln.seq2 == 'WEW', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 15.0, aln

    def test_global2(self):
        alns = aligner('WW', 'WEW', method=self.method, matrix=BLOSUM62,
                       gap_open=-100, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'W-W', aln
        assert aln.seq2 == 'WEW', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == -78.0, aln

    def test_global3(self):
        aln, = aligner('A', 'A', method=self.method, matrix=BLOSUM62,
                       gap_open=-7, max_hits=None)
        assert aln.seq1 == 'A', aln
        assert aln.seq2 == 'A', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 4.0, aln

    def test_global4(self):
        aln, = aligner('R', 'K', method=self.method, matrix=BLOSUM62,
                       gap_open=-7, max_hits=None)
        assert aln.seq1 == 'R', aln
        assert aln.seq2 == 'K', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 2.0, aln

    def test_global5(self):
        alns = aligner('R', 'AR', method=self.method, matrix=BLOSUM62,
                       gap_open=-7, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == '-R', aln
        assert aln.seq2 == 'AR', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == -2.0, aln

    def test_global6(self):
        alns = aligner('AR', 'R', method=self.method, matrix=BLOSUM62,
                       gap_open=-7, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'AR', aln
        assert aln.seq2 == '-R', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == -2.0, aln

    def test_global7(self):
        alns = aligner('AR', 'RA', method=self.method, matrix=BLOSUM62,
                       gap_open=-7, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'AR', aln
        assert aln.seq2 == 'RA', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == -2.0, aln

    def test_global8(self):
        alns = aligner('AR', 'RA', method=self.method, matrix=BLOSUM62,
                       gap_open=-3, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'AR-', aln
        assert aln.seq2 == '-RA', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == -1.0, aln

    def test_global9(self):
        alns = aligner('RAR', 'RR', method=self.method, matrix=BLOSUM62,
                       gap_open=-3, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'RAR', aln
        assert aln.seq2 == 'R-R', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 7.0, aln

    def test_global10(self):
        alns = aligner('RAR', 'RR', method=self.method, matrix=BLOSUM62,
                       gap_open=-10, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'RAR', aln
        assert aln.seq2 == 'R-R', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 0.0, aln

    def test_global11(self):
        alns = aligner('RAAR', 'RR', method=self.method, matrix=BLOSUM62,
                       gap_open=-5, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'RAAR', aln
        assert aln.seq2 == 'R--R', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 0.0, aln

    def test_global12(self):
        alns = aligner('RLR', 'RER', method=self.method, matrix=BLOSUM62,
                       gap_open=-9, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'RLR', aln
        assert aln.seq2 == 'RER', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 7., aln

    def test_global13(self):
        alns = aligner('RLR', 'RER', method=self.method, matrix=BLOSUM62,
                       gap_open=-1, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'RL-R', aln
        assert aln.seq2 == 'R-ER', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 8.0, aln

    def test_global14(self):
        alns = aligner('RLR', 'REER', method=self.method, matrix=BLOSUM62,
                       gap_open=-1, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'RL--R', aln
        assert aln.seq2 == 'R-EER', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 7.0, aln

    def test_global15(self):
        alns = aligner('AGEBAM', 'AGEBAMAM', method=self.method,
                       matrix=BLOSUM62, gap_open=-6, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'AGEBAM--', aln
        assert aln.seq2 == 'AGEBAMAM', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 16.0, aln

    def test_global16(self):
        alns = aligner('CPELIRKNCANTH', 'PREKRLICAN', method=self.method,
                       matrix=BLOSUM62, gap_open=-0.5, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'CP-E--LIRKNCANTH', aln
        assert aln.seq2 == '-PREKRLI---CAN--', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 34.5, aln

    def test_global17(self):
        alns = aligner('CPEL', 'PREK', method=self.method, matrix=BLOSUM62,
                       gap_open=-5, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'CP-EL', aln
        assert aln.seq2 == '-PREK', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 0.0, aln

    def test_global18(self):
        alns = aligner('RLRR', 'RRER', method=self.method, matrix=BLOSUM62,
                       gap_open=-1, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'RLR-R', aln
        assert aln.seq2 == 'R-RER', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 13.0, aln

    def test_global19(self):
        alns = aligner('TAAT', 'TAATTC', method=self.method, matrix=DNAFULL,
                       max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'TAAT--', aln
        assert aln.seq2 == 'TAATTC', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 6.0, aln

    def test_global20(self):
        alns = aligner('WR', 'WRR', method=self.method, matrix=BLOSUM62,
                       max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'WR-', aln
        assert aln.seq2 == 'WRR', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 9.0, aln

    def test_global21(self):
        alns = aligner('AIP', 'AP', method=self.method, matrix=BLOSUM62,
                       max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'AIP', aln
        assert aln.seq2 == 'A-P', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 4.0, aln

    def test_global22(self):
        alns = aligner('PAA', 'PA', method=self.method, matrix=BLOSUM62,
                       max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'PAA', aln
        assert aln.seq2 == 'PA-', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 4.0, aln

    def test_global23(self):
        alns = aligner('TAATTC', 'TAAT', method=self.method, matrix=DNAFULL,
                       gap_open=-10, gap_extend=-1, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'TAATTC', aln
        assert aln.seq2 == 'TAAT--', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 9.0, aln

    def test_global24(self):
        alns = aligner('CELECANTH', 'PELICAN', method=self.method,
                       matrix=BLOSUM62, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'CELECANTH', aln
        assert aln.seq2 == 'PELICAN--', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 8.0, aln

    def test_global25(self):
        alns = aligner('PELICAN', 'CELECANTH', method=self.method,
                       matrix=BLOSUM62, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'PELICAN--', aln
        assert aln.seq2 == 'CELECANTH', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 8.0, aln

    def test_global26(self):
        alns = aligner('AGEBANAN', 'ACEBAN', method=self.method,
                       matrix=BLOSUM62, gap_open=-2, gap_extend=-1,
                       max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'AGEBANAN', aln
        assert aln.seq2 == 'ACEBAN--', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 17.0, aln


class TestGlobalCFE(unittest.TestCase):

    method = 'global_cfe'

    def test_global_cfe1(self):
        alns = aligner('TAAT', 'TAATTC', method=self.method, matrix=DNAFULL,
                       max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'TAAT--', aln
        assert aln.seq2 == 'TAATTC', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 20.0, aln

    def test_global_cfe2(self):
        alns = aligner('TCTAAT', 'TAAT', method=self.method, matrix=DNAFULL,
                       max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'TCTAAT', aln
        assert aln.seq2 == '--TAAT', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 20.0, aln

    def test_global_cfe3(self):
        alns = aligner('PAA', 'PA', method=self.method, matrix=BLOSUM62,
                       max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'PAA', aln
        assert aln.seq2 == 'PA-', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 11.0, aln

    def test_global_cfe4(self):
        alns = aligner('AATGAA', 'AATGAATGAA', method=self.method,
                       matrix=DNAFULL, max_hits=None)
        assert len(alns) == 2, alns
        aln1 = AlignmentResult(
            seq1='AATGAA----', seq2='AATGAATGAA', pos1=0, pos2=0, score=30.0)
        aln2 = AlignmentResult(
            seq1='----AATGAA', seq2='AATGAATGAA', pos1=0, pos2=0, score=30.0)
        assert aln1 in alns, alns
        assert aln2 in alns, alns

    def test_global_cfe5(self):
        alns = aligner('AAAAAAAAAAAAACCTGCGCCCCAAAAAAAAAAAAAAAAAAAA',
                       'CCTGCGCACCCC', method='global_cfe', matrix=DNAFULL,
                       max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'AAAAAAAAAAAAACCTGCGC-CCCAAAAAAAAAAAAAAAAAAAA', aln
        assert aln.seq2 == '-------------CCTGCGCACCCC-------------------', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 39.0, aln


class TestLocal(unittest.TestCase):

    method = 'local'

    def test_local1(self):
        alns = aligner('TCTAAT', 'TAAT', method=self.method, matrix=DNAFULL,
                       max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'TAAT', aln
        assert aln.seq2 == 'TAAT', aln
        assert aln.pos1 == 2, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 20.0, aln

    def test_local2(self):
        alns = aligner('TCTAAT', 'TAATCT', method=self.method, matrix=DNAFULL,
                       max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'TAAT', aln
        assert aln.seq2 == 'TAAT', aln
        assert aln.pos1 == 2, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 20.0, aln

    def test_local3(self):
        alns = aligner('A', 'A', method=self.method, matrix=BLOSUM62)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'A', aln
        assert aln.seq2 == 'A', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 4.0, aln

    def test_local4(self):
        alns = aligner('RA', 'AR', method=self.method, matrix=BLOSUM62,
                       max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'R', aln
        assert aln.seq2 == 'R', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 1, aln
        assert aln.score == 5.0, aln

    def test_local5(self):
        alns = aligner('RRR', 'RR', method=self.method, matrix=BLOSUM62,
                       max_hits=None)
        assert len(alns) == 2, alns
        aln1 = AlignmentResult(seq1='RR', seq2='RR', pos1=0, pos2=0,
                               score=10.0)
        aln2 = AlignmentResult(seq1='RR', seq2='RR', pos1=1, pos2=0,
                               score=10.0)
        assert aln1 in alns, alns
        assert aln2 in alns, alns

    def test_local6(self):
        alns = aligner('PYNCHAN', 'YNCH', method=self.method, matrix=BLOSUM62,
                       max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'YNCH', aln
        assert aln.seq2 == 'YNCH', aln
        assert aln.pos1 == 1, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 30.0, aln

    def test_local7(self):
        alns = aligner('AIP', 'AP', method=self.method, matrix=BLOSUM62,
                       max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'P', aln
        assert aln.seq2 == 'P', aln
        assert aln.pos1 == 2, aln
        assert aln.pos2 == 1, aln
        assert aln.score == 7.0, aln

    def test_local8(self):
        alns = aligner('PAA', 'PA', method=self.method, matrix=BLOSUM62,
                       max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'PA', aln
        assert aln.seq2 == 'PA', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 11.0, aln


class TestGlocal(unittest.TestCase):

    method = 'glocal'

    def test_glocal1(self):
        alns = aligner('AAATAATAAA', 'TAAT', method=self.method,
                       matrix=DNAFULL, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'TAAT', aln
        assert aln.seq2 == 'TAAT', aln
        assert aln.pos1 == 3, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 20.0, aln

    def test_glocal2(self):
        alns = aligner('AAATAATAAA', 'TATAT', method=self.method,
                       matrix=DNAFULL, gap_open=-1, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'TA-AT', aln
        assert aln.seq2 == 'TATAT', aln
        assert aln.pos1 == 3, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 19.0, aln

    def test_glocal3(self):
        alns = aligner('TATATAAA', 'CCTATAT', method=self.method,
                       matrix=DNAFULL, gap_open=-8, gap_double=-8,
                       max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == '--TATAT', aln
        assert aln.seq2 == 'CCTATAT', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 10.0, aln

    def test_glocal4(self):
        alns = aligner('CCTATAT', 'TATATAAA', method=self.method,
                       matrix=DNAFULL, gap_open=-8, gap_double=-8,
                       max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == 'CCTATAT', aln
        assert aln.seq2 == '--TATAT', aln
        assert aln.pos1 == 0, aln
        assert aln.pos2 == 0, aln
        assert aln.score == 10.0, aln


if __name__ == '__main__':
    unittest.main()
