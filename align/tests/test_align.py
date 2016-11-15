# -*- coding: utf-8 -*-

import numpy
import pyximport
pyximport.install(setup_args={
    "include_dirs": numpy.get_include(),
})  # noqa

import unittest

from align import AlignmentResult
from align.align import aligner as pyaligner
from align.calign import aligner as caligner
from align.matrix import DNAFULL, BLOSUM62


class AlignmentTests(unittest.TestCase):

    @property
    def f(self):
        return pyaligner


class TestGlobalPy(AlignmentTests):

    method = 'global'

    def test_global1(self):
        alns = self.f('WW', 'WEW', method=self.method, matrix=BLOSUM62,
                      max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'W-W', aln
        assert aln.seq2 == b'WEW', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 2, aln
        assert aln.end2 == 3, aln
        assert aln.n_gaps1 == 1, aln
        assert aln.n_gaps2 == 0, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 15.0, aln

    def test_global2(self):
        alns = self.f('WW', 'WEW', method=self.method, matrix=BLOSUM62,
                      gap_open=-100, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'W-W', aln
        assert aln.seq2 == b'WEW', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 2, aln
        assert aln.end2 == 3, aln
        assert aln.n_gaps1 == 1, aln
        assert aln.n_gaps2 == 0, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == -78.0, aln

    def test_global3(self):
        aln, = self.f('A', 'A', method=self.method, matrix=BLOSUM62,
                      gap_open=-7, max_hits=None)
        assert aln.seq1 == b'A', aln
        assert aln.seq2 == b'A', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 1, aln
        assert aln.end2 == 1, aln
        assert aln.n_gaps1 == 0, aln
        assert aln.n_gaps2 == 0, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 4.0, aln

    def test_global4(self):
        aln, = self.f('R', 'K', method=self.method, matrix=BLOSUM62,
                      gap_open=-7, max_hits=None)
        assert aln.seq1 == b'R', aln
        assert aln.seq2 == b'K', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 1, aln
        assert aln.end2 == 1, aln
        assert aln.n_gaps1 == 0, aln
        assert aln.n_gaps2 == 0, aln
        assert aln.n_mismatches == 1, aln
        assert aln.score == 2.0, aln

    def test_global5(self):
        alns = self.f('R', 'AR', method=self.method, matrix=BLOSUM62,
                      gap_open=-7, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'-R', aln
        assert aln.seq2 == b'AR', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 1, aln
        assert aln.end2 == 2, aln
        assert aln.n_gaps1 == 1, aln
        assert aln.n_gaps2 == 0, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == -2.0, aln

    def test_global6(self):
        alns = self.f('AR', 'R', method=self.method, matrix=BLOSUM62,
                      gap_open=-7, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'AR', aln
        assert aln.seq2 == b'-R', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 2, aln
        assert aln.end2 == 1, aln
        assert aln.n_gaps1 == 0, aln
        assert aln.n_gaps2 == 1, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == -2.0, aln

    def test_global7(self):
        alns = self.f('AR', 'RA', method=self.method, matrix=BLOSUM62,
                      gap_open=-7, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'AR', aln
        assert aln.seq2 == b'RA', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 2, aln
        assert aln.end2 == 2, aln
        assert aln.n_gaps1 == 0, aln
        assert aln.n_gaps2 == 0, aln
        assert aln.n_mismatches == 2, aln
        assert aln.score == -2.0, aln

    def test_global8(self):
        alns = self.f('AR', 'RA', method=self.method, matrix=BLOSUM62,
                      gap_open=-3, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'AR-', aln
        assert aln.seq2 == b'-RA', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 2, aln
        assert aln.end2 == 2, aln
        assert aln.n_gaps1 == 1, aln
        assert aln.n_gaps2 == 1, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == -1.0, aln

    def test_global9(self):
        alns = self.f('RAR', 'RR', method=self.method, matrix=BLOSUM62,
                      gap_open=-3, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'RAR', aln
        assert aln.seq2 == b'R-R', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 3, aln
        assert aln.end2 == 2, aln
        assert aln.n_gaps1 == 0, aln
        assert aln.n_gaps2 == 1, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 7.0, aln

    def test_global10(self):
        alns = self.f('RAR', 'RR', method=self.method, matrix=BLOSUM62,
                      gap_open=-10, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'RAR', aln
        assert aln.seq2 == b'R-R', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 3, aln
        assert aln.end2 == 2, aln
        assert aln.n_gaps1 == 0, aln
        assert aln.n_gaps2 == 1, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 0.0, aln

    def test_global11(self):
        alns = self.f('RAAR', 'RR', method=self.method, matrix=BLOSUM62,
                      gap_open=-5, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'RAAR', aln
        assert aln.seq2 == b'R--R', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 4, aln
        assert aln.end2 == 2, aln
        assert aln.n_gaps1 == 0, aln
        assert aln.n_gaps2 == 1, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 0.0, aln

    def test_global12(self):
        alns = self.f('RLR', 'RER', method=self.method, matrix=BLOSUM62,
                      gap_open=-9, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'RLR', aln
        assert aln.seq2 == b'RER', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 3, aln
        assert aln.end2 == 3, aln
        assert aln.n_gaps1 == 0, aln
        assert aln.n_gaps2 == 0, aln
        assert aln.n_mismatches == 1, aln
        assert aln.score == 7., aln

    def test_global13(self):
        alns = self.f('RLR', 'RER', method=self.method, matrix=BLOSUM62,
                      gap_open=-1, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'RL-R', aln
        assert aln.seq2 == b'R-ER', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 3, aln
        assert aln.end2 == 3, aln
        assert aln.n_gaps1 == 1, aln
        assert aln.n_gaps2 == 1, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 8.0, aln

    def test_global14(self):
        alns = self.f('RLR', 'REER', method=self.method, matrix=BLOSUM62,
                      gap_open=-1, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'RL--R', aln
        assert aln.seq2 == b'R-EER', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 3, aln
        assert aln.end2 == 4, aln
        assert aln.n_gaps1 == 1, aln
        assert aln.n_gaps2 == 1, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 7.0, aln

    def test_global15(self):
        alns = self.f('AGEBAM', 'AGEBAMAM', method=self.method,
                      matrix=BLOSUM62, gap_open=-6, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'AGEBAM--', aln
        assert aln.seq2 == b'AGEBAMAM', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 6, aln
        assert aln.end2 == 8, aln
        assert aln.n_gaps1 == 1, aln
        assert aln.n_gaps2 == 0, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 16.0, aln

    def test_global16(self):
        alns = self.f('CPELIRKNCANTH', 'PREKRLICAN', method=self.method,
                      matrix=BLOSUM62, gap_open=-0.5, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'CP-E--LIRKNCANTH', aln
        assert aln.seq2 == b'-PREKRLI---CAN--', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 13, aln
        assert aln.end2 == 10, aln
        assert aln.n_gaps1 == 2, aln
        assert aln.n_gaps2 == 3, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 34.5, aln

    def test_global17(self):
        alns = self.f('CPEL', 'PREK', method=self.method, matrix=BLOSUM62,
                      gap_open=-5, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'CP-EL', aln
        assert aln.seq2 == b'-PREK', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 4, aln
        assert aln.end2 == 4, aln
        assert aln.n_gaps1 == 1, aln
        assert aln.n_gaps2 == 1, aln
        assert aln.n_mismatches == 1, aln
        assert aln.score == 0.0, aln

    def test_global18(self):
        alns = self.f('RLRR', 'RRER', method=self.method, matrix=BLOSUM62,
                      gap_open=-1, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'RLR-R', aln
        assert aln.seq2 == b'R-RER', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 4, aln
        assert aln.end2 == 4, aln
        assert aln.n_gaps1 == 1, aln
        assert aln.n_gaps2 == 1, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 13.0, aln

    def test_global19(self):
        alns = self.f('TAAT', 'TAATTC', method=self.method, matrix=DNAFULL,
                      max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'TAAT--', aln
        assert aln.seq2 == b'TAATTC', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 4, aln
        assert aln.end2 == 6, aln
        assert aln.n_gaps1 == 1, aln
        assert aln.n_gaps2 == 0, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 6.0, aln

    def test_global20(self):
        alns = self.f('WR', 'WRR', method=self.method, matrix=BLOSUM62,
                      max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'WR-', aln
        assert aln.seq2 == b'WRR', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 2, aln
        assert aln.end2 == 3, aln
        assert aln.n_gaps1 == 1, aln
        assert aln.n_gaps2 == 0, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 9.0, aln

    def test_global21(self):
        alns = self.f('AIP', 'AP', method=self.method, matrix=BLOSUM62,
                      max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'AIP', aln
        assert aln.seq2 == b'A-P', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 3, aln
        assert aln.end2 == 2, aln
        assert aln.n_gaps1 == 0, aln
        assert aln.n_gaps2 == 1, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 4.0, aln

    def test_global22(self):
        alns = self.f('PAA', 'PA', method=self.method, matrix=BLOSUM62,
                      max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'PAA', aln
        assert aln.seq2 == b'PA-', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 3, aln
        assert aln.end2 == 2, aln
        assert aln.n_gaps1 == 0, aln
        assert aln.n_gaps2 == 1, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 4.0, aln

    def test_global23(self):
        alns = self.f('TAATTC', 'TAAT', method=self.method, matrix=DNAFULL,
                      gap_open=-10, gap_extend=-1, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'TAATTC', aln
        assert aln.seq2 == b'TAAT--', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 6, aln
        assert aln.end2 == 4, aln
        assert aln.n_gaps1 == 0, aln
        assert aln.n_gaps2 == 1, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 9.0, aln

    def test_global24(self):
        alns = self.f('CELECANTH', 'PELICAN', method=self.method,
                      matrix=BLOSUM62, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'CELECANTH', aln
        assert aln.seq2 == b'PELICAN--', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 9, aln
        assert aln.end2 == 7, aln
        assert aln.n_gaps1 == 0, aln
        assert aln.n_gaps2 == 1, aln
        assert aln.n_mismatches == 2, aln
        assert aln.score == 8.0, aln

    def test_global25(self):
        alns = self.f('PELICAN', 'CELECANTH', method=self.method,
                      matrix=BLOSUM62, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'PELICAN--', aln
        assert aln.seq2 == b'CELECANTH', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 7, aln
        assert aln.end2 == 9, aln
        assert aln.n_gaps1 == 1, aln
        assert aln.n_gaps2 == 0, aln
        assert aln.n_mismatches == 2, aln
        assert aln.score == 8.0, aln

    def test_global26(self):
        alns = self.f('AGEBANAN', 'ACEBAN', method=self.method,
                      matrix=BLOSUM62, gap_open=-2, gap_extend=-1,
                      max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'AGEBANAN', aln
        assert aln.seq2 == b'ACEBAN--', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 8, aln
        assert aln.end2 == 6, aln
        assert aln.n_gaps1 == 0, aln
        assert aln.n_gaps2 == 1, aln
        assert aln.n_mismatches == 1, aln
        assert aln.score == 17.0, aln

    def test_global27(self):
        alns = self.f('AATCAAG', 'AATGAATGAGTCAATG', method=self.method,
                      matrix=DNAFULL, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'AAT--------CAA-G', aln
        assert aln.seq2 == b'AATGAATGAGTCAATG', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 7, aln
        assert aln.end2 == 16, aln
        assert aln.n_gaps1 == 2, aln
        assert aln.n_gaps2 == 0, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == -28.0, aln


class TestGlobalCFEPy(AlignmentTests):

    method = 'global_cfe'

    def test_global_cfe1(self):
        alns = self.f('TAAT', 'TAATTC', method=self.method, matrix=DNAFULL,
                      max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'TAAT--', aln
        assert aln.seq2 == b'TAATTC', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 4, aln
        assert aln.end2 == 6, aln
        assert aln.n_gaps1 == 1, aln
        assert aln.n_gaps2 == 0, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 20.0, aln

    def test_global_cfe2(self):
        alns = self.f('TCTAAT', 'TAAT', method=self.method, matrix=DNAFULL,
                      max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'TCTAAT', aln
        assert aln.seq2 == b'--TAAT', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 6, aln
        assert aln.end2 == 4, aln
        assert aln.n_gaps1 == 0, aln
        assert aln.n_gaps2 == 1, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 20.0, aln

    def test_global_cfe3(self):
        alns = self.f('PAA', 'PA', method=self.method, matrix=BLOSUM62,
                      max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'PAA', aln
        assert aln.seq2 == b'PA-', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 3, aln
        assert aln.end2 == 2, aln
        assert aln.n_gaps1 == 0, aln
        assert aln.n_gaps2 == 1, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 11.0, aln

    def test_global_cfe4(self):
        alns = self.f('AATGAA', 'AATGAATGAA', method=self.method,
                      matrix=DNAFULL, max_hits=None)
        assert len(alns) == 2, alns
        aln1 = AlignmentResult(
            seq1=b'AATGAA----', seq2=b'AATGAATGAA', start1=0, start2=0,
            end1=6, end2=10, n_gaps1=1, n_gaps2=0,
            n_mismatches=0, score=30.0)
        aln2 = AlignmentResult(
            seq1=b'----AATGAA', seq2=b'AATGAATGAA', start1=0, start2=0,
            end1=6, end2=10, n_gaps1=1, n_gaps2=0,
            n_mismatches=0, score=30.0)
        assert aln1 in alns, alns
        assert aln2 in alns, alns

    def test_global_cfe5(self):
        alns = self.f('AAAAAAAAAAAAACCTGCGCCCCAAAAAAAAAAAAAAAAAAAA',
                      'CCTGCGCACCCC', method='global_cfe', matrix=DNAFULL,
                      max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'AAAAAAAAAAAAACCTGCGC-CCCAAAAAAAAAAAAAAAAAAAA', aln
        assert aln.seq2 == b'-------------CCTGCGCACCCC-------------------', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 43, aln
        assert aln.end2 == 12, aln
        assert aln.n_gaps1 == 1, aln
        assert aln.n_gaps2 == 2, aln
        assert aln.n_mismatches == 1, aln
        assert aln.score == 39.0, aln

    def test_global_cfe6(self):
        alns = self.f('AATCAAG', 'AATGAATGAGTCAATG', method=self.method,
                      matrix=DNAFULL, max_hits=None)
        assert len(alns) == 2, alns
        aln1 = AlignmentResult(seq1=b'AATCAA-G--------',
                               seq2=b'AATGAATGAGTCAATG',
                               start1=0, start2=0, end1=7, end2=16,
                               n_gaps1=2, n_gaps2=0, n_mismatches=1,
                               score=19.0)
        aln2 = AlignmentResult(seq1=b'--------AATCAA-G',
                               seq2=b'AATGAATGAGTCAATG',
                               start1=0, start2=0, end1=7, end2=16,
                               n_gaps1=2, n_gaps2=0, n_mismatches=1,
                               score=19.0)
        assert aln1 in alns, alns
        assert aln2 in alns, alns


class TestLocalPy(AlignmentTests):

    method = 'local'

    def test_local1(self):
        alns = self.f('TCTAAT', 'TAAT', method=self.method, matrix=DNAFULL,
                      max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'TAAT', aln
        assert aln.seq2 == b'TAAT', aln
        assert aln.start1 == 2, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 6, aln
        assert aln.end2 == 4, aln
        assert aln.n_gaps1 == 0, aln
        assert aln.n_gaps2 == 0, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 20.0, aln

    def test_local2(self):
        alns = self.f('TCTAAT', 'TAATCT', method=self.method, matrix=DNAFULL,
                      max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'TAAT', aln
        assert aln.seq2 == b'TAAT', aln
        assert aln.start1 == 2, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 6, aln
        assert aln.end2 == 4, aln
        assert aln.n_gaps1 == 0, aln
        assert aln.n_gaps2 == 0, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 20.0, aln

    def test_local3(self):
        alns = self.f('A', 'A', method=self.method, matrix=BLOSUM62)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'A', aln
        assert aln.seq2 == b'A', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 1, aln
        assert aln.end2 == 1, aln
        assert aln.n_gaps1 == 0, aln
        assert aln.n_gaps2 == 0, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 4.0, aln

    def test_local4(self):
        alns = self.f('RA', 'AR', method=self.method, matrix=BLOSUM62,
                      max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'R', aln
        assert aln.seq2 == b'R', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 1, aln
        assert aln.end1 == 1, aln
        assert aln.end2 == 2, aln
        assert aln.n_gaps1 == 0, aln
        assert aln.n_gaps2 == 0, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 5.0, aln

    def test_local5(self):
        alns = self.f('RRR', 'RR', method=self.method, matrix=BLOSUM62,
                      max_hits=None)
        assert len(alns) == 2, alns
        aln1 = AlignmentResult(seq1=b'RR', seq2=b'RR', start1=0, start2=0,
                               end1=2, end2=2, n_gaps1=0, n_gaps2=0,
                               n_mismatches=0, score=10.0)
        aln2 = AlignmentResult(seq1=b'RR', seq2=b'RR', start1=1, start2=0,
                               end1=3, end2=2, n_gaps1=0, n_gaps2=0,
                               n_mismatches=0, score=10.0)
        assert aln1 in alns, alns
        assert aln2 in alns, alns

    def test_local6(self):
        alns = self.f('PYNCHAN', 'YNCH', method=self.method, matrix=BLOSUM62,
                      max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'YNCH', aln
        assert aln.seq2 == b'YNCH', aln
        assert aln.start1 == 1, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 5, aln
        assert aln.end2 == 4, aln
        assert aln.n_gaps1 == 0, aln
        assert aln.n_gaps2 == 0, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 30.0, aln

    def test_local7(self):
        alns = self.f('AIP', 'AP', method=self.method, matrix=BLOSUM62,
                      max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'P', aln
        assert aln.seq2 == b'P', aln
        assert aln.start1 == 2, aln
        assert aln.start2 == 1, aln
        assert aln.end1 == 3, aln
        assert aln.end2 == 2, aln
        assert aln.n_gaps1 == 0, aln
        assert aln.n_gaps2 == 0, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 7.0, aln

    def test_local8(self):
        alns = self.f('PAA', 'PA', method=self.method, matrix=BLOSUM62,
                      max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'PA', aln
        assert aln.seq2 == b'PA', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 2, aln
        assert aln.end2 == 2, aln
        assert aln.n_gaps1 == 0, aln
        assert aln.n_gaps2 == 0, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 11.0, aln

    def test_local9(self):
        alns = self.f('AATCAAG', 'AATGAATGAGTCAATG', method=self.method,
                      matrix=DNAFULL, max_hits=None)
        assert len(alns) == 2, alns
        aln1 = AlignmentResult(seq1=b'AATCAA', seq2=b'AATGAA', start1=0,
                               start2=0, end1=6, end2=6,
                               n_gaps1=0, n_gaps2=0, n_mismatches=1,
                               score=21.0)
        aln2 = AlignmentResult(seq1=b'AATCAA', seq2=b'AGTCAA', start1=0,
                               start2=8, end1=6, end2=14,
                               n_gaps1=0, n_gaps2=0, n_mismatches=1,
                               score=21.0)
        assert aln1 in alns, alns
        assert aln2 in alns, alns

    def test_local10(self):
        alns = self.f('AT', 'ATCATCATC', method=self.method,
                      matrix=DNAFULL, max_hits=None)
        assert len(alns) == 3, alns
        aln1 = AlignmentResult(seq1=b'AT', seq2=b'AT', start1=0, start2=0,
                               end1=2, end2=2, n_gaps1=0, n_gaps2=0,
                               n_mismatches=0, score=10.0)
        aln2 = AlignmentResult(seq1=b'AT', seq2=b'AT', start1=0, start2=3,
                               end1=2, end2=5, n_gaps1=0, n_gaps2=0,
                               n_mismatches=0, score=10.0)
        aln3 = AlignmentResult(seq1=b'AT', seq2=b'AT', start1=0, start2=6,
                               end1=2, end2=8, n_gaps1=0, n_gaps2=0,
                               n_mismatches=0, score=10.0)
        assert aln1 in alns, alns
        assert aln2 in alns, alns
        assert aln3 in alns, alns

    def test_local11(self):
        alns = self.f('AT', 'ATCATCATC', method=self.method,
                      matrix=DNAFULL, max_hits=1)
        assert len(alns) == 1, alns
        aln = alns.pop()
        aln1 = AlignmentResult(seq1=b'AT', seq2=b'AT', start1=0, start2=0,
                               end1=2, end2=2, n_gaps1=0, n_gaps2=0,
                               n_mismatches=0, score=10.0)
        aln2 = AlignmentResult(seq1=b'AT', seq2=b'AT', start1=0, start2=3,
                               end1=2, end2=5, n_gaps1=0, n_gaps2=0,
                               n_mismatches=0, score=10.0)
        aln3 = AlignmentResult(seq1=b'AT', seq2=b'AT', start1=0, start2=6,
                               end1=2, end2=8, n_gaps1=0, n_gaps2=0,
                               n_mismatches=0, score=10.0)
        assert aln == aln1 or aln == aln2 or aln == aln3, alns

    def test_local12(self):
        alns = self.f('AT', 'ATCATCATC', method=self.method,
                      matrix=DNAFULL, max_hits=2)
        assert len(alns) == 2, alns
        aln1 = AlignmentResult(seq1=b'AT', seq2=b'AT', start1=0, start2=0,
                               end1=2, end2=2, n_gaps1=0, n_gaps2=0,
                               n_mismatches=0, score=10.0)
        aln2 = AlignmentResult(seq1=b'AT', seq2=b'AT', start1=0, start2=3,
                               end1=2, end2=5, n_gaps1=0, n_gaps2=0,
                               n_mismatches=0, score=10.0)
        aln3 = AlignmentResult(seq1=b'AT', seq2=b'AT', start1=0, start2=6,
                               end1=2, end2=8, n_gaps1=0, n_gaps2=0,
                               n_mismatches=0, score=10.0)
        assert len(set([aln1, aln2, aln3]).intersection(alns)) == 2, alns


class TestGlocalPy(AlignmentTests):

    method = 'glocal'

    def test_glocal1(self):
        alns = self.f('AAATAATAAA', 'TAAT', method=self.method,
                      matrix=DNAFULL, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'TAAT', aln
        assert aln.seq2 == b'TAAT', aln
        assert aln.start1 == 3, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 7, aln
        assert aln.end2 == 4, aln
        assert aln.n_gaps1 == 0, aln
        assert aln.n_gaps2 == 0, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 20.0, aln

    def test_glocal2(self):
        alns = self.f('AAATAATAAA', 'TATAT', method=self.method,
                      matrix=DNAFULL, gap_open=-1, max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'TA-AT', aln
        assert aln.seq2 == b'TATAT', aln
        assert aln.start1 == 3, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 7, aln
        assert aln.end2 == 5, aln
        assert aln.n_gaps1 == 1, aln
        assert aln.n_gaps2 == 0, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 19.0, aln

    def test_glocal3(self):
        alns = self.f('TATATAAA', 'CCTATAT', method=self.method,
                      matrix=DNAFULL, gap_open=-8, gap_double=-8,
                      max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'--TATAT', aln
        assert aln.seq2 == b'CCTATAT', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 5, aln
        assert aln.end2 == 7, aln
        assert aln.n_gaps1 == 1, aln
        assert aln.n_gaps2 == 0, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 10.0, aln

    def test_glocal4(self):
        alns = self.f('CCTATAT', 'TATATAAA', method=self.method,
                      matrix=DNAFULL, gap_open=-8, gap_double=-8,
                      max_hits=None)
        assert len(alns) == 1, alns
        aln = alns.pop()
        assert aln.seq1 == b'CCTATAT', aln
        assert aln.seq2 == b'--TATAT', aln
        assert aln.start1 == 0, aln
        assert aln.start2 == 0, aln
        assert aln.end1 == 7, aln
        assert aln.end2 == 5, aln
        assert aln.n_gaps1 == 0, aln
        assert aln.n_gaps2 == 1, aln
        assert aln.n_mismatches == 0, aln
        assert aln.score == 10.0, aln

    def test_glocal5(self):
        alns = self.f('AATCAAG', 'AATGAATGAGTCAATG', method=self.method,
                      matrix=DNAFULL, max_hits=None)
        assert len(alns) == 2, alns
        aln1 = AlignmentResult(seq1=b'AATCAA-G', seq2=b'AATGAATG', start1=0,
                               start2=0, end1=7, end2=8,
                               n_gaps1=1, n_gaps2=0, n_mismatches=1,
                               score=19.0)
        aln2 = AlignmentResult(seq1=b'AATCAA-G', seq2=b'AGTCAATG', start1=0,
                               start2=8, end1=7, end2=16,
                               n_gaps1=1, n_gaps2=0, n_mismatches=1,
                               score=19.0)
        assert aln1 in alns, alns
        assert aln2 in alns, alns


class TestGlobalC(TestGlobalPy):

    @property
    def f(self):
        return caligner


class TestLocalC(TestLocalPy):

    @property
    def f(self):
        return caligner


class TestGlocalC(TestGlocalPy):

    @property
    def f(self):
        return caligner


class TestGlobalCFEC(TestGlobalCFEPy):

    @property
    def f(self):
        return caligner


if __name__ == '__main__':
    unittest.main()
