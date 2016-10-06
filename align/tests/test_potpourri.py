# -*- coding: utf-8 -*-

import unittest

from align import aligner
from align.matrix import DNAFULL


class TestPotpourri(unittest.TestCase):

    def test_global1(self):
        aln = aligner('WW', 'WEW', method='global')
        assert list(aln.seq1) == ['W', '-', 'W']
        assert list(aln.seq2) == ['W', 'E', 'W']

    def test_global2(self):
        aln = aligner('WW', 'WEW', method='global', gap_open=-100)
        assert list(aln.seq1) == ['W', '-', 'W']
        assert list(aln.seq2) == ['W', 'E', 'W']

    def test_global3(self):
        aln = aligner('A', 'A', method='global', gap_open=-7)
        assert list(aln.seq1) == ['A']
        assert list(aln.seq2) == ['A']

    def test_global4(self):
        aln = aligner('R', 'K', method='global', gap_open=-7)
        assert list(aln.seq1) == ['R']
        assert list(aln.seq2) == ['K']

    def test_global5(self):
        aln = aligner('R', 'AR', method='global', gap_open=-7)
        assert list(aln.seq1) == ['-', 'R'], (aln)
        assert list(aln.seq2) == ['A', 'R']

    def test_global6(self):
        aln = aligner('AR', 'R', method='global', gap_open=-7)
        assert list(aln.seq1) == ['A', 'R']
        assert list(aln.seq2) == ['-', 'R']

    def test_global7(self):
        aln = aligner('AR', 'RA', method='global', gap_open=-7)
        assert list(aln.seq1) == ['A', 'R']
        assert list(aln.seq2) == ['R', 'A']

    def test_global8(self):
        aln = aligner('AR', 'RA', method='global', gap_open=-3)
        assert list(aln.seq1) == ['A', 'R', '-']
        assert list(aln.seq2) == ['-', 'R', 'A']

    def test_global9(self):
        aln = aligner('RAR', 'RR', method='global', gap_open=-3)
        assert list(aln.seq2) == ['R', '-', 'R']
        assert list(aln.seq1) == ['R', 'A', 'R']

    def test_global10(self):
        aln = aligner('RAR', 'RR', method='global', gap_open=-10)
        assert list(aln.seq2) == ['R', '-', 'R']
        assert list(aln.seq1) == ['R', 'A', 'R']

    def test_global11(self):
        aln = aligner('RAAR', 'RR', method='global', gap_open=-5)
        assert list(aln.seq1) == ['R', 'A', 'A', 'R']
        assert list(aln.seq2) == ['R', '-', '-', 'R']

    def test_global12(self):
        aln = aligner('RLR', 'RER', method='global', gap_open=-9)
        assert list(aln.seq1) == ['R', 'L', 'R']
        assert list(aln.seq2) == ['R', 'E', 'R']

    def test_global13(self):
        aln = aligner('RLR', 'RER', method='global', gap_open=-1)
        assert list(aln.seq1) == ['R', 'L', '-', 'R']
        assert list(aln.seq2) == ['R', '-', 'E', 'R']

    def test_global14(self):
        aln = aligner('RLR', 'REER', method='global', gap_open=-1)
        assert list(aln.seq1) == ['R', 'L', '-', '-', 'R']
        assert list(aln.seq2) == ['R', '-', 'E', 'E', 'R']

    def test_global15(self):
        aln = aligner('AGEBAM', 'AGEBAMAM', method='global', gap_open=-6)
        assert list(aln.seq1) == ['A', 'G', 'E', 'B', 'A', 'M', '-', '-']
        assert list(aln.seq2) == ['A', 'G', 'E', 'B', 'A', 'M', 'A', 'M']

    def test_global16(self):
        aln = aligner('CPELIRKNCANTH', 'PREKRLICAN', method='global',
                      gap_open=-0.5)
        assert list(aln.seq1) == ['C', 'P', '-', 'E', '-', '-', 'L', 'I', 'R',
                                  'K', 'N', 'C', 'A', 'N', 'T', 'H']
        assert list(aln.seq2) == ['-', 'P', 'R', 'E', 'K', 'R', 'L', 'I', '-',
                                  '-', '-', 'C', 'A', 'N', '-', '-']

    def test_global17(self):
        aln = aligner('CPEL', 'PREK', method='global', gap_open=-5)
        assert list(aln.seq1) == ['C', 'P', '-', 'E', 'L']
        assert list(aln.seq2) == ['-', 'P', 'R', 'E', 'K']

    def test_global18(self):
        aln = aligner('RLRR', 'RRER', method='global', gap_open=-1)
        assert list(aln.seq1) == ['R', 'L', 'R', '-', 'R']
        assert list(aln.seq2) == ['R', '-', 'R', 'E', 'R']

    def test_global19(self):
        aln = aligner('TAAT', 'TAATTC', method='global', matrix=DNAFULL)
        assert list(aln.seq1) == ['T', 'A', 'A', 'T', '-', '-']
        assert list(aln.seq2) == ['T', 'A', 'A', 'T', 'T', 'C']

    def test_global20(self):
        aln = aligner('WR', 'WRR', method='global')
        assert list(aln.seq2) == ['W', 'R', 'R']
        assert list(aln.seq1) == ['W', 'R', '-']

    def test_global21(self):
        aln = aligner('AIP', 'AP', method='global')
        assert list(aln.seq1) == ['A', 'I', 'P']
        assert list(aln.seq2) == ['A', '-', 'P']

    def test_global22(self):
        aln = aligner('PAA', 'PA', method='global')
        assert list(aln.seq1) == ['P', 'A', 'A']
        assert list(aln.seq2) == ['P', 'A', '-']

    def test_global23(self):
        aln = aligner('TAATTC', 'TAAT', method='global', matrix=DNAFULL,
                      gap_open=-10, gap_extend=-1)
        assert list(aln.seq1) == ['T', 'A', 'A', 'T', 'T', 'C']
        assert list(aln.seq2) == ['T', 'A', 'A', 'T', '-', '-']

        # global_cfe
    def test_global_cfe1(self):
        aln = aligner('TAAT', 'TAATTC', method='global_cfe', matrix=DNAFULL)
        assert list(aln.seq1) == ['T', 'A', 'A', 'T', '-', '-']
        assert list(aln.seq2) == ['T', 'A', 'A', 'T', 'T', 'C']

    def test_global_cfe2(self):
        aln = aligner('TCTAAT', 'TAAT', method='global_cfe', matrix=DNAFULL)
        assert list(aln.seq2) == ['-', '-', 'T', 'A', 'A', 'T']
        assert list(aln.seq1) == ['T', 'C', 'T', 'A', 'A', 'T']

    def test_global_cfe3(self):
        aln = aligner('PAA', 'PA', method='global_cfe')
        assert list(aln.seq1) == ['P', 'A', 'A']
        assert list(aln.seq2) == ['P', 'A', '-']

        # local
    def test_local1(self):
        aln = aligner('TCTAAT', 'TAAT', method='local', matrix=DNAFULL)
        assert list(aln.seq2) == ['T', 'A', 'A', 'T']
        assert list(aln.seq1) == ['T', 'A', 'A', 'T']

    def test_local2(self):
        aln = aligner('TCTAAT', 'TAATCT', method='local', matrix=DNAFULL)
        assert list(aln.seq2) == ['T', 'A', 'A', 'T']
        assert list(aln.seq1) == ['T', 'A', 'A', 'T']

    def test_local3(self):
        aln = aligner('A', 'A', method='local')
        assert list(aln.seq1) == ['A']
        assert list(aln.seq2) == ['A']

    def test_local4(self):
        aln = aligner('RA', 'AR', method='local')
        assert list(aln.seq1) == ['R']
        assert list(aln.seq2) == ['R']

    def test_local5(self):
        aln = aligner('RRR', 'RR', method='local')
        assert list(aln.seq1) == ['R', 'R']
        assert list(aln.seq2) == ['R', 'R']

    def test_local6(self):
        aln = aligner('PYNCHAN', 'YNCH', method='local')
        assert list(aln.seq1) == ['Y', 'N', 'C', 'H']
        assert list(aln.seq2) == ['Y', 'N', 'C', 'H']

    def test_local7(self):
        aln = aligner('AIP', 'AP', method='local')
        assert list(aln.seq1) == ['P']
        assert list(aln.seq2) == ['P']

    def test_local8(self):
        aln = aligner('PAA', 'PA', method='local')
        assert list(aln.seq1) == ['P', 'A']
        assert list(aln.seq2) == ['P', 'A']

    def test_glocal1(self):
        aln = aligner('AAATAATAAA', 'TAAT', method='glocal', matrix=DNAFULL)
        assert list(aln.seq2) == ['T', 'A', 'A', 'T']
        assert list(aln.seq1) == ['T', 'A', 'A', 'T']

    def test_glocal2(self):
        aln = aligner('AAATAATAAA', 'TATAT', method='glocal', gap_open=-1,
                      matrix=DNAFULL)
        assert list(aln.seq1) == ['T', 'A', '-', 'A', 'T']
        assert list(aln.seq2) == ['T', 'A', 'T', 'A', 'T']

    def test_glocal3(self):
        aln = aligner('TATATAAA', 'CCTATAT', method='glocal', gap_open=-1,
                      matrix=DNAFULL)
        assert (aln.seq1, aln.seq2) == ('--TATAT', 'CCTATAT'), aln

    def test_glocal4(self):
        aln = aligner('CCTATAT', 'TATATAAA', method='glocal', gap_open=-1,
                      matrix=DNAFULL)
        assert list(aln.seq2) == ['-', '-', 'T', 'A', 'T', 'A', 'T']
        assert list(aln.seq1) == ['C', 'C', 'T', 'A', 'T', 'A', 'T']


if __name__ == '__main__':
    unittest.main()
