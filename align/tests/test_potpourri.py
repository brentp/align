# -*- coding: utf-8 -*-

import unittest

from align import aligner
from align.matrix import DNAFULL


class TestPotpourri(unittest.TestCase):

    def test_all(self):
        # global
        aln = aligner('WW', 'WEW', method='global')
        assert list(aln.seq1) == ['W', '-', 'W']
        assert list(aln.seq2) == ['W', 'E', 'W']
        aln = aligner('WW', 'WEW', method='global', gap_open=-100)
        assert list(aln.seq1) == ['W', '-', 'W']
        assert list(aln.seq2) == ['W', 'E', 'W']
        aln = aligner('A', 'A', method='global', gap_open=-7)
        assert list(aln.seq1) == ['A']
        assert list(aln.seq2) == ['A']
        aln = aligner('R', 'K', method='global', gap_open=-7)
        assert list(aln.seq1) == ['R']
        assert list(aln.seq2) == ['K']
        aln = aligner('R', 'AR', method='global', gap_open=-7)
        assert list(aln.seq1) == ['-', 'R'], (aln)
        assert list(aln.seq2) == ['A', 'R']
        aln = aligner('AR', 'R', method='global', gap_open=-7)
        assert list(aln.seq1) == ['A', 'R']
        assert list(aln.seq2) == ['-', 'R']
        aln = aligner('AR', 'RA', method='global', gap_open=-7)
        assert list(aln.seq1) == ['A', 'R']
        assert list(aln.seq2) == ['R', 'A']
        aln = aligner('AR', 'RA', method='global', gap_open=-3)
        assert list(aln.seq1) == ['A', 'R', '-']
        assert list(aln.seq2) == ['-', 'R', 'A']
        aln = aligner('RAR', 'RR', method='global', gap_open=-3)
        assert list(aln.seq2) == ['R', '-', 'R']
        assert list(aln.seq1) == ['R', 'A', 'R']
        aln = aligner('RAR', 'RR', method='global', gap_open=-10)
        assert list(aln.seq2) == ['R', '-', 'R']
        assert list(aln.seq1) == ['R', 'A', 'R']
        aln = aligner('RAAR', 'RR', method='global', gap_open=-5)
        assert list(aln.seq1) == ['R', 'A', 'A', 'R']
        assert list(aln.seq2) == ['R', '-', '-', 'R']
        aln = aligner('RLR', 'RER', method='global', gap_open=-9)
        assert list(aln.seq1) == ['R', 'L', 'R']
        assert list(aln.seq2) == ['R', 'E', 'R']
        aln = aligner('RLR', 'RER', method='global', gap_open=-1)
        assert list(aln.seq1) == ['R', 'L', '-', 'R']
        assert list(aln.seq2) == ['R', '-', 'E', 'R']
        aln = aligner('RLR', 'REER', method='global', gap_open=-1)
        assert list(aln.seq1) == ['R', 'L', '-', '-', 'R']
        assert list(aln.seq2) == ['R', '-', 'E', 'E', 'R']
        aln = aligner('AGEBAM', 'AGEBAMAM', method='global', gap_open=-6)
        assert list(aln.seq1) == ['A', 'G', 'E', 'B', 'A', 'M', '-', '-']
        assert list(aln.seq2) == ['A', 'G', 'E', 'B', 'A', 'M', 'A', 'M']
        aln = aligner('CPELIRKNCANTH', 'PREKRLICAN', method='global',
                      gap_open=-0.5)
        assert list(aln.seq1) == ['C', 'P', '-', 'E', '-', '-', 'L', 'I', 'R',
                                  'K', 'N', 'C', 'A', 'N', 'T', 'H']
        assert list(aln.seq2) == ['-', 'P', 'R', 'E', 'K', 'R', 'L', 'I', '-',
                                  '-', '-', 'C', 'A', 'N', '-', '-']
        aln = aligner('CPEL', 'PREK', method='global', gap_open=-5)
        assert list(aln.seq1) == ['C', 'P', '-', 'E', 'L']
        assert list(aln.seq2) == ['-', 'P', 'R', 'E', 'K']
        aln = aligner('RLRR', 'RRER', method='global', gap_open=-1)
        assert list(aln.seq1) == ['R', 'L', 'R', '-', 'R']
        assert list(aln.seq2) == ['R', '-', 'R', 'E', 'R']
        aln = aligner('TAAT', 'TAATTC', method='global', matrix=DNAFULL)
        assert list(aln.seq1) == ['T', 'A', 'A', 'T', '-', '-']
        assert list(aln.seq2) == ['T', 'A', 'A', 'T', 'T', 'C']
        # global_cfe
        aln = aligner('TAAT', 'TAATTC', method='global_cfe', matrix=DNAFULL)
        assert list(aln.seq1) == ['T', 'A', 'A', 'T', '-', '-']
        assert list(aln.seq2) == ['T', 'A', 'A', 'T', 'T', 'C']
        aln = aligner('TCTAAT', 'TAAT', method='global_cfe', matrix=DNAFULL)
        assert list(aln.seq2) == ['-', '-', 'T', 'A', 'A', 'T']
        assert list(aln.seq1) == ['T', 'C', 'T', 'A', 'A', 'T']
        # local
        aln = aligner('TCTAAT', 'TAAT', method='local', matrix=DNAFULL)
        assert list(aln.seq2) == ['T', 'A', 'A', 'T']
        assert list(aln.seq1) == ['T', 'A', 'A', 'T']
        aln = aligner('TCTAAT', 'TAATCT', method='local', matrix=DNAFULL)
        assert list(aln.seq2) == ['T', 'A', 'A', 'T']
        assert list(aln.seq1) == ['T', 'A', 'A', 'T']
        # glocal
        aln = aligner('AAATAATAAA', 'TAAT', method='glocal', matrix=DNAFULL)
        assert list(aln.seq2) == ['T', 'A', 'A', 'T']
        assert list(aln.seq1) == ['T', 'A', 'A', 'T']
        aln = aligner('AAATAATAAA', 'TATAT', method='glocal', gap_open=-1,
                      matrix=DNAFULL)
        assert list(aln.seq1) == ['T', 'A', '-', 'A', 'T']
        assert list(aln.seq2) == ['T', 'A', 'T', 'A', 'T']

        aln = aligner('TATATAAA', 'CCTATAT', method='glocal', gap_open=-1,
                      matrix=DNAFULL)
        assert (aln.seq1, aln.seq2) == ('--TATAT', 'CCTATAT'), aln

        aln = aligner('CCTATAT', 'TATATAAA', method='glocal', gap_open=-1,
                      matrix=DNAFULL)
        assert list(aln.seq2) == ['-', '-', 'T', 'A', 'T', 'A', 'T']
        assert list(aln.seq1) == ['C', 'C', 'T', 'A', 'T', 'A', 'T']
        # old
        aln = aligner('A', 'A', method='local')
        assert list(aln.seq1) == ['A']
        assert list(aln.seq2) == ['A']
        aln = aligner('RA', 'AR', method='local')
        assert list(aln.seq1) == ['R']
        assert list(aln.seq2) == ['R']
        aln = aligner('RRR', 'RR', method='local')
        assert list(aln.seq1) == ['R', 'R']
        assert list(aln.seq2) == ['R', 'R']
        aln = aligner('WR', 'WRR', method='global')
        assert list(aln.seq2) == ['W', 'R', 'R']
        assert list(aln.seq1) == ['W', 'R', '-']
        aln = aligner('PYNCHAN', 'YNCH', method='local')
        assert list(aln.seq1) == ['Y', 'N', 'C', 'H']
        assert list(aln.seq2) == ['Y', 'N', 'C', 'H']
        aln = aligner('AIP', 'AP', method='local')
        assert list(aln.seq1) == ['P']
        assert list(aln.seq2) == ['P']
        aln = aligner('AIP', 'AP', method='global')
        assert list(aln.seq1) == ['A', 'I', 'P']
        assert list(aln.seq2) == ['A', '-', 'P']
        aln = aligner('PAA', 'PA', method='local')
        assert list(aln.seq1) == ['P', 'A']
        assert list(aln.seq2) == ['P', 'A']
        aln = aligner('PAA', 'PA', method='global')
        assert list(aln.seq1) == ['P', 'A', 'A']
        assert list(aln.seq2) == ['P', 'A', '-']
        aln = aligner('PAA', 'PA', method='global_cfe')
        assert list(aln.seq1) == ['P', 'A', 'A']
        assert list(aln.seq2) == ['P', 'A', '-']
        aln = aligner('TAATTC', 'TAAT', method='global', matrix=DNAFULL,
                      gap_open=-10, gap_extend=-1)
        assert list(aln.seq1) == ['T', 'A', 'A', 'T', 'T', 'C']
        assert list(aln.seq2) == ['T', 'A', 'A', 'T', '-', '-']


if __name__ == '__main__':
    unittest.main()
