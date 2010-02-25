import unittest
from align.matrix import DNAFULL
from align import aligner

class TestPotpourri(unittest.TestCase):

    def test_all(self):
        # global
        a, b = aligner('WW','WEW', method= 'global')
        assert list(a) == ['W', '-', 'W']
        assert list(b) == ['W', 'E', 'W']
        a, b = aligner('WW','WEW', method= 'global', gap_open=-100)
        assert list(a) == ['W', '-', 'W']
        assert list(b) == ['W', 'E', 'W']
        a,b = aligner('A', 'A', method='global', gap_open=-7)
        assert list(a) == ['A']
        assert list(b) == ['A']
        a,b = aligner('R','K', method='global', gap_open=-7)
        assert list(a) == ['R']
        assert list(b) == ['K']
        a,b = aligner('R','AR', method='global', gap_open=-7)
        assert list(a) == ['-','R'], (a, b)
        assert list(b) == ['A','R']
        a,b = aligner('AR','R', method='global', gap_open=-7)
        assert list(a) == ['A','R']
        assert list(b) == ['-','R']
        a,b = aligner('AR','RA', method='global', gap_open=-7)
        assert list(a) == ['A','R']
        assert list(b) == ['R','A']
        a,b = aligner('AR','RA', method='global', gap_open=-3)
        assert list(a) == ['A', 'R', '-']
        assert list(b) == ['-', 'R', 'A']
        a,b = aligner('RAR','RR', method='global', gap_open=-3)
        assert list(b) == ['R', '-', 'R']
        assert list(a) == ['R', 'A', 'R']
        a,b = aligner('RAR','RR', method='global', gap_open=-10)
        assert list(b) == ['R', '-', 'R']
        assert list(a) == ['R', 'A', 'R']
        a,b = aligner('RAAR','RR', method='global', gap_open=-5)
        assert list(a) == ['R', 'A', 'A', 'R']
        assert list(b) == ['R', '-', '-', 'R']
        a,b = aligner('RLR','RER', method='global', gap_open=-9)
        assert list(a) == ['R', 'L', 'R']
        assert list(b) == ['R', 'E', 'R']
        a,b = aligner('RLR','RER', method='global', gap_open=-1)
        assert list(a) == ['R', 'L', '-', 'R']
        assert list(b) == ['R', '-', 'E', 'R']
        a,b = aligner('RLR','REER', method='global', gap_open=-1)
        assert list(a) == ['R', 'L', '-', '-', 'R']
        assert list(b) == ['R', '-', 'E', 'E', 'R']
        a, b = aligner('AGEBAM', 'AGEBAMAM', method='global', gap_open=-6)
        assert list(a) == ['A', 'G', 'E', 'B', 'A', 'M', '-', '-']
        assert list(b) == ['A', 'G', 'E', 'B', 'A', 'M', 'A', 'M']
        a, b= aligner('CPELIRKNCANTH', 'PREKRLICAN', method='global', gap_open=-0.5)
        assert list(a) == ['C', 'P', '-', 'E', '-', '-', 'L', 'I', 'R', 'K', 'N', 'C', 'A', 'N', 'T', 'H']
        assert list(b) == ['-', 'P', 'R', 'E', 'K', 'R', 'L', 'I', '-', '-', '-', 'C', 'A', 'N', '-', '-']
        a, b= aligner('CPEL', 'PREK', method='global', gap_open=-6.)
        # assert list(a) == ['C', 'P', 'E', 'L']
        # assert list(b) == ['P', 'R', 'E', 'K']
        a, b= aligner('CPEL', 'PREK', method='global', gap_open=-5)
        assert list(a) == ['C', 'P', '-','E', 'L']
        assert list(b) == ['-','P', 'R', 'E', 'K']
        a,b = aligner('RLRR','RRER', method='global', gap_open=-1)
        assert list(a) == ['R', 'L', 'R', '-', 'R']
        assert list(b) == ['R', '-', 'R', 'E', 'R']
        a, b = aligner('TAAT', 'TAATTC', method='global', matrix=DNAFULL)
        assert list(a) == ['T', 'A', 'A', 'T', '-', '-']
        assert list(b) == ['T', 'A', 'A', 'T', 'T', 'C']
        # global_cfe
        a, b = aligner('TAAT', 'TAATTC', method='global_cfe', matrix=DNAFULL)
        assert list(a) == ['T', 'A', 'A', 'T', '-', '-']
        assert list(b) == ['T', 'A', 'A', 'T', 'T', 'C']
        a, b = aligner('TCTAAT', 'TAAT', method='global_cfe', matrix=DNAFULL)
        assert list(b) == ['-', '-','T', 'A', 'A', 'T' ]
        assert list(a) == ['T', 'C','T', 'A', 'A', 'T' ]
        # local
        a, b = aligner('TCTAAT', 'TAAT', method='local', matrix=DNAFULL)
        assert list(b) == ['T', 'A', 'A', 'T' ]
        assert list(a) == ['T', 'A', 'A', 'T' ]
        a, b = aligner('TCTAAT', 'TAATCT', method='local', matrix=DNAFULL)
        assert list(b) == ['T', 'A', 'A', 'T' ]
        assert list(a) == ['T', 'A', 'A', 'T' ]
        # glocal
        a, b = aligner('AAATAATAAA', 'TAAT', method='glocal', matrix=DNAFULL)
        assert list(b) == ['T', 'A', 'A', 'T' ]
        assert list(a) == ['T', 'A', 'A', 'T' ]
        a, b = aligner('AAATAATAAA', 'TATAT', method='glocal', gap_open=-1, matrix=DNAFULL)
        assert list(a) == ['T', 'A', '-', 'A', 'T' ]
        assert list(b) == ['T', 'A', 'T', 'A', 'T' ]
        a, b = aligner('TATATAAA', 'CCTATAT', method='glocal', gap_open=-1, matrix=DNAFULL)
        assert list(a) == ['-', '-', 'T', 'A', 'T', 'A', 'T' ]
        assert list(b) == ['C', 'C', 'T', 'A', 'T', 'A', 'T' ]
        a, b = aligner('CCTATAT', 'TATATAAA',method='glocal', gap_open=-1, matrix=DNAFULL)
        assert list(b) == ['-', '-', 'T', 'A', 'T', 'A', 'T' ]
        assert list(a) == ['C', 'C', 'T', 'A', 'T', 'A', 'T' ]
        # old
        a, b = aligner('A', 'A', method ='local')
        assert list(a) == ['A']
        assert list(b) == ['A']
        a, b = aligner('RA', 'AR', method ='local')
        assert list(a) == ['R']
        assert list(b) == ['R']
        a, b = aligner('RRR', 'RR', method ='local')
        assert list(a) == ['R', 'R']
        assert list(b) == ['R', 'R']
        a, b = aligner('WR', 'WRR', method ='global')
        assert list(b) == ['W', 'R', 'R']
        assert list(a) == ['W', 'R', '-']
        a,b = aligner('PYNCHAN', 'YNCH', method='local')
        assert list(a) == ['Y', 'N', 'C', 'H']
        assert list(b) == ['Y', 'N', 'C', 'H']
        a, b = aligner('AIP', 'AP', method='local')
        assert list(a) == ['P']
        assert list(b) == ['P']
        a, b = aligner('AIP', 'AP', method='global')
        assert list(a) == ['A','I','P']
        assert list(b) == ['A','-','P']
        a, b = aligner('PAA', 'PA', method='local')
        assert list(a) == ['P','A']
        assert list(b) == ['P','A']
        a, b = aligner('PAA', 'PA', method='global')
        assert list(a) == ['P','A','A']
        assert list(b) == ['P','A','-']
        a, b = aligner('PAA', 'PA', method='global_cfe')
        assert list(a) == ['P','A','A']
        assert list(b) == ['P','A','-']
        a, b = aligner('TAATTC', 'TAAT', method='global', matrix=DNAFULL, gap_open=-10, gap_extend=-1)
        assert list(a) == ['T', 'A', 'A', 'T', 'T', 'C']
        assert list(b) == ['T', 'A', 'A', 'T', '-', '-']

if __name__ == '__main__':
    unittest.main()
