# -*- coding: utf-8 -*-

import unittest

from align.matrix import BLOSUM62
from align import aligner


class TestAlignGlobal(unittest.TestCase):

    method = "global"

    def test_all(self):
        r = aligner("CELECANTH", "PELICAN", method=self.method)
        assert (r.seq1, r.seq2) == ('CELECANTH', 'PELICAN--'), r
        r = aligner("PELICAN", "CELECANTH", method=self.method)
        assert (r.seq1, r.seq2) == ('PELICAN--', 'CELECANTH')

    def test_score(self):
        s0 = "AGEBANAN"
        s1 = "ACEBAN"
        r = aligner(s0, s1, gap_extend=-1, gap_open=-2, matrix=BLOSUM62,
                    method=self.method)
        assert (r.seq1, r.seq2) == ('AGEBANAN', 'ACEBAN--')


class TestGlobalCFE(unittest.TestCase):

    method = 'global_cfe'

    def test_it(self):
        r = aligner(
            'AAAAAAAAAAAAACCTGCGCCCCAAAAAAAAAAAAAAAAAAAA',
            'CCTGCGCACCCC', method=self.method)
        assert (r.seq1, r.seq2) == \
            ('AAAAAAAAAAAAACCTGCGC-CCCAAAAAAAAAAAAAAAAAAAA',
             '-------------CCTGCGCACCCC-------------------')


if __name__ == '__main__':
    unittest.main()
