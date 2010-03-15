import unittest
#from align.matrix import DNAFULL, BLOSUM62
from align import aligner

class TestAlignGlobal(unittest.TestCase):
    method = "global"

    def test_all(self):
        r = aligner("CELECANTH", "PELICAN", method=self.method)
        assert r ==  ('CELECANTH', 'PELICAN--'), r
        r = aligner("PELICAN", "CELECANTH", method=self.method)
        assert r == ('PELICAN--', 'CELECANTH')

    def test_score(self):
        s0 = "AGEBANAN"
        s1 = "ACEBAN"
        r = aligner(s0, s1, gap_extend=-1, gap_open=-2, matrix="BLOSUM62",
                    method=self.method)
        assert r == ('AGEBANAN', 'ACEBAN--')



if __name__ == '__main__':
    unittest.main()
