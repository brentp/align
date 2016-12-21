++++++++++++++++++++++++++++++++++++++++
Align: polite, proper sequence alignment
++++++++++++++++++++++++++++++++++++++++

    :Authors: Marcin Cieślik, Brent Pedersen (brentp), Wibowo Arindrarto (bow)
    :Email: (marcin), bpederse@gmail.com, bow@bow.web.id
    :License: BSD

.. contents ::


About
=====
Various packages implement sequence alignment algorithms with various levels of
success. This package is an attempt to handle the sequence alignment properly,
including edge-cases.


Usage
=====

usage will change. currently ::

    >>> from align.matrix import DNAFULL
    >>> from align import aligner

    >>> aligner('WW','WEW', method= 'global')
    ('W-W', 'WEW')

    >>> aligner('WW','WEWWEW', method= 'glocal')
    ('WW', 'WW')


    >>> aligner('TAATTC', 'TAAT', method='global', matrix=DNAFULL, gap_open=-10, gap_extend=-1)
    ('TAATTC', 'TAAT--')

    >>> aligner('PYNCHAN', 'YNCH', method='local')
    ('YNCH', 'YNCH')

    >>> a, b = aligner('AAAAAAAAAAAAACCTGCGCCCCAAAAAAAAAAAAAAAAAAAA', 'CCTGCGCACCCC', method='global_cfe')
    >>> print "%s\n%s" % (a, b)
    AAAAAAAAAAAAACCTGCGC-CCCAAAAAAAAAAAAAAAAAAAA
    -------------CCTGCGCACCCC-------------------
