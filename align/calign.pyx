# -*- coding: utf-8 -*-

from collections import namedtuple

import numpy as np
import os.path as op
import sys

cimport numpy as np
from libc.string cimport strlen
from cpython.mem cimport PyMem_Malloc, PyMem_Free

from .matrix import BLOSUM62, DNAFULL


METHODS = {"global": 0, "local": 1, "glocal": 2, "global_cfe": 3}

ctypedef np.int_t DTYPE_INT
ctypedef np.uint_t DTYPE_UINT
ctypedef np.float32_t DTYPE_FLOAT


cdef inline DTYPE_FLOAT max3(DTYPE_FLOAT a, DTYPE_FLOAT b, DTYPE_FLOAT c):
    """Returns a maximum of three numpy 32-bit floats."""
    if c > b:
        return c if c > a else a
    return b if b > a else a


cdef inline DTYPE_FLOAT max2(DTYPE_FLOAT a, DTYPE_FLOAT b):
    """Returns a maximum of two numpy 32-bit floats."""
    return b if b > a else a


cdef DTYPE_FLOAT[:, :] make_cmatrix(dict pymatrix):
    """Transforms the given dictionary scoring matrix into a numpy matrix.

    so here, we read a matrix in the NCBI format and put
    it into a numpy array. so the score for a 'C' changing
    to an 'A' is stored in the matrix as:
        mat[ord('C'), ord('A')] = score
    as such, it's a direct array lookup from each pair in the alignment
    to a score. this makes it very fast. the cost is in terms of space.
    though it's usually less than 100*100.

    """
    cdef:
        DTYPE_INT size = sorted([ord(c) for c in pymatrix.keys()]).pop() + 1
        DTYPE_FLOAT score
        np.int8_t c1, c2
        np.ndarray[DTYPE_FLOAT, ndim=2] cmatrix = np.zeros((size, size),
                                                           dtype=np.float32)

    for char1 in pymatrix.keys():
        for char2, score in pymatrix[char1].items():
            c1 = ord(char1)
            c2 = ord(char2)
            cmatrix[c1, c2] = score

    return cmatrix


cdef:
    DTYPE_FLOAT[:, :] m_BLOSUM62 = make_cmatrix(BLOSUM62)
    DTYPE_FLOAT[:, :] m_DNAFULL = make_cmatrix(DNAFULL)
    unsigned char GAP_CHAR = '-'

cdef enum:
    NONE = 0, LEFT = 1, UP = 2, DIAG = 3


cdef tuple caligner(
    const unsigned char* seqj, const unsigned char* seqi, const int imethod,
    const DTYPE_FLOAT gap_open, const DTYPE_FLOAT gap_extend, const DTYPE_FLOAT gap_double,
    const DTYPE_FLOAT[:, :] amatrix, const bint flipped):

    cdef:
        unsigned char* align_j
        unsigned char* align_i
        unsigned char ci, cj
        size_t max_j = strlen(<char *>seqj)
        size_t max_i = strlen(<char *>seqi)
        size_t align_counter = 0
        size_t i = 1, j = 1
        DTYPE_UINT p, matrix_max, row_max, col_max
        DTYPE_FLOAT diag_score, left_score, up_score, max_score
        np.ndarray[DTYPE_FLOAT, ndim=2] agap_i = np.empty((max_i + 1, max_j + 1), dtype=np.float32)
        np.ndarray[DTYPE_FLOAT, ndim=2] agap_j = np.empty((max_i + 1, max_j + 1), dtype=np.float32)
        np.ndarray[DTYPE_FLOAT, ndim=2] score = np.zeros((max_i + 1, max_j + 1), dtype=np.float32)
        np.ndarray[DTYPE_UINT, ndim=2] pointer = np.zeros((max_i + 1, max_j + 1), dtype=np.uint)

    agap_i.fill(-np.inf)
    agap_j.fill(-np.inf)

    # START HERE:
    if imethod == 0:
        pointer[0, 1:] = LEFT
        pointer[1:, 0] = UP
        score[0, 1:] = gap_open + gap_extend * np.arange(0, max_j, dtype=np.float32)
        score[1:, 0] = gap_open + gap_extend * np.arange(0, max_i, dtype=np.float32)
    elif imethod == 3:
        pointer[0, 1:] = LEFT
        pointer[1:, 0] = UP
    elif imethod == 2:
        pointer[0, 1:] = LEFT
        score[0, 1:] = gap_open + gap_extend * np.arange(0, max_j, dtype=np.float32)

    for i in range(1, max_i + 1):
        ci = seqi[i - 1]
        for j in range(1, max_j + 1):
            cj = seqj[j - 1]
            # agap_i
            agap_i[i,j] = max3(
                         score[i, j - 1] + gap_open,
                         agap_i[i, j - 1] + gap_extend,
                         agap_j[i, j - 1] + gap_double)
            # agap_j
            agap_j[i,j] = max3(
                         score[i - 1, j] + gap_open,
                         agap_j[i - 1, j] + gap_extend,
                         agap_i[i - 1, j] + gap_double)
            # score
            diag_score = score[i - 1, j - 1] + amatrix[ci, cj]
            left_score = agap_i[i, j]
            up_score   = agap_j[i, j]
            max_score = max3(diag_score, up_score, left_score)

            score[i, j] = max2(0, max_score) if imethod == 1 else max_score

            if imethod == 1:
                if score[i,j] == 0:
                    pass # point[i,j] = NONE
                elif max_score == diag_score:
                    pointer[i,j] = DIAG
                elif max_score == up_score:
                    pointer[i,j] = UP
                elif max_score == left_score:
                    pointer[i,j] = LEFT
            elif imethod == 2:
                # In a semi-global alignment we want to consume as much as
                # possible of the longer sequence.
                if max_score == up_score:
                    pointer[i,j] = UP
                elif max_score == diag_score:
                    pointer[i,j] = DIAG
                elif max_score == left_score:
                    pointer[i,j] = LEFT
            else:
                # global
                if max_score == up_score:
                    pointer[i,j] = UP
                elif max_score == left_score:
                    pointer[i,j] = LEFT
                else:
                    pointer[i,j] = DIAG

    if imethod == 1:  # local
        # max anywhere
        matrix_max = score.argmax()
        i, j = np.unravel_index(matrix_max, (score.shape[0], score.shape[1]))
    elif imethod == 2:  # glocal
        # max in last col
        i, j = (score[:,-1].argmax(), max_j)
    elif imethod == 3:  # global_cfe
        # from i,j to max(max(last row), max(last col)) for free
        row_max, col_idx = score[-1].max(), score[-1].argmax()
        col_max, row_idx = score[:, -1].max(), score[:, -1].argmax()
        if row_max > col_max:
            pointer[-1,col_idx+1:] = LEFT
        else:
            pointer[row_idx+1:,-1] = UP

    seqlen = max_i + max_j
    align_i = <unsigned char *>PyMem_Malloc(seqlen * sizeof(unsigned char))
    align_j = <unsigned char *>PyMem_Malloc(seqlen * sizeof(unsigned char))

    p = pointer[i, j]
    while p != NONE:
        if p == DIAG:
            i -= 1
            j -= 1
            align_j[align_counter] = seqj[j]
            align_i[align_counter] = seqi[i]
        elif p == LEFT:
            j -= 1
            align_j[align_counter] = seqj[j]
            align_i[align_counter] = GAP_CHAR
        elif p == UP:
            i -= 1
            align_j[align_counter] = GAP_CHAR
            align_i[align_counter] = seqi[i]
        else:
            raise Exception('wtf!:pointer: %i', p)
        align_counter += 1
        p = pointer[i, j]

    if flipped:
        seq1 = bytes(align_i[:align_counter][::-1])
        seq2 = bytes(align_j[:align_counter][::-1])
    else:
        seq1 = bytes(align_j[:align_counter][::-1])
        seq2 = bytes(align_i[:align_counter][::-1])

    PyMem_Free(align_i)
    PyMem_Free(align_j)

    return seq1, seq2


def aligner(seqj, seqi, method='global', gap_open=-7, gap_extend=-7,
            gap_double=-7, matrix="BLOSUM62"):
    '''Calculates the alignment of two sequences.

    The supported 'methods' are:
        * 'global' for a global Needleman-Wunsh algorithm
        * 'local' for a local Smith-Waterman alignment
        * 'global_cfe' for a global alignment with cost-free ends
        * 'glocal' for an alignment which is 'global' only with respect to
          the shorter sequence (also known as a 'semi-global' alignment)

    Returns the aligned (sub)sequences as character arrays.

    Gotoh, O. (1982). J. Mol. Biol. 162, 705-708.
    Needleman, S. & Wunsch, C. (1970). J. Mol. Biol. 48(3), 443-53.
    Smith, T.F. & Waterman M.S. (1981). J. Mol. Biol. 147, 195-197.

    Arguments:

        - seqj (``sequence``) First aligned iterable object of symbols.
        - seqi (``sequence``) Second aligned iterable object of symbols.
        - method (``str``) Type of alignment: 'global', 'global_cfe', 'local',
          'glocal'.
        - gap_open (``float``) The gap-opening cost.
        - gap_extend (``float``) The cost of extending an open gap.
        - gap_double (``float``) The gap-opening cost if a gap is already open
          in the other sequence.
        - matrix (``dict``) A score matrix dictionary.
        - max_hits (``int``) The maximum number of results to return in
          case multiple alignments with the same score are found. If set to 1,
          a single ``AlignmentResult`` object is returned. If set to values
          larger than 1, a list containing ``AlignmentResult`` objects are
          returned. If set to `None`, all alignments with the maximum score
          are returned.
    '''
    assert gap_extend <= 0, "gap_extend penalty must be <= 0"
    assert gap_open <= 0, "gap_open must be <= 0"

    max_j = len(seqj)
    max_i = len(seqi)
    seq1 = seqj if isinstance(seqj, bytes) else bytes(seqj, 'ascii')
    seq2 = seqi if isinstance(seqi, bytes) else bytes(seqi, 'ascii')

    if max_j > max_i:
        flipped = 1
        seq1, seq2 = seq2, seq1
        max_i, max_j = max_j, max_i
    else:
        flipped = 0

    if isinstance(matrix, basestring):
        if matrix == "BLOSUM62":
            score_matrix = m_BLOSUM62
        elif matrix == "DNAFULL":
            score_matrix = m_DNAFULL
        else:
            raise ValueError("Invalid premade scoring matrix:"
                             " {0}.".format(matrix))
    else:
        score_matrix = make_cmatrix(matrix)

    return caligner(seq1, seq2, METHODS[method],
                    gap_open, gap_extend, gap_double,
                    score_matrix, flipped)
