# -*- coding: utf-8 -*-

from collections import namedtuple
from itertools import islice

import numpy as np
import os.path as op
import sys

cimport numpy as np
from libc.string cimport strlen
from cpython.mem cimport PyMem_Malloc, PyMem_Free

from .align import AlignmentResult
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


cdef list caligner(
    const unsigned char* seqj, const unsigned char* seqi, const int imethod,
    const DTYPE_FLOAT gap_open, const DTYPE_FLOAT gap_extend, const DTYPE_FLOAT gap_double,
    const DTYPE_FLOAT[:, :] amatrix, const bint flipped, max_hits):

    cdef:
        unsigned char* align_j
        unsigned char* align_i
        unsigned char ci, cj
        size_t max_j = strlen(<char *>seqj)
        size_t max_i = strlen(<char *>seqi)
        size_t align_counter = 0
        size_t i = 1, j = 1
        int ncol_idces, nrow_idces, idx
        int end_i, end_j, n_gaps_i, n_gaps_j, n_mmatch
        DTYPE_UINT p
        DTYPE_FLOAT diag_score, left_score, up_score, max_score, aln_score
        np.ndarray[DTYPE_FLOAT, ndim=2] agap_i = np.empty((max_i + 1, max_j + 1), dtype=np.float32)
        np.ndarray[DTYPE_FLOAT, ndim=2] agap_j = np.empty((max_i + 1, max_j + 1), dtype=np.float32)
        np.ndarray[DTYPE_FLOAT, ndim=2] score = np.zeros((max_i + 1, max_j + 1), dtype=np.float32)
        np.ndarray[DTYPE_UINT, ndim=2] pointer = np.zeros((max_i + 1, max_j + 1), dtype=np.uint)
        list indices = [], row_idces = [], col_idces = []
        list results = [], tracestart_coords = []

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

    cdef DTYPE_FLOAT matrix_max = 0, row_max = score[-1, 0], col_max = score[0, -1]

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
            if imethod == 1:
                max_score = max2(0, max_score)

            score[i, j] = max_score

            if imethod == 1:
                if score[i,j] == 0:
                    pass # point[i,j] = NONE
                elif max_score == diag_score:
                    pointer[i,j] = DIAG
                elif max_score == up_score:
                    pointer[i,j] = UP
                elif max_score == left_score:
                    pointer[i,j] = LEFT

                # Manual tracking of [i, j] coordinates where score is max
                if max_score > matrix_max:
                    matrix_max = max_score
                    indices = [(i, j)]
                elif max_score == matrix_max:
                    indices.append((i, j))

            elif imethod == 2:
                # In a semi-global alignment we want to consume as much as
                # possible of the longer sequence.
                if max_score == up_score:
                    pointer[i,j] = UP
                elif max_score == diag_score:
                    pointer[i,j] = DIAG
                elif max_score == left_score:
                    pointer[i,j] = LEFT

                # Manual tracking of [i, j] coordinates where col score is max
                if j == max_j:
                    if max_score > col_max:
                        col_max = max_score
                        col_idces = [(i, j)]
                    elif max_score == col_max:
                        col_idces.append((i, j))

            else:
                # global
                if max_score == up_score:
                    pointer[i,j] = UP
                elif max_score == left_score:
                    pointer[i,j] = LEFT
                else:
                    pointer[i,j] = DIAG

                if imethod == 3:
                    # Manual tracking of [i, j] coordinates where col score is max
                    if j == max_j:
                        if max_score > col_max:
                            col_max = max_score
                            col_idces = [(i, j)]
                        elif max_score == col_max:
                            col_idces.append((i, j))
                    if i == max_i:
                        if max_score > row_max:
                            row_max = max_score
                            row_idces = [(i, j)]
                        elif max_score == row_max:
                            row_idces.append((i, j))

    if imethod == 1:
        for index in islice(indices, max_hits):
            tracestart_coords.append(index)

    elif imethod == 2:
        for index in islice(col_idces, max_hits):
            tracestart_coords.append(index)

    elif imethod == 3:
        # from i,j to max(max(last row), max(last col)) for free
        # expecting max to exist on either last column or last row
        if row_max > col_max:
            for index in islice(row_idces, max_hits):
                tracestart_coords.append(index)
        elif row_max < col_max:
            for index in islice(col_idces, max_hits):
                tracestart_coords.append(index)
        # special case: max is on last row, last col
        elif row_max == col_max == score[max_i, max_j]:
            ncol_idces = len(col_idces)
            nrow_idces = len(row_idces)
            # tiebreaker between row/col is whichever has more max scores
            if ncol_idces > nrow_idces:
                for index in islice(col_idces, max_hits):
                    tracestart_coords.append(index)
            elif ncol_idces < nrow_idces:
                for index in islice(row_idces, max_hits):
                    tracestart_coords.append(index)
            elif ncol_idces == nrow_idces == 1:
                tracestart_coords.append((max_i, max_j))
            else:
                raise RuntimeError('Unexpected multiple maximum global_cfe'
                                   ' scores.')
        else:
            raise RuntimeError('Unexpected global_cfe scenario.')
    else:
        # method must be global at this point
        tracestart_coords.append((max_i, max_j))

    seqlen = max_i + max_j
    for i, j in tracestart_coords:
        if imethod == 0 or imethod == 3:
            end_i, end_j = max_i, max_j
        else:
            end_i, end_j = i, j
        aln_score = score[i, j]
        p = pointer[i, j]
        n_gaps_i, n_gaps_j, n_mmatch = 0, 0, 0
        align_counter = 0
        align_i = <unsigned char *>PyMem_Malloc(seqlen * sizeof(unsigned char))
        align_j = <unsigned char *>PyMem_Malloc(seqlen * sizeof(unsigned char))

        # special case for global_cfe ~ one cell may contain multiple pointer
        # directions
        if imethod == 3:
            if i < max_i:
                n_gaps_j += 1
                align_counter = max_i - i
                for idx in range(align_counter):
                    align_j[idx] = GAP_CHAR
                    align_i[idx] = seqi[max_i - 1 * idx - 1]
            elif j < max_j:
                n_gaps_i += 1
                align_counter = max_j - j
                for idx in range(align_counter):
                    align_i[idx] = GAP_CHAR
                    align_j[idx] = seqj[max_j - 1 * idx - 1]

        while p != NONE:
            if p == DIAG:
                i -= 1
                j -= 1
                if seqi[i] != seqj[j]:
                    n_mmatch += 1
                align_j[align_counter] = seqj[j]
                align_i[align_counter] = seqi[i]
            elif p == LEFT:
                j -= 1
                align_j[align_counter] = seqj[j]
                if align_i[align_counter - 1] != GAP_CHAR or align_counter == 0:
                    n_gaps_i += 1
                align_i[align_counter] = GAP_CHAR
            elif p == UP:
                i -= 1
                align_i[align_counter] = seqi[i]
                if align_j[align_counter - 1] != GAP_CHAR or align_counter == 0:
                    n_gaps_j += 1
                align_j[align_counter] = GAP_CHAR
            else:
                raise Exception('wtf!')
            p = pointer[i, j]
            align_counter += 1

        alns_i = bytes(align_i[:align_counter][::-1])
        alns_j = bytes(align_j[:align_counter][::-1])

        PyMem_Free(align_i)
        PyMem_Free(align_j)

        aln = (AlignmentResult(alns_i, alns_j, i, j, end_i, end_j,
                            n_gaps_i, n_gaps_j, n_mmatch, aln_score)
            if flipped else
            AlignmentResult(alns_j, alns_i, j, i, end_j, end_i,
                            n_gaps_j, n_gaps_i, n_mmatch, aln_score))

        results.append(aln)

    return results


def aligner(seqj, seqi, method='global', gap_open=-7, gap_extend=-7,
            gap_double=-7, matrix="BLOSUM62", max_hits=1):
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
    assert max_hits is None or max_hits > 0

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

    if isinstance(matrix, dict):
        score_matrix = make_cmatrix(matrix)
    elif matrix == "BLOSUM62":
        score_matrix = m_BLOSUM62
    elif matrix == "DNAFULL":
        score_matrix = m_DNAFULL
    else:
        raise ValueError("Invalid premade scoring matrix:"
                            " {0}.".format(matrix))

    return caligner(seq1, seq2, METHODS[method],
                    gap_open, gap_extend, gap_double,
                    score_matrix, flipped, max_hits)
