# -*- coding: utf-8 -*-

from collections import namedtuple

import numpy as np

from .matrix import BLOSUM62


# Container for alignment result
AlignmentResult = namedtuple(
    'AlignmentResult',
    ['seq1', 'seq2', 'start1', 'start2', 'score'])


def aligner(seqj, seqi, method='global', gap_open=-7, gap_extend=-7,
            gap_double=-7, matrix=BLOSUM62, max_hits=1):
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
    assert max_hits is None or max_hits > 0
    NONE, LEFT, UP, DIAG = range(4)  # NONE is 0
    max_j = len(seqj)
    max_i = len(seqi)

    if max_j > max_i:
        flip = 1
        seqi, seqj = seqj, seqi
        max_i, max_j = max_j, max_i
    else:
        flip = 0

    F = np.zeros((max_i + 1, max_j + 1), dtype=np.float32)
    I = np.ndarray((max_i + 1, max_j + 1), dtype=np.float32)
    I.fill(-np.inf)
    J = np.ndarray((max_i + 1, max_j + 1), dtype=np.float32)
    J.fill(-np.inf)
    pointer = np.zeros((max_i + 1, max_j + 1), dtype=np.uint)  # NONE

    if method == 'global':
        pointer[0, 1:] = LEFT
        pointer[1:, 0] = UP
        F[0, 1:] = gap_open + gap_extend * \
            np.arange(0, max_j, dtype=np.float32)
        F[1:, 0] = gap_open + gap_extend * \
            np.arange(0, max_i, dtype=np.float32)
    elif method == 'global_cfe':
        pointer[0, 1:] = LEFT
        pointer[1:, 0] = UP
    elif method == 'glocal':
        pointer[0, 1:] = LEFT
        F[0, 1:] = gap_open + gap_extend * \
            np.arange(0, max_j, dtype=np.float32)

    for i in range(1, max_i + 1):
        ci = seqi[i - 1]
        for j in range(1, max_j + 1):
            cj = seqj[j - 1]
            # I
            I[i, j] = max(
                         F[i, j - 1] + gap_open,
                         I[i, j - 1] + gap_extend,
                         J[i, j - 1] + gap_double)
            # J
            J[i, j] = max(
                         F[i - 1, j] + gap_open,
                         J[i - 1, j] + gap_extend,
                         I[i - 1, j] + gap_double)
            # F
            diag_score = F[i - 1, j - 1] + matrix[cj][ci]
            left_score = I[i, j]
            up_score = J[i, j]
            max_score = max(diag_score, up_score, left_score)

            F[i, j] = max(0, max_score) if method == 'local' else max_score

            if method == 'local':
                if F[i, j] == 0:
                    pass  # point[i,j] = NONE
                elif max_score == diag_score:
                    pointer[i, j] = DIAG
                elif max_score == up_score:
                    pointer[i, j] = UP
                elif max_score == left_score:
                    pointer[i, j] = LEFT
            elif method == 'glocal':
                # In a semi-global alignment we want to consume as much as
                # possible of the longer sequence.
                if max_score == up_score:
                    pointer[i, j] = UP
                elif max_score == diag_score:
                    pointer[i, j] = DIAG
                elif max_score == left_score:
                    pointer[i, j] = LEFT
            else:
                # global
                if max_score == up_score:
                    pointer[i, j] = UP
                elif max_score == left_score:
                    pointer[i, j] = LEFT
                else:
                    pointer[i, j] = DIAG

    # container for traceback coordinates
    ij_pairs = []
    if method == 'local':
        # max anywhere
        maxv_indices = np.argwhere(F == F.max())[:max_hits]
        for index in maxv_indices:
            ij_pairs.append(index)
    elif method == 'glocal':
        # max in last col
        max_score = F[:, -1].max()
        maxi_indices = np.argwhere(F[:, -1] == F[:, -1].max())\
            .flatten()[:max_hits]
        for i in maxi_indices:
            ij_pairs.append((i, max_j))
    elif method == 'global_cfe':
        # from i,j to max(max(last row), max(last col)) for free
        row_max = F[-1].max()
        col_max = F[:, -1].max()
        # expecting max to exist on either last column or last row
        if row_max > col_max:
            col_idces = np.argwhere(F[-1] == row_max).flatten()
            for cid in col_idces:
                ij_pairs.append((i, cid))
        elif row_max < col_max:
            row_idces = np.argwhere(F[:, -1] == col_max).flatten()
            for rid in row_idces:
                ij_pairs.append((rid, j))
        # special case: max is on last row, last col
        elif row_max == col_max == F[i, j]:
            # check if max score also exist on other cells in last row
            # or last col. we expect only one of the case.
            col_idces = np.argwhere(F[-1] == row_max).flatten()
            row_idces = np.argwhere(F[:, -1] == col_max).flatten()
            ncol_idces = len(col_idces)
            nrow_idces = len(row_idces)

            # tiebreaker between row/col is whichever has more max scores
            if ncol_idces > nrow_idces:
                for cid in col_idces:
                    ij_pairs.append((i, cid))
            elif ncol_idces < nrow_idces:
                for rid in row_idces:
                    ij_pairs.append((rid, j))
            elif ncol_idces == nrow_idces == 1:
                ij_pairs.append((i, j))
            else:
                raise RuntimeError('Unexpected multiple maximum global_cfe'
                                   ' scores.')
        else:
            raise RuntimeError('Unexpected global_cfe scenario.')
    else:
        # method must be global at this point
        ij_pairs.append((i, j))

    results = []
    for cur_i, (i, j) in enumerate(ij_pairs):
        align_j = []
        align_i = []
        score = F[i, j]
        p = pointer[i, j]

        # special case for global_cfe ~ one cell may contain multiple pointer
        # directions
        if method == 'global_cfe':
            if i < max_i:
                align_i.extend([c for c in seqi[i:][::-1]])
                align_j.extend(['-'] * (max_i - i))
            elif j < max_j:
                align_i.extend(['-'] * (max_j - j))
                align_j.extend([c for c in seqj[j:][::-1]])

        while p != NONE:
            if p == DIAG:
                i -= 1
                j -= 1
                align_j.append(seqj[j])
                align_i.append(seqi[i])
            elif p == LEFT:
                j -= 1
                align_j.append(seqj[j])
                align_i.append('-')
            elif p == UP:
                i -= 1
                align_j.append('-')
                align_i.append(seqi[i])
            else:
                raise Exception('wtf!')
            p = pointer[i, j]
        align_i = ''.join(align_i[::-1])
        align_j = ''.join(align_j[::-1])
        aln = (AlignmentResult(align_i, align_j, i, j, score)
               if flip else
               AlignmentResult(align_j, align_i, j, i, score))

        results.append(aln)
        if max_hits == 1:
            break

    return results
