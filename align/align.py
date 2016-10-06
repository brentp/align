import numpy as np
from .matrix import BLOSUM62, DNAFULL
from collections import namedtuple

AlignmentResult = namedtuple(
    'AlignmentResult',
    ['seq1', 'seq2', 'pos1', 'pos2', 'score'])

def max_index(array):
    """
    """
    max_value = array.argmax()
    idx = np.unravel_index(max_value, array.shape)
    return idx

def aligner(seqj, seqi, method='global', gap_open=-7, gap_extend=-7, \
            gap_double=-7, matrix=BLOSUM62, n_max_return=1):
    """
    Calculates the alignment of two sequences. The supported "methods" are
    "global" for a global Needleman-Wunsh algorithm, "local" for a local
    Smith-Waterman alignment, "global_cfe" for a global alignment with cost-free
    ends and "glocal" for an alignment which is "global" only with respect to
    the shorter sequence, this is also known as a "semi-global" alignment."
    Returns the aligned (sub)sequences as character arrays.

    Gotoh, O. (1982). J. Mol. Biol. 162, 705-708.
    Needleman, S. & Wunsch, C. (1970). J. Mol. Biol. 48(3), 443-53.
    Smith, T.F. & Waterman M.S. (1981). J. Mol. Biol. 147, 195-197.

    Arguments:

        - seqj (``sequence``) First aligned iterable object of symbols.
        - seqi (``sequence``) Second aligned iterable object of symbols.
        - method (``str``) Type of alignment: "global", "global_cfe", "local",
          "glocal".
        - gap_open (``float``) The gap-opening cost.
        - gap_extend (``float``) The cost of extending an open gap.
        - gap_double (``float``) The gap-opening cost if a gap is already open
          in the other sequence.
        - matrix (``dict``) A score matrix dictionary.
        - n_max_return (``int``) The maximum number of results to return in
          case multiple alignments with the same score are found. If set to 1,
          a single ``AlignmentResult`` object is returned. If set to values
          larger than 1, a list containing ``AlignmentResult`` objects are
          returned.
    """
    assert n_max_return > 0
    NONE, LEFT, UP, DIAG = range(4) # NONE is 0
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
    pointer = np.zeros((max_i + 1, max_j + 1), dtype=np.uint) # NONE

    if method == 'global':
        pointer[0, 1:] = LEFT
        pointer[1:, 0] = UP
        F[0, 1:] = gap_open + gap_extend * np.arange(0, max_j, dtype=np.float32)
        F[1:, 0] = gap_open + gap_extend * np.arange(0, max_i, dtype=np.float32)
    elif method == 'global_cfe':
        pointer[0, 1:] = LEFT
        pointer[1:, 0] = UP
    elif method == 'glocal':
        pointer[0, 1:] = LEFT
        F[0, 1:] = gap_open + gap_extend * np.arange(0, max_j, dtype=np.float32)

    for i in range(1, max_i + 1):
        ci = seqi[i - 1]
        for j in range(1, max_j + 1):
            cj = seqj[j - 1]
            # I
            I[i,j] = max(
                         F[i, j - 1] + gap_open,
                         I[i, j - 1] + gap_extend,
                         J[i, j - 1] + gap_double)
            # J
            J[i,j] = max(
                         F[i - 1, j] + gap_open,
                         J[i - 1, j] + gap_extend,
                         I[i - 1, j] + gap_double)
            # F
            diag_score = F[i - 1, j - 1] + matrix[cj][ci]
            left_score = I[i, j]
            up_score   = J[i, j]
            max_score = max(diag_score, up_score, left_score)

            F[i,j] = max(0, max_score) if method == 'local' else max_score

            if method == 'local':
                if F[i,j] == 0:
                    pass # point[i,j] = NONE
                elif max_score == diag_score:
                    pointer[i,j] = DIAG
                elif max_score == up_score:
                    pointer[i,j] = UP
                elif max_score == left_score:
                    pointer[i,j] = LEFT
            elif method == 'glocal':
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

    # TODO: break this big if-block apart
    if n_max_return == 1:
        if method == 'local':
            # max anywhere
            i, j = max_index(F)
        elif method == 'glocal':
            # max in last col
            i, j = (F[:,-1].argmax(), max_j)
        elif method == 'global_cfe':
            # from i,j to max(max(last row), max(last col)) for free
            row_max, col_idx = F[-1].max(), F[-1].argmax()
            col_max, row_idx = F[:, -1].max(), F[:, -1].argmax()
            if row_max > col_max:
                pointer[-1,col_idx+1:] = LEFT
            else:
                pointer[row_idx+1:,-1] = UP

        align_j = []
        align_i = []
        p = pointer[i, j]
        while p != NONE:
            if p == DIAG:
                i -= 1
                j -= 1
                align_j.append(seqj[j])
                align_i.append(seqi[i])
            elif p == LEFT:
                j -= 1
                align_j.append(seqj[j])
                align_i.append("-")
            elif p == UP:
                i -= 1
                align_j.append("-")
                align_i.append(seqi[i])
            else:
                raise Exception('wtf!')
            p = pointer[i, j]
        align_i = "".join(align_i[::-1])
        align_j = "".join(align_j[::-1])
        #np.array(align_i.reverse())
        return (AlignmentResult(align_i, align_j, None, None, None)
                if flip else AlignmentResult(align_j, align_i, None, None, None))
    else:
        ijs = []
        if method == "glocal":
            max_score = F[:, -1].max()
            maxv_indices = np.argwhere(F[:, -1] == max_score)[0]
            ijs = [(i, max_j, max_score) for i in maxv_indices][:n_max_return]
        else:
            raise NotImplementedError

        results = []

        for i, j, score in ijs:
            align_j = []
            align_i = []
            p = pointer[i, j]
            while p != NONE:
                if p == DIAG:
                    i -= 1
                    j -= 1
                    align_j.append(seqj[j])
                    align_i.append(seqi[i])
                elif p == LEFT:
                    j -= 1
                    align_j.append(seqj[j])
                    align_i.append("-")
                elif p == UP:
                    i -= 1
                    align_j.append("-")
                    align_i.append(seqi[i])
                else:
                    raise Exception('wtf!')
                p = pointer[i, j]
            align_i = "".join(align_i[::-1])
            align_j = "".join(align_j[::-1])
            aln = (AlignmentResult(align_i, align_j, i, j, score)
                   if flip else AlignmentResult(align_j, align_i, j, i, score))
            results.append(aln)

        return results


if __name__ == '__main__':
    # global
    a, b = aligner('WW','WEW', method= 'global')
    assert a == 'W-W'
    assert b == 'WEW'
