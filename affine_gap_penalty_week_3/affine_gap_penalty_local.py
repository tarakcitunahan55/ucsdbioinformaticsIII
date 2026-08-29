"""
Solve the local alignment with affine gap penalties problem.

this is the same three-layer affine gap dynamic programming as the global
version, but adapted for local alignment (smith-waterman style):

    - scores are never allowed to drop below 0 (a fresh alignment can
      restart anywhere), so 0 acts as a "give up and start over" option
    - the borders of the table are all initialized to 0, instead of
      charging gap penalties for leading gaps
    - the best score is the maximum value found anywhere in the middle
      matrix, not just the value in the bottom-right corner
    - traceback starts at the cell with that maximum value and stops as
      soon as it reaches a cell with score 0 (or hits row/column 0)

the three layers again represent:
    middle[i][j] : best score for some alignment of a suffix of v[0..i)
                   and a suffix of w[0..j), ending in a match/mismatch
    lower[i][j]  : best score ending in a gap in w (extra v symbol)
    upper[i][j]  : best score ending in a gap in v (extra w symbol)

Caveat: Returns only one alignment among others (if possible)
"""

NEG_INF = float('-inf')


def local_affine_gap_alignment(match_reward, mismatch_penalty, gap_open, gap_extend, v, w):
    """
    compute the optimal local alignment score and alignment of v and w
    under an affine gap penalty scheme.

    parameters:
        match_reward     -- reward added for a matching pair of symbols (positive number)
        mismatch_penalty -- penalty subtracted for a mismatching pair (positive number)
        gap_open         -- penalty subtracted for opening a new gap (positive number)
        gap_extend       -- penalty subtracted for extending an existing gap (positive number)
        v, w             -- the two strings to align

    returns:
        (best_score, aligned_v, aligned_w)
    """
    n, m = len(v), len(w)

    middle = [[0] * (m + 1) for _ in range(n + 1)]
    lower = [[0] * (m + 1) for _ in range(n + 1)]
    upper = [[0] * (m + 1) for _ in range(n + 1)]

    # local alignment: every border cell starts at 0, since an alignment
    # can begin fresh at any position rather than being forced to consume
    # a leading prefix as gaps
    for i in range(n + 1):
        middle[i][0] = 0
        lower[i][0] = 0
        upper[i][0] = NEG_INF
    for j in range(m + 1):
        middle[0][j] = 0
        upper[0][j] = 0
        lower[0][j] = NEG_INF

    best_score = 0
    best_i, best_j = 0, 0

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            lower[i][j] = max(lower[i - 1][j] - gap_extend,
                               middle[i - 1][j] - gap_open)

            upper[i][j] = max(upper[i][j - 1] - gap_extend,
                               middle[i][j - 1] - gap_open)

            if v[i - 1] == w[j - 1]:
                diag = middle[i - 1][j - 1] + match_reward
            else:
                diag = middle[i - 1][j - 1] - mismatch_penalty

            # the key local-alignment twist: never let the score go below 0,
            # since a new local alignment can always start here instead
            middle[i][j] = max(0, diag, lower[i][j], upper[i][j])

            if middle[i][j] > best_score:
                best_score = middle[i][j]
                best_i, best_j = i, j

    aligned_v, aligned_w = _traceback_local(v, w, middle, lower, upper,
                                             match_reward, mismatch_penalty,
                                             gap_open, gap_extend,
                                             best_i, best_j)

    return best_score, aligned_v, aligned_w


def _traceback_local(v, w, middle, lower, upper, match_reward, mismatch_penalty,
                      gap_open, gap_extend, start_i, start_j):
    """
    walk backwards from the cell holding the maximum score, staying inside
    whichever layer currently holds the optimal value, and stop as soon as
    a cell with score 0 is reached (that marks the start of the local
    alignment) or the edge of the table is hit.
    """
    i, j = start_i, start_j

    if middle[i][j] == 0:
        return '', ''

    if middle[i][j] == lower[i][j]:
        layer = 'lower'
    elif middle[i][j] == upper[i][j]:
        layer = 'upper'
    else:
        layer = 'middle'

    aligned_v = []
    aligned_w = []

    while i > 0 or j > 0:
        # stop the local alignment once we hit a cell whose middle score is 0
        if layer == 'middle' and middle[i][j] == 0:
            break

        if layer == 'lower':
            aligned_v.append(v[i - 1])
            aligned_w.append('-')
            if i > 1 and lower[i][j] == lower[i - 1][j] - gap_extend:
                layer = 'lower'
            else:
                layer = 'middle'
            i -= 1

        elif layer == 'upper':
            aligned_v.append('-')
            aligned_w.append(w[j - 1])
            if j > 1 and upper[i][j] == upper[i][j - 1] - gap_extend:
                layer = 'upper'
            else:
                layer = 'middle'
            j -= 1

        else:  # layer == 'middle'
            if i > 0 and j > 0:
                score = match_reward if v[i - 1] == w[j - 1] else -mismatch_penalty
                if middle[i][j] == middle[i - 1][j - 1] + score:
                    aligned_v.append(v[i - 1])
                    aligned_w.append(w[j - 1])
                    i -= 1
                    j -= 1
                    layer = 'middle'
                    continue
            if i > 0 and middle[i][j] == lower[i][j]:
                layer = 'lower'
            elif j > 0 and middle[i][j] == upper[i][j]:
                layer = 'upper'
            else:
                # score must have been 0 here; loop condition above will
                # catch this and stop on the next iteration
                break

    aligned_v.reverse()
    aligned_w.reverse()
    return ''.join(aligned_v), ''.join(aligned_w)

if __name__ == "__main__":
    match_reward = 1
    mismatch_penalty = 3
    gap_open = 2
    gap_extend = 1
    v = "GA"
    w = "GTTA"

    score, aligned_v, aligned_w = local_affine_gap_alignment(
        match_reward, mismatch_penalty, gap_open, gap_extend, v, w
    )

    print(score)
    print(aligned_v)
    print(aligned_w)
