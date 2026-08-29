"""
Solve the global alignment with affine gap penalties problem.

Given a match reward, a mismatch penalty, a gap opening penalty, a gap
extension penalty, and two nucleotide strings v and w, this script computes
the maximum global alignment score between v and w under the affine gap
penalty model, and reconstructs "one optimal alignment" achieving that score.

Affine gap model:
    - opening a gap (the first '-' in a run) costs "gap_open"
    - every additional '-' in the same run costs "gap_extend"
    - so a gap of length L costs gap_open + (L - 1) * gap_extend

The classic way to handle this with dynamic programming is to keep three
score matrices instead of one:

    middle[i][j] : best score for v[0..i), w[0..j) ending in a match/mismatch
    lower[i][j]  : best score ending in a gap in w (i.e. v[i-1] aligned to '-')
    upper[i][j]  : best score ending in a gap in v (i.e. w[j-1] aligned to '-')

"lower" is called that because in the classic 3-layer alignment graph it
sits "below" the middle layer (extra vertical edges, consuming v only), and
"upper" sits "above" it (extra horizontal edges, consuming w only).

Caveat: Returns only one alignment among others (if possible)

"""

NEG_INF = float('-inf')


def affine_gap_alignment(match_reward, mismatch_penalty, gap_open, gap_extend, v, w):
    """
    compute the optimal global alignment score and alignment of v and w
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

    # three score matrices, each of size (n+1) x (m+1)
    middle = [[0] * (m + 1) for _ in range(n + 1)]
    lower = [[0] * (m + 1) for _ in range(n + 1)]   # gap in w, consumes a v symbol
    upper = [[0] * (m + 1) for _ in range(n + 1)]   # gap in v, consumes a w symbol

    # base case: aligning empty prefixes
    middle[0][0] = 0
    lower[0][0] = NEG_INF
    upper[0][0] = NEG_INF

    # base case: v prefix against empty w -> all gaps in w (vertical moves only)
    for i in range(1, n + 1):
        lower[i][0] = -gap_open - gap_extend * (i - 1)
        middle[i][0] = lower[i][0]
        upper[i][0] = NEG_INF

    # base case: empty v against w prefix -> all gaps in v (horizontal moves only)
    for j in range(1, m + 1):
        upper[0][j] = -gap_open - gap_extend * (j - 1)
        middle[0][j] = upper[0][j]
        lower[0][j] = NEG_INF

    # fill in the rest of the table
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # extend an existing vertical gap, or open a fresh one from middle
            lower[i][j] = max(lower[i - 1][j] - gap_extend,
                               middle[i - 1][j] - gap_open)

            # extend an existing horizontal gap, or open a fresh one from middle
            upper[i][j] = max(upper[i][j - 1] - gap_extend,
                               middle[i][j - 1] - gap_open)

            # diagonal move: match or mismatch v[i-1] with w[j-1]
            if v[i - 1] == w[j - 1]:
                diag = middle[i - 1][j - 1] + match_reward
            else:
                diag = middle[i - 1][j - 1] - mismatch_penalty

            # middle takes the best of a diagonal step or closing either gap layer
            middle[i][j] = max(diag, lower[i][j], upper[i][j])

    best_score = middle[n][m]

    aligned_v, aligned_w = _traceback(v, w, middle, lower, upper,
                                       match_reward, mismatch_penalty,
                                       gap_open, gap_extend)

    return best_score, aligned_v, aligned_w


def _traceback(v, w, middle, lower, upper, match_reward, mismatch_penalty, gap_open, gap_extend):
    """
    walk backwards from (n, m) to (0, 0), always staying inside whichever
    layer ("middle", "lower", "upper") currently holds the optimal score,
    and prepend the corresponding symbols to the two output strings.
    """
    n, m = len(v), len(w)
    i, j = n, m

    # figure out which layer produced the optimum at (n, m)
    if middle[i][j] == lower[i][j]:
        layer = 'lower'
    elif middle[i][j] == upper[i][j]:
        layer = 'upper'
    else:
        layer = 'middle'

    aligned_v = []
    aligned_w = []

    while i > 0 or j > 0:
        if layer == 'lower':
            # a vertical move: v[i-1] aligned against a gap
            aligned_v.append(v[i - 1])
            aligned_w.append('-')
            # decide whether this gap symbol just opened or extends the one above
            if i > 1 and lower[i][j] == lower[i - 1][j] - gap_extend:
                layer = 'lower'
            else:
                layer = 'middle'
            i -= 1

        elif layer == 'upper':
            # a horizontal move: gap aligned against w[j-1]
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
                # fallback for edge (should only trigger at the borders)
                if i > 0:
                    layer = 'lower'
                else:
                    layer = 'upper'

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

    score, aligned_v, aligned_w = affine_gap_alignment(
        match_reward, mismatch_penalty, gap_open, gap_extend, v, w
    )

    print(score)
    print(aligned_v)
    print(aligned_w)

