"""
Solve the global alignment with affine gap penalties problem, returning all
optimal alignments instead of just one.

Warning: the number of optimal alignments can grow exponentially with
sequence length whenever there are many ties in the dp recurrence, so this
is only practical for short sequences (length <= 10). For longer sequences, use the single-path
version.
"""

NEG_INF = float('-inf')


def affine_gap_alignment_all(match_reward, mismatch_penalty, gap_open, gap_extend, v, w):
    """
    compute the optimal global alignment score under an affine gap penalty
    scheme, and return every distinct optimal alignment achieving it.

    """
    n, m = len(v), len(w)

    middle = [[0] * (m + 1) for _ in range(n + 1)]
    lower = [[0] * (m + 1) for _ in range(n + 1)]
    upper = [[0] * (m + 1) for _ in range(n + 1)]

    middle[0][0] = 0
    lower[0][0] = NEG_INF
    upper[0][0] = NEG_INF

    for i in range(1, n + 1):
        lower[i][0] = -gap_open - gap_extend * (i - 1)
        middle[i][0] = lower[i][0]
        upper[i][0] = NEG_INF

    for j in range(1, m + 1):
        upper[0][j] = -gap_open - gap_extend * (j - 1)
        middle[0][j] = upper[0][j]
        lower[0][j] = NEG_INF

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
            middle[i][j] = max(diag, lower[i][j], upper[i][j])

    best_score = middle[n][m]

    results = set()  # use a set to dedupe alignments reached via different branches
    aligned_v, aligned_w = [], []

    def recurse(i, j, layer):
        """
        explore every branch of the traceback that matches the optimal
        score at the current cell, in the current layer. base case: reached
        (0, 0) with no cell left to explain, so record the alignment built
        so far.
        """
        if i == 0 and j == 0:
            results.add((''.join(reversed(aligned_v)), ''.join(reversed(aligned_w))))
            return

        if layer == 'lower':
            aligned_v.append(v[i - 1])
            aligned_w.append('-')
            # branch 1: this gap symbol extends a gap that was already open
            if i > 1 and lower[i][j] == lower[i - 1][j] - gap_extend:
                recurse(i - 1, j, 'lower')
            # branch 2: this gap symbol is the first one of a freshly opened gap
            if lower[i][j] == middle[i - 1][j] - gap_open:
                recurse(i - 1, j, 'middle')
            aligned_v.pop()
            aligned_w.pop()

        elif layer == 'upper':
            aligned_v.append('-')
            aligned_w.append(w[j - 1])
            if j > 1 and upper[i][j] == upper[i][j - 1] - gap_extend:
                recurse(i, j - 1, 'upper')
            if upper[i][j] == middle[i][j - 1] - gap_open:
                recurse(i, j - 1, 'middle')
            aligned_v.pop()
            aligned_w.pop()

        else:  # layer == 'middle'
            # branch A: reached here via a diagonal match/mismatch step
            if i > 0 and j > 0:
                score = match_reward if v[i - 1] == w[j - 1] else -mismatch_penalty
                if middle[i][j] == middle[i - 1][j - 1] + score:
                    aligned_v.append(v[i - 1])
                    aligned_w.append(w[j - 1])
                    recurse(i - 1, j - 1, 'middle')
                    aligned_v.pop()
                    aligned_w.pop()
            # branch B: reached here by closing a gap that was open in 'lower'
            if i > 0 and middle[i][j] == lower[i][j]:
                recurse(i, j, 'lower')
            # branch C: reached here by closing a gap that was open in 'upper'
            if j > 0 and middle[i][j] == upper[i][j]:
                recurse(i, j, 'upper')

    # start the search from whichever layer(s) at (n, m) actually hold the
    # optimal score -- there could be more than one
    if middle[n][m] == middle[n - 1][m - 1] + (
        match_reward if n > 0 and m > 0 and v[n - 1] == w[m - 1] else -mismatch_penalty
    ) if n > 0 and m > 0 else False:
        pass  # this check is redundant with the 'middle' branch inside recurse; just start there

    recurse(n, m, 'middle')
    if middle[n][m] == lower[n][m]:
        recurse(n, m, 'lower')
    if middle[n][m] == upper[n][m]:
        recurse(n, m, 'upper')

    return best_score, list(results)


match_reward = 1
mismatch_penalty = 3
gap_open = 2
gap_extend = 1
v = "GA"
w = "GTTA"

score, all_alignments = affine_gap_alignment_all(match_reward, mismatch_penalty, gap_open, gap_extend, v, w)

print(score)
print(f"number of optimal alignments: {len(all_alignments)}")
for av, aw in all_alignments:
    print(av)
    print(aw)
    print()