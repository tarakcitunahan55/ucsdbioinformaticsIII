"""
Global alignment of three sequences with constant (linear) indel penalty,
using sum-of-pairs column scoring.

for each column of the alignment:
    - for every pair of sequences that both contribute a real character
      in that column: add match_reward if they're equal, subtract
      mismatch_penalty if they differ
    - for every sequence that has a gap in that column: subtract a flat
      indel_penalty (same cost regardless of how long the gap run is --
      no gap-open vs. gap-extend distinction, unlike the affine version)

this is a direct 3-string generalization of ordinary global alignment:
same dp structure as the multiple lcs code, but with real match/mismatch/
indel scores plugged into each of the 7 possible moves instead of the
simple 0/1 "all three agree" rule.

Caveat: it is computationally faster than multiple_align_affine_gap_penalty code (linear indel is easier to run at 4+ sequences than affine, but with realistic datasets with length > 15 both are impractical)
but biologoically less accurate since the indel penalty is linear/constant
"""


def three_way_linear_gap_alignment(match_reward, mismatch_penalty, indel_penalty, v, w, u):
    """
    compute the optimal global alignment score of three strings v, w, u
    under sum-of-pairs scoring with a constant indel penalty, and
    reconstruct one optimal 3-way alignment achieving that score.
    """
    n1, n2, n3 = len(v), len(w), len(u)

    dp = [[[0] * (n3 + 1) for _ in range(n2 + 1)] for _ in range(n1 + 1)]
    backpointer = [[[None] * (n3 + 1) for _ in range(n2 + 1)] for _ in range(n1 + 1)]

    # base cases: aligning a prefix of one string against nothing but gaps
    # in the other two, so every real character costs indel_penalty twice
    # (once against each of the other two all-gap sequences)
    for i in range(1, n1 + 1):
        dp[i][0][0] = dp[i - 1][0][0] - 2 * indel_penalty
        backpointer[i][0][0] = (1, 0, 0)
    for j in range(1, n2 + 1):
        dp[0][j][0] = dp[0][j - 1][0] - 2 * indel_penalty
        backpointer[0][j][0] = (0, 1, 0)
    for k in range(1, n3 + 1):
        dp[0][0][k] = dp[0][0][k - 1] - 2 * indel_penalty
        backpointer[0][0][k] = (0, 0, 1)

    for i in range(n1 + 1):
        for j in range(n2 + 1):
            for k in range(n3 + 1):
                if i == 0 and j == 0 and k == 0:
                    continue
                if (i, j, 0) == (i, j, 0) and j > 0 and i > 0 and k == 0:
                    pass  # placeholder, real 2d base cases handled below

    # handle the three 2d "face" base cases (one sequence empty, other two aligned pairwise)
    for i in range(1, n1 + 1):
        for j in range(1, n2 + 1):
            best_value = NEG_INF = float('-inf')
            best_move = None
            pair_score = match_reward if v[i-1] == w[j-1] else -mismatch_penalty
            candidate = dp[i-1][j-1][0] + pair_score - 2 * indel_penalty
            # note: when k=0, u contributes nothing, so this column has v,w real + u gap
            # sum-of-pairs: (v,w) pair scored normally, (v,u) and (w,u) both charged indel
            if candidate > best_value:
                best_value, best_move = candidate, (1, 1, 0)
            candidate = dp[i-1][j][0] - 2 * indel_penalty
            if candidate > best_value:
                best_value, best_move = candidate, (1, 0, 0)
            candidate = dp[i][j-1][0] - 2 * indel_penalty
            if candidate > best_value:
                best_value, best_move = candidate, (0, 1, 0)
            dp[i][j][0] = best_value
            backpointer[i][j][0] = best_move

    for i in range(1, n1 + 1):
        for k in range(1, n3 + 1):
            best_value = float('-inf')
            best_move = None
            pair_score = match_reward if v[i-1] == u[k-1] else -mismatch_penalty
            candidate = dp[i-1][0][k-1] + pair_score - 2 * indel_penalty
            if candidate > best_value:
                best_value, best_move = candidate, (1, 0, 1)
            candidate = dp[i-1][0][k] - 2 * indel_penalty
            if candidate > best_value:
                best_value, best_move = candidate, (1, 0, 0)
            candidate = dp[i][0][k-1] - 2 * indel_penalty
            if candidate > best_value:
                best_value, best_move = candidate, (0, 0, 1)
            dp[i][0][k] = best_value
            backpointer[i][0][k] = best_move

    for j in range(1, n2 + 1):
        for k in range(1, n3 + 1):
            best_value = float('-inf')
            best_move = None
            pair_score = match_reward if w[j-1] == u[k-1] else -mismatch_penalty
            candidate = dp[0][j-1][k-1] + pair_score - 2 * indel_penalty
            if candidate > best_value:
                best_value, best_move = candidate, (0, 1, 1)
            candidate = dp[0][j-1][k] - 2 * indel_penalty
            if candidate > best_value:
                best_value, best_move = candidate, (0, 1, 0)
            candidate = dp[0][j][k-1] - 2 * indel_penalty
            if candidate > best_value:
                best_value, best_move = candidate, (0, 0, 1)
            dp[0][j][k] = best_value
            backpointer[0][j][k] = best_move

    # main interior: all three indices positive
    for i in range(1, n1 + 1):
        for j in range(1, n2 + 1):
            for k in range(1, n3 + 1):
                best_value = float('-inf')
                best_move = None

                def pair(a, b):
                    return match_reward if a == b else -mismatch_penalty

                # move (1,1,1): all three contribute real characters -- sum-of-pairs over all 3 pairs
                col_score = pair(v[i-1], w[j-1]) + pair(v[i-1], u[k-1]) + pair(w[j-1], u[k-1])
                candidate = dp[i-1][j-1][k-1] + col_score
                if candidate > best_value:
                    best_value, best_move = candidate, (1, 1, 1)

                # moves with exactly two real characters, one gap -- the gapped
                # sequence is charged indel_penalty against each of the other two
                candidate = dp[i-1][j-1][k] + pair(v[i-1], w[j-1]) - 2 * indel_penalty
                if candidate > best_value:
                    best_value, best_move = candidate, (1, 1, 0)
                candidate = dp[i-1][j][k-1] + pair(v[i-1], u[k-1]) - 2 * indel_penalty
                if candidate > best_value:
                    best_value, best_move = candidate, (1, 0, 1)
                candidate = dp[i][j-1][k-1] + pair(w[j-1], u[k-1]) - 2 * indel_penalty
                if candidate > best_value:
                    best_value, best_move = candidate, (0, 1, 1)

                # moves with exactly one real character, two gaps -- charged
                # indel_penalty against each of the other two gapped sequences,
                # plus indel_penalty between the two gaps themselves (by convention,
                # gap-vs-gap columns are usually scored 0; adjust if your rubric differs)
                candidate = dp[i-1][j][k] - 2 * indel_penalty
                if candidate > best_value:
                    best_value, best_move = candidate, (1, 0, 0)
                candidate = dp[i][j-1][k] - 2 * indel_penalty
                if candidate > best_value:
                    best_value, best_move = candidate, (0, 1, 0)
                candidate = dp[i][j][k-1] - 2 * indel_penalty
                if candidate > best_value:
                    best_value, best_move = candidate, (0, 0, 1)

                dp[i][j][k] = best_value
                backpointer[i][j][k] = best_move

    best_score = dp[n1][n2][n3]
    aligned_v, aligned_w, aligned_u = _traceback(v, w, u, backpointer, n1, n2, n3)
    return best_score, aligned_v, aligned_w, aligned_u


def _traceback(v, w, u, backpointer, i, j, k):
    av, aw, au = [], [], []
    while i > 0 or j > 0 or k > 0:
        dv, dw, du = backpointer[i][j][k]
        av.append(v[i-1] if dv else '-')
        aw.append(w[j-1] if dw else '-')
        au.append(u[k-1] if du else '-')
        i -= dv; j -= dw; k -= du
    av.reverse(); aw.reverse(); au.reverse()
    return ''.join(av), ''.join(aw), ''.join(au)


match_reward = 1
mismatch_penalty = 1
indel_penalty = 2
v = "GATTACA"
w = "GCATGCU"
u = "GATTGCA"

score, av, aw, au = three_way_linear_gap_alignment(match_reward, mismatch_penalty, indel_penalty, v, w, u)
print(score)
print(av)
print(aw)
print(au)