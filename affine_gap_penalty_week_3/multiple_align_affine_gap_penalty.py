"""
Global alignment of three sequences with affine gap penalties, using
sum-of-pairs column scoring.

For each column of the alignment:
    - for every pair of sequences in that column where both have real
      characters: score += match_reward or -mismatch_penalty
    - for every sequence that has a gap in that column: subtract gap_open
      if this is the first gap in a new run, or gap_extend if the gap run
      is continuing from the previous column

State space: dp[i][j][k][mask], where mask is a 3-bit flag indicating,
for each of v, w, u, whether that sequence was in the middle of a gap run
at the previous column (needed to know whether a new gap opens or an
existing one extends).

Caveat: this is only feasible for very short sequences (roughly
length <= 15-20) since the table has n1*n2*n3*8 states and 7 possible
moves per cell.

For 4+ sequences, exact affine MSA becomes computationally impractical 
this is precisely why the earlier "Multiple_align_LCS" problem sidesteps the issue entirely by using the trivial 0/1 column score (no gap penalty at all), 
and why real MSA software (ClustalW, MUSCLE, T-Coffee) uses heuristics instead of exact DP: progressive alignment (align pairs first, then merge into a growing profile), or iterative refinement — trading optimality guarantees for tractable runtime.
"""

NEG_INF = float('-inf')


def three_way_affine_alignment(match_reward, mismatch_penalty, gap_open, gap_extend, v, w, u):
    n1, n2, n3 = len(v), len(w), len(u)

    # dp[i][j][k][mask] -> best score; mask bit 0 = v in gap, bit 1 = w in gap, bit 2 = u in gap
    dp = {}
    backptr = {}

    def get(i, j, k, mask):
        return dp.get((i, j, k, mask), NEG_INF)

    dp[(0, 0, 0, 0)] = 0

    # the 7 non-empty subsets of {v, w, u} that could advance in one move
    moves = [(1, 0, 0), (0, 1, 0), (0, 0, 1),
             (1, 1, 0), (1, 0, 1), (0, 1, 1),
             (1, 1, 1)]

    for i in range(n1 + 1):
        for j in range(n2 + 1):
            for k in range(n3 + 1):
                for mask in range(8):
                    cur = get(i, j, k, mask)
                    if cur == NEG_INF and (i, j, k, mask) != (0, 0, 0, 0):
                        continue
                    for dv, dw, du in moves:
                        ni, nj, nk = i + dv, j + dw, k + du
                        if ni > n1 or nj > n2 or nk > n3:
                            continue

                        # column score: sum-of-pairs over the (at most 3) real chars,
                        # plus gap open/extend penalties for the sequences left behind
                        col_score = 0
                        chars = {}
                        if dv:
                            chars['v'] = v[i]
                        if dw:
                            chars['w'] = w[j]
                        if du:
                            chars['u'] = u[k]

                        keys = list(chars.keys())
                        for a in range(len(keys)):
                            for b in range(a + 1, len(keys)):
                                if chars[keys[a]] == chars[keys[b]]:
                                    col_score += match_reward
                                else:
                                    col_score -= mismatch_penalty

                        new_mask = 0
                        if not dv:
                            was_gap = mask & 1
                            col_score -= gap_extend if was_gap else gap_open
                            new_mask |= 1
                        if not dw:
                            was_gap = mask & 2
                            col_score -= gap_extend if was_gap else gap_open
                            new_mask |= 2
                        if not du:
                            was_gap = mask & 4
                            col_score -= gap_extend if was_gap else gap_open
                            new_mask |= 4

                        candidate = cur + col_score
                        if candidate > get(ni, nj, nk, new_mask):
                            dp[(ni, nj, nk, new_mask)] = candidate
                            backptr[(ni, nj, nk, new_mask)] = (i, j, k, mask, dv, dw, du)

    best_score = NEG_INF
    best_mask = None
    for mask in range(8):
        val = get(n1, n2, n3, mask)
        if val > best_score:
            best_score = val
            best_mask = mask

    aligned_v, aligned_w, aligned_u = _traceback(v, w, u, backptr, n1, n2, n3, best_mask)
    return best_score, aligned_v, aligned_w, aligned_u


def _traceback(v, w, u, backptr, i, j, k, mask):
    av, aw, au = [], [], []
    while (i, j, k, mask) != (0, 0, 0, 0):
        pi, pj, pk, pmask, dv, dw, du = backptr[(i, j, k, mask)]
        av.append(v[pi] if dv else '-')
        aw.append(w[pj] if dw else '-')
        au.append(u[pk] if du else '-')
        i, j, k, mask = pi, pj, pk, pmask
    av.reverse(); aw.reverse(); au.reverse()
    return ''.join(av), ''.join(aw), ''.join(au)


match_reward = 1
mismatch_penalty = 1
gap_open = 2
gap_extend = 1
v = "GATTACA"
w = "GCATGCU"
u = "GATTGCA"

score, av, aw, au = three_way_affine_alignment(match_reward, mismatch_penalty, gap_open, gap_extend, v, w, u)
print(score)
print(av)
print(aw)
print(au)