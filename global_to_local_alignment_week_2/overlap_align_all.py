"""Solve the overlap alignment problem for two nucleotide strings, but
find every overlap alignment that achieves the best possible score,
not just one.

An overlap alignment aligns some suffix v' of v against some prefix
w' of w. Just like earlier "all alignments" versions, this records
every tied direction at every cell (instead of picking one), finds
every column where the best score could end (instead of picking the
first one), and then uses memoized recursion to explore every
combination of ties and collect every distinct resulting alignment.
"""

import sys


def build_score_table(v, w, match_reward, mismatch_penalty, indel_penalty):
    """Fill in the scoring table and a matching backtrack table for
    overlap alignment.

    score[i][j] holds the best possible score of an alignment between
    some suffix of v ending at v[i-1], and the first j characters of
    w (a genuine prefix of w, since w' is never allowed to skip its
    own beginning).

    Unlike a single-path version, backtrack[i][j] stores a list of
    every direction that reaches the best score at that cell, not
    just one, so ties can carry two or even three recorded directions.

    Column zero gets a free ride, fixed at zero for every row, since
    the alignment is allowed to begin anywhere along v. Row zero is
    forced to pay an indel penalty for every step, since w' must
    start from the very first character of w.
    """
    n, m = len(v), len(w)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    backtrack = [[[] for _ in range(m + 1)] for _ in range(n + 1)]

    for i in range(n + 1):
        score[i][0] = 0
        backtrack[i][0] = ["start"]

    for j in range(1, m + 1):
        score[0][j] = score[0][j - 1] - indel_penalty
        backtrack[0][j] = ["right"]

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if v[i - 1] == w[j - 1]:
                diagonal_score = score[i - 1][j - 1] + match_reward
            else:
                diagonal_score = score[i - 1][j - 1] - mismatch_penalty

            down_score = score[i - 1][j] - indel_penalty
            right_score = score[i][j - 1] - indel_penalty

            best = max(down_score, right_score, diagonal_score)
            score[i][j] = best

            # every option is checked on its own, so every tied
            # direction gets recorded together instead of only
            # keeping the first one found
            directions = []
            if best == down_score:
                directions.append("down")
            if best == right_score:
                directions.append("right")
            if best == diagonal_score:
                directions.append("diagonal")

            backtrack[i][j] = directions

    return score, backtrack


def find_best_ending_columns(score, n):
    """Find every column in the last row that achieves the overall
    best score.

    v' must be a suffix of v, so the alignment always ends at row n,
    but w' only needs to be a prefix of w, so it can end at any
    column of that last row -- and more than one column can share the
    same top score. Every one of those tied columns is a valid ending
    point.
    """
    m = len(score[0]) - 1
    best_score = score[n][0]

    for j in range(m + 1):
        if score[n][j] > best_score:
            best_score = score[n][j]

    best_columns = [j for j in range(m + 1) if score[n][j] == best_score]

    return best_score, best_columns


def all_alignments(backtrack, v, w, i, j, memo=None):
    """Recursively collect every alignment achieving the best score
    that ends at cell (i, j), working back toward wherever each one
    truly begins along v.

    Returns a set of (aligned_v, aligned_w) pairs. Whenever "start" is
    one of the recorded directions at a cell, that branch stops right
    there and contributes the pair of empty strings, since that marks
    where v was trimmed and the alignment begins. Memoization keyed
    on (i, j) avoids repeating work when different paths land on the
    same cell.
    """
    if memo is None:
        memo = {}

    if (i, j) in memo:
        return memo[(i, j)]

    results = set()

    for direction in backtrack[i][j]:
        if direction == "start":
            results.add(("", ""))
        elif direction == "down":
            for prefix_v, prefix_w in all_alignments(backtrack, v, w, i - 1, j, memo):
                results.add((prefix_v + v[i - 1], prefix_w + "-"))
        elif direction == "right":
            for prefix_v, prefix_w in all_alignments(backtrack, v, w, i, j - 1, memo):
                results.add((prefix_v + "-", prefix_w + w[j - 1]))
        else:  # diagonal
            for prefix_v, prefix_w in all_alignments(backtrack, v, w, i - 1, j - 1, memo):
                results.add((prefix_v + v[i - 1], prefix_w + w[j - 1]))

    memo[(i, j)] = results
    return results


def all_overlap_alignments(v, w, match_reward, mismatch_penalty, indel_penalty):
    """Return the best overlap alignment score along with every
    distinct alignment that achieves it.

    Builds the scoring table with every tied direction recorded per
    cell, finds every column in the last row sharing the overall best
    score since more than one ending column can tie, and then collects
    every alignment reachable by exploring every tied direction from
    each of those ending columns.
    """
    sys.setrecursionlimit(10000)  # protects against deep recursion on long strings

    score, backtrack = build_score_table(v, w, match_reward, mismatch_penalty, indel_penalty)
    best_score, best_columns = find_best_ending_columns(score, len(v))

    memo = {}
    alignments = set()
    for j in best_columns:
        alignments |= all_alignments(backtrack, v, w, len(v), j, memo)

    return best_score, alignments


if __name__ == "__main__":
    match_reward = 1
    mismatch_penalty = 1
    indel_penalty = 2

    v = "GAGA" # contributes a suffix
    w = "GAT" # contributes a prefix
    # the order you put v,w or w,v is important
    
    best_score, alignments = all_overlap_alignments(
        v, w, match_reward, mismatch_penalty, indel_penalty
    )

    print(f"best score: {best_score}")
    print(f"number of distinct optimal overlap alignments: {len(alignments)}")
    for aligned_v, aligned_w in sorted(alignments):
        print(aligned_v)
        print(aligned_w)
        print()