"""Solve the global alignment problem for two nucleotide strings, but
find every alignment that achieves the best possible score, not just
one.

Given a match reward, a mismatch penalty, an indel penalty, and two
strings, this builds a scoring table using dynamic programming, then
uses memoized recursion to explore every tied optimal direction and
collect every distinct alignment that reaches the best score.

No affine gap penalties (opening and extending the gap/indels get the same penalty)
"""

import sys


def build_score_table(v, w, match_reward, mismatch_penalty, indel_penalty):
    """Fill in the scoring table and a matching backtrack table.

    score[i][j] holds the best possible alignment score between the
    first i characters of v and the first j characters of w.

    Unlike a single-path version, backtrack[i][j] stores a list of
    every direction that reaches the best score at that cell, not
    just one. A cell can end up with one, two, or all three of
    "down", "right", and "diagonal" recorded, whenever they tie.
    """
    n, m = len(v), len(w)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    backtrack = [[[] for _ in range(m + 1)] for _ in range(n + 1)]

    for i in range(1, n + 1):
        score[i][0] = score[i - 1][0] - indel_penalty
        backtrack[i][0] = ["down"]

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

            # every option is checked on its own, so ties are all
            # recorded together instead of only keeping the first one
            directions = []
            if best == down_score:
                directions.append("down")
            if best == right_score:
                directions.append("right")
            if best == diagonal_score:
                directions.append("diagonal")

            backtrack[i][j] = directions

    return score, backtrack


def all_alignments(backtrack, v, w, i, j, memo=None):
    """
    A recursive, top-down dynamic programming algorithm using memoization.
    Recursively collect every alignment achieving the best score,
    starting from cell (i, j) and working back toward the top left
    corner.

    Returns a set of (aligned_v, aligned_w) pairs. Memoization keyed
    on (i, j) avoids repeating work when different paths through the
    table land on the same cell, which happens often once several
    directions are tied at various points.
    """
    if memo is None:
        memo = {}

    if i == 0 and j == 0:
        return {("", "")} #every path recurses to 0,0 and two strings get created -> every distinct path through the tied backtrack pointers produces its own pair of aligned strings, and all of them end up collected into the final results set, with duplicates (identical string pairs arising from genuinely different paths) automatically merged away since it's a set.

    if (i, j) in memo:
        return memo[(i, j)]

    results = set()

    for direction in backtrack[i][j]:
        if direction == "down":
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


def all_global_alignments(v, w, match_reward, mismatch_penalty, indel_penalty):
    """Return the best alignment score along with every distinct
    alignment that achieves it.
    """
    score, backtrack = build_score_table(
        v, w, match_reward, mismatch_penalty, indel_penalty
    )
    best_score = score[len(v)][len(w)]
    alignments = all_alignments(backtrack, v, w, len(v), len(w))
    return best_score, alignments


if __name__ == "__main__":
    sys.setrecursionlimit(10000)  # protects against deep recursion on long strings

    match_reward = 1
    mismatch_penalty = 1
    indel_penalty = 2

    v = "GAGA"
    w = "GAT"

    best_score, alignments = all_global_alignments(
        v, w, match_reward, mismatch_penalty, indel_penalty
    )

    print(f"best score: {best_score}")
    print(f"number of distinct optimal alignments: {len(alignments)}")
    for aligned_v, aligned_w in sorted(alignments):
        print(aligned_v)
        print(aligned_w)
        print()