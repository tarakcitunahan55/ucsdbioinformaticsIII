"""Solve the local alignment problem for two nucleotide strings, but
find every local alignment that achieves the best possible score, not
just one.

Given a match reward, a mismatch penalty, an indel penalty, and two
strings, this builds a scoring table using dynamic programming, then
uses memoized recursion to explore every tied optimal direction --
including every tied place an alignment could end and every tied
place one could start -- and collects every distinct alignment that
reaches the best score.

No affine gap penalties (opening and extending the gap/indels get the same penalty)
"""

import sys


def build_score_table(v, w, match_reward, mismatch_penalty, indel_penalty):
    """Fill in the scoring table and a matching backtrack table.

    score[i][j] holds the best possible score of a local alignment
    ending at v[i-1] paired with w[j-1] -- or, thanks to the free
    ride option, possibly just a score of zero if starting fresh here
    beats continuing whatever came before.

    Unlike a single-path version, backtrack[i][j] stores a list of
    every source that reaches the best score at that cell, not just
    one. A cell can end up with any combination of "down", "right",
    "diagonal", and "start" recorded, whenever they tie.
    """
    n, m = len(v), len(w)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    backtrack = [[[] for _ in range(m + 1)] for _ in range(n + 1)]

    # the first row and column (index zero) have nowhere to backtrack to, so every cell along them is simply treated as a starting point
    #Since indel penalties only ever subtract from the score, dragging in leading gaps to reach row 0 or column 0 can never beat just starting fresh with a score of zero — so every cell in 1st row/column (index zero) there is always a valid "start," and score zero, never a forced gap step.
    for i in range(n + 1):
        backtrack[i][0] = ["start"]
    for j in range(m + 1):
        backtrack[0][j] = ["start"]

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if v[i - 1] == w[j - 1]:
                diagonal_score = score[i - 1][j - 1] + match_reward
            else:
                diagonal_score = score[i - 1][j - 1] - mismatch_penalty

            down_score = score[i - 1][j] - indel_penalty
            right_score = score[i][j - 1] - indel_penalty
            free_ride_score = 0

            best = max(down_score, right_score, diagonal_score, free_ride_score)
            score[i][j] = best

            # every option is checked on its own, so every tied
            # source gets recorded together instead of only keeping
            # the first one found
            sources = []
            if best == down_score:
                sources.append("down")
            if best == right_score:
                sources.append("right")
            if best == diagonal_score:
                sources.append("diagonal")
            if best == free_ride_score:
                sources.append("start")

            backtrack[i][j] = sources

    return score, backtrack


def find_best_ending_cells(score):
    """Find every cell that achieves the overall best score.

    A local alignment can end anywhere in the table, and more than
    one cell can share the same top score, so every one of those
    tied cells needs to be treated as a valid ending point.
    """
    n = len(score) - 1
    m = len(score[0]) - 1

    best_score = score[0][0]
    for i in range(n + 1):
        for j in range(m + 1):
            if score[i][j] > best_score:
                best_score = score[i][j]

    best_cells = [
        (i, j)
        for i in range(n + 1)
        for j in range(m + 1)
        if score[i][j] == best_score
    ]

    return best_score, best_cells


def all_alignments(backtrack, v, w, i, j, memo=None):
    """Recursively collect every alignment achieving the best score
    that ends at cell (i, j), working back toward wherever each one
    actually began.

    Returns a set of (aligned_v, aligned_w) pairs. Whenever "start"
    is one of the recorded sources at a cell, that branch stops right
    there and contributes the pair of empty strings, since that is
    where a fresh alignment begins. Memoization keyed on (i, j) avoids
    repeating work when different paths land on the same cell, which
    happens often once ties multiply across the table.
    """
    if memo is None:
        memo = {}

    if (i, j) in memo:
        return memo[(i, j)]

    results = set()

    for source in backtrack[i][j]:
        if source == "start":
            results.add(("", ""))
        elif source == "down":
            for prefix_v, prefix_w in all_alignments(backtrack, v, w, i - 1, j, memo):
                results.add((prefix_v + v[i - 1], prefix_w + "-"))
        elif source == "right":
            for prefix_v, prefix_w in all_alignments(backtrack, v, w, i, j - 1, memo):
                results.add((prefix_v + "-", prefix_w + w[j - 1]))
        else:  # diagonal
            for prefix_v, prefix_w in all_alignments(backtrack, v, w, i - 1, j - 1, memo):
                results.add((prefix_v + v[i - 1], prefix_w + w[j - 1]))

    memo[(i, j)] = results
    return results


def all_local_alignments(v, w, match_reward, mismatch_penalty, indel_penalty):
    """
    A recursive, top-down dynamic programming algorithm using memoization.
    Return the best local alignment score along with every
    distinct alignment that achieves it.

    Builds the scoring table (with the free ride option, and every
    tied source recorded per cell), finds every cell sharing the
    overall best score since more than one ending point can tie, and
    then collects every alignment reachable by exploring every tied
    direction from each of those ending points.
    """
    score, backtrack = build_score_table(v, w, match_reward, mismatch_penalty, indel_penalty)
    best_score, best_cells = find_best_ending_cells(score)

    memo = {}
    alignments = set()
    for i, j in best_cells:
        alignments |= all_alignments(backtrack, v, w, i, j, memo)

    return best_score, alignments


if __name__ == "__main__":
    sys.setrecursionlimit(10000)  # protects against deep recursion on long strings

    match_reward = 1
    mismatch_penalty = 1
    indel_penalty = 2

    v = "GAGA"
    w = "GAT"

    best_score, alignments = all_local_alignments(
        v, w, match_reward, mismatch_penalty, indel_penalty
    )

    print(f"best score: {best_score}")
    print(f"number of distinct optimal local alignments: {len(alignments)}")
    for aligned_v, aligned_w in sorted(alignments):
        print(aligned_v)
        print(aligned_w)
        print()