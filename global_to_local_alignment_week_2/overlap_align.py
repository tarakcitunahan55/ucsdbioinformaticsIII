"""Solve the overlap alignment problem for two nucleotide strings.

An overlap alignment looks for the best way to align some suffix v'
of v against some prefix w' of w. This models, for example, checking
whether the tail end of one sequencing read lines up with the front
end of another. It sits between global and fitting alignment: v is
allowed to be trimmed, but only at its start (since v' is a suffix,
it must run all the way to the end of v); w is allowed to be
trimmed, but only at its end (since w' is a prefix, it must start
right at the beginning of w).

Caveat: if there are multiple best overlap alignments, this returns only one of them.
"""

import sys


def build_score_table(v, w, match_reward, mismatch_penalty, indel_penalty):
    """Fill in the scoring table and a matching backtrack table for
    overlap alignment.

    score[i][j] holds the best possible score of an alignment between
    some suffix of v ending at v[i-1], and the first j characters of
    w (i.e. a genuine prefix of w, since w' is never allowed to skip
    its own beginning).

    Column zero gets a free ride, fixed at zero for every row, since
    the alignment is allowed to begin anywhere along v -- none of w
    has been used yet at that point, so there's nothing to pay for.
    Row zero, on the other hand, is forced to pay an indel penalty
    for every step, since w' must start from the very first
    character of w, so any characters of w used with zero characters
    of v have to be paid for as gaps, the same as in global alignment.
    """
    n, m = len(v), len(w)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    backtrack = [[None] * (m + 1) for _ in range(n + 1)]

    for i in range(n + 1):
        score[i][0] = 0
        backtrack[i][0] = "start"

    for j in range(1, m + 1):
        score[0][j] = score[0][j - 1] - indel_penalty
        backtrack[0][j] = "right"

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

            # picking one direction that achieves the best value;
            # diagonal is checked last so it's only chosen when it's
            # strictly the best among the three real moves
            if best == down_score:
                backtrack[i][j] = "down"
            elif best == right_score:
                backtrack[i][j] = "right"
            else:
                backtrack[i][j] = "diagonal"

    return score, backtrack


def find_best_ending_column(score, n):
    """Search the entire last row for the highest score.

    v' must be a suffix of v, so the alignment always ends at row n
    (all of v used, from wherever it started). But w' only needs to
    be a prefix of w, so it's free to stop early -- meaning the
    alignment can end at any column of that last row. Every column in
    row n is checked to find where the best score actually sits.
    """
    m = len(score[0]) - 1
    best_score = score[n][0]
    best_j = 0

    for j in range(m + 1):
        if score[n][j] > best_score:
            best_score = score[n][j]
            best_j = j

    return best_score, best_j


def reconstruct_alignment(backtrack, v, w, i, j):
    """Recursively walk the backtrack table from cell (i, j) back
    toward wherever the overlap alignment actually begins along v.

    The recursion stops as soon as it reaches a cell marked "start",
    since that marks the point where v was trimmed and the
    alignment truly begins -- not necessarily the very first
    character of v.
    """
    direction = backtrack[i][j]

    if direction == "start":
        return "", ""

    if direction == "down":
        prefix_v, prefix_w = reconstruct_alignment(backtrack, v, w, i - 1, j)
        return prefix_v + v[i - 1], prefix_w + "-"
    elif direction == "right":
        prefix_v, prefix_w = reconstruct_alignment(backtrack, v, w, i, j - 1)
        return prefix_v + "-", prefix_w + w[j - 1]
    else:  # diagonal
        prefix_v, prefix_w = reconstruct_alignment(backtrack, v, w, i - 1, j - 1)
        return prefix_v + v[i - 1], prefix_w + w[j - 1]


def overlap_alignment(v, w, match_reward, mismatch_penalty, indel_penalty):
    """Return the best overlap alignment score along with one
    alignment that achieves it.

    A suffix of v gets aligned against a prefix of w. The table is
    built with a free ride down column zero (representing where
    along v the alignment may begin) and a forced gap cost along row
    zero (since w' must start from its own beginning). The best score
    is then found by scanning the whole last row, since w' is allowed
    to end early, and reconstruction walks back only until reaching
    the cell marking where the alignment truly started along v.

    
    """
    sys.setrecursionlimit(10000)  # protects against deep recursion on long strings

    score, backtrack = build_score_table(v, w, match_reward, mismatch_penalty, indel_penalty)
    best_score, best_j = find_best_ending_column(score, len(v))
    aligned_v, aligned_w = reconstruct_alignment(backtrack, v, w, len(v), best_j)
    return best_score, aligned_v, aligned_w


if __name__ == "__main__":
    match_reward = 1
    mismatch_penalty = 1
    indel_penalty = 2

    v = "GAGA" # contributes a suffix
    w = "GAT" # contributes a prefix
    # the order you put v,w or w,v is important


    best_score, aligned_v, aligned_w = overlap_alignment(
        v, w, match_reward, mismatch_penalty, indel_penalty
    )

    print(best_score)
    print(aligned_v)
    print(aligned_w)