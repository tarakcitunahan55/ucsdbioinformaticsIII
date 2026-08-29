"""Solve the local alignment problem for two nucleotide strings.

Given a match reward, a mismatch penalty, an indel penalty, and two
strings, this builds a scoring table using dynamic programming, then
uses recursion to walk the table backwards and recover one local
alignment that achieves the best possible score.

Unlike global alignment, a local alignment does not have to use the
entire length of either string -- it can be a shorter, well-matching
stretch found anywhere inside them. This is modeled by giving every
cell in the table a "free ride" back to a score of zero, representing
the option of simply starting a new alignment at that point instead
of extending a poor one.

Caveat: if there are multiple best local alignments, this returns
only one.
No affine gap penalties (opening and extending the gap/indels get the same penalty)
"""

import sys


def build_score_table(v, w, match_reward, mismatch_penalty, indel_penalty):
    """Fill in the scoring table and a matching backtrack table.

    score[i][j] holds the best possible score of a local alignment
    ending at v[i-1] paired with w[j-1] -- or, thanks to the free
    ride option, possibly just a score of zero if starting over here
    is better than continuing whatever came before.

    Four moves are possible when filling each cell:
      - move down: align v's next character with a gap in w
      - move right: align w's next character with a gap in v
      - move diagonally: align v's next character with w's next one,
        earning the match reward if they agree, or paying the
        mismatch penalty if they don't
      - take the free ride: reset the score to zero, treating this
        cell as the start of a fresh alignment
    """
    n, m = len(v), len(w)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    backtrack = [[None] * (m + 1) for _ in range(n + 1)]

    # the first row and column (index zero) have nowhere to backtrack to, so every cell along them is simply treated as a starting point
    #Since indel penalties only ever subtract from the score, dragging in leading gaps to reach row 0 or column 0 can never beat just starting fresh with a score of zero — so every cell in 1st row/column (index zero) there is always a valid "start," and score zero, never a forced gap step.
    for i in range(n + 1):
        backtrack[i][0] = "start" 
    for j in range(m + 1):
        backtrack[0][j] = "start"

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

            # picking one source that achieves the best value; the
            # free ride is checked last, so it only gets chosen when
            # nothing else reaches as high a score as starting fresh
            if best == down_score:
                backtrack[i][j] = "down"
            elif best == right_score:
                backtrack[i][j] = "right"
            elif best == diagonal_score:
                backtrack[i][j] = "diagonal"
            else:
                backtrack[i][j] = "start"

    return score, backtrack


def find_best_ending_cell(score):
    """Search every cell of the table for the highest score.

    A local alignment is allowed to end anywhere, not just at the
    bottom right corner the way a global alignment must, so the whole
    table needs to be scanned to find where the best score actually
    lives.
    """
    n = len(score) - 1
    m = len(score[0]) - 1

    best_score = score[0][0]
    best_i, best_j = 0, 0

    for i in range(n + 1):
        for j in range(m + 1):
            if score[i][j] > best_score:
                best_score = score[i][j]
                best_i, best_j = i, j

    return best_score, best_i, best_j


def reconstruct_alignment(backtrack, v, w, i, j):
    """Recursively walk the backtrack table from cell (i, j) back
    toward wherever the alignment actually began.

    Unlike global alignment, this does not necessarily walk all the
    way back to row zero and column zero. It stops as soon as it
    reaches a cell marked "start", since that is the point where the
    free ride was taken -- the true beginning of this alignment,
    which can be somewhere in the middle of the original strings.
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


def local_alignment(v, w, match_reward, mismatch_penalty, indel_penalty):
    """Return the best local alignment score along with one alignment
    that achieves it.

    Builds the scoring table (with the free ride option baked into
    every cell), searches the whole table for the highest score since
    a local alignment can end anywhere, then reconstructs the
    alignment by walking backward only as far as the cell where that
    best-scoring alignment actually started.
    """
    score, backtrack = build_score_table(v, w, match_reward, mismatch_penalty, indel_penalty)
    best_score, best_i, best_j = find_best_ending_cell(score)
    aligned_v, aligned_w = reconstruct_alignment(backtrack, v, w, best_i, best_j)
    return best_score, aligned_v, aligned_w


if __name__ == "__main__":
    sys.setrecursionlimit(10000)  # protects against deep recursion on long strings

    match_reward = 1
    mismatch_penalty = 1
    indel_penalty = 2

    v = "GAGA"
    w = "GAT"

    best_score, aligned_v, aligned_w = local_alignment(
        v, w, match_reward, mismatch_penalty, indel_penalty
    )

    print(best_score)
    print(aligned_v)
    print(aligned_w)