"""Solve the fitting alignment problem for two amino acid strings.

A fitting alignment finds the best way to align an entire copy of the
shorter string w against some substring of the longer string v -- v is
allowed to be trimmed at both ends, but every character of w must be
used. This sits between global and local alignment: like local
alignment, part of the table gets a shortcut so v doesn't have to be
used from the very start; like global alignment, every character of w
still has to be accounted for.

Scoring for matches and mismatches uses the blosum62 substitution
matrix, which gives an amino-acid-pair-specific reward or penalty
instead of a single flat match reward and mismatch penalty.

Caveat: if there are multiple best alignments, it returns only one.
No affine gap penalties (opening and extending the gap/indels get the same penalty)
"""

import sys


blosum62_text = """
   A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
A  4  0 -2 -1 -2  0 -2 -1 -1 -1 -1 -2 -1 -1 -1  1  0  0 -3 -2
C  0  9 -3 -4 -2 -3 -3 -1 -3 -1 -1 -3 -3 -3 -3 -1 -1 -1 -2 -2
D -2 -3  6  2 -3 -1 -1 -3 -1 -4 -3  1 -1  0 -2  0 -1 -3 -4 -3
E -1 -4  2  5 -3 -2  0 -3  1 -3 -2  0 -1  2  0  0 -1 -2 -3 -2
F -2 -2 -3 -3  6 -3 -1  0 -3  0  0 -3 -4 -3 -3 -2 -2 -1  1  3
G  0 -3 -1 -2 -3  6 -2 -4 -2 -4 -3  0 -2 -2 -2  0 -2 -3 -2 -3
H -2 -3 -1  0 -1 -2  8 -3 -1 -3 -2  1 -2  0  0 -1 -2 -3 -2  2
I -1 -1 -3 -3  0 -4 -3  4 -3  2  1 -3 -3 -3 -3 -2 -1  3 -3 -1
K -1 -3 -1  1 -3 -2 -1 -3  5 -2 -1  0 -1  1  2  0 -1 -2 -3 -2
L -1 -1 -4 -3  0 -4 -3  2 -2  4  2 -3 -3 -2 -2 -2 -1  1 -2 -1
M -1 -1 -3 -2  0 -3 -2  1 -1  2  5 -2 -2  0 -1 -1 -1  1 -1 -1
N -2 -3  1  0 -3  0  1 -3  0 -3 -2  6 -2  0  0  1  0 -3 -4 -2
P -1 -3 -1 -1 -4 -2 -2 -3 -1 -3 -2 -2  7 -1 -2 -1 -1 -2 -4 -3
Q -1 -3  0  2 -3 -2  0 -3  1 -2  0  0 -1  5  1  0 -1 -2 -2 -1
R -1 -3 -2  0 -3 -2  0 -3  2 -2 -1  0 -2  1  5 -1 -1 -3 -3 -2
S  1 -1  0  0 -2  0 -1 -2  0 -2 -1  1 -1  0 -1  4  1 -2 -3 -2
T  0 -1 -1 -1 -2 -2 -2 -1 -1 -1 -1  0 -1 -1 -1  1  5  0 -2 -2
V  0 -1 -3 -2 -1 -3 -3  3 -2  1  1 -3 -2 -2 -3 -2  0  4 -3 -1
W -3 -2 -4 -3  1 -2 -2 -3 -3 -2 -1 -4 -4 -2 -3 -3 -2 -3 11  2
Y -2 -2 -3 -2  3 -3  2 -1 -2 -1 -1 -2 -3 -1 -2 -2 -2 -1  2  7
"""


def parse_blosum62(text):
    """Turn the raw blosum62 text block into a lookup table.

    Returns a dictionary of dictionaries, so that score_table["A"]["C"]
    gives the substitution score between amino acids a and c. Reading
    it this way avoids hand-typing a large nested dictionary directly
    into the code and keeps the matrix easy to double check against
    the original table.
    """
    lines = [line for line in text.strip().splitlines() if line.strip()]
    header = lines[0].split()
    score_table = {amino_acid: {} for amino_acid in header}

    for line in lines[1:]:
        parts = line.split()
        row_amino_acid = parts[0]
        values = [int(value) for value in parts[1:]]
        for column_amino_acid, value in zip(header, values):
            score_table[row_amino_acid][column_amino_acid] = value

    return score_table


blosum62 = parse_blosum62(blosum62_text)


def build_score_table(v, w, score_table, indel_penalty):
    """Fill in the scoring table and a matching backtrack table for
    fitting alignment.

    score[i][j] holds the best possible score of a fitting alignment
    between some suffix-trimmed, prefix-trimmed stretch of v ending at
    v[i-1], and the first j characters of w.

    Column zero gets a free ride, fixed at zero for every row, since a
    fitting alignment is allowed to begin anywhere along v -- none of
    w has been used yet at that point, so there's nothing to pay for.
    Row zero, on the other hand, is forced to pay an indel penalty for
    every step, since w must be used in full, and with zero characters
    of v available, every character of w has to be a gap.
    """
    n, m = len(v), len(w)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    backtrack = [[None] * (m + 1) for _ in range(n + 1)]

    for i in range(n + 1):
        score[i][0] = 0
        backtrack[i][0] = "start" #first column (index zero) gets a free ride and score zero as in local alignment problem

    for j in range(1, m + 1):
        score[0][j] = score[0][j - 1] - indel_penalty
        backtrack[0][j] = "right" #first row (index zero) has to pay for indel penalties as in global alignment problem

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal_score = score[i - 1][j - 1] + score_table[v[i - 1]][w[j - 1]]
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

                #local alignment has 0 as an edge into every cell, because both strings can restart. Fitting alignment only lets v restart, and only along column zero — once you're past column zero, w is locked in from that point forward, so there's no more zero-option to offer.

    return score, backtrack


def find_best_ending_row(score, m):
    """Search the entire last column for the highest score.

    A fitting alignment must use all of w, so it always ends at
    column m, but it's allowed to end at any row of v, since v can be
    trimmed from the end as well as the start. Every row in that last
    column is checked to find where the best score actually sits.
    """
    n = len(score) - 1
    best_score = score[0][m]
    best_i = 0

    for i in range(n + 1):
        if score[i][m] > best_score:
            best_score = score[i][m]
            best_i = i

    return best_score, best_i


def reconstruct_alignment(backtrack, v, w, i, j):
    """Recursively walk the backtrack table from cell (i, j) back
    toward wherever the fitting alignment actually begins along v.

    The recursion stops as soon as it reaches a cell marked "start",
    since that marks the point where v was trimmed and the alignment
    truly begins -- not necessarily the very first character of v.
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


def fitting_alignment(v, w, score_table, indel_penalty=1):
    """Return the best fitting alignment score along with one
    alignment that achieves it.

    Every character of w gets used, aligned against some contiguous
    stretch of v, which can start and end anywhere within v. The
    table is built with a free ride down column zero (representing
    where along v the alignment may begin) and a forced gap cost
    along row zero (since none of w may ever be skipped). The best
    score is then found by scanning the whole last column, since the
    alignment can end at any row of v, and reconstruction walks back
    only until reaching the cell marking where it truly started.
    """
    sys.setrecursionlimit(10000)  # protects against deep recursion on long strings

    score, backtrack = build_score_table(v, w, score_table, indel_penalty)
    best_score, best_i = find_best_ending_row(score, len(w))
    aligned_v, aligned_w = reconstruct_alignment(backtrack, v, w, best_i, len(w))
    return best_score, aligned_v, aligned_w


if __name__ == "__main__":
    v = "DISCREPANTLY" # longer one
    w = "PATENT" # shorter one
    # the order you put v,w or w,v is important

    best_score, aligned_v, aligned_w = fitting_alignment(v, w, blosum62)

    print(best_score)
    print(aligned_v)
    print(aligned_w)