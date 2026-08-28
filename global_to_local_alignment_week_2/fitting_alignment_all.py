"""Solve the fitting alignment problem for two amino acid strings, but
find every fitting alignment that achieves the best possible score,
not just one.

A fitting alignment aligns all of w against some contiguous, trimmed
stretch of v. Just like earlier "all alignments" versions, this
records every tied direction at every cell (instead of picking one),
finds every row where the best score could end (instead of picking
the first one), and then uses memoized recursion to explore every
combination of ties and collect every distinct resulting alignment.

Scoring for matches and mismatches uses the blosum62 substitution
matrix, which gives an amino-acid-pair-specific reward or penalty
instead of a single flat match reward and mismatch penalty.
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
    gives the substitution score between amino acids a and c.
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
    between some trimmed stretch of v ending at v[i-1], and the first
    j characters of w.

    Unlike a single-path version, backtrack[i][j] stores a list of
    every direction that reaches the best score at that cell, not
    just one, so ties can carry two or even three recorded directions.

    Column zero gets a free ride, fixed at zero for every row, since a
    fitting alignment is allowed to begin anywhere along v. Row zero
    is forced to pay an indel penalty for every step, since all of w
    must be used, and with zero characters of v available, every
    character of w has to be a gap.
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
            diagonal_score = score[i - 1][j - 1] + score_table[v[i - 1]][w[j - 1]]
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


def find_best_ending_rows(score, m):
    """Find every row in the last column that achieves the overall
    best score.

    A fitting alignment must use all of w, so it always ends at
    column m, but it can end at any row of v, and more than one row
    can share the same top score -- every one of those tied rows is a
    valid ending point.
    """
    n = len(score) - 1
    best_score = score[0][m]

    for i in range(n + 1):
        if score[i][m] > best_score:
            best_score = score[i][m]

    best_rows = [i for i in range(n + 1) if score[i][m] == best_score]

    return best_score, best_rows


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


def all_fitting_alignments(v, w, score_table, indel_penalty=1):
    """Return the best fitting alignment score along with every
    distinct alignment that achieves it.

    Builds the scoring table with every tied direction recorded per
    cell, finds every row in the last column sharing the overall best
    score since more than one ending row can tie, and then collects
    every alignment reachable by exploring every tied direction from
    each of those ending rows.
    """
    sys.setrecursionlimit(10000)  # protects against deep recursion on long strings

    score, backtrack = build_score_table(v, w, score_table, indel_penalty)
    best_score, best_rows = find_best_ending_rows(score, len(w))

    memo = {}
    alignments = set()
    for i in best_rows:
        alignments |= all_alignments(backtrack, v, w, i, len(w), memo)

    return best_score, alignments


if __name__ == "__main__":
    v = "DISCREPANTLY" # longer one
    w = "PATENT" # shorter one
    # the order you put v,w or w,v is important

    best_score, alignments = all_fitting_alignments(v, w, blosum62)

    print(f"best score: {best_score}")
    print(f"number of distinct optimal fitting alignments: {len(alignments)}")
    for aligned_v, aligned_w in sorted(alignments):
        print(aligned_v)
        print(aligned_w)
        print()