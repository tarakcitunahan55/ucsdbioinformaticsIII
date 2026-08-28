"""Solve the edit distance problem for two strings.

Edit distance counts the "fewest" single-character insertions,
deletions, and substitutions needed to turn one string into the
other. This can be seen as a close cousin of "global alignment": the
same three moves (down, right, diagonal) are available at every cell,
but instead of maximizing a reward/match, this minimizes a cost, where a
mismatch costs one edit, a match costs nothing, and every gap
(insertion or deletion) costs one edit as well.

"""

import sys


def build_distance_table(v, w):
    """Fill in the edit distance table using dynamic programming.

    distance[i][j] holds the fewest edits needed to turn the first i
    characters of v into the first j characters of w.

    Three moves are possible when filling each cell:
      - move down: delete v's next character (costs one edit)
      - move right: insert w's next character (costs one edit)
      - move diagonally: substitute v's next character with w's next
        one if they differ (costs one edit), or simply line them up
        for free if they already agree (costs nothing)
    """
    n, m = len(v), len(w)
    distance = [[0] * (m + 1) for _ in range(n + 1)]

    # turning the first i characters of v into an empty string takes exactly i deletions
    for i in range(n + 1):
        distance[i][0] = i

    # turning an empty string into the first j characters of w takes exactly j insertions
    for j in range(m + 1):
        distance[0][j] = j

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            substitution_cost = 0 if v[i - 1] == w[j - 1] else 1

            deletion = distance[i - 1][j] + 1
            insertion = distance[i][j - 1] + 1
            substitution = distance[i - 1][j - 1] + substitution_cost

            # unlike alignment, which maximizes a reward, this problem
            # minimizes a cost, so the best choice at each cell is the
            # smallest of the three options rather than the largest
            distance[i][j] = min(deletion, insertion, substitution)

    return distance


def edit_distance(v, w):
    """Return the edit distance between two strings, i.e. the fewest
    single-character insertions, deletions, and substitutions needed
    to turn one string into the other.
    """
    distance = build_distance_table(v, w)
    return distance[len(v)][len(w)] # by the time you reach distance[len(v)][len(w)], you've fully transformed all of v, character by character, into all of w.
#That said, edit distance is symmetric — edit_distance(v, w) == edit_distance(w, v) — since turning v into w takes exactly as many edits as turning w into v (an insertion in one direction is just a deletion in the other, and substitutions work the same both ways). So while w is the "target" in the way the code is written, the actual numeric answer wouldn't change if you swapped which string you called v and which you called w



if __name__ == "__main__":
    sys.setrecursionlimit(10000)  # kept for consistency with earlier problems in this set

    v = "GAGA"
    w = "GAT"

    print(edit_distance(v, w))