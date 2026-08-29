""" Given two DNA strings use dynamic programming to the find the longest common subsequence (LCS).
Recursion to get the output.
Caveat: if there are multiple LCSs, it returns only one.
No penalty for gaps and mismatches (since indels would occur heavily we can get very long but biologically not meaningful matches)
"""
import sys

def lcs_backtrack(v, w):
    """
    Bottom-up Dynamic Programming
    Builds both the DP table (s) and the backtrack matrix, following
    the exact same edge-weight logic as ManhattanTourist -- except here,
    the "diagonal" edge only has weight 1 if the characters match
    (weight 0 otherwise), while "down" and "right" edges always have
    weight 0. So s[i][j] is really just tracking the LCS length as we
    take the best of three incoming edges instead of two.
    """
    n, m = len(v), len(w)

    # s[i][j] = length of LCS of v[0:i] and w[0:j]
    s = [[0] * (m + 1) for _ in range(n + 1)]

    # backtrack[i][j] will hold one of "down", "right", "diag"
    backtrack = [[None] * (m + 1) for _ in range(n + 1)]

    # Initialize first column: s[i][0] = 0 for all i (matches "for i <- 0 to |v|: s_i,0 <- 0")
    for i in range(n + 1):
        s[i][0] = 0

    # Initialize first row: s[0][j] = 0 for all j (matches "for j <- 0 to |w|: s_0,j <- 0")
    for j in range(m + 1):
        s[0][j] = 0

    # Main fill loop
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # match = 1 if the characters at this position agree, else 0
            match = 1 if v[i - 1] == w[j - 1] else 0

            # Take the best of three options: come from above, from the
            # left, or diagonally (with a +1 bonus only if characters match)
            s[i][j] = max(
                s[i - 1][j],
                s[i][j - 1],
                s[i - 1][j - 1] + match
            )

            # Tie-break order is down -> right -> diag, so "diag" is only chosen
            # when it's the strict unique max. Since s[i-1][j] and s[i][j-1] are
            # both >= s[i-1][j-1] (monotonicity), a non-match (match=0) can never
            # make diag the strict winner -- so diag is only ever picked on an
            # actual character match. On matches that tie with down/right, this
            # code favors down/right, meaning some valid matches get skipped in
            # the reconstruction. LCS length is unaffected either way; only the
            # particular LCS string returned (of possibly several optimal ones)
            # can differ.

            if s[i][j] == s[i - 1][j]:
                backtrack[i][j] = "down"
            elif s[i][j] == s[i][j - 1]:
                backtrack[i][j] = "right"
            elif s[i][j] == s[i - 1][j - 1] + match:
                backtrack[i][j] = "diag"

    return backtrack


def output_lcs(backtrack, v, i, j):
    """Walks the backtrack matrix to reconstruct the LCS string.
    Recursive Algorithm is still efficient here as output_lcs never branches — it makes exactly one recursive call per invocation
    """
    if i == 0 or j == 0:
        return ""
    if backtrack[i][j] == "down":
        return output_lcs(backtrack, v, i - 1, j)
    elif backtrack[i][j] == "right":
        return output_lcs(backtrack, v, i, j - 1)
    else: # "diag"
        return output_lcs(backtrack, v, i - 1, j - 1) + v[i - 1]


def longest_common_subsequence(v, w):
    backtrack = lcs_backtrack(v, w)
    return output_lcs(backtrack, v, len(v), len(w))


if __name__ == "__main__":
    sys.setrecursionlimit(10000) # otherwise our recursion depth would be insufficient and may not run

    s = "GCGATC"   
    t = "CTGACG"

    result = longest_common_subsequence(s, t)
    print(f"LCS: {result}")