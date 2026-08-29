"""
Solve the multiple longest common subsequence for three strings.

the score of a column in a multiple global alignment is defined as:
    1  if all three symbols in that column are identical
    0  if at least one symbol disagrees (including a symbol being a gap (indel))

Maximizing the total score over an alignment is equivalent to finding the
longest common subsequence (lcs) shared by all three strings, since only
columns where all three strings contribute the same matching character
count toward the score.

this is solved with a 3-dimensional dynamic programming table, generalizing
the classic 2-string lcs recurrence to 3 strings. from any cell (i, j, k)
there are 7 possible "predecessor" moves (every non-empty subset of
{decrement i, decrement j, decrement k}), and only the move that decrements
all three indices together can add 1 to the score, and only when the three
consumed characters are all equal.

Caveat: We don't worry about gaps or mismatches at all (no gap or mismatch penalty) -> simply finds the maximum number of matches bt. three seqs
Since indels would occur heavily, we can get very long but biologically not meaningful matches
"""


def multiple_lcs(v, w, u):
    """
    compute the length of the longest common subsequence shared by three
    strings v, w, u, and reconstruct one optimal 3-way alignment achieving
    that score.

    parameters:
        v, w, u -- three dna strings

    returns:
        (lcs_length, aligned_v, aligned_w, aligned_u)
    """
    n1, n2, n3 = len(v), len(w), len(u)

    # dp[i][j][k] = length of the lcs of v[0:i], w[0:j], u[0:k]
    dp = [[[0] * (n3 + 1) for _ in range(n2 + 1)] for _ in range(n1 + 1)]

    # backpointer[i][j][k] stores which move produced dp[i][j][k], as a
    # tuple (di, dj, dk) of how much each index was decremented to get here
    backpointer = [[[None] * (n3 + 1) for _ in range(n2 + 1)] for _ in range(n1 + 1)]

    for i in range(n1 + 1):
        for j in range(n2 + 1):
            for k in range(n3 + 1):
                if i == 0 and j == 0 and k == 0:
                    continue

                best_value = -1
                best_move = None

                # try the move that consumes one character from all three
                # strings at once -- the only move that can score a point
                if i > 0 and j > 0 and k > 0:
                    bonus = 1 if v[i - 1] == w[j - 1] == u[k - 1] else 0
                    candidate = dp[i - 1][j - 1][k - 1] + bonus
                    if candidate > best_value:
                        best_value = candidate
                        best_move = (1, 1, 1)

                # try the three moves that consume characters from exactly
                # two of the three strings (never adds to the score)
                if i > 0 and j > 0:
                    candidate = dp[i - 1][j - 1][k]
                    if candidate > best_value:
                        best_value = candidate
                        best_move = (1, 1, 0)
                if i > 0 and k > 0:
                    candidate = dp[i - 1][j][k - 1]
                    if candidate > best_value:
                        best_value = candidate
                        best_move = (1, 0, 1)
                if j > 0 and k > 0:
                    candidate = dp[i][j - 1][k - 1]
                    if candidate > best_value:
                        best_value = candidate
                        best_move = (0, 1, 1)

                # try the three moves that consume a character from just
                # one of the three strings (never adds to the score)
                if i > 0:
                    candidate = dp[i - 1][j][k]
                    if candidate > best_value:
                        best_value = candidate
                        best_move = (1, 0, 0)
                if j > 0:
                    candidate = dp[i][j - 1][k]
                    if candidate > best_value:
                        best_value = candidate
                        best_move = (0, 1, 0)
                if k > 0:
                    candidate = dp[i][j][k - 1]
                    if candidate > best_value:
                        best_value = candidate
                        best_move = (0, 0, 1)

                dp[i][j][k] = best_value
                backpointer[i][j][k] = best_move

    lcs_length = dp[n1][n2][n3]

    aligned_v, aligned_w, aligned_u = _traceback(v, w, u, backpointer, n1, n2, n3)

    return lcs_length, aligned_v, aligned_w, aligned_u


def _traceback(v, w, u, backpointer, i, j, k):
    """
    walk backwards from (n1, n2, n3) to (0, 0, 0) using the stored
    backpointers, building the three aligned strings one column at a time.
    a '-' is placed in a string's row whenever that string's index was not
    decremented by the chosen move (i.e. it did not contribute a real
    character to that column).
    """
    aligned_v = []
    aligned_w = []
    aligned_u = []

    while i > 0 or j > 0 or k > 0:
        di, dj, dk = backpointer[i][j][k]

        aligned_v.append(v[i - 1] if di else '-')
        aligned_w.append(w[j - 1] if dj else '-')
        aligned_u.append(u[k - 1] if dk else '-')

        i -= di
        j -= dj
        k -= dk

    aligned_v.reverse()
    aligned_w.reverse()
    aligned_u.reverse()

    return ''.join(aligned_v), ''.join(aligned_w), ''.join(aligned_u)


v = "ATATCGG"
w = "TCCGA"
u = "ATGTACTG"

lcs_length, aligned_v, aligned_w, aligned_u = multiple_lcs(v, w, u)

print(lcs_length)
print(aligned_v)
print(aligned_w)
print(aligned_u)