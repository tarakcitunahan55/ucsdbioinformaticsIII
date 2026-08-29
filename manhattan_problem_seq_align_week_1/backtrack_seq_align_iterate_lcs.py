""" Given two DNA strings use dynamic programming to the find the longest common subsequence.
Iteration to get the output.
Caveat: if there are multiple LCSs, it returns only one.
No penalty for gaps and mismatches (since indels would occur heavily we can get very long but biologically not meaningful matches)
"""
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

            if s[i][j] == s[i - 1][j]:
                backtrack[i][j] = "down"
            elif s[i][j] == s[i][j - 1]:
                backtrack[i][j] = "right"
            elif s[i][j] == s[i - 1][j - 1] + match:
                backtrack[i][j] = "diag"

    return backtrack


def iterative_output_lcs(backtrack, v, w):
    """
    Iterative version of OutputLCS. Instead of recursing down to the
    base case and building the string on the way back up the call
    stack, this walks backward from (i, j) toward (0, 0) directly,
    building the LCS string as it goes.

    Note on indexing: the pseudocode's v[i] refers to the i-th character
    using 1-indexed convention (matching how backtrack[i][j] is indexed).
    In Python (0-indexed strings), that character is v[i - 1].
    """
    lcs = ""              # build the result here
    i = len(v)             # start at the bottom-right corner of the table
    j = len(w)

    # Keep walking backward until we hit either edge of the table
    # (matches the recursive base case: "if i = 0 or j = 0: return")
    while i > 0 and j > 0:
        if backtrack[i][j] == "down":
            # This character of v wasn't part of the LCS -- skip it.
            i -= 1
        elif backtrack[i][j] == "right":
            # This character of w wasn't part of the LCS -- skip it.
            j -= 1
        else:  # "diag"
            # v[i-1] IS part of the LCS. Since we're walking BACKWARD
            # (from the end toward the start), we must place this
            # character at the FRONT of what we've built so far --
            # otherwise the LCS would come out reversed.
            lcs = v[i - 1] + lcs
            i -= 1
            j -= 1

    return lcs


def longest_common_subsequence(v, w):
    backtrack = lcs_backtrack(v, w)
    return iterative_output_lcs(backtrack, v, w)


if __name__ == "__main__":

    s = "GAGGGTATAAAGAGCATGTAAAAAGCTCCAATATCGACTGCCATGTTGACTACTGACTAACTTGGTACGCGTTTGGGGGTATGATCCTTTGCAAGGGACGGGTCAGACTAATCAGGTCTTCACCCCTTAGTTACGGGCTCGCCCAAATAGAGATGATCTTTTCCCTCTGTACCCGAGTATTTCAAAATAAACTCCAGCGCCACTTTCGTTATACCTTATTAGGGTACTTAACACGTAGTTTAGGACCTTCATATGCATCCAGATGCAAGCCTGGTTGAAAAATGTATGATAGGCCTTATAGAGCGTGTGTTGATTTTCCAAGGCCGAACTGGAGTCGTCGCTGGAATTCCGTCTCCCAGGCACCTTAGGACATTTAAAAAATAAGCCCGGGATGACAGAATGCTTCAGTCTCTATGAATGAAACACTAACACCGCAACTGGGAATGATTAAGACAGTTATGGTACTCCCGCACTCCGGCCATACTCCACTGATCAAGCACCCGCCGTTCGGGTCGGTTCCACAGTTCGGCGCTATCCGATTCGTATAAAGAACGTTAATGCTCCGCTTGCGTGGTGGTCGGGCCAAAGCGGATCGCGTTTCCTATCTCAATTGAGGGAACGGCGATCCTAATGAGGTGTTTACTGCTGAAATGATGAGCCGCCCTTTACGGGTACCGCATCCCCGATAATGGAGGAGAATGCGAGGATCCAATCAACAGCAGAAGAGCTCGTCTACGCCAGTAGATGGTTCGGTACGTTAGGGCTGAGAAGTCAATTTAGTGTCGACCTCTGATGAACGAGATCGTTCGCGCACCCAAGGACATGGGGTGTGTCGTATAGGGATCCTAGGGCCAGCGAGATACCTAGGAACCTCGTTCGCTCCAATGCCCGCGGAAGGAACTGTTAGTCGATCTGTGGCATGAGGACGTGAGATACGTCTAGTCTTACGGTAAAACGTCGCTC"
    t = "CTACCATGAGGCACTTACGGCGGGGGTGTTCTACGAGCAACTGTAACCGCATCTCTGCAGAGGCGTGGATGTAATCACGGAGTCGGCGTCTAATCCATAGTCATCTCCCGTCTTCGGCATACACTTCAGCTGAAATCGTGCCAGCTGAGCCGCGCGGTTTTTTTGACTGTAAGAGAACTCCATTGCTGCTATCCCAGAACTCTTGAAAGTGCTTACGAGGACAGGACTGCGTTTCTACCAACAGCTGTCCAATATCATACCGCAGCGGTCGTTTGGCGGGCACCATTATTCAGGGCGTGAAGCCCCTACTGGCGGTCGGCGGGCAAAAAACGACTAAATGAGAGAAGATGTCCATACAACACAACGGATTTTTTGTCTGGACCGCTTTACAGGGAAAGGGGGCTTATTGGGACCCACATCCTCCTTTCTTATTCTTTGCGTTCCTTGACTACAAGCTACTGCCTCATTCTTTCACACAATGTAGTCAACTGGACCAGGGTTTTACAATACATCTTCTTGGCAATCGTAAAAGCCGTGTCAGCAAACGGTACTTGCAATTGCTGACGGAAGACCCTGGAGTCTTTATATTTTCTTTCGAGTAGACGTCTGGCCGCGACGTATTACCTTTGCTGACAGAAGCTCGTTATGTACGTGCGTAGTGCATTGCAGATGCTTGAGGATATTAGCGCAATTACATCTTGAGACGAAGGCGATATGGAATATCTGTCTGTGTTATTAGGATTCGATCACGATGATCATCACTGCACGAAACTTGTGGCGGACCAGCAAGAGCTATTCTTTACGTTGTAGAAAG"

    result = longest_common_subsequence(s, t)
    print(f"LCS: {result}")