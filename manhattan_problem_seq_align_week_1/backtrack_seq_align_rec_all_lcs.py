""" Given two DNA strings use dynamic programming to the find the longest common subsequence(s) (LCS).
Memoization to get the output.
If there are multiple LCSs, it returns all
No penalty for gaps and mismatches (since indels would occur heavily we can get very long but biologically not meaningful matches)
"""

def lcs_backtrack_all(v, w):
    """
    Bottom-up Dynamic Programming
    Modified LCSBackTrack: instead of storing ONE direction per cell,
    store ALL directions that achieve the maximum value at that cell.
    This is what lets us later explore every optimal path, not just one.
    """
    n, m = len(v), len(w)
    s = [[0] * (m + 1) for _ in range(n + 1)]
    backtrack = [[[] for _ in range(m + 1)] for _ in range(n + 1)]

    for i in range(n + 1):
        s[i][0] = 0
    for j in range(m + 1):
        s[0][j] = 0

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = 1 if v[i - 1] == w[j - 1] else 0
            s[i][j] = max(s[i - 1][j], s[i][j - 1], s[i - 1][j - 1] + match)

            directions = []
            # Check EVERY option independently (not if/elif!) --
            # a cell can have 2 or even 3 valid directions tied for the max.
            if s[i][j] == s[i - 1][j]:
                directions.append("down")
            if s[i][j] == s[i][j - 1]:
                directions.append("right")
            # Only count "diag" if it's a REAL match (as we established
            # earlier, an unmatched diagonal can never actually tie for
            # the max -- but we check match explicitly here for clarity
            # and safety, rather than relying on that proof alone).
            if match == 1 and s[i][j] == s[i - 1][j - 1] + 1:
                directions.append("diag")

            backtrack[i][j] = directions

    return backtrack


def all_lcs(backtrack, v, i, j, memo=None):
    """
    A recursive, top-down dynamic programming algorithm using memoization.

    Without the memo dictionary, it would be pure recursion. With memo, it becomes memoized recursion / top-down DP.

    Modified OutputLCS: instead of following a single direction, this
    explores EVERY valid direction at each cell and collects every
    resulting LCS string into a set (which automatically removes
    duplicates, since different paths can sometimes produce the same
    string).

    Uses memoization keyed on (i, j), since many different paths through
    the table can revisit the exact same cell -- without caching, we'd
    redo the same exploration work repeatedly (the same "overlapping
    subproblems" issue we've seen in recursions, making them inefficient).
    Here once you've solved the subproblem (i, j), you store its answer.
    If another recursive path reaches (i, j), you don't calculate it again.
    """
    if memo is None:
        memo = {}

    if i == 0 or j == 0: #Base case: hit the edge of the table, nothing left to build with — return a set containing just the empty string (not an empty set — empty string is needed so merging/adding works correctly one level up).
        return {""}

    if (i, j) in memo: #If we've already solved "all LCSs reachable from (i,j)" before, reuse it instantly instead of recomputing.
        return memo[(i, j)] #dictionary access using a tuple as the key.

    results = set()

    for direction in backtrack[i][j]:
        if direction == "down": #down/right → recurse and merge in the whole set of strings unchanged (no character consumed).
            results |= all_lcs(backtrack, v, i - 1, j, memo) # (|=) Union merges two sets without deleting what's already there or use results.update(any iterable)
        elif direction == "right":
            results |= all_lcs(backtrack, v, i, j - 1, memo)
        else:  # "diag" 
            for suffix in all_lcs(backtrack, v, i - 1, j - 1, memo): #diag → recurse to (i-1, j-1), and for every string that comes back, append v[i-1] (the matched character) before adding it in.
                results.add(suffix + v[i - 1])

    memo[(i, j)] = results
    return results


def all_longest_common_subsequences(v, w):
    """Convenience wrapper: builds the backtrack matrix and finds every LCS."""
    backtrack = lcs_backtrack_all(v, w)
    return all_lcs(backtrack, v, len(v), len(w))


if __name__ == "__main__":
    s = "AACCTTGG"
    t = "ACACTGTGA"

    results = all_longest_common_subsequences(s, t)
    print(f"Number of distinct longest common subsequences: {len(results)}")
    for lcs in sorted(results):
        print(lcs)