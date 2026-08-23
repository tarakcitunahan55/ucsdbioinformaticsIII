def manhattan_tourist(n, m, down, right):
    """
    Compute the length of the longest path from (0,0) to (n,m) in a Manhattan Tourist Problem grid, 
    using dynamic programming.

    Parameters:
        n, m  - grid dimensions (n rows of vertical moves, m columns of horizontal moves)
        down  - n x (m+1) matrix; down[i][j] = weight of the vertical edge
                from node (i, j) DOWN to node (i+1, j)
        right - (n+1) x m matrix; right[i][j] = weight of the horizontal edge
                from node (i, j) RIGHT to node (i, j+1)

    Returns:
        The length of the longest path from (0,0) to (n,m).
    """

    # s[i][j] will store the longest path length from (0,0) to (i,j).
    # Build an (n+1) x (m+1) table, one entry per grid node.
    s = [[0] * (m + 1) for _ in range(n + 1)]

    # s[0][0] = 0: the starting node requires 0 length to "reach itself".
    s[0][0] = 0

    # Fill in the FIRST COLUMN (j = 0): the only way to reach (i, 0) is
    # by moving straight down from (i-1, 0) -- there's no "left" option
    # since j can't go below 0. This mirrors the first for-loop in the
    # pseudocode.
    for i in range(1, n + 1):
        s[i][0] = s[i - 1][0] + down[i - 1][0]

    # Fill in the FIRST ROW (i = 0): similarly, the only way to reach
    # (0, j) is by moving right from (0, j-1) -- no "down" option exists
    # since i can't go below 0.
    for j in range(1, m + 1):
        s[0][j] = s[0][j - 1] + right[0][j - 1]

    # Fill in the REST OF THE TABLE. Every other node (i, j) can be reached
    # from EITHER above (moving down) OR from the left (moving right).
    # We take whichever gives the longer total path, exactly like the
    # "compare and keep the better one" pattern from DPChange.
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            from_above = s[i - 1][j] + down[i - 1][j]
            from_left = s[i][j - 1] + right[i][j - 1]
            s[i][j] = max(from_above, from_left)

    # The answer is the longest path to the final destination node (n, m).
    return s[n][m]


def parse_input(text):
    """
    Parses the problem's input format:
        line 1: "n m"
        next n lines: Down matrix rows (each with m+1 integers)
        a line containing just "-"
        next (n+1) lines: Right matrix rows (each with m integers)
    """
    lines = [line.strip() for line in text.strip().splitlines() if line.strip() != '']

    n, m = map(int, lines[0].split())

    down = []
    idx = 1
    for _ in range(n):
        down.append(list(map(int, lines[idx].split())))
        idx += 1

    # skip the "-" separator line
    idx += 1

    right = []
    for _ in range(n + 1):
        right.append(list(map(int, lines[idx].split())))
        idx += 1

    return n, m, down, right


# ---- Example usage with the sample input ----
if __name__ == "__main__":
    sample_input = """18 16
1 1 0 3 1 2 1 2 0 0 2 2 2 4 3 4 4
0 3 4 3 1 3 4 2 3 2 1 3 3 2 0 3 2
4 3 4 0 0 0 1 2 2 0 4 2 3 3 0 1 3
3 4 0 1 2 4 4 0 3 1 0 0 2 3 0 2 4
0 3 3 2 1 0 2 4 4 4 4 3 1 0 1 3 2
1 0 1 1 2 2 0 1 2 0 4 1 0 4 4 3 4
1 0 3 4 1 0 2 3 3 0 0 0 0 2 2 4 2
0 2 3 4 4 2 2 2 2 2 0 4 1 4 0 2 2
0 4 0 2 2 4 2 0 3 3 2 2 1 0 1 3 4
4 0 4 1 1 3 2 0 2 1 3 3 2 0 1 2 2
1 3 4 3 1 4 4 1 3 4 0 3 0 1 0 1 2
3 4 4 2 1 0 1 3 2 4 2 4 1 0 0 0 0
2 4 1 1 1 3 0 3 0 1 1 3 2 0 2 2 3
2 4 4 4 2 2 0 4 2 1 3 2 2 4 1 0 0
3 1 4 4 2 3 0 3 3 1 2 1 0 1 3 2 4
4 4 1 1 1 1 1 2 3 4 2 2 3 3 2 4 4
0 3 1 3 4 2 1 2 4 4 3 4 1 0 0 2 3
4 1 2 2 4 3 2 3 1 1 3 0 1 3 4 1 1
-
4 2 4 3 0 2 1 3 4 0 3 1 0 1 2 1
2 4 3 2 4 4 1 1 4 1 0 0 1 2 4 4
1 4 0 2 0 3 2 4 4 4 3 0 4 0 3 2
2 4 3 3 0 0 1 2 4 3 1 1 3 2 1 3
2 1 3 4 4 4 1 4 3 4 3 2 2 3 2 1
4 1 1 4 3 0 2 0 4 3 2 1 0 0 1 4
0 2 3 3 1 0 2 1 0 3 3 0 1 3 2 1
2 3 1 4 2 4 3 2 3 2 1 1 4 2 0 4
0 4 3 0 3 0 4 0 4 0 3 0 2 2 3 4
2 2 2 1 0 2 1 1 1 2 3 2 2 1 1 3
0 0 1 3 4 3 0 1 2 0 4 3 2 2 0 2
4 2 3 0 2 2 4 4 0 2 4 1 0 2 2 0
0 0 0 2 0 0 4 4 0 4 1 0 3 4 3 2
3 2 1 4 4 0 3 4 2 4 1 2 0 1 4 1
2 3 2 1 1 0 0 4 4 3 3 4 2 4 3 4
3 2 3 3 3 0 2 0 3 4 1 2 2 2 2 2
0 1 2 4 0 0 3 2 2 3 2 1 0 4 3 1
4 1 0 2 3 2 3 0 1 0 2 1 2 4 1 0
1 4 0 3 4 3 2 3 1 3 4 3 4 1 2 4
"""

    n, m, down, right = parse_input(sample_input)
    result = manhattan_tourist(n, m, down, right)
    print(f"Longest path length: {result}")