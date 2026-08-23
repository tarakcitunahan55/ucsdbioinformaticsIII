def south_or_east(i, j, down, right):
    """
    Recursively compute the length of the longest path from (0,0) to (i,j)
    in a Manhattan Tourist Problem grid.

    Parameters:
        i, j  - the coordinates of the destination node
        down  - a 2D grid of weights for vertical (south) edges;
                down[i][j] = weight of the edge INTO (i,j) coming from (i-1,j)
        right - a 2D grid of weights for horizontal (east) edges;
                right[i][j] = weight of the edge INTO (i,j) coming from (i,j-1)

    Returns:
        The length of the longest path from (0,0) to (i,j).
    """

    # Base case: we've reached the starting corner. No path needed, length 0.
    if i == 0 and j == 0:
        return 0

    # Start both candidate paths as "impossible" (negative infinity),
    # in case one direction isn't reachable (e.g. i=0 means no "south" option).
    x = float('-inf')
    y = float('-inf')

    # Option 1: arrive at (i,j) by moving SOUTH from (i-1, j).
    # Only valid if i > 0 (otherwise there's no row above to come from).
    if i > 0:
        x = south_or_east(i - 1, j, down, right) + down[i][j]

    # Option 2: arrive at (i,j) by moving EAST from (i, j-1).
    # Only valid if j > 0 (otherwise there's no column to the left to come from).
    if j > 0:
        y = south_or_east(i, j - 1, down, right) + right[i][j]

    # The longest path to (i,j) is whichever direction gives the bigger total.
    return max(x, y)


# ---- Example usage ----
if __name__ == "__main__":
    # A small 2x2 example grid (3 rows x 3 cols of nodes, i.e. n=2, m=2).
    # down[i][j]  = weight of vertical edge landing on (i,j)
    # right[i][j] = weight of horizontal edge landing on (i,j)
    # Row/col 0 of these arrays are unused placeholders where not applicable.

    down = [
        [0, 0, 0],   # edges into row 0 (unused, no vertical edges land in row 0)
        [1, 0, 2],   # edges into row 1: down[1][0]=1, down[1][1]=0, down[1][2]=2
        [4, 3, 3],   # edges into row 2
    ]

    right = [
        [0, 3, 2],   # edges into col 1,2 of row 0
        [0, 1, 5],   # edges into col 1,2 of row 1
        [0, 6, 5],   # edges into col 1,2 of row 2
    ]

    n, m = 2, 2   # destination node (bottom-right corner)
    result = south_or_east(n, m, down, right)
    print(f"Longest path from (0,0) to ({n},{m}): {result}")