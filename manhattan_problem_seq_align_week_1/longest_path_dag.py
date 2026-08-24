def parse_input(text):
    """
    Parses the input format:
        line 1: "source sink"
        remaining lines: "from to weight" (one edge per line)
    """
    lines = [line.strip() for line in text.strip().splitlines() if line.strip() != '']

    source, sink = map(int, lines[0].split())

    edges = []
    for line in lines[1:]:
        u, v, w = map(int, line.split())
        edges.append((u, v, w))

    return source, sink, edges


def longest_path_dag(source, sink, edges):
    """
    Finds the length of the longest path from 'source' to 'sink' in a DAG,
    along with the actual path, using dynamic programming over nodes in
    topological order.

    Key idea (same as ManhattanTourist): process nodes in an order where
    every predecessor is already finalized before we need it. Since node
    labels are given in increasing topological order, sorting nodes
    numerically gives us that order for free.
    """

    # Build a lookup: for each node v, which edges point INTO it (u -> v)?
    # This mirrors how ManhattanTourist looked up "down" and "right" edges
    # feeding into each cell -- here it's just a variable number of
    # incoming edges instead of always exactly two.
    incoming = {}
    all_nodes = set([source, sink])
    for u, v, w in edges:
        incoming.setdefault(v, []).append((u, w))
        all_nodes.add(u)
        all_nodes.add(v)

    # s[v] = length of the longest path from source to v.
    # backtrack[v] = which predecessor node achieved that best path.
    # Using a dict (not a list) since node labels might not start at 0
    # or be contiguous.
    s = {source: 0}
    backtrack = {}

    # Process nodes in increasing order -- this IS the topological order,
    # per the problem's guarantee. Skip anything before 'source' (irrelevant)
    # and stop once we've handled 'sink' (nothing beyond it matters).
    for node in sorted(all_nodes):
        if node <= source:
            continue   # source itself is already initialized to 0; anything before it is irrelevant
        if node > sink:
            break       # nothing past the sink matters for this problem

        best = float('-inf')
        best_predecessor = None

        # Look at every edge feeding INTO this node, exactly like checking
        # "down" and "right" edges into a grid cell -- just with a
        # variable-length list of incoming edges instead of always two.
        for (u, w) in incoming.get(node, []):
            if u in s:   # only consider predecessors that are actually reachable from source
                candidate = s[u] + w
                if candidate > best:
                    best = candidate
                    best_predecessor = u

        # Only record a value if at least one reachable incoming edge was found.
        # If a node has no reachable predecessor, it stays absent from s
        # (equivalent to "unreachable" / -infinity).
        if best_predecessor is not None:
            s[node] = best
            backtrack[node] = best_predecessor

    # Reconstruct the path by walking backtrack pointers from sink back to source.
    path = [sink]
    current = sink
    while current != source:
        current = backtrack[current]
        path.append(current)
    path.reverse()

    return s[sink], path


if __name__ == "__main__":
    sample_input = """0 49
43 45 9
3 32 5
39 42 16
16 36 3
46 48 12
9 15 19
15 33 18
7 39 9
34 39 9
21 25 1
39 49 18
26 31 3
37 38 13
33 46 7
21 42 18
7 29 2
6 23 3
38 42 14
31 43 19
43 47 7
31 40 19
32 45 3
41 42 9
36 38 17
26 49 15
20 31 10
29 38 15
13 23 16
46 49 15
2 19 4
3 38 2
9 24 9
27 38 10
39 46 17
46 47 3
18 48 5
28 48 6
33 47 19
36 43 10
13 23 6
36 37 12
13 39 8
28 33 2
13 42 10
35 46 2
40 49 6
42 48 14
4 10 17
14 21 20
7 15 18
38 43 15
23 24 2
1 18 8
9 43 8
13 40 17
35 37 2
2 24 11
14 38 13
30 33 6
13 26 8
47 49 10
0 3 20
11 43 11
27 49 16
23 34 10
37 42 4
17 30 5
7 9 20
0 24 20
34 36 19
8 37 12
26 38 4
31 44 10
41 46 7
17 24 19
40 46 8
26 42 7
45 48 16
7 32 2
34 43 10
7 47 20
11 28 6
1 29 7
11 39 18
16 17 14
29 31 19
5 27 8
5 13 9
16 21 11
23 49 2
0 5 1
21 32 1
20 23 9
12 45 16
0 1 11
1 2 3
2 4 11
4 6 1
6 7 20
7 8 19
8 11 17
11 12 5
12 14 17
14 16 12
16 20 3
20 35 17
35 41 14
"""

    source, sink, edges = parse_input(sample_input)
    length, path = longest_path_dag(source, sink, edges)

    print(length)
    print(" ".join(map(str, path)))