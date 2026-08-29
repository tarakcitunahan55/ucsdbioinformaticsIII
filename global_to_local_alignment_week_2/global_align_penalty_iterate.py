"""Solve the global alignment problem for two nucleotide strings.

Given a match reward, a mismatch penalty, an indel penalty, and two
strings, this builds a scoring table using dynamic programming, then
walks the table backwards "iteratively "to recover one alignment that achieves the
best possible score.
Caveat: if there are multiple best alignments, it returns only one.
No affine gap penalties (opening and extending the gap/indels get the same penalty)
"""


def build_score_table(v, w, match_reward, mismatch_penalty, indel_penalty):
    """Fill in the scoring table and a matching backtrack table.

    score[i][j] holds the best possible alignment score between the
    first i characters of v and the first j characters of w.

    Three moves are possible when filling each cell:
      - move down: align v's next character with a gap in w
      - move right: align w's next character with a gap in v
      - move diagonally: align v's next character with w's next one,
        earning the match reward if they agree, or paying the
        mismatch penalty if they don't
    """
    n, m = len(v), len(w)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    backtrack = [[None] * (m + 1) for _ in range(n + 1)]

    # filling the first column: aligning some prefix of v with nothing
    # from w costs an indel penalty for every character used
    for i in range(1, n + 1):
        score[i][0] = score[i - 1][0] - indel_penalty
        backtrack[i][0] = "down"

    # filling the first row: same idea, but for prefixes of w
    for j in range(1, m + 1):
        score[0][j] = score[0][j - 1] - indel_penalty
        backtrack[0][j] = "right"

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if v[i - 1] == w[j - 1]:
                diagonal_score = score[i - 1][j - 1] + match_reward
            else:
                diagonal_score = score[i - 1][j - 1] - mismatch_penalty

            down_score = score[i - 1][j] - indel_penalty
            right_score = score[i][j - 1] - indel_penalty

            best = max(down_score, right_score, diagonal_score)
            score[i][j] = best

            # picking one direction that achieves the best value;
            # diagonal is checked last so it's only chosen when it's
            # strictly the best, avoiding an unearned free alignment
            if best == down_score:
                backtrack[i][j] = "down"
            elif best == right_score:
                backtrack[i][j] = "right"
            else:
                backtrack[i][j] = "diagonal"

    return score, backtrack


def reconstruct_alignment(backtrack, v, w, i, j):
    """Walk the backtrack table from the bottom right corner back to
    the top left corner, building the two aligned strings along the
    way.

    A dash character represents a gap introduced into one of the
    strings so that both aligned strings end up the same length.
    """
    aligned_v = []
    aligned_w = []

    while i > 0 or j > 0:
        direction = backtrack[i][j]
        if direction == "down":
            aligned_v.append(v[i - 1])
            aligned_w.append("-")
            i -= 1
        elif direction == "right":
            aligned_v.append("-")
            aligned_w.append(w[j - 1])
            j -= 1
        else:  # diagonal
            aligned_v.append(v[i - 1])
            aligned_w.append(w[j - 1])
            i -= 1
            j -= 1

    # characters were appended while walking backwards, so the
    # collected pieces need to be reversed to read in the right order
    aligned_v.reverse()
    aligned_w.reverse()

    return "".join(aligned_v), "".join(aligned_w)


def global_alignment(v, w, match_reward, mismatch_penalty, indel_penalty):
    """Return the best alignment score along with one alignment that
    achieves it.
    """
    score, backtrack = build_score_table(
        v, w, match_reward, mismatch_penalty, indel_penalty
    )
    best_score = score[len(v)][len(w)]
    aligned_v, aligned_w = reconstruct_alignment(backtrack, v, w, len(v), len(w))
    return best_score, aligned_v, aligned_w


if __name__ == "__main__":

    match_reward = 1
    mismatch_penalty = 1
    indel_penalty = 2

    v = "GAGA"
    w = "GAT"

    best_score, aligned_v, aligned_w = global_alignment(
        v, w, match_reward, mismatch_penalty, indel_penalty
    )

    print(best_score)
    print(aligned_v)
    print(aligned_w)