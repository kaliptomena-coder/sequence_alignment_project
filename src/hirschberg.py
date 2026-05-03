def last_line_nw(seq1, seq2, match=1, mismatch=-1, gap=-2):

    n, m = len(seq1), len(seq2)

    # We only keep two rows at a time: the previous and the current one
    prev = [i * gap for i in range(m + 1)]
    curr = [0] * (m + 1)

    for i in range(1, n + 1):
        curr[0] = i * gap
        for j in range(1, m + 1):
            s       = match if seq1[i-1] == seq2[j-1] else mismatch
            curr[j] = max(prev[j-1] + s, prev[j] + gap, curr[j-1] + gap)
        prev = curr[:]

    return prev


def nw_small(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """
       Full Needleman-Wunsch for small base cases (when n==1 or m==1).
       Returns sequences and score.
       """

    n, m = len(seq1), len(seq2)
    dp = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(n + 1): dp[i][0] = i * gap
    for j in range(m + 1): dp[0][j] = j * gap

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s       = match if seq1[i-1] == seq2[j-1] else mismatch
            dp[i][j] = max(dp[i-1][j-1] + s, dp[i-1][j] + gap, dp[i][j-1] + gap)

    # Traceback
    a1, a2, i, j = "", "", n, m
    while i > 0 or j > 0:
        s = match if (i > 0 and j > 0 and seq1[i-1] == seq2[j-1]) else mismatch
        if i > 0 and j > 0 and dp[i][j] == dp[i-1][j-1] + s:
            a1 += seq1[i-1]; a2 += seq2[j-1]; i -= 1; j -= 1
        elif i > 0 and dp[i][j] == dp[i-1][j] + gap:
            a1 += seq1[i-1]; a2 += '-'; i -= 1
        else:
            a1 += '-'; a2 += seq2[j-1]; j -= 1

    return a1[::-1], a2[::-1], dp[n][m]


def hirschberg(seq1, seq2, match=1, mismatch=-1, gap=-2):
    n, m = len(seq1), len(seq2)

    # Base cases: one sequence is empty → fill entirely with gaps
    if n == 0:
        return '-' * m, seq2, gap*m
    if m == 0:
        return seq1, '-' * n, gap*n

    # Base case: one sequence is a single character → use full NW
    if n == 1 or m == 1:
        a1, a2, score= nw_small(seq1, seq2, match, mismatch, gap)
        return a1, a2, score

    # Divide seq1 at its midpoint
    mid1 = n // 2

    # Forward pass on the top half of seq1
    score_left  = last_line_nw(seq1[:mid1],       seq2,       match, mismatch, gap)
    # Backward pass on the bottom half (reversed sequences)
    score_right = last_line_nw(seq1[mid1:][::-1], seq2[::-1], match, mismatch, gap)

    # Finding the split point in seq2 that maximizes the combined score
    partition = [score_left[j] + score_right[m - j] for j in range(m + 1)]
    mid2      = partition.index(max(partition))

    # Recursing on the left and right sub-problems
    left1,  left2, left_score  = hirschberg(seq1[:mid1],  seq2[:mid2],  match, mismatch, gap)
    right1, right2, right_score = hirschberg(seq1[mid1:],  seq2[mid2:],  match, mismatch, gap)

    return left1 + right1, left2 + right2, left_score+right_score


if __name__ == "__main__":
    s1 = "AGTAACG"
    s2 = "ACATAG"
    r1, r2, score = hirschberg(s1, s2)
    print("=== Hirschberg Alignment ===")
    print(f"Seq1: {r1}")
    print(f"Seq2: {r2}")
    print(f"Score: {score}")