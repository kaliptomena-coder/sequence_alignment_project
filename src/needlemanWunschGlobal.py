def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):

    n, m = len(seq1), len(seq2)

    # Allocating the DP and traceback matrices
    score_matrix = [[0] * (m + 1) for _ in range(n + 1)]
    trace_matrix = [[None] * (m + 1) for _ in range(n + 1)]

    # Initializing the first column (gaps in seq2)
    for i in range(1, n + 1):
        score_matrix[i][0] = i * gap
        trace_matrix[i][0] = 'up'

    # Initializing the first row (gaps in seq1)
    for j in range(1, m + 1):
        score_matrix[0][j] = j * gap
        trace_matrix[0][j] = 'left'

    # Filling the rest of the matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = score_matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            up   = score_matrix[i-1][j] + gap
            left = score_matrix[i][j-1] + gap

            score_matrix[i][j], trace_matrix[i][j] = max(
                (diag, 'diag'), (up, 'up'), (left, 'left')
            )

    # Traceback from bottom-right to top-left
    align1, align2 = [], []
    i, j = n, m

    while i > 0 or j > 0:
        direction = trace_matrix[i][j]
        if direction == 'diag':
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            i -= 1; j -= 1
        elif direction == 'up':
            align1.append(seq1[i-1])
            align2.append('-')
            i -= 1
        else:
            align1.append('-')
            align2.append(seq2[j-1])
            j -= 1

    return ''.join(reversed(align1)), ''.join(reversed(align2)), score_matrix[n][m]

if __name__ == "__main__":
    # Parameter sensitivity experiment: same sequences, different gap penalties
    seq_a = "GATTACAACTTG"
    seq_b = "GATCCAGTTCAAA"

    print("=== Parameter Sensitivity Experiment ===")
    for gap_val in [-1, -2, -5, -10]:
        a1, a2, sc = needleman_wunsch(seq_a, seq_b, gap=gap_val)
        print(f"Gap={gap_val:3d}  Score={sc:4d}  {a1}  {a2}")
