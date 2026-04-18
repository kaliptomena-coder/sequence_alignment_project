def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """
    Performing global alignment while allowing for parameter sensitivity testing.
    """
    n, m = len(seq1), len(seq2)
    score_matrix = [[0]*(m+1) for _ in range(n + 1)]
    # Traceback matrix: 'D' (diag), 'U' (up), 'L' (left)
    trace_matrix = [[None]*(m+1) for _ in range(n+1)]

    # Initializing boundaries using the provided gap penalty
    for i in range(n + 1):
        score_matrix[i][0] = i * gap
        trace_matrix[i][0] = 'U'
    for j in range(m + 1):
        score_matrix[0][j] = j * gap
        trace_matrix[0][j] = 'L'

    # Filling the matrix based on dynamic programming rules
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = score_matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            up   = score_matrix[i-1][j] + gap
            left = score_matrix[i][j-1] + gap

            best = max(diag, up, left)
            score_matrix[i][j] = best

            # deterministic tie-breaking
            if best == diag:
                trace_matrix[i][j] = 'D'
            elif best == up:
                trace_matrix[i][j] = 'U'
            else:
                trace_matrix[i][j] = 'L'

    # Executing traceback to reconstruct the alignment
    align1, align2 = "", ""
    i, j = n, m
    while i > 0 or j > 0:
        t = trace_matrix[i][j]
        if t == 'D':
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            i -= 1;
            j -= 1
        elif t=='U':
            align1 += seq1[i-1]
            align2 += "-"
            i -= 1
        else:  # 'L'
            align1 += "-"
            align2 += seq2[j-1]
            j -= 1

    return align1[::-1], align2[::-1], score_matrix[n][m]

if __name__ == "__main__":
    # Defining the sequences to test
    seq_a = "GATTACAACTTG"
    seq_b = "GATCCAGTTCA"

    # Testing sensitivity to different gap penalties
    gap_penalties = [-1, -2, -5, -10]

    print(" Parameter Sensitivity Experiment ")
    for p in gap_penalties:
        # Calling the function with the specific gap penalty 'p'
        a1, a2, final_score = needleman_wunsch(seq_a, seq_b, gap=p)
        print(f"Gap Penalty {p} -> Score: {final_score}, Alignment: {a1} / {a2}")