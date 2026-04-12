def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """
    Performing global alignment while allowing for parameter sensitivity testing.
    """
    n, m = len(seq1), len(seq2)
    score_matrix = [[0 for _ in range(m + 1)] for _ in range(n + 1)]

    # Initializing boundaries using the provided gap penalty
    for i in range(n + 1): score_matrix[i][0] = i * gap
    for j in range(m + 1): score_matrix[0][j] = j * gap

    # Filling the matrix based on dynamic programming rules
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if seq1[i-1] == seq2[j-1] else mismatch
            score_matrix[i][j] = max(score_matrix[i-1][j-1] + s,
                                     score_matrix[i-1][j] + gap,
                                     score_matrix[i][j-1] + gap)

    # Executing traceback to reconstruct the alignment
    align1, align2 = "", ""
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and score_matrix[i][j] == score_matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
            align1 += seq1[i-1]; align2 += seq2[j-1]; i -= 1; j -= 1
        elif i > 0 and score_matrix[i][j] == score_matrix[i-1][j] + gap:
            align1 += seq1[i-1]; align2 += "-"; i -= 1
        else:
            align1 += "-"; align2 += seq2[j-1]; j -= 1

    return align1[::-1], align2[::-1], score_matrix[n][m]

if __name__ == "__main__":
    # Defining the sequences to test
    seq_a = "GATTACA"
    seq_b = "GATCA"

    # Testing sensitivity to different gap penalties
    gap_penalties = [-1, -2, -5, -10]

    print(" Parameter Sensitivity Experiment ")
    for p in gap_penalties:
        # Calling the function with the specific gap penalty 'p'
        a1, a2, score = needleman_wunsch(seq_a, seq_b, gap=p)
        print(f"Gap Penalty {p} -> Score: {score}, Alignment: {a1} / {a2}")