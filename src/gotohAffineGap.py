def gotoh_alignment(seq1, seq2, match=2, mismatch=-1, gap_open=-5, gap_extend=-1):
    """
    Performing global alignment using Gotoh's algorithm with affine gap penalties.
    """
    n, m = len(seq1), len(seq2)
    NEG = float('-inf')

    # Initializing matrices for matching (M), gapping seq1 (P), and gapping seq2 (Q)
    M = [[0] * (m + 1) for _ in range(n + 1)]
    P = [[NEG] * (m + 1) for _ in range(n + 1)]
    Q = [[NEG] * (m + 1) for _ in range(n + 1)]
    trace = [[''] * (m + 1) for _ in range(n + 1)]

    # Setting initial gap penalties for the vertical boundary
    for i in range(1, n + 1):
        M[i][0] = gap_open + (i - 1) * gap_extend
        P[i][0] = M[i][0]

    # Setting initial gap penalties for the horizontal boundary
    for j in range(1, m + 1):
        M[0][j] = gap_open + (j - 1) * gap_extend
        Q[0][j] = M[0][j]

    # Filling the scoring matrices using dynamic programming
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # Calculating the cost for opening or extending a gap in seq2
            P[i][j] = max(M[i - 1][j] + gap_open, P[i - 1][j] + gap_extend)

            # Calculating the cost for opening or extending a gap in seq1
            Q[i][j] = max(M[i][j - 1] + gap_open, Q[i][j - 1] + gap_extend)

            # Determining the score for matching or mismatching characters
            s = match if seq1[i - 1] == seq2[j - 1] else mismatch
            diag = M[i - 1][j - 1] + s

            # Identifying the best score among the three possible moves
            best = max(diag, P[i][j], Q[i][j])
            M[i][j] = best

            # Recording the chosen path for performing the traceback
            if best == diag:
                trace[i][j] = 'D'  # Diagonal move
            elif best == P[i][j]:
                trace[i][j] = 'P'  # Upward move
            else:
                trace[i][j] = 'Q'  # Leftward move

    # Starting the traceback from the bottom-right cell
    align1, align2 = '', ''
    i, j = n, m

    # Iterating until both sequences are fully processed
    while i > 0 or j > 0:
        # Handling the top boundary by moving left
        if i == 0:
            align1 += '-'
            align2 += seq2[j - 1]
            j -= 1
        # Handling the left boundary by moving up
        elif j == 0:
            align1 += seq1[i - 1]
            align2 += '-'
            i -= 1
        else:
            t = trace[i][j]
            # Processing a diagonal match or mismatch
            if t == 'D':
                align1 += seq1[i - 1]
                align2 += seq2[j - 1]
                i -= 1
                j -= 1
            # Processing a gap in sequence 2
            elif t == 'P':
                align1 += seq1[i - 1]
                align2 += '-'
                i -= 1
            # Processing a gap in sequence 1
            else: # t == 'Q'
                align1 += '-'
                align2 += seq2[j - 1]
                j -= 1

    # Reversing the strings and returning the final alignment results
    return align1[::-1], align2[::-1], M[n][m]

if __name__ == "__main__":
    # Testing the Gotoh implementation with sample sequences
    s1, s2 = "GATTACA", "GATCA"
    res1, res2, score = gotoh_alignment(s1, s2)
    print(f"Final Alignment 1: {res1}")
    print(f"Final Alignment 2: {res2}")
    print(f"Final Score: {score}")