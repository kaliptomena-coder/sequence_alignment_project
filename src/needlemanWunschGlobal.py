def needleman_wunsch(seq1, seq2):
    # Setting scoring parameters
    match = 1
    mismatch = -1
    gap = -2

    # Determining matrix dimensions based on sequence lengths
    n = len(seq1)
    m = len(seq2)

    # Initializing the scoring matrix with zeros
    score_matrix = [[0 for _ in range(m + 1)] for _ in range(n + 1)]

    # Filling the first column with cumulative gap penalties
    for i in range(n + 1):
        score_matrix[i][0] = i * gap

    # Filling the first row with cumulative gap penalties
    for j in range(m + 1):
        score_matrix[0][j] = j * gap

    # Filling the matrix using dynamic programming
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # Calculating score for match or mismatch
            if seq1[i-1] == seq2[j-1]:
                score = match
            else:
                score = mismatch

            # Computing potential scores from diagonal, upward, and leftward moves
            diagonal = score_matrix[i-1][j-1] + score
            up = score_matrix[i-1][j] + gap
            left = score_matrix[i][j-1] + gap

            # Storing the highest score in the current cell
            score_matrix[i][j] = max(diagonal, up, left)

    # Starting traceback from the bottom-right corner to find the optimal alignment
    align1 = ""
    align2 = ""
    i = n
    j = m

    while i > 0 or j > 0:
        current_score = score_matrix[i][j]

        # Checking if the current path originated from a diagonal move
        if i > 0 and j > 0 and (current_score == score_matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)):
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            i -= 1
            j -= 1

        # Checking if the current path originated from an upward move (gap in seq2)
        elif i > 0 and current_score == score_matrix[i-1][j] + gap:
            align1 += seq1[i-1]
            align2 += "-"
            i -= 1

        # Identifying a leftward move (gap in seq1) as the remaining possibility
        else:
            align1 += "-"
            align2 += seq2[j-1]
            j -= 1

    # Printing the alignment results and final score
    print(f"\nAlignment Result:")
    print(align1[::-1])
    print(align2[::-1])
    print(f"Final Score: {score_matrix[n][m]}")

    # Reversing the strings and returning the results
    return align1[::-1], align2[::-1], score_matrix[n][m]

if __name__ == "__main__":
    # Testing the algorithm with example sequences
    result = needleman_wunsch("GATTACA", "GATCA")

    # Displaying the return values
    print("\nFunction Return Values (Align1, Align2, Score):")
    print(result)