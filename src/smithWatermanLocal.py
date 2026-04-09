def smith_waterman(seq1, seq2):
    # Scoring system
    match = 1
    mismatch = -1
    gap = -2

    n, m = len(seq1), len(seq2)

    # Initialize matrix with zeros (no gap penalties for first row/column in local)
    score_matrix = [[0 for _ in range(m + 1)] for _ in range(n + 1)]

    # Fill the matrix
    max_score = 0
    max_pos = (0, 0)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score = match if seq1[i-1] == seq2[j-1] else mismatch

            diagonal = score_matrix[i-1][j-1] + score
            up = score_matrix[i-1][j] + gap
            left = score_matrix[i][j-1] + gap

            # Local alignment rule: value cannot be negative
            score_matrix[i][j] = max(0, diagonal, up, left)

            # Track where the highest score is for traceback
            if score_matrix[i][j] >= max_score:
                max_score = score_matrix[i][j]
                max_pos = (i, j)

    # Traceback
    align1, align2 = "", ""
    i, j = max_pos

    # Stop when we hit a cell with score 0
    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        current = score_matrix[i][j]
        diag_score = match if seq1[i-1] == seq2[j-1] else mismatch

        if current == score_matrix[i-1][j-1] + diag_score:
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            i -= 1
            j -= 1
        elif i > 0 and current == score_matrix[i-1][j] + gap:
            align1 += seq1[i-1]
            align2 += "-"
            i -= 1
        else:
            align1 += "-"
            align2 += seq2[j-1]
            j -= 1

    print(f"\nLocal Alignment Result:")
    print(align1[::-1])
    print(align2[::-1])
    print(f"Max Alignment Score: {max_score}")

    return score_matrix

if __name__ == "__main__":
    result = smith_waterman("GATTACA", "GATCA")

    print("\nLocal Score Matrix:")
    for row in result:
        print(row)