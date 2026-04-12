def smith_waterman(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """
    Finding the best local alignment by zeroing out negative scores.
    """
    n, m = len(seq1), len(seq2)
    score_matrix = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    max_score, max_pos = 0, (0, 0)

    # Filling the matrix and tracking the peak score position
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if seq1[i-1] == seq2[j-1] else mismatch
            score_matrix[i][j] = max(0, score_matrix[i-1][j-1] + s,
                                     score_matrix[i-1][j] + gap,
                                     score_matrix[i][j-1] + gap)
            if score_matrix[i][j] >= max_score:
                max_score = score_matrix[i][j]
                max_pos = (i, j)

    # Tracing back from the peak until hitting a zero-score cell
    align1, align2 = "", ""
    i, j = max_pos
    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        s = match if seq1[i-1] == seq2[j-1] else mismatch
        if score_matrix[i][j] == score_matrix[i-1][j-1] + s:
            align1 += seq1[i-1]; align2 += seq2[j-1]; i -= 1; j -= 1
        elif i > 0 and score_matrix[i][j] == score_matrix[i-1][j] + gap:
            align1 += seq1[i-1]; align2 += "-"; i -= 1
        else:
            align1 += "-"; align2 += seq2[j-1]; j -= 1

    return align1[::-1], align2[::-1], max_score

if __name__ == "__main__":
    # Testing the local alignment
    aln1, aln2, score = smith_waterman("GATTACA", "GATCA")

    print(f"\nFinal Aligned Sequences:")
    print(f"Seq 1: {aln1}")
    print(f"Seq 2: {aln2}")
    print(f"Score: {score}")