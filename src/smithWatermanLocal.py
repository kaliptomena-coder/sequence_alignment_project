# Finds the highest-scoring local match between two sequences.
# Unlike Needleman-Wunsch, this method can ignore poorly-matching ends,
# making it ideal for finding conserved motifs or domains.
#
# Key difference from NW:
#   - Negative scores are reset to 0 (prevents extending bad alignments)
#   - Traceback starts from the maximum score cell, not the bottom-right
#   - Traceback stops when a zero-score cell is reached
#
# Time complexity:  O(n * m)
# Space complexity: O(n * m)

def smith_waterman(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """
    Finding the best local alignment between seq1 and seq2.

    Parameters
    ----------
    seq1, seq2 : str  - the two sequences to align
    match      : int  - reward for matching characters  (default +1)
    mismatch   : int  - penalty for mismatching chars   (default -1)
    gap        : int  - gap penalty per character       (default -2)

    Returns
    -------
    align1, align2 : str  - locally aligned strings with '-' for gaps
    max_score      : int  - score of the best local alignment found
    """
    n, m = len(seq1), len(seq2)
    score_matrix = [[0] * (m + 1) for _ in range(n + 1)]
    max_score, max_pos = 0, (0, 0)

    # Filling the scoring matrix; zero is the floor (no negative values)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if seq1[i-1] == seq2[j-1] else mismatch
            score_matrix[i][j] = max(
                0,
                score_matrix[i-1][j-1] + s,
                score_matrix[i-1][j]   + gap,
                score_matrix[i][j-1]   + gap
            )
            # Tracking where the overall maximum score is
            if score_matrix[i][j] >= max_score:
                max_score = score_matrix[i][j]
                max_pos   = (i, j)

    # Tracing back from the maximum score position until we hit a zero cell
    align1, align2 = "", ""
    i, j = max_pos

    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        s = match if seq1[i-1] == seq2[j-1] else mismatch
        if score_matrix[i][j] == score_matrix[i-1][j-1] + s:
            align1 += seq1[i-1]; align2 += seq2[j-1]; i -= 1; j -= 1
        elif score_matrix[i][j] == score_matrix[i-1][j] + gap:
            align1 += seq1[i-1]; align2 += '-'; i -= 1
        else:
            align1 += '-'; align2 += seq2[j-1]; j -= 1

    return align1[::-1], align2[::-1], max_score


if __name__ == "__main__":
    a1, a2, score = smith_waterman("GATTACA", "GATCA")
    print("=== Smith-Waterman Local Alignment ===")
    print(f"Seq1: {a1}")
    print(f"Seq2: {a2}")
    print(f"Score: {score}")