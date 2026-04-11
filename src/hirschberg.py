def last_line_nw(seq1, seq2, match=2, mismatch=-1, gap=-1):
    """Calculating only the last row of the Needleman-Wunsch matrix to save memory."""
    n, m = len(seq1), len(seq2)
    # Initializing the first row with cumulative gap penalties
    prev = [i * gap for i in range(m + 1)]
    curr = [0] * (m + 1)

    for i in range(1, n + 1):
        # Setting the first cell of the current row as a vertical gap
        curr[0] = i * gap
        for j in range(1, m + 1):
            # Determining the score for the current character pair
            score = match if seq1[i-1] == seq2[j-1] else mismatch
            # Selecting the maximum score from diagonal, vertical, and horizontal moves
            curr[j] = max(prev[j-1] + score,
                          prev[j] + gap,
                          curr[j-1] + gap)
        # Copying the current row to the previous row for the next iteration
        prev = curr[:]
    return prev

def nw_small(seq1, seq2, match, mismatch, gap):
    """Handling base cases for minimal string lengths during recursion."""
    # Returning all gaps if the first sequence is empty
    if not seq1: return "-" * len(seq2), seq2
    # Returning all gaps if the second sequence is empty
    if not seq2: return seq1, "-" * len(seq1)

    # Returning exact matches if sequences are identical
    if seq1 == seq2:
        return seq1, seq2
    else:
        # Padding the shorter string with gaps to match lengths
        return seq1 + "-" * (max(0, len(seq2)-len(seq1))), seq2 + "-" * (max(0, len(seq1)-len(seq2)))

def hirschberg(seq1, seq2, match=2, mismatch=-1, gap=-1):
    """Executing the divide and conquer alignment to achieve linear space complexity."""
    n, m = len(seq1), len(seq2)

    # Returning gaps for empty sequence scenarios
    if n == 0:
        return "-" * m, seq2
    elif m == 0:
        return seq1, "-" * n

    # Resolving the alignment using the helper function for single-character sequences
    if n == 1 or m == 1:
        return nw_small(seq1, seq2, match, mismatch, gap)

    # Splitting the first sequence into two equal halves
    mid1 = n // 2

    # Processing the forward pass on the first half of the sequences
    score_left = last_line_nw(seq1[:mid1], seq2, match, mismatch, gap)
    # Processing the backward pass on the reversed second half of the sequences
    score_right = last_line_nw(seq1[mid1:][::-1], seq2[::-1], match, mismatch, gap)

    # Summing scores to identify the optimal split point in the second sequence
    partition = [score_left[j] + score_right[m-j] for j in range(m + 1)]
    mid2 = partition.index(max(partition))

    # Triggering recursive calls for the resulting left and right sub-problems
    aln1_left, aln2_left = hirschberg(seq1[:mid1], seq2[:mid2], match, mismatch, gap)
    aln1_right, aln2_right = hirschberg(seq1[mid1:], seq2[mid2:], match, mismatch, gap)

    # Concatenating the sub-alignments into the final resulting sequences
    return aln1_left + aln1_right, aln2_left + aln2_right

if __name__ == "__main__":
    # Testing the algorithm with DNA sample strings
    s1 = "AGTAACG"
    s2 = "ACATAG"
    res1, res2 = hirschberg(s1, s2)
    print("--- Hirschberg Result ---")
    print(f"Seq1: {res1}")
    print(f"Seq2: {res2}")