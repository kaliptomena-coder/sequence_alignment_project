# Classic dynamic programming algorithm for optimal global sequence alignment.
# Aligns two sequences end-to-end using a scoring scheme for matches,
# mismatches, and linear gap penalties.
#
# Algorithm overview:
#   1. Initialize a (n+1) x (m+1) DP matrix with gap penalties along edges
#   2. Fill each cell with the best of: diagonal (match/mismatch), up (gap), left (gap)
#   3. Traceback from bottom-right corner to reconstruct the alignment
#
# Time complexity:  O(n * m)
# Space complexity: O(n * m)
#
# Tie-breaking: diagonal > up > left (deterministic)

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """
    Running global alignment on seq1 and seq2 with configurable scoring.

    Parameters
    ----------
    seq1, seq2 : str   - the two sequences to align
    match      : int   - reward for matching characters  (default +1)
    mismatch   : int   - penalty for mismatching chars   (default -1)
    gap        : int   - penalty per gap character       (default -2)

    Returns
    -------
    align1, align2 : str  - aligned strings with '-' representing gaps
    score          : int  - total alignment score
    """
    n, m = len(seq1), len(seq2)

    # Allocating the DP and traceback matrices
    score_matrix = [[0] * (m + 1) for _ in range(n + 1)]
    trace_matrix = [[None] * (m + 1) for _ in range(n + 1)]

    # Initializing the first column (gaps in seq2)
    for i in range(1, n + 1):
        score_matrix[i][0] = i * gap
        trace_matrix[i][0] = 'U'  # Up

    # Initializing the first row (gaps in seq1)
    for j in range(1, m + 1):
        score_matrix[0][j] = j * gap
        trace_matrix[0][j] = 'L'  # Left

    # Filling the rest of the matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = score_matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            up   = score_matrix[i-1][j] + gap
            left = score_matrix[i][j-1] + gap

            best = max(diag, up, left)
            score_matrix[i][j] = best

            # Recording where we came from (deterministic tie-breaking)
            if best == diag:
                trace_matrix[i][j] = 'D'
            elif best == up:
                trace_matrix[i][j] = 'U'
            else:
                trace_matrix[i][j] = 'L'

    # Traceback from bottom-right to top-left
    align1, align2 = "", ""
    i, j = n, m

    while i > 0 or j > 0:
        direction = trace_matrix[i][j]
        if direction == 'D':
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            i -= 1; j -= 1
        elif direction == 'U':
            align1 += seq1[i-1]
            align2 += '-'
            i -= 1
        else:  # 'L'
            align1 += '-'
            align2 += seq2[j-1]
            j -= 1

    # Reversing because we built the strings backwards during traceback
    return align1[::-1], align2[::-1], score_matrix[n][m]

if __name__ == "__main__":
    # Parameter sensitivity experiment: same sequences, different gap penalties
    seq_a = "GATTACAACTTG"
    seq_b = "GATCCAGTTCAAA"

    print("=== Parameter Sensitivity Experiment ===")
    for gap_val in [-1, -2, -5, -10]:
        a1, a2, sc = needleman_wunsch(seq_a, seq_b, gap=gap_val)
        print(f"Gap={gap_val:3d}  Score={sc:4d}  {a1}  {a2}")
