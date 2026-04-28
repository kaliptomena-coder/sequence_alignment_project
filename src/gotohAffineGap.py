# Extension of Needleman-Wunsch that uses a more biologically realistic
# gap cost model: gap_cost = gap_open + k * gap_extend
#
# Opening a gap is expensive (reflects the rare event of an insertion/deletion),
# but extending an existing gap is cheap (consecutive indels are common).
#
# Three matrices are maintained simultaneously:
#   M[i][j] — best score ending with a match or mismatch at (i, j)
#   P[i][j] — best score ending with a gap in seq2 (vertical gap / deletion)
#   Q[i][j] — best score ending with a gap in seq1 (horizontal gap / insertion)
#
# Traceback pointers store (previous_matrix, previous_i, previous_j).
#
# Time complexity:  O(n * m)
# Space complexity: O(n * m)

def gotoh_affine_gap(seq1, seq2, match=2, mismatch=-1, gap_open=-5, gap_extend=-1):
    """
    Running global alignment with affine gap penalties using the Gotoh algorithm.

    Parameters
    ----------
    seq1, seq2  : str   - the two sequences to align
    match       : int   - reward for matching characters    (default +2)
    mismatch    : int   - penalty for mismatching chars     (default -1)
    gap_open    : int   - penalty for opening a gap         (default -5)
    gap_extend  : int   - penalty for extending a gap by 1  (default -1)

    Returns
    -------
    align1, align2 : str   - aligned strings with '-' for gaps
    score          : float - best alignment score
    """
    n, m = len(seq1), len(seq2)
    NEG_INF = float('-inf')

    # --- Initializing the three DP matrices ---
    M = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    P = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    Q = [[NEG_INF] * (m + 1) for _ in range(n + 1)]

    trace_M = [[None] * (m + 1) for _ in range(n + 1)]
    trace_P = [[None] * (m + 1) for _ in range(n + 1)]
    trace_Q = [[None] * (m + 1) for _ in range(n + 1)]

    # --- Boundary conditions ---
    M[0][0] = 0

    # First column: consuming seq1 with deletions (gaps in seq2)
    for i in range(1, n + 1):
        P[i][0]       = gap_open + (i - 1) * gap_extend
        M[i][0]       = P[i][0]
        trace_P[i][0] = ('P', i - 1, 0) if i > 1 else ('M', i - 1, 0)
        trace_M[i][0] = ('P', i, 0)

    # First row: consuming seq2 with insertions (gaps in seq1)
    for j in range(1, m + 1):
        Q[0][j]       = gap_open + (j - 1) * gap_extend
        M[0][j]       = Q[0][j]
        trace_Q[0][j] = ('Q', 0, j - 1) if j > 1 else ('M', 0, j - 1)
        trace_M[0][j] = ('Q', 0, j)

    # --- Filling the three matrices ---
    for i in range(1, n + 1):
        for j in range(1, m + 1):

            # P: gap in seq2 (vertical move)
            open_gap   = M[i-1][j] + gap_open
            extend_gap = P[i-1][j] + gap_extend
            if open_gap >= extend_gap:
                P[i][j]       = open_gap
                trace_P[i][j] = ('M', i - 1, j)
            else:
                P[i][j]       = extend_gap
                trace_P[i][j] = ('P', i - 1, j)

            # Q: gap in seq1 (horizontal move)
            open_gap   = M[i][j-1] + gap_open
            extend_gap = Q[i][j-1] + gap_extend
            if open_gap >= extend_gap:
                Q[i][j]       = open_gap
                trace_Q[i][j] = ('M', i, j - 1)
            else:
                Q[i][j]       = extend_gap
                trace_Q[i][j] = ('Q', i, j - 1)

            # M: match or mismatch (diagonal move)
            s    = match if seq1[i-1] == seq2[j-1] else mismatch
            diag = M[i-1][j-1] + s
            best = max(diag, P[i][j], Q[i][j])
            M[i][j] = best

            if best == diag:
                trace_M[i][j] = ('M', i - 1, j - 1)
            elif best == P[i][j]:
                trace_M[i][j] = ('P', i, j)
            else:
                trace_M[i][j] = ('Q', i, j)

    # --- Choosing the best terminal matrix ---
    final_score = max(M[n][m], P[n][m], Q[n][m])
    if M[n][m] >= P[n][m] and M[n][m] >= Q[n][m]:
        matrix = 'M'
    elif P[n][m] >= Q[n][m]:
        matrix = 'P'
    else:
        matrix = 'Q'

    # --- Traceback ---
    align1, align2 = "", ""
    i, j = n, m

    while i > 0 or j > 0:
        if matrix == 'M':
            prev_matrix, ni, nj = trace_M[i][j]
            if ni == i - 1 and nj == j - 1:
                align1 += seq1[i-1]
                align2 += seq2[j-1]
            else:
                # Transitioning into a gap state — don't emit characters yet
                matrix = prev_matrix
                continue
        elif matrix == 'P':
            prev_matrix, ni, nj = trace_P[i][j]
            align1 += seq1[i-1]
            align2 += '-'
        else:  # 'Q'
            prev_matrix, ni, nj = trace_Q[i][j]
            align1 += '-'
            align2 += seq2[j-1]

        i, j, matrix = ni, nj, prev_matrix

    return align1[::-1], align2[::-1], final_score


def test_gotoh():
    s1, s2 = "GATTACA", "GATCA"
    a1, a2, score = gotoh_affine_gap(s1, s2)
    print("=== Gotoh Affine Gap Alignment ===")
    print(f"Seq1: {a1}\nSeq2: {a2}\nScore: {score}")


if __name__ == "__main__":
    test_gotoh()