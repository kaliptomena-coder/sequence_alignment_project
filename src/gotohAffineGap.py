def gotoh_affine_gap(seq1, seq2, match=2, mismatch=-1, gap_open=-5, gap_extend=-1):
    """
    Global alignment using Gotoh's algorithm with affine gap penalties.

    Affine gap model:
        gap cost = gap_open + k * gap_extend

    We use three matrices:
        M[i][j] = best score ending with a match/mismatch
        P[i][j] = best score ending with a gap in seq2 (vertical gap)
        Q[i][j] = best score ending with a gap in seq1 (horizontal gap)

    Traceback matrices store:
        (previous_matrix, previous_i, previous_j)
    """

    n, m = len(seq1), len(seq2)
    NEG_INF = float('-inf')

    # -------------------------
    # 1. Initialize matrices
    # -------------------------
    M = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    P = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    Q = [[NEG_INF] * (m + 1) for _ in range(n + 1)]

    # Traceback pointers
    trace_M = [[None] * (m + 1) for _ in range(n + 1)]
    trace_P = [[None] * (m + 1) for _ in range(n + 1)]
    trace_Q = [[None] * (m + 1) for _ in range(n + 1)]

    # -------------------------
    # 2. Boundary conditions
    # -------------------------
    M[0][0] = 0

    # First column (gaps in seq2)
    for i in range(1, n + 1):
        P[i][0] = gap_open + (i - 1) * gap_extend
        M[i][0] = P[i][0]
        trace_P[i][0] = ('P', i - 1, 0) if i > 1 else ('M', i - 1, 0)
        trace_M[i][0] = ('P', i, 0)

    # First row (gaps in seq1)
    for j in range(1, m + 1):
        Q[0][j] = gap_open + (j - 1) * gap_extend
        M[0][j] = Q[0][j]
        trace_Q[0][j] = ('Q', 0, j - 1) if j > 1 else ('M', 0, j - 1)
        trace_M[0][j] = ('Q', 0, j)

    # -------------------------
    # 3. Fill matrices
    # -------------------------
    for i in range(1, n + 1):
        for j in range(1, m + 1):

            # P (gap in seq2)
            open_gap = M[i - 1][j] + gap_open
            extend_gap = P[i - 1][j] + gap_extend

            if open_gap >= extend_gap:
                P[i][j] = open_gap
                trace_P[i][j] = ('M', i - 1, j)
            else:
                P[i][j] = extend_gap
                trace_P[i][j] = ('P', i - 1, j)

            # Q (gap in seq1)
            open_gap = M[i][j - 1] + gap_open
            extend_gap = Q[i][j - 1] + gap_extend

            if open_gap >= extend_gap:
                Q[i][j] = open_gap
                trace_Q[i][j] = ('M', i, j - 1)
            else:
                Q[i][j] = extend_gap
                trace_Q[i][j] = ('Q', i, j - 1)

            # M (match/mismatch)
            score = match if seq1[i - 1] == seq2[j - 1] else mismatch
            diag = M[i - 1][j - 1] + score

            best = max(diag, P[i][j], Q[i][j])
            M[i][j] = best

            if best == diag:
                trace_M[i][j] = ('M', i - 1, j - 1)
            elif best == P[i][j]:
                trace_M[i][j] = ('P', i, j)
            else:
                trace_M[i][j] = ('Q', i, j)

    # -------------------------
    # 4. Choose best end
    # -------------------------
    if M[n][m] >= P[n][m] and M[n][m] >= Q[n][m]:
        matrix = 'M'
    elif P[n][m] >= Q[n][m]:
        matrix = 'P'
    else:
        matrix = 'Q'

    # -------------------------
    # 5. Traceback
    # -------------------------
    align1, align2 = "", ""
    i, j = n, m

    while i > 0 or j > 0:

        if matrix == 'M':
            prev_matrix, ni, nj = trace_M[i][j]

            if ni == i - 1 and nj == j - 1:
                # diagonal
                align1 += seq1[i - 1]
                align2 += seq2[j - 1]
            else:
                # transition into gap state → DO NOT emit characters here
                matrix = prev_matrix
                continue

        elif matrix == 'P':
            prev_matrix, ni, nj = trace_P[i][j]
            align1 += seq1[i - 1]
            align2 += '-'

        elif matrix == 'Q':
            prev_matrix, ni, nj = trace_Q[i][j]
            align1 += '-'
            align2 += seq2[j - 1]

        i, j, matrix = ni, nj, prev_matrix

    return align1[::-1], align2[::-1], max(M[n][m], P[n][m], Q[n][m])

def test_gotoh():
    # Тест: вставка одного гэпа (ожидается штраф открытия)
    s1, s2 = "GATTACA", "GATCA"
    a1, a2, score = gotoh_affine_gap(s1, s2)
    print(f"Alignment: {a1}\n           {a2}\nScore: {score}")
    # Проверьте, совпадает ли скор с ручным расчетом

if __name__ == "__main__":
    test_gotoh()