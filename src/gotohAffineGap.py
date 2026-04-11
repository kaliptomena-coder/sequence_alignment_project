def gotoh_alignment(seq1, seq2, match=2, mismatch=-1, gap_open=-5, gap_extend=-1):
    n, m = len(seq1), len(seq2)

    # Инициализация матриц
    # M - основная матрица, P - для разрывов в первой строке, Q - во второй
    M = [[0] * (m + 1) for _ in range(n + 1)]
    P = [[float('-inf')] * (m + 1) for _ in range(n + 1)]
    Q = [[float('-inf')] * (m + 1) for _ in range(n + 1)]

    # Базовые условия
    M[0][0] = 0
    for i in range(1, n + 1):
        M[i][0] = gap_open + (i - 1) * gap_extend
        P[i][0] = M[i][0]
    for j in range(1, m + 1):
        M[0][j] = gap_open + (j - 1) * gap_extend
        Q[0][j] = M[0][j]

    # Заполнение матриц
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # Матрица P (разрыв в seq1)
            P[i][j] = max(M[i-1][j] + gap_open, P[i-1][j] + gap_extend)
            # Матрица Q (разрыв в seq2)
            Q[i][j] = max(M[i][j-1] + gap_open, Q[i][j-1] + gap_extend)
            # Основная матрица M
            score = match if seq1[i-1] == seq2[j-1] else mismatch
            M[i][j] = max(M[i-1][j-1] + score, P[i][j], Q[i][j])

    return M[n][m]


s1 = "GATTACAGG"
s2 = "GCATGCU"
print(f"Gotoh Score: {gotoh_alignment(s1, s2)}")