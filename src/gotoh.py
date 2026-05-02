import math

def gotoh(seq1, seq2, match=1, mismatch=-1, gap_open=-3, gap_extend=-1):
    n, m = len(seq1), len(seq2)
    M = [[-math.inf]*(m+1) for _ in range(n+1)]
    Ix = [[-math.inf]*(m+1) for _ in range(n+1)]
    Iy = [[-math.inf]*(m+1) for _ in range(n+1)]

    trace_M = [[None]*(m+1) for _ in range(n+1)]
    trace_Ix = [[None]*(m+1) for _ in range(n+1)]
    trace_Iy = [[None]*(m+1) for _ in range(n+1)]

    M[0][0] = 0

    for i in range(1, n+1):
        Ix[i][0] = gap_open + (i-1)*gap_extend
        trace_Ix[i][0] = ('Ix', i-1, 0) if i > 1 else ('M', i-1, 0)

    for j in range(1, m+1):
        Iy[0][j] = gap_open + (j-1)*gap_extend
        trace_Iy[0][j] = ('Iy', 0, j-1) if j > 1 else ('M', 0, j-1)

    for i in range(1, n+1):
        for j in range(1, m+1):
            s = match if seq1[i-1] == seq2[j-1] else mismatch

            diag_options = [M[i-1][j-1], Ix[i-1][j-1], Iy[i-1][j-1]]
            M[i][j] = max(diag_options) + s
            if M[i][j] == diag_options[0] + s:
                trace_M[i][j] = ('M', i-1, j-1)
            elif M[i][j] == diag_options[1] + s:
                trace_M[i][j] = ('Ix', i-1, j-1)
            else:
                trace_M[i][j] = ('Iy', i-1, j-1)

            ix_open = M[i-1][j] + gap_open
            ix_extend = Ix[i-1][j] + gap_extend
            if ix_open >= ix_extend:
                Ix[i][j] = ix_open
                trace_Ix[i][j] = ('M', i-1, j)
            else:
                Ix[i][j] = ix_extend
                trace_Ix[i][j] = ('Ix', i-1, j)

            iy_open = M[i][j-1] + gap_open
            iy_extend = Iy[i][j-1] + gap_extend
            if iy_open >= iy_extend:
                Iy[i][j] = iy_open
                trace_Iy[i][j] = ('M', i, j-1)
            else:
                Iy[i][j] = iy_extend
                trace_Iy[i][j] = ('Iy', i, j-1)

    final_score = max(M[n][m], Ix[n][m], Iy[n][m])
    if final_score == M[n][m]:
        matrix, i, j = 'M', n, m
    elif final_score == Ix[n][m]:
        matrix, i, j = 'Ix', n, m
    else:
        matrix, i, j = 'Iy', n, m

    align1 = []
    align2 = []

    while i > 0 or j > 0:
        if matrix == 'M':
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            matrix, i, j = trace_M[i][j]
        elif matrix == 'Ix':
            align1.append(seq1[i-1])
            align2.append('-')
            matrix, i, j = trace_Ix[i][j]
        else:
            align1.append('-')
            align2.append(seq2[j-1])
            matrix, i, j = trace_Iy[i][j]

    align1 = ''.join(reversed(align1))
    align2 = ''.join(reversed(align2))

    return align1, align2, final_score


def test_gotoh():
    s1, s2 = "GATTACA", "GATCA"
    a1, a2, score = gotoh(s1, s2)
    print("=== Gotoh Affine Gap Alignment ===")
    print(f"Seq1: {a1}")
    print(f"Seq2: {a2}")
    print(f"Score: {score}")


if __name__ == "__main__":
    test_gotoh()