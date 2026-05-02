def smith_waterman(seq1, seq2, match=2, mismatch=-1, gap=-2):
    n, m = len(seq1), len(seq2)
    score = [[0]*(m+1) for _ in range(n+1)]
    traceback = [[None]*(m+1) for _ in range(n+1)]
    max_score, max_pos = 0, (0,0)

    for i in range(1, n+1):
        for j in range(1, m+1):
            diag = score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            up = score[i-1][j] + gap
            left = score[i][j-1] + gap
            score[i][j] = max(0, diag, up, left)
            if score[i][j] == diag:
                traceback[i][j] = 'diag'
            elif score[i][j] == up:
                traceback[i][j] = 'up'
            elif score[i][j] == left:
                traceback[i][j] = 'left'
            if score[i][j] > max_score:
                max_score = score[i][j]
                max_pos = (i, j)

    align1, align2 = [], []
    i, j = max_pos
    while score[i][j] > 0:
        tb = traceback[i][j]
        if tb == 'diag':
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            i -= 1; j -= 1
        elif tb == 'up':
            align1.append(seq1[i-1])
            align2.append('-')
            i -= 1
        else:
            align1.append('-')
            align2.append(seq2[j-1])
            j -= 1

    return ''.join(reversed(align1)), ''.join(reversed(align2)), max_score


if __name__ == "__main__":
    a1, a2, score = smith_waterman("GATTACA", "GATCA")
    print("=== Smith-Waterman Local Alignment ===")
    print(f"Seq1: {a1}")
    print(f"Seq2: {a2}")
    print(f"Score: {score}")