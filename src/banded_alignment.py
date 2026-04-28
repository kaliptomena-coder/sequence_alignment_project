# For sequences that are very similar (e.g., >80% identical), the optimal
# alignment will stay close to the main diagonal of the DP matrix.
# Banded alignment exploits this by only filling cells where |i - j| <= k,
# reducing complexity from O(n*m) to O(k*n).
#
# When to use:
#   - Sequences that differ by only a few insertions/deletions (small k)
#   - Aligning millions of short reads where speed is critical
#
# When NOT to use:
#   - Distantly related sequences that may have large indels
#   - When a guaranteed optimal alignment is required at all costs
#
# Cells outside the band are set to NEG_INF and cannot influence the result.
# If the band is too narrow to reach the bottom-right corner, the algorithm
# automatically falls back to full Needleman-Wunsch.

NEG_INF = float('-inf')


def banded_nw(seq1, seq2, k=10, match=2, mismatch=-1, gap=-2):
    """
    Running Needleman-Wunsch restricted to a diagonal band of width 2k+1.

    Parameters
    ----------
    seq1, seq2 : str  - the two sequences to align
    k          : int  - half-bandwidth; only cells where |i-j| <= k are filled
    match      : int  - match reward     (default +2)
    mismatch   : int  - mismatch penalty (default -1)
    gap        : int  - gap penalty      (default -2)

    Returns
    -------
    align1, align2 : str  - aligned strings (with '-' for gaps)
    score          : int  - alignment score (or fallback NW score)
    """
    n = len(seq1)
    m = len(seq2)

    # Checking if the band is wide enough to even reach the corner
    if abs(n - m) > k:
        print(f"WARNING: |len(seq1)-len(seq2)| = {abs(n-m)} > k={k}. Band too narrow.")
        print("Falling back to full Needleman-Wunsch...")
        return _full_nw_fallback(seq1, seq2, match, mismatch, gap)

    # Allocating DP and traceback matrices (full size but mostly unused)
    dp    = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    trace = [[None]    * (m + 1) for _ in range(n + 1)]

    # Base case: (0,0) is the starting point
    dp[0][0] = 0

    # Initializing the boundary cells that fall within the band
    for j in range(1, min(k + 1, m + 1)):
        dp[0][j]    = j * gap
        trace[0][j] = 'L'

    for i in range(1, min(k + 1, n + 1)):
        dp[i][0]    = i * gap
        trace[i][0] = 'U'

    # Filling only cells within the band
    for i in range(1, n + 1):
        j_start = max(1, i - k)
        j_end   = min(m, i + k)

        for j in range(j_start, j_end + 1):

            # Diagonal: match or mismatch
            if dp[i-1][j-1] != NEG_INF:
                s          = match if seq1[i-1] == seq2[j-1] else mismatch
                diag_score = dp[i-1][j-1] + s
            else:
                diag_score = NEG_INF

            # Up: gap in seq2
            up_score   = dp[i-1][j] + gap if dp[i-1][j]   != NEG_INF else NEG_INF

            # Left: gap in seq1
            left_score = dp[i][j-1] + gap if dp[i][j-1]   != NEG_INF else NEG_INF

            best      = max(diag_score, up_score, left_score)
            dp[i][j]  = best

            if best == NEG_INF:
                trace[i][j] = None
            elif best == diag_score:
                trace[i][j] = 'D'
            elif best == up_score:
                trace[i][j] = 'U'
            else:
                trace[i][j] = 'L'

    # Checking whether we reached the corner
    if dp[n][m] == NEG_INF:
        print("ERROR: Band too narrow, could not reach dp[n][m]. Falling back...")
        return _full_nw_fallback(seq1, seq2, match, mismatch, gap)

    # Traceback from the bottom-right corner
    align1, align2 = "", ""
    i, j = n, m

    while i > 0 or j > 0:
        direction = trace[i][j]
        if direction == 'D':
            align1 += seq1[i-1]; align2 += seq2[j-1]; i -= 1; j -= 1
        elif direction == 'U':
            align1 += seq1[i-1]; align2 += '-'; i -= 1
        elif direction == 'L':
            align1 += '-'; align2 += seq2[j-1]; j -= 1
        else:
            print(f"ERROR: No traceback direction at ({i},{j}). Stopping.")
            break

    return align1[::-1], align2[::-1], dp[n][m]


def _full_nw_fallback(seq1, seq2, match, mismatch, gap):
    """Standard O(n*m) Needleman-Wunsch used when the band is too narrow."""
    n, m = len(seq1), len(seq2)
    dp    = [[0]    * (m + 1) for _ in range(n + 1)]
    trace = [[None] * (m + 1) for _ in range(n + 1)]

    for i in range(n + 1): dp[i][0] = i * gap; trace[i][0] = 'U'
    for j in range(m + 1): dp[0][j] = j * gap; trace[0][j] = 'L'
    dp[0][0] = 0; trace[0][0] = None

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s    = match if seq1[i-1] == seq2[j-1] else mismatch
            best = max(dp[i-1][j-1] + s, dp[i-1][j] + gap, dp[i][j-1] + gap)
            dp[i][j] = best
            if   best == dp[i-1][j-1] + s:  trace[i][j] = 'D'
            elif best == dp[i-1][j]   + gap: trace[i][j] = 'U'
            else:                            trace[i][j] = 'L'

    a1, a2, i, j = "", "", n, m
    while i > 0 or j > 0:
        d = trace[i][j]
        if   d == 'D': a1 += seq1[i-1]; a2 += seq2[j-1]; i -= 1; j -= 1
        elif d == 'U': a1 += seq1[i-1]; a2 += '-';        i -= 1
        else:          a1 += '-';        a2 += seq2[j-1];          j -= 1

    return a1[::-1], a2[::-1], dp[n][m]


if __name__ == "__main__":
    print("=== Banded Needleman-Wunsch Demo ===")

    # Test 1: Identical sequences — a small band should work fine
    s1, s2 = "GATTACAACTTG", "GATTACAACTTG"
    a1, a2, score = banded_nw(s1, s2, k=3)
    print(f"\nTest 1 (identical, k=3): Score={score}\n  {a1}\n  {a2}")

    # Test 2: One mismatch
    s1, s2 = "GATTACAACTTG", "GATCACAACTTG"
    a1, a2, score = banded_nw(s1, s2, k=3)
    print(f"\nTest 2 (one mismatch, k=3): Score={score}\n  {a1}\n  {a2}")

    # Test 3: One small gap
    s1, s2 = "GATTACA", "GATACA"
    a1, a2, score = banded_nw(s1, s2, k=2)
    print(f"\nTest 3 (one gap, k=2): Score={score}\n  {a1}\n  {a2}")

    # Test 4: Band too narrow — expect automatic fallback
    s1, s2 = "AAAAAAAAAAA", "TTTTTT"
    a1, a2, score = banded_nw(s1, s2, k=2)
    print(f"\nTest 4 (fallback, k=2): Score={score}\n  {a1}\n  {a2}")