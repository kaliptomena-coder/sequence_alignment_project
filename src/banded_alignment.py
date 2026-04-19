# =============================================================================
# FILE: banded_alignment.py
# PURPOSE: Banded Needleman-Wunsch alignment (space and time efficient)
#
# WHAT IS BANDED ALIGNMENT?
#   Normal NW fills every cell in an (n+1) x (m+1) matrix → O(n*m) time.
#   Banded alignment says: "I only care about alignments that stay close
#   to the main diagonal."  We only fill cells where |i - j| <= k.
#   This makes it O(k * n) time, which is MUCH faster when sequences
#   are similar and k can be small.
#
# WHEN TO USE IT:
#   - Sequences that are very similar (e.g., >80% identical)
#   - You know the sequences won't have huge insertions/deletions
#   - Speed is important (e.g., aligning millions of short reads)
#
# WHEN NOT TO USE IT:
#   - Distantly related sequences (need a larger band or full NW)
#   - You need a guaranteed optimal alignment always
# =============================================================================

# We use a sentinel value for "unreachable" cells (outside the band)
NEG_INF = float('-inf')


def banded_nw(seq1, seq2, k=10, match=2, mismatch=-1, gap=-2):
    """
    Banded Needleman-Wunsch global alignment.

    Parameters
    ----------
    seq1   : str   – first sequence  (length n)
    seq2   : str   – second sequence (length m)
    k      : int   – half-bandwidth; only cells where |i-j| <= k are filled
    match  : int   – score for matching characters  (positive, e.g. +2)
    mismatch: int  – score for mismatching chars    (negative, e.g. -1)
    gap    : int   – penalty per gap character      (negative, e.g. -2)

    Returns
    -------
    align1 : str   – aligned version of seq1 (may contain '-')
    align2 : str   – aligned version of seq2 (may contain '-')
    score  : int   – alignment score (or None if band was too narrow)
    """

    n = len(seq1)   # number of rows  (seq1 along rows)
    m = len(seq2)   # number of columns (seq2 along columns)

    # -------------------------------------------------------------------------
    # SAFETY CHECK: if sequences are very different in length and k is too
    # small, the band won't even reach the bottom-right corner.
    # -------------------------------------------------------------------------
    if abs(n - m) > k:
        print(f"WARNING: |len(seq1) - len(seq2)| = {abs(n-m)} > k={k}.")
        print("The band is too narrow to reach the end. Increase k.")
        print("Falling back to full NW...")
        return _full_nw_fallback(seq1, seq2, match, mismatch, gap)

    # -------------------------------------------------------------------------
    # STEP 1: Create DP matrix and traceback matrix.
    #         We create a full (n+1)×(m+1) matrix but will only WRITE to
    #         cells inside the band.  Cells outside get NEG_INF so they
    #         never contribute to the max() calls.
    # -------------------------------------------------------------------------
    # dp[i][j] = best alignment score for seq1[0..i-1] vs seq2[0..j-1]
    dp    = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    trace = [[None]    * (m + 1) for _ in range(n + 1)]

    # -------------------------------------------------------------------------
    # STEP 2: Initialise the starting cell and the boundary rows/columns.
    #         Only initialise cells that are INSIDE the band (|i-j| <= k).
    # -------------------------------------------------------------------------
    dp[0][0] = 0   # no sequences aligned yet → score 0

    # First row: gaps in seq1 (we're consuming seq2 with deletions)
    for j in range(1, min(k + 1, m + 1)):   # only up to the band edge
        dp[0][j]    = j * gap
        trace[0][j] = 'L'   # came from the Left

    # First column: gaps in seq2 (we're consuming seq1 with insertions)
    for i in range(1, min(k + 1, n + 1)):   # only up to the band edge
        dp[i][0]    = i * gap
        trace[i][0] = 'U'   # came from Up

    # -------------------------------------------------------------------------
    # STEP 3: Fill the DP matrix — but ONLY within the band.
    #         For each row i, j can range from max(1, i-k) to min(m, i+k).
    # -------------------------------------------------------------------------
    for i in range(1, n + 1):

        # Calculate the column range for this row that stays inside the band
        j_start = max(1, i - k)   # don't go further left than column 1
        j_end   = min(m, i + k)   # don't go further right than column m

        for j in range(j_start, j_end + 1):

            # --- Option A: DIAGONAL move (match or mismatch) ---
            # We can only come from diagonal if dp[i-1][j-1] is reachable
            if dp[i-1][j-1] != NEG_INF:
                s = match if seq1[i-1] == seq2[j-1] else mismatch
                diag_score = dp[i-1][j-1] + s
            else:
                diag_score = NEG_INF

            # --- Option B: UP move (gap in seq2) ---
            # Coming from the cell above means we align seq1[i-1] with a gap
            if dp[i-1][j] != NEG_INF:
                up_score = dp[i-1][j] + gap
            else:
                up_score = NEG_INF

            # --- Option C: LEFT move (gap in seq1) ---
            # Coming from the left means we align seq2[j-1] with a gap
            if dp[i][j-1] != NEG_INF:
                left_score = dp[i][j-1] + gap
            else:
                left_score = NEG_INF

            # --- Choose the best option ---
            best = max(diag_score, up_score, left_score)
            dp[i][j] = best

            # Record where we came from (for traceback later)
            if best == NEG_INF:
                trace[i][j] = None   # unreachable cell
            elif best == diag_score:
                trace[i][j] = 'D'   # Diagonal
            elif best == up_score:
                trace[i][j] = 'U'   # Up
            else:
                trace[i][j] = 'L'   # Left

    # -------------------------------------------------------------------------
    # STEP 4: Check if we actually reached the bottom-right corner.
    #         If not, the band was too narrow.
    # -------------------------------------------------------------------------
    final_score = dp[n][m]
    if final_score == NEG_INF:
        print("ERROR: Band too narrow, couldn't reach dp[n][m]. Increase k.")
        return _full_nw_fallback(seq1, seq2, match, mismatch, gap)

    # -------------------------------------------------------------------------
    # STEP 5: TRACEBACK — walk back from dp[n][m] to dp[0][0]
    #         following the recorded directions to reconstruct the alignment.
    # -------------------------------------------------------------------------
    align1 = ""   # will hold the aligned version of seq1
    align2 = ""   # will hold the aligned version of seq2
    i, j = n, m  # start at the bottom-right corner

    while i > 0 or j > 0:
        direction = trace[i][j]

        if direction == 'D':
            # Diagonal: both sequences contributed a character
            align1 += seq1[i - 1]
            align2 += seq2[j - 1]
            i -= 1
            j -= 1

        elif direction == 'U':
            # Up: seq1 contributed a character, seq2 gets a gap
            align1 += seq1[i - 1]
            align2 += '-'
            i -= 1

        elif direction == 'L':
            # Left: seq2 contributed a character, seq1 gets a gap
            align1 += '-'
            align2 += seq2[j - 1]
            j -= 1

        else:
            # Shouldn't happen if the band was wide enough
            print(f"ERROR: No traceback direction at ({i},{j}). Stopping.")
            break

    # The alignment was built backwards → reverse both strings
    return align1[::-1], align2[::-1], final_score


# =============================================================================
# HELPER: Full NW fallback (used when band is too narrow)
#         This is just a standard NW without band restriction.
# =============================================================================
def _full_nw_fallback(seq1, seq2, match, mismatch, gap):
    """Standard O(nm) Needleman-Wunsch — used as fallback."""
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
            if   best == dp[i-1][j-1] + s:    trace[i][j] = 'D'
            elif best == dp[i-1][j]   + gap:   trace[i][j] = 'U'
            else:                              trace[i][j] = 'L'

    align1, align2, i, j = "", "", n, m
    while i > 0 or j > 0:
        d = trace[i][j]
        if   d == 'D': align1 += seq1[i-1]; align2 += seq2[j-1]; i -= 1; j -= 1
        elif d == 'U': align1 += seq1[i-1]; align2 += '-';        i -= 1
        else:          align1 += '-';        align2 += seq2[j-1];          j -= 1
    return align1[::-1], align2[::-1], dp[n][m]


# =============================================================================
# Tests
# =============================================================================
if __name__ == "__main__":
    print("=" * 60)
    print("BANDED NEEDLEMAN-WUNSCH DEMO")
    print("=" * 60)

    # Test 1: Very similar sequences → small band should work fine
    s1 = "GATTACAACTTG"
    s2 = "GATTACAACTTG"   # identical
    a1, a2, score = banded_nw(s1, s2, k=3)
    print(f"\nTest 1 (identical sequences, k=3):")
    print(f"  Seq1: {a1}")
    print(f"  Seq2: {a2}")
    print(f"  Score: {score}")

    # Test 2: Similar but with one mismatch
    s1 = "GATTACAACTTG"
    s2 = "GATCACAACTTG"   # one mismatch at position 3
    a1, a2, score = banded_nw(s1, s2, k=3)
    print(f"\nTest 2 (one mismatch, k=3):")
    print(f"  Seq1: {a1}")
    print(f"  Seq2: {a2}")
    print(f"  Score: {score}")

    # Test 3: One small gap  → need k >= 1
    s1 = "GATTACA"
    s2 = "GATACA"   # one character shorter
    a1, a2, score = banded_nw(s1, s2, k=2)
    print(f"\nTest 3 (one gap, k=2):")
    print(f"  Seq1: {a1}")
    print(f"  Seq2: {a2}")
    print(f"  Score: {score}")

    # Test 4: Band too narrow → should auto-fallback
    s1 = "AAAAAAAAAAA"   # length 11
    s2 = "TTTTTT"        # length 6 — difference of 5
    a1, a2, score = banded_nw(s1, s2, k=2)   # k=2 < 5 → will warn & fallback
    print(f"\nTest 4 (band too narrow, k=2):")
    print(f"  Seq1: {a1}")
    print(f"  Seq2: {a2}")
    print(f"  Score: {score}")

    print("\nAll tests complete!")
