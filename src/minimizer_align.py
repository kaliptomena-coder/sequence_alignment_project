# =============================================================================
# Minimizer-based approximate alignment with anchor chaining:
#   fill the gaps BETWEEN chained anchors using standard NW alignment.
#   This gives us a complete alignment string with '-' characters.
# =============================================================================


# We need NW for gap-filling between anchors
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Try to import NW implementation; fall back to a local copy if not found
try:
    from needlemanWunschGlobal import needleman_wunsch
except ImportError:
    def needleman_wunsch(s1, s2, match=2, mismatch=-1, gap=-2):
        """Minimal NW fallback."""
        n, m = len(s1), len(s2)
        dp = [[0]*(m+1) for _ in range(n+1)]
        for i in range(n+1): dp[i][0] = i*gap
        for j in range(m+1): dp[0][j] = j*gap
        for i in range(1,n+1):
            for j in range(1,m+1):
                s = match if s1[i-1]==s2[j-1] else mismatch
                dp[i][j] = max(dp[i-1][j-1]+s, dp[i-1][j]+gap, dp[i][j-1]+gap)
        a1,a2,i,j="","",n,m
        while i>0 or j>0:
            if i>0 and j>0 and dp[i][j]==dp[i-1][j-1]+(match if s1[i-1]==s2[j-1] else mismatch):
                a1+=s1[i-1]; a2+=s2[j-1]; i-=1; j-=1
            elif i>0 and dp[i][j]==dp[i-1][j]+gap:
                a1+=s1[i-1]; a2+='-'; i-=1
            else:
                a1+='-'; a2+=s2[j-1]; j-=1
        return a1[::-1], a2[::-1], dp[n][m]


# =============================================================================
# 1. Minimizer sketching (same as before, slightly cleaned up)
# =============================================================================

def get_minimizers(sequence, k, w):
    """
    Select the lexicographically smallest k-mer in each window of size w.

    A minimizer is a compact representation of a sequence region.
    Two sequences that share a minimizer likely share that region.

    Parameters
    ----------
    sequence : str  – the input sequence
    k        : int  – k-mer length
    w        : int  – window size (w >= k)

    Returns
    -------
    list of (kmer, position) tuples — one per window
    """
    minimizers = []
    seen = set()   # avoid exact duplicates

    for i in range(len(sequence) - w + 1):
        window   = sequence[i : i + w]
        # All k-mers in this window
        kmers_in_window = [(window[j : j + k], i + j)
                           for j in range(w - k + 1)]
        if not kmers_in_window:
            continue
        # Pick the lexicographically smallest k-mer
        min_kmer, min_pos = min(kmers_in_window, key=lambda x: x[0])

        key = (min_kmer, min_pos)
        if key not in seen:
            minimizers.append(key)
            seen.add(key)

    return minimizers


# =============================================================================
# 2. Find shared anchors (exact k-mer matches between query and target)
# =============================================================================

def find_anchors(query, target, k, w):
    """
    Find positions where query and target share a minimizer k-mer.

    Each anchor is a (q_start, t_start, length) tuple meaning:
      query[q_start : q_start+k] == target[t_start : t_start+k]

    Parameters
    ----------
    query, target : str  – the two sequences to compare
    k, w          : int  – minimizer parameters

    Returns
    -------
    list of (q_start, t_start, k) tuples, sorted by q_start
    """
    query_mins  = get_minimizers(query,  k, w)
    target_mins = get_minimizers(target, k, w)

    # Build a lookup: kmer → list of target positions
    target_lookup = {}
    for kmer, t_pos in target_mins:
        target_lookup.setdefault(kmer, []).append(t_pos)

    anchors = []
    for kmer, q_pos in query_mins:
        if kmer in target_lookup:
            for t_pos in target_lookup[kmer]:
                anchors.append((q_pos, t_pos, k))   # (q_start, t_start, length)

    # Sort by query start position
    anchors.sort(key=lambda x: x[0])
    return anchors


# =============================================================================
# 3. CO-LINEAR CHAINING
#
# We have a set of anchors (q_pos, t_pos, length).
#          We want to pick the subset that:
#            - is ordered: anchor i comes before anchor j means
#              q_pos[i] < q_pos[j]  AND  t_pos[i] < t_pos[j]
#            - maximises total matched bases (sum of anchor lengths)
#
# This is solved with dynamic programming — similar to Longest Increasing
# Subsequence (LIS), but we maximise weight instead of count.
# =============================================================================

def chain_anchors(anchors):
    """
    Find the best co-linear chain of anchors using dynamic programming.

    An anchor B can follow anchor A only if:
        B.q_start >= A.q_start + A.length   (no overlap in query)
        B.t_start >= A.t_start + A.length   (no overlap in target)

    We maximise the total matched bases (sum of anchor lengths in the chain).

    Parameters
    ----------
    anchors : list of (q_start, t_start, length) tuples

    Returns
    -------
    list of (q_start, t_start, length) – the best chain, in order
    """
    if not anchors:
        return []

    n = len(anchors)

    # dp[i] = best total matched bases for a chain ending at anchor i
    dp      = [anchors[i][2] for i in range(n)]   # initialise with own length
    parent  = [-1] * n                             # backpointer

    for i in range(1, n):
        q_i, t_i, len_i = anchors[i]

        for j in range(i):
            q_j, t_j, len_j = anchors[j]

            # Can anchor i follow anchor j?
            if q_i >= q_j + len_j and t_i >= t_j + len_j:
                candidate = dp[j] + len_i
                if candidate > dp[i]:
                    dp[i]     = candidate
                    parent[i] = j

    # Find the chain with the maximum score
    best_end = max(range(n), key=lambda i: dp[i])

    # Traceback to reconstruct the chain
    chain   = []
    idx     = best_end
    while idx != -1:
        chain.append(anchors[idx])
        idx = parent[idx]

    chain.reverse()   # built backwards → reverse
    return chain


# =============================================================================
# 4. Fill gaps between chained anchors with NW alignment
#
# After chaining, we have a list of anchor regions where query and target
# already agree.  Between consecutive anchors, we need to align the
# remaining characters using NW.
# =============================================================================

def fill_gaps_with_nw(query, target, chain):
    """
    Build a complete alignment string by:
      1. For each gap between chained anchors, run NW alignment
      2. For each anchor region, emit the matching characters directly

    Parameters
    ----------
    query, target : str  – original (un-gapped) sequences
    chain         : list of (q_start, t_start, length) tuples

    Returns
    -------
    align_q, align_t : str – full alignment strings (with '-' gaps)
    score            : int – approximate alignment score (sum of NW sub-scores)
    """
    align_q = ""
    align_t = ""
    total_score = 0

    # Pointers to where we are in each sequence
    q_cursor = 0
    t_cursor = 0

    for (q_start, t_start, length) in chain:

        # -----------------------------------------------------------------
        # Gap BEFORE this anchor: align the region between cursor and anchor
        # -----------------------------------------------------------------
        q_gap_region = query [q_cursor : q_start]
        t_gap_region = target[t_cursor : t_start]

        if q_gap_region or t_gap_region:
            # Use NW to align the gap region
            g1, g2, gscore = needleman_wunsch(q_gap_region, t_gap_region)
            align_q     += g1
            align_t     += g2
            total_score += gscore

        # -----------------------------------------------------------------
        # ANCHOR REGION: characters match directly, no gaps needed
        # -----------------------------------------------------------------
        anchor_seq = query[q_start : q_start + length]
        align_q    += anchor_seq
        align_t    += anchor_seq   # same characters (it's a match)
        total_score += length * 2  # +2 per matching character (default match score)

        # Advance cursors to just after this anchor
        q_cursor = q_start + length
        t_cursor = t_start + length

    # -----------------------------------------------------------------
    # Gap AFTER the last anchor: align remaining tails
    # -----------------------------------------------------------------
    q_tail = query [q_cursor:]
    t_tail = target[t_cursor:]

    if q_tail or t_tail:
        g1, g2, gscore = needleman_wunsch(q_tail, t_tail)
        align_q     += g1
        align_t     += g2
        total_score += gscore

    return align_q, align_t, total_score


# =============================================================================
# 5. Main function — the complete pipeline
# =============================================================================

def minimizer_align(query, target, k=4, w=8):
    """
    Full minimizer-based approximate alignment pipeline:
        1. Sketch both sequences with minimizers
        2. Find shared anchors (exact k-mer matches)
        3. Chain anchors co-linearly
        4. Fill gaps between anchors with NW
        5. Return the complete alignment strings

    Parameters
    ----------
    query, target : str  – input sequences
    k             : int  – k-mer size   (larger = more specific anchors)
    w             : int  – window size  (larger = sparser sketch)

    Returns
    -------
    align_q : str  – aligned query  (may contain '-')
    align_t : str  – aligned target (may contain '-')
    score   : int  – alignment score
    chain   : list – the selected anchor chain (for inspection)
    """

    # Step 1 & 2: Find shared anchors
    anchors = find_anchors(query, target, k, w)

    if not anchors:
        # No shared minimizers → fall back to full NW
        print("  No shared minimizers found. Falling back to full NW.")
        a1, a2, score = needleman_wunsch(query, target)
        return a1, a2, score, []

    # Step 3: Chain the anchors
    chain = chain_anchors(anchors)

    # Step 4: Fill gaps between chained anchors
    align_q, align_t, score = fill_gaps_with_nw(query, target, chain)

    return align_q, align_t, score, chain


# =============================================================================
# DEMO
# =============================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("MINIMIZER ALIGNMENT WITH CHAINING DEMO")
    print("=" * 60)

    # Test 1: Sequences with one small mutation
    q = "GATTACAGATTACA"
    t = "GATTAGAGATTACA"

    print(f"\nQuery:  {q}")
    print(f"Target: {t}")

    a_q, a_t, score, chain = minimizer_align(q, t, k=3, w=5)

    print(f"\nChained anchors ({len(chain)}):")
    for q_s, t_s, ln in chain:
        print(f"  q[{q_s}:{q_s+ln}] = t[{t_s}:{t_s+ln}]  '{q[q_s:q_s+ln]}'")

    print(f"\nAlignment:")
    print(f"  Query:  {a_q}")
    print(f"  Target: {a_t}")
    print(f"  Score:  {score}")

    # Test 2: Nearly identical sequences → chain should cover most positions
    q2 = "ACGTACGTACGTACGT"
    t2 = "ACGTACGTACGTACGT"

    print(f"\n\nTest 2 (identical):")
    print(f"  Query:  {q2}")
    a_q2, a_t2, score2, chain2 = minimizer_align(q2, t2, k=4, w=6)
    print(f"  Align:  {a_q2}")
    print(f"  Score:  {score2}")
    print(f"  Anchors used: {len(chain2)}")

    # Test 3: No common k-mers → should fall back to NW
    print(f"\n\nTest 3 (no common minimizers → fallback to NW):")
    a_q3, a_t3, score3, _ = minimizer_align("AAAAAAA", "TTTTTTT", k=4, w=5)
    print(f"  Query:  {a_q3}")
    print(f"  Target: {a_t3}")
    print(f"  Score:  {score3}")

    print("\nAll tests complete!")
