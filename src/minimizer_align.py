import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from needlemanWunschGlobal import needleman_wunsch

# Step 1: Minimizer Sketching

def greedy_similarity(seq1, seq2):
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return matches / min(len(seq1), len(seq2))

def get_minimizers(sequence, k, w):
    """
    Selecting the lexicographically smallest k-mer in each sliding window of size w.

    Parameters
    ----------
    sequence : str  - the input sequence
    k        : int  - k-mer length (how specific each anchor is)
    w        : int  - window size (how dense the sketch is; must be >= k)

    Returns
    -------
    list of (kmer, position) — one entry per unique (kmer, position) pair
    """
    minimizers = []
    seen = set()

    for i in range(len(sequence) - w + 1):
        window      = sequence[i:i+w]
        kmers_in_window = [(window[j:j+k], i + j) for j in range(w - k + 1)]
        if not kmers_in_window:
            continue
        min_kmer, min_pos = min(kmers_in_window, key=lambda x: x[0])
        key = (min_kmer, min_pos)
        if key not in seen:
            minimizers.append(key)
            seen.add(key)

    return minimizers

# Step 2: Finding Shared Anchors

def find_anchors(query, target, k, w):
    """
    Finding positions where query and target share the same minimizer k-mer.
    Each anchor is a (q_start, t_start, k) tuple meaning those k characters match.

    Returns a list of anchors sorted by query position.
    """
    query_mins  = get_minimizers(query,  k, w)
    target_mins = get_minimizers(target, k, w)

    # Building a lookup from k-mer to all target positions where it appears
    target_lookup = {}
    for kmer, t_pos in target_mins:
        target_lookup.setdefault(kmer, []).append(t_pos)

    anchors = []
    for kmer, q_pos in query_mins:
        if kmer in target_lookup:
            for t_pos in target_lookup[kmer]:
                anchors.append((q_pos, t_pos, k))

    anchors.sort(key=lambda x: x[0])
    return anchors

# Step 3: Co-Linear Chaining (DP similar to Longest Increasing Subsequence)


def chain_anchors(anchors):
    """
    Selecting the best subset of anchors that appear in the same order
    in both sequences (co-linear) and maximizes total matched bases.

    Anchor B can follow anchor A only if:
        B.q_start >= A.q_start + A.length  (no query overlap)
        B.t_start >= A.t_start + A.length  (no target overlap)

    Returns the best chain as a list of (q_start, t_start, length) tuples.
    """
    if not anchors:
        return []

    n = len(anchors)
    dp     = [anchors[i][2] for i in range(n)]  # initializing with own anchor length
    parent = [-1] * n

    for i in range(1, n):
        q_i, t_i, len_i = anchors[i]
        for j in range(i):
            q_j, t_j, len_j = anchors[j]
            # Checking the co-linearity condition
            if q_i >= q_j + len_j and t_i >= t_j + len_j:
                candidate = dp[j] + len_i
                if candidate > dp[i]:
                    dp[i]     = candidate
                    parent[i] = j

    # Reconstructing the chain by following backpointers
    best_end = max(range(n), key=lambda i: dp[i])
    chain, idx = [], best_end
    while idx != -1:
        chain.append(anchors[idx])
        idx = parent[idx]

    chain.reverse()
    return chain

# Step 4: Filling Gaps Between Anchors with NW

def fill_gaps_with_nw(query, target, chain):
    """
    Building the full alignment by:
      - Using NW to align the sequence regions between consecutive anchors
      - Copying anchor regions directly (they are already exact matches)

    Parameters
    ----------
    query, target : str   - original un-gapped sequences
    chain         : list  - co-linear anchor chain from chain_anchors()

    Returns
    -------
    align_q, align_t : str  - complete aligned strings
    total_score      : int  - approximate alignment score
    """
    align_q, align_t = "", ""
    total_score = 0
    q_cursor, t_cursor = 0, 0

    for (q_start, t_start, length) in chain:
        # Aligning the gap region before this anchor
        q_gap = query [q_cursor:q_start]
        t_gap = target[t_cursor:t_start]
        if q_gap or t_gap:
            g1, g2, gscore = needleman_wunsch(q_gap, t_gap)
            align_q     += g1
            align_t     += g2
            total_score += gscore

        # Copying the anchor region directly (it's an exact match)
        anchor_seq   = query[q_start:q_start+length]
        align_q     += anchor_seq
        align_t     += anchor_seq
        total_score += length * 2  # +2 per matching character (default match reward)

        q_cursor = q_start + length
        t_cursor = t_start + length

    # Aligning any remaining tail after the last anchor
    q_tail = query [q_cursor:]
    t_tail = target[t_cursor:]
    if q_tail or t_tail:
        g1, g2, gscore = needleman_wunsch(q_tail, t_tail)
        align_q     += g1
        align_t     += g2
        total_score += gscore

    return align_q, align_t, total_score


# Step 5: Full Pipeline

def minimizer_align(query, target, k=4, w=8):

    anchors = find_anchors(query, target, k, w)

    if not anchors:

        a1, a2, score = needleman_wunsch(query, target)
        return a1, a2, score, []

    chain = chain_anchors(anchors)
    align_q, align_t, score = fill_gaps_with_nw(query, target, chain)

    return align_q, align_t, score, chain


if __name__ == "__main__":

    q = "GATTACAGATTACA"
    t = "GATTAGAGATTACA"
    a_q, a_t, score, chain = minimizer_align(q, t, k=3, w=5)

    print(f"\nQuery:  {q}")
    print(f"Target: {t}")
    print(f"Anchors used: {len(chain)}")
    print(f"Alignment:\n  {a_q}\n  {a_t}\n  Score: {score}")

    similarity = greedy_similarity(a_q, a_t)
    print(f"Percent identity: {similarity:.2%}")

    # Test 2: Identical sequences
    q2 = "ACGTACGTACGTACGT"
    a_q2, a_t2, score2, chain2 = minimizer_align(q2, q2, k=4, w=6)
    print(f"\nIdentical sequences — Score: {score2}, Anchors: {len(chain2)}")

    similarity = greedy_similarity(a_q2, a_t2)
    print(f"Percent identity: {similarity:.2%}")

    # Test 3: No common k-mers → fallback expected
    print("\nNo common minimizers:")
    a_q3, a_t3, score3, chain3=minimizer_align("AAAAAAA", "TTTTTTT", k=4, w=5)

    similarity = greedy_similarity(a_q3, a_t3)
    print(f"Percent identity: {similarity:.2%}")
