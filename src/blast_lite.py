# =============================================================================
#   BLAST-Lite: Heuristic Local Alignment (Seed + X-drop Extension)
#
#   get_kmers    → Index all k-mers for O(1) lookup
#   extend_left  → Walk left from seed, stop on X-drop condition
#   extend_right → Walk right from seed, same X-drop logic
#   blast_lite   → Orchestrate: seed → extend → filter → sort by score
#
#   Time: O((n+m)*k) | Space: O(n) | Trade-off: speed vs optimality
# =============================================================================

def get_kmers(sequence, k):
    """Generating all possible k-mers and their starting positions from a sequence."""
    kmers = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if kmer not in kmers:
            kmers[kmer] = []
        kmers[kmer].append(i)
    return kmers

def extend_left(seq1, seq2, i, j, match_score=2, mismatch_penalty=-1, threshold=5):
    """Extending a seed match to the left using X-dropoff logic."""
    l_score, max_l_score = 0, 0
    l_ext1, l_ext2 = "", ""
    # Starting from the character immediately before the seed
    curr_i, curr_j = i - 1, j - 1

    while curr_i >= 0 and curr_j >= 0:
        # Calculating score for the current pair
        score_change = match_score if seq1[curr_i] == seq2[curr_j] else mismatch_penalty
        l_score += score_change

        if l_score > max_l_score:
            max_l_score = l_score
        elif l_score < max_l_score - threshold:
            break

        l_ext1 += seq1[curr_i]
        l_ext2 += seq2[curr_j]
        curr_i -= 1
        curr_j -= 1

    # Reversing strings since they were collected walking backwards
    return l_ext1[::-1], l_ext2[::-1], max_l_score

def extend_right(seq1, seq2, i, j, match_score=2, mismatch_penalty=-1, threshold=5):
    """Extending a seed match to the right using X-dropoff logic."""
    r_score, max_r_score = 0, 0
    r_ext1, r_ext2 = "", ""
    curr_i, curr_j = i, j

    while curr_i < len(seq1) and curr_j < len(seq2):
        score_change = match_score if seq1[curr_i] == seq2[curr_j] else mismatch_penalty
        r_score += score_change

        if r_score > max_r_score:
            max_r_score = r_score
        elif r_score < max_r_score - threshold:
            break

        r_ext1 += seq1[curr_i]
        r_ext2 += seq2[curr_j]
        curr_i += 1
        curr_j += 1

    return r_ext1, r_ext2, max_r_score

def blast_lite(query, target, k=3, threshold=5):
    """Performing a heuristic bidirectional seed-and-extend alignment."""
    target_kmers = get_kmers(target, k)
    hsp_results = []
    match_val = 2 # Standard match score for calculation

    for i in range(len(query) - k + 1):
        query_kmer = query[i:i+k]

        if query_kmer in target_kmers:
            for j in target_kmers[query_kmer]:
                # 1. Getting the left extension (starts before index i and j)
                l_ext1, l_ext2, l_score = extend_left(query, target, i, j, threshold=threshold)

                # 2. Getting the right extension (starts after the k-mer)
                r_ext1, r_ext2, r_score = extend_right(query, target, i + k, j + k, threshold=threshold)

                # 3. Combining left + seed + right
                full_ext1 = l_ext1 + query_kmer + r_ext1
                full_ext2 = l_ext2 + query_kmer + r_ext2

                # Calculating total score (Left + Seed + Right)
                seed_score = k * match_val
                total_score = l_score + seed_score + r_score

                if total_score > k * 1.5:
                    hsp_results.append({
                        "query_pos": i - len(l_ext1),
                        "target_pos": j - len(l_ext2),
                        "alignment": (full_ext1, full_ext2),
                        "score": total_score
                    })

    return sorted(hsp_results, key=lambda x: x['score'], reverse=True)

if __name__ == "__main__":
    q = "GGAGTCAG"
    t = "GAAGTCGG"

    results = blast_lite(q, t, k=3)

    print(f"--- BLAST-style Bidirectional Results ---")
    seen_alignments = set()
    for res in results:
        if res['alignment'] not in seen_alignments:
            print(f"Score: {res['score']} | Query Start: {res['query_pos']}, Target Start: {res['target_pos']}")
            print(f"Seq1: {res['alignment'][0]}")
            print(f"Seq2: {res['alignment'][1]}\n")
            seen_alignments.add(res['alignment'])