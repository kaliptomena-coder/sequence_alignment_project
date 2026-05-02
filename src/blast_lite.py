def get_kmers(sequence, k):

    kmers = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if kmer not in kmers:
            kmers[kmer] = []
        kmers[kmer].append(i)
    return kmers


def extend_left(seq1, seq2, i, j, match_score=2, mismatch_penalty=-1, threshold=5):
    """
    Extending a seed match to the left using the X-drop stopping rule.
    Stops when the score drops more than `threshold` below the running max.

    Parameters
    ----------
    seq1, seq2        : str  - the full query and target sequences
    i, j              : int  - positions immediately before the seed starts
    match_score       : int  - reward for a match
    mismatch_penalty  : int  - penalty for a mismatch
    threshold         : int  - X-drop value (stop when score falls this much)

    Returns
    -------
    l_ext1, l_ext2 : str  - left extensions for each sequence (reversed)
    max_l_score    : int  - best score accumulated during extension
    """
    l_score, max_l_score = 0, 0
    l_ext1, l_ext2 = "", ""
    curr_i, curr_j = i - 1, j - 1

    while curr_i >= 0 and curr_j >= 0:
        score_change = match_score if seq1[curr_i] == seq2[curr_j] else mismatch_penalty
        l_score += score_change

        if l_score > max_l_score:
            max_l_score = l_score
        elif l_score < max_l_score - threshold:
            break  # X-drop condition met — stop extending

        l_ext1 += seq1[curr_i]
        l_ext2 += seq2[curr_j]
        curr_i -= 1
        curr_j -= 1

    # The extensions were collected walking backwards, so we reverse them
    return l_ext1[::-1], l_ext2[::-1], max_l_score


def extend_right(seq1, seq2, i, j, match_score=2, mismatch_penalty=-1, threshold=5):
    """
    Extending a seed match to the right using the X-drop stopping rule.

    Parameters
    ----------
    seq1, seq2        : str  - the full query and target sequences
    i, j              : int  - positions right after the seed ends
    match_score       : int  - reward for a match
    mismatch_penalty  : int  - penalty for a mismatch
    threshold         : int  - X-drop value

    Returns
    -------
    r_ext1, r_ext2 : str  - right extensions for each sequence
    max_r_score    : int  - best score accumulated during extension
    """
    r_score, max_r_score = 0, 0
    r_ext1, r_ext2 = "", ""
    curr_i, curr_j = i, j

    while curr_i < len(seq1) and curr_j < len(seq2):
        score_change = match_score if seq1[curr_i] == seq2[curr_j] else mismatch_penalty
        r_score += score_change

        if r_score > max_r_score:
            max_r_score = r_score
        elif r_score < max_r_score - threshold:
            break  # X-drop condition met

        r_ext1 += seq1[curr_i]
        r_ext2 += seq2[curr_j]
        curr_i += 1
        curr_j += 1

    return r_ext1, r_ext2, max_r_score

def blast_lite(query, target, k=3, threshold=5):

    target_kmers = get_kmers(target, k)
    hsp_results  = []
    match_val    = 2  # Match reward used for seed scoring

    for i in range(len(query) - k + 1):
        query_kmer = query[i:i+k]

        if query_kmer in target_kmers:
            for j in target_kmers[query_kmer]:

                # Left extension starts immediately before the seed
                l_ext1, l_ext2, l_score = extend_left(
                    query, target, i, j, threshold=threshold)

                # Right extension starts immediately after the seed
                r_ext1, r_ext2, r_score = extend_right(
                    query, target, i + k, j + k, threshold=threshold)

                # Combining: left + seed + right
                full_ext1  = l_ext1 + query_kmer + r_ext1
                full_ext2  = l_ext2 + query_kmer + r_ext2
                seed_score = k * match_val
                total_score = l_score + seed_score + r_score

                # Only keeping HSPs that score better than the seed alone by a margin
                if total_score > k * 1.5:
                    hsp_results.append({
                        "query_pos":  i - len(l_ext1),
                        "target_pos": j - len(l_ext2),
                        "alignment":  (full_ext1, full_ext2),
                        "score":      total_score
                    })

    return sorted(hsp_results, key=lambda x: x['score'], reverse=True)


if __name__ == "__main__":
    q = "GGAGTCAG"
    t = "GAAGTCGG"

    results = blast_lite(q, t, k=3)

    print("=== BLAST-style Bidirectional Results ===")
    seen = set()
    for res in results:
        if res['alignment'] not in seen:
            print(f"Score: {res['score']} | Query Start: {res['query_pos']}, "
                  f"Target Start: {res['target_pos']}")
            print(f"  Query:  {res['alignment'][0]}")
            print(f"  Target: {res['alignment'][1]}\n")
            seen.add(res['alignment'])