def get_kmers(sequence, k):
    """Generating all possible k-mers and their starting positions from a sequence."""
    kmers = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if kmer not in kmers:
            kmers[kmer] = []
        kmers[kmer].append(i)
    return kmers

def extend_match(seq1, seq2, i, j, match_score=2, mismatch_penalty=-1, threshold=5):
    """Extending a seed match to the right until the score drops significantly."""
    r_score = 0
    max_r_score = 0
    r_ext1, r_ext2 = "", ""
    curr_i, curr_j = i, j

    while curr_i < len(seq1) and curr_j < len(seq2):
        # Checking if characters match or mismatch
        score_change = match_score if seq1[curr_i] == seq2[curr_j] else mismatch_penalty
        r_score += score_change

        # Updating the peak score if the current path is better
        if r_score > max_r_score:
            max_r_score = r_score
        # Breaking the loop if the score falls too far below the peak (X-dropoff)
        elif r_score < max_r_score - threshold:
            break

        r_ext1 += seq1[curr_i]
        r_ext2 += seq2[curr_j]
        curr_i += 1
        curr_j += 1

    return r_ext1, r_ext2, max_r_score

def blast_lite(query, target, k=3, threshold=5):
    """Performing a heuristic seed-and-extend alignment between two sequences."""
    # Indexing the target sequence for fast lookup
    target_kmers = get_kmers(target, k)
    hsp_results = []

    # Iterating through the query to find potential seeds
    for i in range(len(query) - k + 1):
        query_kmer = query[i:i+k]

        if query_kmer in target_kmers:
            # Checking every location where this seed exists in the target
            for j in target_kmers[query_kmer]:
                # Sending the coordinates to the extension function
                ext1, ext2, score = extend_match(query, target, i, j, match_score=2, threshold=threshold)

                # Filtering for results that are stronger than just a random seed match
                if score > k * 1.5:
                    hsp_results.append({
                        "query_pos": i,
                        "target_pos": j,
                        "alignment": (ext1, ext2),
                        "score": score
                    })

    # Returning the best hits sorted by their total score
    return sorted(hsp_results, key=lambda x: x['score'], reverse=True)

if __name__ == "__main__":
    # Defining sample sequences for verification
    q = "GGAGTCAG"
    t = "GAAGTCGG"

    results = blast_lite(q, t, k=3)

    print(f"--- BLAST-style Heuristic Results ---")
    # Displaying the top unique hits found
    seen_alignments = set()
    for res in results:
        if res['alignment'] not in seen_alignments:
            print(f"Score: {res['score']} | Query Start: {res['query_pos']}, Target Start: {res['target_pos']}")
            print(f"Seq1: {res['alignment'][0]}")
            print(f"Seq2: {res['alignment'][1]}\n")
            seen_alignments.add(res['alignment'])