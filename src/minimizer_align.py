def get_minimizers(sequence, k, w):
    """Picking the lexicographically smallest k-mer within every window of size w."""
    minimizers = set()
    for i in range(len(sequence) - w + 1):
        window = sequence[i : i + w]
        kmers_in_window = [window[j : j + k] for j in range(len(window) - k + 1)]
        if kmers_in_window:
            m = min(kmers_in_window)
            pos = i + window.find(m)
            minimizers.add((m, pos))
    return minimizers

def align_with_minimizers(query, target, k=3, w=5):
    """Matching minimizers between two sequences to find alignment anchors."""
    # Sketching both sequences
    query_sketch = get_minimizers(query, k, w)
    target_sketch = get_minimizers(target, k, w)

    # Converting target sketch to a dictionary for fast lookup
    target_lookup = {}
    for kmer, pos in target_sketch:
        if kmer not in target_lookup:
            target_lookup[kmer] = []
        target_lookup[kmer].append(pos)

    anchors = []
    # Finding shared minimizers (anchors)
    for kmer, q_pos in query_sketch:
        if kmer in target_lookup:
            for t_pos in target_lookup[kmer]:
                # Calculating the relative offset (helps identify diagonal matches)
                offset = t_pos - q_pos
                anchors.append((q_pos, t_pos, kmer, offset))

    # Sorting anchors by query position
    return sorted(anchors, key=lambda x: x[0])

if __name__ == "__main__":
    # Defining two sequences with a small mutation
    seq_q = "GATTACAGATTACA"
    seq_t = "GATTAGAGATTACA"

    results = align_with_minimizers(seq_q, seq_t, k=3, w=5)

    print(f"--- Minimizer-based Anchors ---")
    for q_pos, t_pos, kmer, offset in results:
        status = "Aligned" if offset == 0 else "Shifted"
        print(f"K-mer: '{kmer}' | Q: {q_pos} -> T: {t_pos} | [{status}]")