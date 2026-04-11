def get_minimizers(sequence, k, w):
    """Picking the lexicographically smallest k-mer within every window of size w."""
    minimizers = set()
    # Sliding a window of size w across the sequence
    for i in range(len(sequence) - w + 1):
        window = sequence[i : i + w]
        # Finding all k-mers in this window
        kmers_in_window = [window[j : j + k] for j in range(len(window) - k + 1)]
        # Selecting the 'smallest' k-mer (alphabetical order)
        if kmers_in_window:
            m = min(kmers_in_window)
            # Storing the minimizer and its position
            # Finding the first occurrence of this minimizer in the current window
            pos = i + window.find(m)
            minimizers.add((m, pos))
    return minimizers

if __name__ == "__main__":
    seq = "GATTACAGATTACA"
    # Using k=3 and window=5
    mins = get_minimizers(seq, 3, 5)
    print(f"--- Minimizer Sketch for {seq} ---")
    for m, pos in sorted(list(mins), key=lambda x: x[1]):
        print(f"Minimizer: {m} at position {pos}")