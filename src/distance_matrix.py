# Used as input to the UPGMA guide tree builder for progressive MSA.
# Distance is calculated from alignment identity:
#   distance = 1.0 - (matches / alignment_length)
#
# A distance of 0.0 means the sequences are identical.
# A distance of 1.0 means no residues matched at all.

import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.append(current_dir)

from needlemanWunschGlobal import needleman_wunsch


def calculate_distance(seq1, seq2):
    """
    Estimating evolutionary distance between two sequences using alignment identity.
    Identity = matched positions / total alignment length.
    Distance = 1 - identity.
    """
    result = needleman_wunsch(seq1, seq2)

    # Handling both tuple return (aligned strings + score) and plain matrix fallback
    if isinstance(result, tuple):
        aln1, aln2 = result[0], result[1]
    else:
        aln1, aln2 = seq1, seq2

    matches = sum(1 for a, b in zip(aln1, aln2) if a == b)
    length  = max(len(aln1), len(aln2))

    identity = matches / length if length > 0 else 0
    return round(1.0 - identity, 4)


def generate_matrix(sequence_dict):
    """
    Computing an all-vs-all pairwise distance matrix for a dictionary of sequences.
    The result is a symmetric matrix where matrix[i][j] = distance(seq_i, seq_j).

    Returns
    -------
    names  : list of str        - sequence identifiers
    matrix : list of list[float] - symmetric distance matrix
    """
    names  = list(sequence_dict.keys())
    n      = len(names)
    matrix = [[0.0] * n for _ in range(n)]

    print(f"Computing pairwise distances for {n} sequences...")

    for i in range(n):
        for j in range(i + 1, n):
            dist          = calculate_distance(sequence_dict[names[i]],
                                               sequence_dict[names[j]])
            matrix[i][j]  = dist
            matrix[j][i]  = dist  # Matrix is symmetric

    return names, matrix


def print_matrix(names, matrix):
    """Printing the distance matrix in a tab-separated format."""
    print("\nDistance Matrix (identity-based)")
    print("\t" + "\t".join(names))
    for i, row in enumerate(matrix):
        print(f"{names[i]}\t" + "\t".join(map(str, row)))


if __name__ == "__main__":
    from data_loader import load_fasta

    real_data = load_fasta("globins.fasta")

    if real_data:
        labels, dist_matrix = generate_matrix(real_data)
        print_matrix(labels, dist_matrix)
    else:
        print("Failed to load data.")
