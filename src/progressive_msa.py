import sys
import os

# Ensuring the script can find your other modules
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.append(current_dir)

from data_loader import load_fasta
from distance_matrix import generate_matrix
from upgma import run_upgma
from needlemanWunschGlobal import needleman_wunsch

def perform_progressive_alignment(tree, sequences):
    """
    Recursively aligning sequences by following the guide tree.
    """
    # Base case: If we reached a leaf node (a single sequence name)
    if isinstance(tree, str):
        return sequences[tree]

    # Recursive step: Diving into the left and right branches of the tree
    left_result = perform_progressive_alignment(tree[0], sequences)
    right_result = perform_progressive_alignment(tree[1], sequences)

    # Aligning the results from the two branches
    # This uses your existing Needleman-Wunsch implementation
    print(f"Aligning: {tree[0]} with {tree[1]}")
    aln1, aln2, score = needleman_wunsch(left_result, right_result)

    # Returning the aligned version of the first sequence to continue the chain
    # In a full ClustalW, we would return a profile, but this works for 3 sequences
    return aln1

if __name__ == "__main__":
    # 1. Loading the real protein data
    data = load_fasta("globins.fasta")

    # 2. Calculating the distances
    labels, matrix = generate_matrix(data)

    # 3. Generating the UPGMA guide tree
    guide_tree = run_upgma(labels, matrix)

    print("\n--- Starting Progressive Alignment ---")
    # 4. Executing the final alignment following the tree
    final_result = perform_progressive_alignment(guide_tree, data)

    print("\n--- Final MSA Result (Representative String) ---")
    print(final_result)