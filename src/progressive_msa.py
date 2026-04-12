from data_loader import load_fasta
from distance_matrix import generate_matrix
from upgma import run_upgma
from needlemanWunschGlobal import needleman_wunsch

def apply_gaps(seq, ref_ungapped, ref_gapped):
    """
    Inserting gaps into a sequence to match the pattern of a reference alignment.
    """
    result, si = '', 0
    for rc in ref_gapped:
        if rc == '-':
            # Adding a gap if the reference has a gap
            result += '-'
        else:
            # Adding the character from the sequence or padding if the sequence ends
            result += seq[si] if si < len(seq) else '-'
            si += 1
    return result

def perform_progressive_alignment(tree, sequences):
    """
    Executing the recursive alignment process by following the guide tree.
    """
    # Returning a list with one sequence when hitting a leaf node
    if isinstance(tree, str):
        return [sequences[tree]]

    # Recursively fetching aligned sequences from the left and right branches
    left_seqs = perform_progressive_alignment(tree[0], sequences)
    right_seqs = perform_progressive_alignment(tree[1], sequences)

    # Selecting the first sequence of each group as the representative for aligning
    rep_left = left_seqs[0].replace('-', '')
    rep_right = right_seqs[0].replace('-', '')

    print(f"Aligning: {tree[0]} with {tree[1]}")
    aln_left, aln_right, score = needleman_wunsch(rep_left, rep_right)

    # Applying the new gap pattern to every sequence in the left group
    aligned_left = [apply_gaps(s.replace('-', ''), rep_left, aln_left)
                    for s in left_seqs]

    # Applying the new gap pattern to every sequence in the right group
    aligned_right = [apply_gaps(s.replace('-', ''), rep_right, aln_right)
                     for s in right_seqs]

    # Returning the combined list of all aligned sequences
    return aligned_left + aligned_right

if __name__ == "__main__":
    # Loading the sequences from the fasta file
    data = load_fasta("globins.fasta")

    # Generating the distance matrix and creating the UPGMA guide tree
    labels, matrix = generate_matrix(data)
    guide_tree = run_upgma(labels, matrix)

    print("\n--- Starting Progressive Alignment ---")

    # Running the final progressive alignment
    msa_results = perform_progressive_alignment(guide_tree, data)

    # Printing the final Multiple Sequence Alignment (MSA)
    print("\n--- Final Multiple Sequence Alignment ---")
    for i, seq in enumerate(msa_results):
        print(f"Sequence {i+1}: {seq}")