# Progressive alignment builds an MSA by following a guide tree:
#   1. Compute pairwise distances (identity-based)
#   2. Build a guide tree with UPGMA clustering
#   3. Align sequences bottom-up along the tree:
#      - Leaf nodes return their single sequence
#      - Internal nodes align the two child groups using representative sequences
#      - Gaps introduced during alignment are propagated to all group members
#
# Known weakness: errors introduced early in the tree cannot be corrected later.
# This is addressed by iterative refinement (see iterative_refinement.py).
#
# Time complexity: O(n² * L²) where n=number of sequences, L=alignment length
# Space complexity: O(n * L)

from data_loader       import load_fasta
from distance_matrix   import generate_matrix
from upgma             import run_upgma
from needlemanWunschGlobal import needleman_wunsch


def apply_gaps(seq, ref_ungapped, ref_gapped):
    """
    Propagating a new gap pattern from a re-aligned representative to a sequence.

    For each character in the newly gapped representative:
      - If it is a '-', we insert a gap in the sequence too
      - Otherwise, we take the next character from the original sequence

    Parameters
    ----------
    seq           : str  - the original aligned sequence (may already have gaps)
    ref_ungapped  : str  - the representative with all gaps removed
    ref_gapped    : str  - the representative after re-alignment (with new gaps)

    Returns
    -------
    str - the sequence with the new gap pattern applied
    """
    result, si = '', 0
    for rc in ref_gapped:
        if rc == '-':
            result += '-'
        else:
            result += seq[si] if si < len(seq) else '-'
            si += 1
    return result


def perform_progressive_alignment(tree, sequences):
    """
    Recursively aligning sequences by following the guide tree structure.

    Parameters
    ----------
    tree      : str or tuple  - leaf (sequence name) or internal node (left, right)
    sequences : dict          - {name: sequence} mapping from the FASTA file

    Returns
    -------
    list of str - all aligned sequences from this subtree
    """
    # Base case: leaf node — just return the single raw sequence in a list
    if isinstance(tree, str):
        return [sequences[tree]]

    # Recursing down both branches of the tree
    left_seqs  = perform_progressive_alignment(tree[0], sequences)
    right_seqs = perform_progressive_alignment(tree[1], sequences)

    # Using the first sequence from each group as the representative
    # We strip gaps before alignment, then propagate back afterwards
    rep_left  = left_seqs[0].replace('-', '')
    rep_right = right_seqs[0].replace('-', '')

    print(f"  Aligning group '{tree[0]}' with group '{tree[1]}'")
    aln_left, aln_right, _ = needleman_wunsch(rep_left, rep_right)

    # Propagating the new gap pattern to every member of each group
    aligned_left  = [apply_gaps(s.replace('-', ''), rep_left,  aln_left)  for s in left_seqs]
    aligned_right = [apply_gaps(s.replace('-', ''), rep_right, aln_right) for s in right_seqs]

    return aligned_left + aligned_right


if __name__ == "__main__":
    data = load_fasta("globins.fasta")
    if not data:
        print("Could not load globins.fasta — check that the file is in the data/ folder.")
        exit()

    labels, matrix   = generate_matrix(data)
    guide_tree        = run_upgma(labels, matrix)

    print("\n=== Progressive Alignment ===")
    msa_results = perform_progressive_alignment(guide_tree, data)

    print("\n=== Final Multiple Sequence Alignment ===")
    for i, seq in enumerate(msa_results):
        print(f"Sequence {i+1}: {seq}")