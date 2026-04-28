# Measures how much a Global (Needleman-Wunsch) and Local (Smith-Waterman)
# alignment agree with each other. The idea comes from the T-Coffee algorithm,
# which builds a consistency library by combining multiple pairwise alignments
# before running progressive MSA.
#
# A higher consistency score indicates that both methods agree on which
# residues should be paired — a sign of a reliable alignment.
#
# Consistency is measured column-by-column on the overlapping portion of
# both alignments, counting positions where both align the same non-gap char.

from needlemanWunschGlobal import needleman_wunsch
from smithWatermanLocal    import smith_waterman


def check_consistency(nw_res, sw_res):
    """
    Comparing the first aligned string from each method column by column.

    Parameters
    ----------
    nw_res : tuple  - return value from needleman_wunsch() i.e. (align1, align2, score)
    sw_res : tuple  - return value from smith_waterman()

    Returns
    -------
    float - percentage of overlapping columns where both methods agree (non-gap match)
    """
    nw_aln = nw_res[0]
    sw_aln = sw_res[0]

    # Only comparing up to the length of the shorter alignment
    len_to_check = min(len(nw_aln), len(sw_aln))
    matches = 0

    for i in range(len_to_check):
        if nw_aln[i] == sw_aln[i] and nw_aln[i] != '-':
            matches += 1

    return (matches / len_to_check) * 100 if len_to_check > 0 else 0.0


if __name__ == "__main__":
    s1, s2 = "GATTACA", "GATCA"

    nw = needleman_wunsch(s1, s2)
    sw = smith_waterman(s1, s2)

    consistency = check_consistency(nw, sw)

    print("=== T-Coffee Consistency Check ===")
    print(f"NW alignment: {nw[0]} / {nw[1]}  (score={nw[2]})")
    print(f"SW alignment: {sw[0]} / {sw[1]}  (score={sw[2]})")
    print(f"Consistency Score: {consistency:.2f}%")
