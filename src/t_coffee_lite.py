# =============================================================================
# T-Coffee Style Consistency Checker
#
# Measures agreement between Global (Needleman-Wunsch) and Local (Smith-Waterman)
# alignments by comparing aligned residues (ignoring gaps). Higher consistency
# indicates the two methods agree on which residues should align.
#
# FUNCTION:
# - check_consistency(nw_res, sw_res) → percentage of matching non-gap positions
#
# USE CASE: Building consistency libraries for progressive MSA (T-Coffee)
# =============================================================================

from needlemanWunschGlobal import needleman_wunsch
from smithWatermanLocal import smith_waterman

def check_consistency(nw_res, sw_res):
    """Evaluating how much the Global and Local alignments agree."""
    # nw_res[0] and sw_res[0] are the first aligned sequences
    matches = 0
    # We compare the characters in the aligned strings
    len_to_check = min(len(nw_res[0]), len(sw_res[0]))

    for i in range(len_to_check):
        if nw_res[0][i] == sw_res[0][i] and nw_res[0][i] != '-':
            matches += 1

    return (matches / len_to_check) * 100

if __name__ == "__main__":
    s1, s2 = "GATTACA", "GATCA"

    # Running both algorithms
    nw = needleman_wunsch(s1, s2)
    sw = smith_waterman(s1, s2)

    # Calculating the score
    consistency = check_consistency(nw, sw)

    print("T-Coffee Consistency Library")
    print(f"Consistency Score: {consistency:.2f}%")