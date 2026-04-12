from needlemanWunschGlobal import needleman_wunsch

def sum_of_pairs(msa):
    """Calculating the total Sum-of-Pairs (SP) score for the entire alignment."""
    total_score = 0
    num_seqs = len(msa)
    for i in range(num_seqs):
        for j in range(i + 1, num_seqs):
            for char1, char2 in zip(msa[i], msa[j]):
                if char1 != '-' and char2 != '-':
                    total_score += 1 if char1 == char2 else -1
                else:
                    total_score += -2
    return total_score

def apply_group_gaps(group, old_rep, new_rep):
    """Applying the new gap pattern from a re-aligned representative to the whole group."""
    new_group = []
    for seq in group:
        new_seq = ""
        idx = 0
        for char in new_rep:
            if char == "-":
                new_seq += "-"
            else:
                # Taking the character (or gap) from the original MSA position
                new_seq += seq[idx]
                idx += 1
        new_group.append(new_seq)
    return new_group

def refine_once(msa, nw_fn):
    """Performing one round of partition-realign-accept to improve SP score."""
    best_score = sum_of_pairs(msa)
    best_msa = list(msa)
    n = len(msa)

    # Trying every possible "split" point in the MSA
    for pivot in range(1, n):
        group_a = msa[:pivot]
        group_b = msa[pivot:]

        # Creating representatives by removing internal gaps
        rep_a = group_a[0].replace('-', '')
        rep_b = group_b[0].replace('-', '')

        # Re-aligning the two representatives
        a1, a2, _ = nw_fn(rep_a, rep_b)

        # Propagating the new gaps to the entire groups
        new_msa = apply_group_gaps(group_a, rep_a, a1) + \
                  apply_group_gaps(group_b, rep_b, a2)

        new_score = sum_of_pairs(new_msa)

        # Keeping the new alignment only if it's better
        if new_score > best_score:
            best_score = new_score
            best_msa = new_msa

    return best_msa, best_score

if __name__ == "__main__":
    # Test set from previous MSA
    test_msa = [
        "-VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFP-T-TK-TYFPH---FDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPV-NFKLLSHCLLVTLAAH---LPAEFTPA-VHA-SL-DKFLASVSTVLTSKYR",
        "-VLSAADKTNVKAAWSKVGGHAGEYGAEALERMFLGFP-T-TK-TYFPH---FDLSHGSAQVKAHGKKVGDALTLAVGHLDDLPGALSNLSDLHAHKLRVDPV-NFKLLSHCLLSTLAVH---LPNDFTPA-VHA-SL-DKFLSSVSTVLTSKYR",
        "-VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPW-T-QR-FFESF---GDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKL-HVDPENFRLLGNVLVC---VLAHHFGK-EFT-PP-VQAAYQKVVAGVANAL"
    ]

    print("Running Iterative Refinement")
    initial_score = sum_of_pairs(test_msa)
    print(f"Initial Sum-of-Pairs Score: {initial_score}")

    refined_msa, final_score = refine_once(test_msa, needleman_wunsch)
    print(f"Refined Sum-of-Pairs Score: {final_score}")

    if final_score > initial_score:
        print(f"Improvement found. New score: {final_score}")
    else:
        print("No better alignment found in this round.")