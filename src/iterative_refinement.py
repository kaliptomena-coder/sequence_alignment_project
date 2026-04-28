# Progressive alignment can make mistakes early in the tree that propagate
# forward. Iterative refinement tries to fix these by:
#   1. Splitting the current MSA into two groups at every possible partition point
#   2. Re-aligning the group representatives with Needleman-Wunsch
#   3. Propagating the new gaps to all sequences in each group
#   4. Accepting the change only if the Sum-of-Pairs (SP) score improves
#   5. Repeating until no improvement is found (convergence)
#
# Sum-of-Pairs score:
#   For every column, for every pair of sequences:
#     +1  for a matching non-gap character
#     -1  for a mismatching non-gap character
#     -2  if either sequence has a gap
#
# Convergence: SP score stops increasing after one full round of all partitions.
#
# Complexity: O(iterations * n² * L) where n=sequences, L=alignment length

def sum_of_pairs(msa):
    """
    Computing the Sum-of-Pairs (SP) score for an MSA.
    Higher SP means a better overall alignment.

    Parameters
    ----------
    msa : list of str  - all aligned sequences (must be the same length)

    Returns
    -------
    int - the total SP score
    """
    total    = 0
    num_seqs = len(msa)
    aln_len  = len(msa[0]) if msa else 0

    for i in range(num_seqs):
        for j in range(i + 1, num_seqs):
            for col in range(aln_len):
                c1, c2 = msa[i][col], msa[j][col]
                if c1 != '-' and c2 != '-':
                    total += 1 if c1 == c2 else -1
                else:
                    total -= 2  # gap penalty

    return total


def apply_group_gaps(group, old_rep, new_rep):
    """
    Propagating a new gap pattern to every sequence in a group.

    After re-aligning two group representatives, the new gaps in each
    representative must be inserted at the same positions in all other
    group members.

    Parameters
    ----------
    group   : list of str  - original aligned sequences in this group
    old_rep : str          - the representative with all gaps removed
    new_rep : str          - the representative after re-alignment (with new gaps)

    Returns
    -------
    list of str - group sequences updated to the new gap pattern
    """
    new_group = []
    for seq in group:
        new_seq = ""
        idx     = 0
        for char in new_rep:
            if char == '-':
                new_seq += '-'        # inserting gap where the representative has one
            else:
                new_seq += seq[idx] if idx < len(seq) else '-'
                idx += 1
        new_group.append(new_seq)
    return new_group


def refine_once(msa, nw_fn):
    """
    Running one round of partition-realign-accept.

    Tries every split point (1 sequence vs the rest, 2 vs the rest, etc.)
    and keeps the partition that gives the biggest SP improvement.

    Parameters
    ----------
    msa   : list of str  - the current MSA
    nw_fn : callable     - the Needleman-Wunsch function to use for re-alignment

    Returns
    -------
    best_msa   : list of str  - the (possibly) improved MSA
    best_score : int          - SP score of the best MSA found this round
    """
    best_score = sum_of_pairs(msa)
    best_msa   = list(msa)

    for pivot in range(1, len(msa)):
        group_a = msa[:pivot]
        group_b = msa[pivot:]

        # Using gap-free representatives for re-alignment
        rep_a = group_a[0].replace('-', '')
        rep_b = group_b[0].replace('-', '')

        a1, a2, _ = nw_fn(rep_a, rep_b)

        new_msa = (apply_group_gaps(group_a, rep_a, a1) +
                   apply_group_gaps(group_b, rep_b, a2))

        new_score = sum_of_pairs(new_msa)
        if new_score > best_score:
            best_score = new_score
            best_msa   = new_msa

    return best_msa, best_score


def iterative_refinement(msa, nw_fn, max_iterations=50):
    """
    Running iterative refinement until the SP score converges.

    Parameters
    ----------
    msa            : list of str  - the initial MSA (e.g., from progressive alignment)
    nw_fn          : callable     - Needleman-Wunsch function
    max_iterations : int          - safety limit to prevent infinite loops

    Returns
    -------
    best_msa    : list of str  - the refined MSA
    final_score : int          - final SP score
    history     : list of int  - SP score recorded after each iteration
    """
    current_msa   = list(msa)
    current_score = sum_of_pairs(current_msa)
    history       = [current_score]

    print(f"  Starting iterative refinement. Initial SP score: {current_score}")

    for iteration in range(1, max_iterations + 1):
        new_msa, new_score = refine_once(current_msa, nw_fn)

        status = 'improved' if new_score > current_score else 'no change'
        print(f"  Round {iteration:3d}: SP = {new_score}  ({status})")

        history.append(new_score)

        if new_score <= current_score:
            print(f"  Converged after {iteration} round(s).")
            break

        current_msa   = new_msa
        current_score = new_score
    else:
        print(f"  WARNING: Reached max_iterations={max_iterations}. May not have converged.")

    return current_msa, current_score, history


if __name__ == "__main__":
    try:
        from needlemanWunschGlobal import needleman_wunsch
    except ImportError:
        # Minimal NW fallback for standalone testing
        def needleman_wunsch(s1, s2, match=1, mismatch=-1, gap=-2):
            n, m = len(s1), len(s2)
            dp = [[0]*(m+1) for _ in range(n+1)]
            for i in range(n+1): dp[i][0] = i*gap
            for j in range(m+1): dp[0][j] = j*gap
            for i in range(1,n+1):
                for j in range(1,m+1):
                    s = match if s1[i-1]==s2[j-1] else mismatch
                    dp[i][j] = max(dp[i-1][j-1]+s, dp[i-1][j]+gap, dp[i][j-1]+gap)
            a1,a2,i,j="","",n,m
            while i>0 or j>0:
                s = match if i>0 and j>0 and s1[i-1]==s2[j-1] else mismatch
                if i>0 and j>0 and dp[i][j]==dp[i-1][j-1]+s:
                    a1+=s1[i-1]; a2+=s2[j-1]; i-=1; j-=1
                elif i>0 and dp[i][j]==dp[i-1][j]+gap:
                    a1+=s1[i-1]; a2+='-'; i-=1
                else:
                    a1+='-'; a2+=s2[j-1]; j-=1
            return a1[::-1], a2[::-1], dp[n][m]

    print("=== Iterative Refinement Demo ===")

    initial_msa = [
        "-VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFP-T-TK-TYFPH---FDLSHGSAQVK",
        "-VLSAADKTNVKAAWSKVGGHAGEYGAEALERMFLGFP-T-TK-TYFPH---FDLSHGSAQVK",
        "-VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPW-T-QR-FFESF---GDLSTPDAVMG",
    ]

    initial_score = sum_of_pairs(initial_msa)
    print(f"\nInitial SP Score: {initial_score}")

    refined_msa, final_score, history = iterative_refinement(
        initial_msa, needleman_wunsch, max_iterations=20
    )

    print(f"\nFinal SP Score:  {final_score}")
    print(f"Improvement:     {final_score - initial_score:+d}")
    print(f"Score history:   {history}")
    print("\nFinal MSA:")
    for seq in refined_msa:
        print(f"  {seq}")
