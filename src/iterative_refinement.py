# =============================================================================
# MUSCLE-style Iterative Refinement for Multiple Sequence Alignment (MSA)
#
# Core idea: Repeatedly split MSA into two groups, re-align their gap-free
# representatives via NW, propagate new gaps to all sequences, accept if
# Sum-of-Pairs (SP) score improves. Converges when no improvement found.
#
# FUNCTIONS:
# - sum_of_pairs(msa)     → SP score (+1 match, -1 mismatch, -2 gap)
# - apply_group_gaps()    → transfers new gap pattern to group members
# - refine_once(msa, nw)  → tries all splits, returns best improvement
# - iterative_refinement()→ convergence loop around refine_once()
#
# COMPLEXITY: O(iterations * n² * L) where n=sequences, L=alignment length
# =============================================================================

def sum_of_pairs(msa):
    """
    Calculate Sum-of-Pairs (SP) score for an MSA.

    For every PAIR of sequences (i, j), for every COLUMN:
      +1 if both have the same non-gap character (match)
      -1 if they have different non-gap characters (mismatch)
      -2 if either has a gap (gap penalty)

    A higher SP score means a better alignment.
    """
    total     = 0
    num_seqs  = len(msa)
    aln_len   = len(msa[0]) if msa else 0

    for i in range(num_seqs):
        for j in range(i + 1, num_seqs):
            for col in range(aln_len):
                c1 = msa[i][col]
                c2 = msa[j][col]
                if c1 != '-' and c2 != '-':
                    total += 1 if c1 == c2 else -1
                else:
                    total += -2   # gap penalty

    return total


def apply_group_gaps(group, old_rep, new_rep):
    """
    Propagate a new gap pattern from a re-aligned representative
    to every sequence in the group.

    Parameters
    ----------
    group   : list of str – the original aligned sequences
    old_rep : str         – the representative WITHOUT gaps (was used for re-alignment)
    new_rep : str         – the representative WITH NEW gaps (output of NW)

    Returns
    -------
    list of str – all sequences with the new gap pattern applied
    """
    new_group = []
    for seq in group:
        new_seq = ""
        idx     = 0   # position in the un-gapped sequence
        for char in new_rep:
            if char == '-':
                new_seq += '-'        # insert a gap where new_rep has a gap
            else:
                # Take the character from the original aligned sequence
                new_seq += seq[idx] if idx < len(seq) else '-'
                idx += 1
        new_group.append(new_seq)
    return new_group


def refine_once(msa, nw_fn):
    """
    One round of partition-realign-accept.
    Tries every possible split of the MSA into two groups,
    realigns group representatives, and keeps the improvement if SP increases.

    Returns
    -------
    best_msa   : list of str – (possibly) improved MSA
    best_score : int         – SP score of the best MSA found
    """
    best_score = sum_of_pairs(msa)
    best_msa   = list(msa)
    n          = len(msa)

    for pivot in range(1, n):
        group_a = msa[:pivot]
        group_b = msa[pivot:]

        # Use first sequence of each group as representative (remove its gaps)
        rep_a = group_a[0].replace('-', '')
        rep_b = group_b[0].replace('-', '')

        # Re-align the two representatives
        a1, a2, _ = nw_fn(rep_a, rep_b)

        # Propagate the new gap pattern to all sequences in each group
        new_msa = (apply_group_gaps(group_a, rep_a, a1) +
                   apply_group_gaps(group_b, rep_b, a2))

        new_score = sum_of_pairs(new_msa)

        if new_score > best_score:
            best_score = new_score
            best_msa   = new_msa

    return best_msa, best_score


# =============================================================================
# iterative_refinement() — wraps refine_once() in a convergence loop
# =============================================================================

def iterative_refinement(msa, nw_fn, max_iterations=50):
    """
    Run iterative refinement until the SP score stops improving.

    This is the MUSCLE-style approach:
      1. Start with an initial MSA (e.g., from progressive alignment)
      2. Try to improve it by splitting into groups and re-aligning
      3. Accept improvements; repeat
      4. Stop when no improvement is found in a full round

    Parameters
    ----------
    msa            : list of str – the initial MSA to refine
    nw_fn          : callable    – your needleman_wunsch function
    max_iterations : int         – safety limit (stops even if still improving)
                                   prevents infinite loops on tricky inputs

    Returns
    -------
    best_msa    : list of str – the refined MSA
    final_score : int         – SP score of the refined MSA
    history     : list of int – SP score after each iteration (for plotting)
    """
    current_msa   = list(msa)
    current_score = sum_of_pairs(current_msa)
    history       = [current_score]   # track score per iteration

    print(f"  Iterative refinement starting. Initial SP score: {current_score}")

    for iteration in range(1, max_iterations + 1):

        # Run one round of refinement
        new_msa, new_score = refine_once(current_msa, nw_fn)

        print(f"  Round {iteration:3d}: SP score = {new_score}"
              f"  ({'improved' if new_score > current_score else 'no change'})")

        history.append(new_score)

        if new_score <= current_score:
            # No improvement this round → CONVERGENCE — stop the loop
            print(f"  Converged after {iteration} round(s).")
            break

        # Accept the improvement and continue
        current_msa   = new_msa
        current_score = new_score

    else:
        # Loop finished without break → hit max_iterations
        print(f"  WARNING: Reached max_iterations={max_iterations}. "
              f"May not have fully converged.")

    return current_msa, current_score, history


# =============================================================================
# DEMO
# =============================================================================

if __name__ == "__main__":
    try:
        from needlemanWunschGlobal import needleman_wunsch
    except ImportError:
        # Minimal fallback NW for demo
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
                if i>0 and j>0 and dp[i][j]==dp[i-1][j-1]+(match if s1[i-1]==s2[j-1] else mismatch):
                    a1+=s1[i-1]; a2+=s2[j-1]; i-=1; j-=1
                elif i>0 and dp[i][j]==dp[i-1][j]+gap:
                    a1+=s1[i-1]; a2+='-'; i-=1
                else:
                    a1+='-'; a2+=s2[j-1]; j-=1
            return a1[::-1], a2[::-1], dp[n][m]

    print("=" * 60)
    print("ITERATIVE REFINEMENT DEMO")
    print("=" * 60)

    # Initial MSA from progressive alignment output
    initial_msa = [
        "-VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFP-T-TK-TYFPH---FDLSHGSAQVK",
        "-VLSAADKTNVKAAWSKVGGHAGEYGAEALERMFLGFP-T-TK-TYFPH---FDLSHGSAQVK",
        "-VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPW-T-QR-FFESF---GDLSTPDAVMG",
    ]

    initial_score = sum_of_pairs(initial_msa)
    print(f"\nInitial SP Score: {initial_score}")

    # Run full iterative refinement
    refined_msa, final_score, history = iterative_refinement(
        initial_msa, needleman_wunsch, max_iterations=20
    )

    print(f"\nFinal SP Score: {final_score}")
    print(f"Improvement: {final_score - initial_score:+d} points")
    print(f"\nScore history: {history}")

    print(f"\nFinal MSA:")
    for seq in refined_msa:
        print(f"  {seq}")
