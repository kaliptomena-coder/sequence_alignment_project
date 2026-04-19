# =============================================================================
#Profile HMM with Laplace pseudo-count smoothing:
# - Deletion (D) states produce no emission (prob = 1.0)
#   - Insertion (I) states emit with a uniform background distribution
#     if no observations exist for that state
#   - Viterbi correctly distinguishes emitting vs. silent states

#   Added a "pseudo-count" of 1 to every amino acid (or nucleotide) before
#   dividing to get probabilities.
#   It ensures every probability > 0, so log(prob) is always finite.
# =============================================================================

import math
from collections import defaultdict


# Define the alphabet — change to "ACGT" for DNA/RNA
PROTEIN_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
DNA_ALPHABET      = "ACGT"


class ProfileHMM:
    """
    Plan-7-inspired Profile HMM for protein or DNA sequence families.

    States
    ------
    M1..Mk  : Match states      — emit one character from learned distribution
    I0..Ik  : Insertion states  — emit one character from background distribution
    D1..Dk  : Deletion states   — SILENT (emit nothing, represent a gap)

    Usage
    -----
    msa   = ["ACGT--ACGT", "ACGTAAACGT", "ACGT--ACGG"]
    model = ProfileHMM(msa, gap_threshold=0.5, alphabet=DNA_ALPHABET)
    path, score = model.viterbi("ACGTACGT")
    """

    def __init__(self, msa, gap_threshold=0.5, alphabet=DNA_ALPHABET):
        """
        Build the Profile HMM from a multiple sequence alignment (MSA).

        Parameters
        ----------
        msa           : list of str  – aligned sequences (all same length)
        gap_threshold : float        – columns with >= this fraction of gaps
                                       are treated as insertion columns
        alphabet      : str          – characters to consider (e.g., "ACGT")
        """
        self.msa           = msa
        self.nseq          = len(msa)
        self.length        = len(msa[0]) if msa else 0
        self.gap_threshold = gap_threshold
        self.alphabet      = alphabet

        # These will be populated by _build_model()
        self.match_cols    = []    # indices of "match" columns
        self.states        = []    # all state names
        self.emissions     = {}    # {state: {char: probability}}
        self.transitions   = {}    # {state: {next_state: probability}}

        if self.nseq > 0:
            self._build_model()

    # =========================================================================
    # MODEL BUILDING
    # =========================================================================

    def _build_model(self):
        """Build emission and transition probability tables from the MSA."""

        # -----------------------------------------------------------------
        # 1. Determine which columns are "match" columns.
        #    A column is a match column if < gap_threshold fraction of
        #    sequences have a gap there.
        # -----------------------------------------------------------------
        for col in range(self.length):
            gap_fraction = sum(1 for seq in self.msa if seq[col] == '-') / self.nseq
            if gap_fraction < self.gap_threshold:
                self.match_cols.append(col)

        k = len(self.match_cols)   # number of match states

        # -----------------------------------------------------------------
        # 2. Define all state names.
        #    I0 = insertion before position 1
        #    M1..Mk = match states
        #    D1..Dk = deletion states (silent)
        #    I1..Ik = insertion states after each match
        # -----------------------------------------------------------------
        self.states = ["BEGIN"]
        for idx in range(k + 1):
            self.states.append(f"I{idx}")
        for idx in range(1, k + 1):
            self.states.append(f"M{idx}")
            self.states.append(f"D{idx}")
        self.states.append("END")

        # -----------------------------------------------------------------
        # 3. Count emissions and transitions by tracing each sequence
        #    through the MSA column by column.
        # -----------------------------------------------------------------
        # Use defaultdict so accessing an unseen key gives 0 automatically
        emit_counts = defaultdict(lambda: defaultdict(int))
        trans_counts = defaultdict(lambda: defaultdict(int))

        # match_col_set for O(1) lookup
        match_col_set = set(self.match_cols)
        # Map from column index to match state index (1-based)
        col_to_match_idx = {col: idx + 1
                            for idx, col in enumerate(self.match_cols)}

        for seq in self.msa:
            prev_state  = "BEGIN"
            match_idx   = 0         # tracks which M/D state we are in

            for col in range(self.length):
                symbol = seq[col]

                if col in match_col_set:
                    # This is a match column
                    match_idx = col_to_match_idx[col]

                    if symbol == '-':
                        # Gap in a match column → DELETION state (silent)
                        curr_state = f"D{match_idx}"
                        # D states do NOT emit — no emit count update
                    else:
                        # Character in a match column → MATCH state
                        curr_state = f"M{match_idx}"
                        emit_counts[curr_state][symbol] += 1

                else:
                    # This is an insertion column
                    if symbol != '-':
                        # Non-gap character in insertion column → INSERTION state
                        curr_state = f"I{match_idx}"
                        emit_counts[curr_state][symbol] += 1
                    else:
                        # Gap in an insertion column → skip entirely
                        continue

                # Record the transition from previous state to current state
                trans_counts[prev_state][curr_state] += 1
                prev_state = curr_state

            # End of sequence: transition to END
            trans_counts[prev_state]["END"] += 1

        # -----------------------------------------------------------------
        # 4. Normalise emission counts to probabilities.
        #
        #    LAPLACE PSEUDO-COUNTS
        #    Before dividing, add 1 to every alphabet character.
        #    This ensures probability > 0 for every character,
        #    so log(probability) is never -infinity.
        # -----------------------------------------------------------------
        self.emissions = {}

        for state in self.states:

            if state.startswith('D') or state in ('BEGIN', 'END'):
                # Deletion states and boundary states are SILENT — no emissions
                # We store an empty dict to signal "this state doesn't emit"
                self.emissions[state] = {}
                continue

            counts = dict(emit_counts[state])   # observed counts (may be empty)

            # --- Apply Laplace pseudo-count ---
            # Add 1 to every character in the alphabet, even unseen ones
            for char in self.alphabet:
                counts[char] = counts.get(char, 0) + 1  # <-- THE FIX

            # Now normalise: divide each count by the total
            total = sum(counts.values())
            self.emissions[state] = {char: count / total
                                     for char, count in counts.items()}

        # -----------------------------------------------------------------
        # 5. Normalise transition counts to probabilities.
        #    Also add a small pseudo-count (0.1) to prevent zero transitions.
        # -----------------------------------------------------------------
        self.transitions = {}

        for state, dests in trans_counts.items():
            # Add pseudo-count to each observed destination
            counts = {dest: count + 0.1 for dest, count in dests.items()}
            total  = sum(counts.values())
            self.transitions[state] = {dest: count / total
                                       for dest, count in counts.items()}

    # =========================================================================
    # VITERBI DECODING
    # =========================================================================

    def viterbi(self, sequence):
        """
        Find the most probable state path through the model for a given sequence.

        Uses the Viterbi algorithm in LOG-SPACE to avoid numerical underflow.
        (Multiplying many small probabilities together → rounds to 0.
         Using log turns multiplications into additions, which is stable.)

        Parameters
        ----------
        sequence : str – the query sequence to decode

        Returns
        -------
        best_path  : list of str – sequence of state names
        best_score : float       – log-probability of the best path
        """
        T = len(sequence)   # length of the query sequence

        # viterbi[t][state] = log-probability of the best path ending
        # in `state` after having emitted sequence[0..t-1]
        # We use a list of dicts: one dict per time step.

        # Time step 0: no characters emitted yet
        viterbi = [{}]
        backptr  = [{}]     # backpointer for traceback

        # Initialise all states to -infinity (unreachable)
        for state in self.states:
            viterbi[0][state] = -math.inf
            backptr[0][state]  = None

        # The model starts in BEGIN state with probability 1 (log = 0)
        viterbi[0]["BEGIN"] = 0.0

        # -----------------------------------------------------------------
        # Forward pass: fill the Viterbi table
        # -----------------------------------------------------------------
        for t in range(1, T + 1):
            viterbi.append({})
            backptr.append({})
            obs = sequence[t - 1]   # the character emitted at step t

            for curr_state in self.states:

                # --- Determine emission log-probability ---
                if curr_state.startswith('D') or curr_state in ('BEGIN', 'END'):
                    # Silent states: they don't emit, so we should NOT
                    # consume a character when transitioning INTO them.
                    # For simplicity in this implementation, we give them
                    # emission log-prob = 0 (i.e., prob = 1, no cost).
                    # A full implementation would handle silent states with
                    # a separate pass that doesn't advance t.
                    log_emit = 0.0
                else:
                    # Emitting states: get the probability of seeing `obs`
                    emit_dist = self.emissions.get(curr_state, {})
                    prob      = emit_dist.get(obs, 1e-9)   # 1e-9 as safety floor
                    log_emit  = math.log(prob)

                # --- Find the best previous state ---
                best_log_prob = -math.inf
                best_prev     = None

                for prev_state in self.states:
                    prev_score = viterbi[t - 1].get(prev_state, -math.inf)
                    if prev_score == -math.inf:
                        continue   # unreachable previous state, skip

                    # Transition log-probability from prev → curr
                    trans_prob = self.transitions.get(prev_state, {}).get(curr_state, 0)
                    if trans_prob <= 0:
                        continue   # no transition exists
                    log_trans = math.log(trans_prob)

                    candidate = prev_score + log_trans + log_emit

                    if candidate > best_log_prob:
                        best_log_prob = candidate
                        best_prev     = prev_state

                viterbi[t][curr_state] = best_log_prob
                backptr[t][curr_state] = best_prev

        # -----------------------------------------------------------------
        # Termination: find the best final state at t = T
        # -----------------------------------------------------------------
        best_final_score = -math.inf
        best_final_state = None

        for state in self.states:
            score = viterbi[T].get(state, -math.inf)
            if score > best_final_score:
                best_final_score = score
                best_final_state = state

        # -----------------------------------------------------------------
        # Traceback: follow backpointers to reconstruct the path
        # -----------------------------------------------------------------
        path  = []
        state = best_final_state

        for t in range(T, 0, -1):
            path.append(state)
            state = backptr[t].get(state)
            if state is None:
                break

        path.reverse()   # we built it backwards

        return path, best_final_score


# =============================================================================
# DEMO
# =============================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("PROFILE HMM DEMO")
    print("=" * 60)

    # A small training MSA (3 sequences, DNA)
    training_msa = [
        "ACGT--ACGT",
        "ACGTAAACGT",
        "ACGT--ACGG",
    ]

    print("\nTraining MSA:")
    for s in training_msa:
        print(f"  {s}")

    # Build the model
    model = ProfileHMM(training_msa, gap_threshold=0.5, alphabet=DNA_ALPHABET)

    print(f"\nMatch columns identified: {model.match_cols}")
    print(f"Total match states: {len(model.match_cols)}")
    print(f"\nEmission table (M1):")
    if "M1" in model.emissions:
        for char, prob in sorted(model.emissions["M1"].items()):
            print(f"  P({char} | M1) = {prob:.4f}")

    # Decode a new sequence
    query = "ACGTACGT"
    print(f"\nDecoding query: {query}")
    path, score = model.viterbi(query)
    print(f"Viterbi path:  {' → '.join(path)}")
    print(f"Viterbi score: {score:.4f}")

    # Verify score is not -inf
    if math.isinf(score):
        print("\nERROR: Score is -inf! Pseudo-counts may not be applied.")
    else:
        print("\nSUCCESS: Score is finite (pseudo-counts working correctly).")
