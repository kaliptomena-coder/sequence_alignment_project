# A Profile HMM captures the statistical properties of a protein/DNA family
# from a multiple sequence alignment (MSA). It has three types of states:
#
#   M (Match)     — emits one character from a learned distribution
#   I (Insertion) — emits one character from a background distribution
#   D (Deletion)  — SILENT state (no emission; represents a skipped column)
#
# Building the model:
#   - Columns with fewer than gap_threshold fraction of gaps → Match columns
#   - Characters in Match columns → train M state emission probabilities
#   - Characters in Insertion columns → train I state emission probabilities
#   - Transitions counted from each sequence's path through the columns
#
# Laplace pseudo-count smoothing:
#   Adding 1 to every character's count before normalizing ensures every
#   emission probability is > 0. Without this, unseen characters would give
#   log(0) = -inf and break the Viterbi score.
#
# Viterbi decoding (log-space):
#   Finding the single most probable state path for a query sequence.
#   Log-space arithmetic avoids numerical underflow from multiplying many
#   small probabilities together.

import math
from collections import defaultdict

# Alphabet choices — change to PROTEIN_ALPHABET for amino acid sequences
PROTEIN_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
DNA_ALPHABET      = "ACGT"

class ProfileHMM:
    """
    Plan-7-inspired Profile HMM for protein or DNA sequence families.

    Usage:
        msa   = ["ACGT--ACGT", "ACGTAAACGT", "ACGT--ACGG"]
        model = ProfileHMM(msa, gap_threshold=0.5, alphabet=DNA_ALPHABET)
        path, score = model.viterbi("ACGTACGT")
    """

    def __init__(self, msa, gap_threshold=0.5, alphabet=DNA_ALPHABET):
        """
        Building the Profile HMM from a multiple sequence alignment.

        Parameters
        ----------
        msa           : list of str  - aligned sequences (all the same length)
        gap_threshold : float        - columns with >= this fraction of gaps
                                       become insertion columns
        alphabet      : str          - characters to model (e.g., "ACGT")
        """
        self.msa           = msa
        self.nseq          = len(msa)
        self.length        = len(msa[0]) if msa else 0
        self.gap_threshold = gap_threshold
        self.alphabet      = alphabet

        self.match_cols  = []
        self.states      = []
        self.emissions   = {}
        self.transitions = {}

        if self.nseq > 0:
            self._build_model()

    # Model Building

    def _build_model(self):
        """Computing emission and transition probability tables from the MSA."""

        # Step 1: Identifying which columns are Match vs Insertion columns
        for col in range(self.length):
            gap_frac = sum(1 for seq in self.msa if seq[col] == '-') / self.nseq
            if gap_frac < self.gap_threshold:
                self.match_cols.append(col)

        k = len(self.match_cols)  # total number of match states

        # Step 2: Naming all states
        self.states = ["BEGIN"]
        for idx in range(k + 1):
            self.states.append(f"I{idx}")
        for idx in range(1, k + 1):
            self.states.append(f"M{idx}")
            self.states.append(f"D{idx}")
        self.states.append("END")

        # Step 3: Counting emissions and transitions by tracing each sequence
        emit_counts  = defaultdict(lambda: defaultdict(int))
        trans_counts = defaultdict(lambda: defaultdict(int))

        match_col_set    = set(self.match_cols)
        col_to_match_idx = {col: idx + 1 for idx, col in enumerate(self.match_cols)}

        for seq in self.msa:
            prev_state = "BEGIN"
            match_idx  = 0

            for col in range(self.length):
                symbol = seq[col]

                if col in match_col_set:
                    match_idx = col_to_match_idx[col]
                    if symbol == '-':
                        curr_state = f"D{match_idx}"  # gap in a match column → deletion
                    else:
                        curr_state = f"M{match_idx}"
                        emit_counts[curr_state][symbol] += 1
                else:
                    if symbol != '-':
                        curr_state = f"I{match_idx}"
                        emit_counts[curr_state][symbol] += 1
                    else:
                        continue  # gap in insertion column → skip

                trans_counts[prev_state][curr_state] += 1
                prev_state = curr_state

            trans_counts[prev_state]["END"] += 1

        # Step 4: Normalizing emission counts to probabilities with Laplace smoothing
        # Adding 1 to every alphabet character ensures no probability is exactly 0
        self.emissions = {}
        for state in self.states:
            if state.startswith('D') or state in ('BEGIN', 'END'):
                self.emissions[state] = {}  # silent states have no emissions
                continue

            counts = dict(emit_counts[state])
            for char in self.alphabet:
                counts[char] = counts.get(char, 0) + 1  # Laplace pseudo-count

            total = sum(counts.values())
            self.emissions[state] = {c: count / total for c, count in counts.items()}

        # Step 5: Normalizing transition counts (small pseudo-count to prevent zeros)
        self.transitions = {}
        for state, dests in trans_counts.items():
            counts = {dest: count + 0.1 for dest, count in dests.items()}
            total  = sum(counts.values())
            self.transitions[state] = {dest: count / total for dest, count in counts.items()}

    # Viterbi Decoding

    def viterbi(self, sequence):
        """
        Finding the most probable state path through the model for a query sequence.

        Uses log-space arithmetic to avoid numerical underflow:
        Instead of multiplying many small probabilities (which rounds to 0),
        we add their logarithms (which remain finite).

        Parameters
        ----------
        sequence : str  - the query sequence to decode

        Returns
        -------
        best_path  : list of str  - the sequence of state names
        best_score : float        - log-probability of the best path
        """
        T = len(sequence)

        # viterbi[t][state] = best log-prob path to `state` after emitting t characters
        viterbi = [{}]
        backptr  = [{}]

        for state in self.states:
            viterbi[0][state] = -math.inf
            backptr[0][state]  = None

        viterbi[0]["BEGIN"] = 0.0  # always start in BEGIN with log-prob = 0

        # Forward pass
        for t in range(1, T + 1):
            viterbi.append({})
            backptr.append({})
            obs = sequence[t - 1]

            for curr_state in self.states:
                # Silent states don't consume characters; log-emission cost = 0
                if curr_state.startswith('D') or curr_state in ('BEGIN', 'END'):
                    log_emit = 0.0
                else:
                    prob     = self.emissions.get(curr_state, {}).get(obs, 1e-9)
                    log_emit = math.log(prob)

                best_log_prob = -math.inf
                best_prev     = None

                for prev_state in self.states:
                    prev_score = viterbi[t-1].get(prev_state, -math.inf)
                    if prev_score == -math.inf:
                        continue
                    trans_prob = self.transitions.get(prev_state, {}).get(curr_state, 0)
                    if trans_prob <= 0:
                        continue
                    candidate = prev_score + math.log(trans_prob) + log_emit
                    if candidate > best_log_prob:
                        best_log_prob = candidate
                        best_prev     = prev_state

                viterbi[t][curr_state] = best_log_prob
                backptr[t][curr_state] = best_prev

        # Termination: find the highest-scoring final state
        best_final_score = -math.inf
        best_final_state = None
        for state in self.states:
            s = viterbi[T].get(state, -math.inf)
            if s > best_final_score:
                best_final_score = s
                best_final_state = state

        # Traceback: follow backpointers to reconstruct the path
        path, state = [], best_final_state
        for t in range(T, 0, -1):
            path.append(state)
            state = backptr[t].get(state)
            if state is None:
                break

        path.reverse()
        return path, best_final_score


if __name__ == "__main__":
    print("=== Profile HMM Demo ===")

    training_msa = [
        "ACGT--ACGT",
        "ACGTAAACGT",
        "ACGT--ACGG",
    ]

    print("\nTraining MSA:")
    for s in training_msa:
        print(f"  {s}")

    model = ProfileHMM(training_msa, gap_threshold=0.5, alphabet=DNA_ALPHABET)

    print(f"\nMatch columns: {model.match_cols}")
    print(f"\nM1 emission probabilities:")
    if "M1" in model.emissions:
        for char, prob in sorted(model.emissions["M1"].items()):
            print(f"  P({char} | M1) = {prob:.4f}")

    query = "ACGTACGT"
    print(f"\nDecoding: {query}")
    path, score = model.viterbi(query)
    print(f"Viterbi path:  {' -> '.join(path)}")
    print(f"Viterbi score: {score:.4f}")

    if math.isinf(score):
        print("\nERROR: Score is -inf — pseudo-counts may not be applied correctly.")
    else:
        print("\nSUCCESS: Score is finite (pseudo-count smoothing is working).")
