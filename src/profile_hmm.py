import math
from collections import defaultdict

class ProfileHMM:
    def __init__(self, msa, gap_threshold=0.5):
        self.msa = msa
        self.nseq = len(msa)
        self.length = len(msa[0])
        self.gap_threshold = gap_threshold

        self.match_cols = []
        self.states = []
        self.emissions = {}
        self.transitions = {}

        self._build_model()

    def _build_model(self):
        # 1. Determine match columns
        for col in range(self.length):
            gaps = sum(1 for seq in self.msa if seq[col] == '-')
            if gaps / self.nseq < self.gap_threshold:
                self.match_cols.append(col)

        k = len(self.match_cols)

        # States: M1..Mk, I0..Ik, D1..Dk
        for i in range(k + 1):
            self.states.append(f"I{i}")
        for i in range(1, k + 1):
            self.states.append(f"M{i}")
            self.states.append(f"D{i}")

        # Initialize counts
        emit_counts = defaultdict(lambda: defaultdict(int))
        trans_counts = defaultdict(lambda: defaultdict(int))

        # 2. Count emissions and transitions
        for seq in self.msa:
            prev_state = "M0"

            match_idx = 0

            for col in range(self.length):
                symbol = seq[col]

                if col in self.match_cols:
                    match_idx += 1

                    if symbol == '-':
                        state = f"D{match_idx}"
                    else:
                        state = f"M{match_idx}"
                        emit_counts[state][symbol] += 1
                else:
                    if symbol != '-':
                        state = f"I{match_idx}"
                        emit_counts[state][symbol] += 1
                    else:
                        continue  # skip pure gap insertions

                trans_counts[prev_state][state] += 1
                prev_state = state

        # 3. Normalize emissions
        self.emissions = {}
        for state, counts in emit_counts.items():
            total = sum(counts.values())
            self.emissions[state] = {
                k: v / total for k, v in counts.items()
            }

        # 4. Normalize transitions
        self.transitions = {}
        for s1, dests in trans_counts.items():
            total = sum(dests.values())
            self.transitions[s1] = {
                s2: v / total for s2, v in dests.items()
            }

    # ---------------------------
    # VITERBI DECODING
    # ---------------------------
    def viterbi(self, sequence):
        V = [{}]
        path = {}

        # init
        for state in self.states:
            V[0][state] = -math.inf
            path[state] = []

        V[0]["M0"] = 0

        # dynamic programming
        for t in range(len(sequence)):
            V.append({})
            new_path = {}

            for curr in self.states:
                max_prob = -math.inf
                best_prev = None

                for prev in self.transitions:
                    if curr in self.transitions.get(prev, {}):
                        trans_p = math.log(self.transitions[prev][curr])
                    else:
                        continue

                    emit_p = 0
                    if curr in self.emissions:
                        emit_p = math.log(
                            self.emissions[curr].get(sequence[t], 1e-6)
                        )

                    prob = V[t].get(prev, -math.inf) + trans_p + emit_p

                    if prob > max_prob:
                        max_prob = prob
                        best_prev = prev

                V[t + 1][curr] = max_prob
                new_path[curr] = path.get(best_prev, []) + [curr]

            path = new_path

        # termination
        final_state = max(V[-1], key=V[-1].get)
        return path[final_state], V[-1][final_state]