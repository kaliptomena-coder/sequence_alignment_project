# Sequence Alignment — Unit Tests
#
# Coverage: all 10 required algorithms + UPGMA, progressive MSA,
#           real FASTA identity stats, and HBB long-sequence tests.
#
# Run:
#   python test_all.py          (built-in runner, no dependencies)
#   pytest test_all.py -v       (pytest)

import sys
import os
import math
import itertools

src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src'))
if src_path not in sys.path:
    sys.path.insert(0, src_path)

data_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'data'))

from needlemanWunschGlobal import needleman_wunsch
from smithWatermanLocal   import smith_waterman
from gotoh                import gotoh
from hirschberg           import hirschberg
from banded_alignment     import banded_nw
from blast_lite           import blast_lite
from minimizer_align      import minimizer_align
from iterative_refinement import sum_of_pairs, refine_once
from profile_hmm          import ProfileHMM
from progressive_msa      import perform_progressive_alignment
from distance_matrix      import generate_matrix
from upgma                import run_upgma


# ── Shared fixtures ───────────────────────────────────────────────────────────

SMALL_MSA = [
    "ACGT--ACGT",
    "ACGTAAACGT",
    "ACGT--ACGG",
]


def _score_from_strings(a1, a2, match=1, mismatch=-1, gap=-2):
    """Recompute alignment score directly from two aligned strings."""
    score = 0
    for c1, c2 in zip(a1, a2):
        if c1 == '-' or c2 == '-':
            score += gap
        elif c1 == c2:
            score += match
        else:
            score += mismatch
    return score


def _load_fasta(filename):
    """Parse a FASTA file from the data/ directory; return {name: seq}."""
    path = os.path.join(data_path, filename)
    if not os.path.exists(path):
        return {}
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name:
                    seqs[name] = ''.join(buf)
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line.upper())
    if name:
        seqs[name] = ''.join(buf)
    return seqs


def _pairwise_identities(seqs):
    """Return list of NW pairwise identity % for all pairs in a dict of sequences."""
    ids = []
    for a, b in itertools.combinations(seqs.keys(), 2):
        a1, a2, _ = needleman_wunsch(seqs[a], seqs[b])
        matches = sum(c1 == c2 and c1 != '-' for c1, c2 in zip(a1, a2))
        cols    = sum(1 for c1, c2 in zip(a1, a2) if c1 != '-' and c2 != '-')
        ids.append(100.0 * matches / cols if cols else 0.0)
    return ids


# =============================================================================
# 1. Needleman-Wunsch (Global Alignment)
# =============================================================================

class TestNeedlemanWunsch:

    def test_identical_sequences(self):
        """Identical input → no gaps, positive score, aligned strings equal input."""
        a1, a2, score = needleman_wunsch("ACGT", "ACGT")
        assert a1 == "ACGT" and a2 == "ACGT"
        assert score > 0

    def test_aligned_strings_equal_length(self):
        """Fundamental invariant: both aligned strings must always be the same length."""
        for s1, s2 in [("ACGT", "AGT"), ("GATTACA", "GATCA"), ("A", "ACGT")]:
            a1, a2, _ = needleman_wunsch(s1, s2)
            assert len(a1) == len(a2), f"Unequal lengths for '{s1}' vs '{s2}'"

    def test_gap_introduced_in_classic_example(self):
        """GATTACA vs GATCA: optimal alignment requires at least one gap."""
        a1, a2, score = needleman_wunsch("GATTACA", "GATCA")
        assert '-' in a1 or '-' in a2
        assert score > 0

    def test_harsher_penalty_cannot_improve_score(self):
        """Increasing gap cost must not raise the alignment score."""
        _, _, lenient = needleman_wunsch("GATTACA", "GATCA", gap=-1)
        _, _, harsh   = needleman_wunsch("GATTACA", "GATCA", gap=-10)
        assert lenient >= harsh

    def test_symmetry(self):
        """Swapping the two sequences must produce the same optimal score."""
        _, _, ab = needleman_wunsch("GATTACA", "GATCA")
        _, _, ba = needleman_wunsch("GATCA",   "GATTACA")
        assert ab == ba

    def test_dp_score_matches_traceback_score(self):
        """Score from the DP cell must equal the score recomputed from aligned strings."""
        a1, a2, dp_score = needleman_wunsch("GCATGCU", "GATTACA")
        assert dp_score == _score_from_strings(a1, a2)


# =============================================================================
# 2. Smith-Waterman (Local Alignment)
# =============================================================================

class TestSmithWaterman:

    def test_score_never_negative(self):
        """SW floors cells at zero; score must always be ≥ 0."""
        _, _, score = smith_waterman("AAAA", "TTTT")
        assert score >= 0

    def test_finds_embedded_motif_better_than_nw(self):
        """SW must outscore NW when a shared motif is flanked by dissimilar sequence."""
        _, _, sw = smith_waterman("AAAACGTAAAA", "TTTTCGTTTTT")
        _, _, nw = needleman_wunsch("AAAACGTAAAA", "TTTTCGTTTTT")
        assert sw > nw, "SW should isolate the CGT core; NW penalises mismatched flanks"

    def test_no_common_subsequence_returns_zero(self):
        """With severe penalties and no shared bases, SW must return score = 0."""
        _, _, score = smith_waterman("AAAA", "TTTT", match=1, mismatch=-10, gap=-10)
        assert score == 0


# =============================================================================
# 3. Gotoh (Affine Gap Penalties)
# =============================================================================

class TestGotoh:

    def test_identical_sequences_produce_no_gaps(self):
        a1, a2, _ = gotoh("ACGT", "ACGT")
        assert '-' not in a1 and '-' not in a2

    def test_deletion_introduces_gap(self):
        """One deleted character must produce at least one gap in the alignment."""
        a1, a2, _ = gotoh("GATTACA", "GATACA",
                          match=2, mismatch=-1, gap_open=-5, gap_extend=-1)
        assert len(a1) == len(a2)
        assert '-' in a1 or '-' in a2

    def test_higher_open_penalty_lowers_score(self):
        """A more expensive gap-open cost must not improve the alignment score."""
        s1, s2 = "GATTACAACTTG", "GATCCAGTTCAAA"
        _, _, cheap     = gotoh(s1, s2, gap_open=-3,  gap_extend=-1)
        _, _, expensive = gotoh(s1, s2, gap_open=-10, gap_extend=-1)
        assert cheap >= expensive


# =============================================================================
# 4. Hirschberg (Space-Efficient Global Alignment)
# =============================================================================

class TestHirschberg:

    PAIRS = [
        ("AGTAACG", "ACATAG"),
        ("GATTACA", "GATCA"),
        ("ACGT",    "ACGT"),
        ("AAAA",    "TTTT"),
    ]

    def test_score_matches_needleman_wunsch(self):
        """Hirschberg must produce the same optimal score as NW (space-efficiency only)."""
        for s1, s2 in self.PAIRS:
            hb_a1, hb_a2, _ = hirschberg(s1, s2, match=2, mismatch=-1, gap=-1)
            _, _, nw_score   = needleman_wunsch(s1, s2, match=2, mismatch=-1, gap=-1)
            hb_score = _score_from_strings(hb_a1, hb_a2, match=2, mismatch=-1, gap=-1)
            assert hb_score == nw_score, \
                f"Hirschberg {hb_score} ≠ NW {nw_score} for '{s1}'/'{s2}'"

    def test_empty_sequence_handled(self):
        a1, a2, _ = hirschberg("", "ACG")
        assert len(a1) == len(a2)

    def test_output_only_valid_characters(self):
        """Aligned strings may only contain A/C/G/T and dashes."""
        a1, a2, _ = hirschberg("GATTACA", "GATCACA")
        valid = set("ACGT-")
        assert len(a1) == len(a2)
        assert all(c in valid for c in a1 + a2)


# =============================================================================
# 5. Banded Alignment
# =============================================================================

class TestBandedAlignment:

    def test_score_matches_full_nw_for_identical_sequences(self):
        """For identical sequences, banded NW must return the same score as full NW."""
        s = "GATTACAACTTG"
        _, _, banded = banded_nw(s, s, k=5, match=2, mismatch=-1, gap=-2)
        _, _, full   = needleman_wunsch(s, s, match=2, mismatch=-1, gap=-2)
        assert banded == full

    def test_narrow_band_falls_back_without_crash(self):
        """When band is too narrow the algorithm must fall back gracefully."""
        a1, a2, _ = banded_nw("AAAAAAAAAAA", "TTTTTT", k=1)
        assert a1 is not None and len(a1) == len(a2)


# =============================================================================
# 6. BLAST-lite (Heuristic Seed-and-Extend)
# =============================================================================

class TestBLASTLite:

    def test_identical_sequences_find_at_least_one_hsp(self):
        assert len(blast_lite("GATTACA", "GATTACA", k=3)) > 0

    def test_hsp_has_required_fields(self):
        for r in blast_lite("GGAGTCAG", "GAAGTCGG", k=3):
            assert {'score', 'alignment', 'query_pos', 'target_pos'} <= r.keys()

    def test_no_shared_kmer_returns_empty_list(self):
        assert blast_lite("AAA", "TTT", k=3) == []

    def test_results_sorted_highest_score_first(self):
        results = blast_lite("GATTACAGATTACA", "GATTACAGATTACA", k=3)
        if len(results) > 1:
            scores = [r['score'] for r in results]
            assert scores == sorted(scores, reverse=True)


# =============================================================================
# 7. Minimizer-Based Alignment
# =============================================================================

class TestMinimizerAlign:

    def test_identical_sequences_produce_anchors(self):
        _, _, _, chain = minimizer_align("GATTACAGATTACA", "GATTACAGATTACA", k=3, w=5)
        assert len(chain) > 0

    def test_each_anchor_is_integer_triple(self):
        _, _, _, chain = minimizer_align("GATTACAGATTACA", "GATTACAGATTACA", k=3, w=5)
        for anchor in chain:
            assert len(anchor) == 3 and all(isinstance(x, int) for x in anchor)

    def test_dissimilar_sequences_produce_fewer_anchors(self):
        _, _, _, same = minimizer_align("GATTACAGATTACA", "GATTACAGATTACA", k=3, w=5)
        _, _, _, diff = minimizer_align("GATTACAGATTACA", "CCCCCCCCCCCCC",  k=3, w=5)
        assert len(same) >= len(diff)


# =============================================================================
# 8. Iterative Refinement (MUSCLE-style)
# =============================================================================

class TestIterativeRefinement:

    def test_refinement_never_worsens_sp_score(self):
        """One round of partition-realign-accept must not decrease the SP score."""
        initial = sum_of_pairs(SMALL_MSA)
        _, refined = refine_once(SMALL_MSA, needleman_wunsch)
        assert refined >= initial

    def test_all_sequences_same_length_after_refinement(self):
        """All sequences in the refined MSA must have equal length."""
        msa, _ = refine_once(SMALL_MSA, needleman_wunsch)
        assert len({len(s) for s in msa}) == 1

    def test_perfect_alignment_outscores_mismatched(self):
        """SP score of a perfect alignment must exceed a mismatched one."""
        assert sum_of_pairs(["AAAA", "AAAA"]) > sum_of_pairs(["AAAA", "TTTT"])


# =============================================================================
# 9. Profile HMM (Viterbi Decoding)
# =============================================================================

class TestProfileHMM:

    def test_emission_probabilities_sum_to_one_and_are_positive(self):
        """Laplace smoothing must keep all probabilities > 0 and summing to 1."""
        model = ProfileHMM(SMALL_MSA)
        for state, probs in model.emissions.items():
            if not probs:
                continue
            assert all(p > 0 for p in probs.values())
            assert abs(sum(probs.values()) - 1.0) < 1e-6

    def test_transition_probabilities_sum_to_one(self):
        model = ProfileHMM(SMALL_MSA)
        for state, dests in model.transitions.items():
            assert abs(sum(dests.values()) - 1.0) < 1e-6

    def test_viterbi_score_is_finite(self):
        """With Laplace smoothing the Viterbi log-probability must not be -inf."""
        _, score = ProfileHMM(SMALL_MSA).viterbi("ACGT")
        assert not math.isinf(score)

    def test_viterbi_path_contains_only_valid_states(self):
        model = ProfileHMM(SMALL_MSA)
        path, _ = model.viterbi("ACGTACGT")
        assert len(path) > 0
        assert all(s in model.states for s in path)


# =============================================================================
# 10. Progressive MSA + UPGMA
# =============================================================================

class TestProgressiveMSA:

    SEQS = {"s1": "ACGT", "s2": "AGGT", "s3": "ACGA"}

    def _build(self, seqs):
        labels, matrix = generate_matrix(seqs)
        return perform_progressive_alignment(run_upgma(labels, matrix), seqs)

    def test_output_count_matches_input(self):
        assert len(self._build(self.SEQS)) == len(self.SEQS)

    def test_all_sequences_same_length(self):
        """Core MSA invariant: every aligned sequence must be the same length."""
        assert len({len(s) for s in self._build(self.SEQS)}) == 1


class TestUPGMA:

    def test_all_labels_present_in_guide_tree(self):
        labels = ["A", "B", "C"]
        matrix = [[0, 1, 2], [1, 0, 3], [2, 3, 0]]
        tree   = run_upgma(labels, matrix)
        def leaves(t):
            return [t] if isinstance(t, str) else leaves(t[0]) + leaves(t[1])
        assert sorted(leaves(tree)) == sorted(labels)

    def test_single_sequence_returns_itself(self):
        assert run_upgma(["A"], [[0]]) == "A"


# =============================================================================
# 11. Long-sequence tests (HBB locus) — skipped if hbb.fasta is absent
# =============================================================================

class TestLongSequence:

    def _hbb(self):
        seqs = _load_fasta("hbb.fasta")
        return next(iter(seqs.values())) if seqs else None

    def test_hbb_file_loads(self):
        seq = self._hbb()
        if seq is None:
            print("    SKIP: hbb.fasta not found"); return
        assert len(seq) > 1000

    def test_hirschberg_completes_on_5kb(self):
        """Hirschberg (O(n) space) must finish on a 5 kb window without error."""
        seq = self._hbb()
        if seq is None:
            print("    SKIP: hbb.fasta not found"); return
        a1, a2, score = hirschberg(seq[:5000], seq[100:5100])
        assert len(a1) == len(a2) and score > 0

    def test_blast_is_not_slower_than_nw_on_2kb(self):
        """BLAST-lite must not be dramatically slower than NW (allows 2× tolerance)."""
        import time
        seq = self._hbb()
        if seq is None:
            print("    SKIP: hbb.fasta not found"); return
        s1, s2 = seq[:2000], seq[50:2050]
        t0 = time.perf_counter(); blast_lite(s1, s2);       t_b = time.perf_counter() - t0
        t0 = time.perf_counter(); needleman_wunsch(s1, s2); t_n = time.perf_counter() - t0
        assert t_b <= t_n * 2, f"BLAST {t_b:.3f}s vs NW {t_n:.3f}s"


# =============================================================================
# 12. Real FASTA identity statistics
# =============================================================================

class TestFASTAIdentityStats:
    """
    Computes mean pairwise NW identity over all pairs in each real FASTA dataset
    and asserts the result falls within a biologically plausible range.
    Results are printed so they can be cited in the report.
    """

    # (file, expected_mean_lo, expected_mean_hi)
    DATASETS = [
        ("cytochrome_c.fasta",     80, 100),  # closely related cross-species
        ("serine_proteases.fasta", 20,  80),  # distant homologs
        ("synthetic_snps.fasta",   90, 100),  # controlled ~1% mutation
        ("globins.fasta",          20, 100),  # wide range: closely related to distant
    ]

    def _check(self, fname, lo, hi):
        from statistics import mean, stdev
        seqs = _load_fasta(fname)
        if not seqs:
            print(f"    SKIP: {fname} not found"); return
        ids = _pairwise_identities(seqs)
        m   = mean(ids)
        sd  = stdev(ids) if len(ids) > 1 else 0.0
        print(f"    {fname}: {m:.1f}% ± {sd:.1f}%  (n={len(ids)} pairs)")
        assert lo <= m <= hi, \
            f"{fname}: mean identity {m:.1f}% outside expected range {lo}–{hi}%"

    def test_cytochrome_c_identity(self):
        self._check(*self.DATASETS[0])

    def test_serine_proteases_identity(self):
        self._check(*self.DATASETS[1])

    def test_synthetic_snps_identity(self):
        self._check(*self.DATASETS[2])

    def test_globins_identity(self):
        self._check(*self.DATASETS[3])


# =============================================================================
# Test runner
# =============================================================================

if __name__ == "__main__":
    import traceback

    ALL_CLASSES = [
        TestNeedlemanWunsch,
        TestSmithWaterman,
        TestGotoh,
        TestHirschberg,
        TestBandedAlignment,
        TestBLASTLite,
        TestMinimizerAlign,
        TestIterativeRefinement,
        TestProfileHMM,
        TestProgressiveMSA,
        TestUPGMA,
        TestLongSequence,
        TestFASTAIdentityStats,
    ]

    passed = failed = errors = 0

    for cls in ALL_CLASSES:
        instance = cls()
        methods  = [m for m in dir(instance) if m.startswith("test_")]
        print(f"\n{cls.__name__}  ({len(methods)} tests)")
        for name in methods:
            try:
                getattr(instance, name)()
                print(f"  PASS  {name}")
                passed += 1
            except AssertionError as e:
                print(f"  FAIL  {name}: {e}")
                failed += 1
            except Exception as e:
                print(f"  ERROR {name}: {e}")
                traceback.print_exc()
                errors += 1

    print(f"RESULTS: {passed} passed  |  {failed} failed  |  {errors} errors")