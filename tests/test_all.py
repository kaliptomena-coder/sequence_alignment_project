# Tests are organized by algorithm. Each test class checks:
#   - Correctness on known examples
#   - Edge cases (empty sequences, identical sequences, etc.)
#   - Fundamental invariants (e.g., aligned strings must have equal length)
#   - Score properties (e.g., harsher gap penalty cannot increase the score)
#
# Run with:
#   python test_all.py

import sys
import os

# Adding the src directory to the import path
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src'))
if src_path not in sys.path:
    sys.path.insert(0, src_path)

from needlemanWunschGlobal  import needleman_wunsch
from smithWatermanLocal     import smith_waterman
from gotohAffineGap         import gotoh_affine_gap
from hirschberg             import hirschberg
from banded_alignment       import banded_nw
from blast_lite             import blast_lite
from minimizer_align        import minimizer_align
from iterative_refinement   import sum_of_pairs, refine_once
from profile_hmm            import ProfileHMM


# =============================================================================
# 1. Needleman-Wunsch Tests
# =============================================================================

class TestNeedlemanWunsch:

    def test_identical_sequences_perfect_score(self):
        """Two identical sequences should produce a perfect alignment with no gaps."""
        a1, a2, score = needleman_wunsch("ACGT", "ACGT")
        assert a1 == "ACGT", f"Expected ACGT but got {a1}"
        assert a2 == "ACGT", f"Expected ACGT but got {a2}"
        assert score > 0, f"Perfect match should be positive, got {score}"

    def test_completely_different_sequences(self):
        """Completely different sequences should still return an alignment without crashing."""
        a1, a2, score = needleman_wunsch("AAAA", "TTTT")
        assert a1 is not None
        assert len(a1) == len(a2), f"Aligned strings must be equal length"

    def test_empty_sequence_returns_all_gaps(self):
        """Aligning an empty sequence should produce an all-gap result for seq1."""
        a1, a2, score = needleman_wunsch("", "ACG")
        assert set(a1) <= {'-', ''}, f"Expected all gaps, got: {a1}"

    def test_single_char_match(self):
        """Single matching characters should give a positive score."""
        a1, a2, score = needleman_wunsch("A", "A")
        assert score > 0

    def test_single_char_mismatch(self):
        """A single mismatch should return the mismatch penalty as the score."""
        a1, a2, score = needleman_wunsch("A", "T", match=1, mismatch=-1, gap=-2)
        assert score == -1, f"Expected -1 for one mismatch, got {score}"

    def test_gattaca_alignment(self):
        """Classic GATTACA vs GATCA — optimal alignment introduces one gap."""
        a1, a2, score = needleman_wunsch("GATTACA", "GATCA")
        assert len(a1) == len(a2)
        assert '-' in a1 or '-' in a2, "Should contain at least one gap"
        assert score > 0

    def test_aligned_strings_equal_length(self):
        """Fundamental invariant: both aligned strings must always be equal length."""
        for s1, s2 in [("ACGT", "AGT"), ("HELLO", "HELP"), ("A", "ACGT")]:
            a1, a2, _ = needleman_wunsch(s1, s2)
            assert len(a1) == len(a2), f"Unequal lengths for {s1}/{s2}"

    def test_gap_penalty_affects_score(self):
        """A harsher gap penalty should produce a lower or equal score."""
        _, _, score_lenient = needleman_wunsch("GATTACA", "GATCA", gap=-1)
        _, _, score_harsh   = needleman_wunsch("GATTACA", "GATCA", gap=-10)
        assert score_lenient >= score_harsh

    def test_symmetry(self):
        """Swapping the two sequences should give the same score."""
        _, _, score_ab = needleman_wunsch("GATTACA", "GATCA")
        _, _, score_ba = needleman_wunsch("GATCA", "GATTACA")
        assert score_ab == score_ba


# =============================================================================
# 2. Smith-Waterman Tests
# =============================================================================

class TestSmithWaterman:

    def test_local_score_never_negative(self):
        """SW scores are always >= 0 because the matrix is floored at zero."""
        _, _, score = smith_waterman("AAAA", "TTTT")
        assert score >= 0

    def test_identical_sequences(self):
        """Identical sequences should produce a positive local score."""
        a1, a2, score = smith_waterman("ACGT", "ACGT")
        assert score > 0

    def test_local_finds_substring(self):
        """SW should find a meaningful local alignment on a standard example."""
        a1, a2, score = smith_waterman("ACACACTA", "AGCACACA",
                                       match=2, mismatch=-1, gap=-1)
        assert score > 0
        assert len(a1) > 0 and len(a2) > 0

    def test_no_common_subsequence(self):
        """If no common subsequence exists, the score should be 0."""
        a1, a2, score = smith_waterman("AAAA", "TTTT",
                                       match=1, mismatch=-10, gap=-10)
        assert score == 0

    def test_aligned_lengths_equal(self):
        """Both returned aligned strings must always have equal length."""
        a1, a2, _ = smith_waterman("GATTACA", "GATCA")
        assert len(a1) == len(a2)


# =============================================================================
# 3. Gotoh Algorithm Tests
# =============================================================================

class TestGotoh:

    def test_returns_three_values(self):
        """Gotoh must return exactly three values: (align1, align2, score)."""
        result = gotoh_affine_gap("GATTACA", "GATCA")
        assert len(result) == 3

    def test_single_gap_introduces_dash(self):
        """One deletion should introduce at least one gap character."""
        a1, a2, score = gotoh_affine_gap("GATTACA", "GATACA",
                                         match=2, mismatch=-1,
                                         gap_open=-5, gap_extend=-1)
        assert len(a1) == len(a2)
        assert '-' in a1 or '-' in a2

    def test_identical_sequences_no_gap(self):
        """Identical sequences should produce no gaps at all."""
        a1, a2, score = gotoh_affine_gap("ACGT", "ACGT")
        assert '-' not in a1
        assert '-' not in a2

    def test_aligned_strings_equal_length(self):
        """Both aligned strings must be equal length."""
        a1, a2, _ = gotoh_affine_gap("GATTACAACTTG", "GATCCAGTTCAAA")
        assert len(a1) == len(a2)


# 4. Hirschberg Tests

class TestHirschberg:

    def test_matches_nw_score(self):
        """Hirschberg must produce the same score as standard NW on all test pairs."""
        test_pairs = [
            ("AGTAACG", "ACATAG"),
            ("GATTACA", "GATCA"),
            ("ACGT",    "ACGT"),
            ("AAAA",    "TTTT"),
        ]
        for s1, s2 in test_pairs:
            hb_a1, hb_a2  = hirschberg(s1, s2, match=2, mismatch=-1, gap=-1)
            _, _, nw_score = needleman_wunsch(s1, s2, match=2, mismatch=-1, gap=-1)
            hb_score       = _score_alignment(hb_a1, hb_a2, match=2, mismatch=-1, gap=-1)
            assert hb_score == nw_score, f"Hirschberg {hb_score} != NW {nw_score} for {s1}/{s2}"

    def test_empty_sequence(self):
        """Aligning an empty sequence should produce equal-length aligned strings."""
        a1, a2 = hirschberg("", "ACG")
        assert len(a1) == len(a2)

    def test_single_character(self):
        """Single identical characters should align to themselves."""
        a1, a2 = hirschberg("A", "A")
        assert a1 == "A" and a2 == "A"

    def test_output_only_valid_chars(self):
        """Aligned strings should only contain DNA characters and dashes."""
        a1, a2 = hirschberg("GATTACA", "GATCACA")
        valid = set("ACGT-")
        assert len(a1) == len(a2)
        assert all(c in valid for c in a1)
        assert all(c in valid for c in a2)

# 5. Banded DP Tests

class TestBandedDP:

    def test_identical_sequences_k1(self):
        """Identical sequences with a minimal band width should align correctly."""
        a1, a2, score = banded_nw("ACGT", "ACGT", k=1)
        assert score > 0
        assert len(a1) == len(a2)

    def test_matches_full_nw_for_similar_seqs(self):
        """For identical sequences, banded NW should give the same score as full NW."""
        s1, s2 = "GATTACAACTTG", "GATTACAACTTG"
        _, _, score_banded = banded_nw(s1, s2, k=5, match=2, mismatch=-1, gap=-2)
        _, _, score_nw     = needleman_wunsch(s1, s2, match=2, mismatch=-1, gap=-2)
        assert score_banded == score_nw

    def test_narrow_band_triggers_fallback(self):
        """A band that is too narrow should fall back gracefully without crashing."""
        a1, a2, score = banded_nw("AAAAAAAAAAA", "TTTTTT", k=1)
        assert a1 is not None
        assert len(a1) == len(a2)

# 6. BLAST-Lite Tests

class TestBLASTLite:

    def test_identical_sequences_finds_hsp(self):
        """Identical sequences should produce at least one HSP."""
        results = blast_lite("GATTACA", "GATTACA", k=3)
        assert len(results) > 0

    def test_hsp_has_required_fields(self):
        """Every result dictionary must contain the four expected keys."""
        results = blast_lite("GGAGTCAG", "GAAGTCGG", k=3)
        for res in results:
            assert 'score'      in res
            assert 'alignment'  in res
            assert 'query_pos'  in res
            assert 'target_pos' in res

    def test_results_sorted_descending(self):
        """Results must be sorted highest-score first."""
        results = blast_lite("GATTACAGATTACA", "GATTACAGATTACA", k=3)
        if len(results) > 1:
            scores = [r['score'] for r in results]
            assert scores == sorted(scores, reverse=True)

    def test_no_common_kmer_returns_empty(self):
        """Sequences with no shared k-mer should return an empty list."""
        results = blast_lite("AAA", "TTT", k=3)
        assert results == []

    def test_alignment_strings_equal_length(self):
        """Each HSP alignment pair must have strings of equal length."""
        results = blast_lite("GATTACAGATTACA", "GATTACAGATTACA", k=3)
        for res in results:
            a1, a2 = res['alignment']
            assert len(a1) == len(a2)

# 7. Minimizer Alignment Tests

class TestMinimizerAlign:

    def test_returns_chain_list(self):
        """The fourth return value must be a list."""
        _, _, _, chain = minimizer_align("GATTACA", "GATTACA", k=3, w=5)
        assert isinstance(chain, list)

    def test_identical_sequences_finds_anchors(self):
        """Identical sequences should produce at least one shared anchor."""
        _, _, _, chain = minimizer_align("GATTACAGATTACA", "GATTACAGATTACA", k=3, w=5)
        assert len(chain) > 0

    def test_anchor_is_a_triple(self):
        """Each anchor should be a (q_pos, t_pos, length) triple."""
        _, _, _, chain = minimizer_align("GATTACAGATTACA", "GATTACAGATTACA", k=3, w=5)
        for anchor in chain:
            assert len(anchor) == 3
            assert all(isinstance(x, int) for x in anchor)

    def test_different_sequences_fewer_anchors(self):
        """More dissimilar sequences should produce fewer anchors."""
        _, _, _, same = minimizer_align("GATTACAGATTACA", "GATTACAGATTACA", k=3, w=5)
        _, _, _, diff = minimizer_align("GATTACAGATTACA", "CCCCCCCCCCCCC",  k=3, w=5)
        assert len(same) >= len(diff)


# 8. Iterative Refinement Tests

class TestIterativeRefinement:

    TEST_MSA = [
        "ACGT--ACGT",
        "ACGTAAACGT",
        "ACGT--ACGG",
    ]

    def test_sp_score_is_numeric(self):
        """Sum-of-Pairs score must be a number."""
        score = sum_of_pairs(self.TEST_MSA)
        assert isinstance(score, (int, float))

    def test_identical_msa_has_higher_sp(self):
        """A perfect MSA (all identical) should outscore a mismatched one."""
        perfect   = ["AAAA", "AAAA", "AAAA"]
        imperfect = ["AAAA", "TTTT", "AAAA"]
        assert sum_of_pairs(perfect) > sum_of_pairs(imperfect)

    def test_refinement_does_not_worsen_score(self):
        """One round of refinement must not decrease the SP score."""
        initial_score = sum_of_pairs(self.TEST_MSA)
        _, refined_score = refine_once(self.TEST_MSA, needleman_wunsch)
        assert refined_score >= initial_score

    def test_sequences_same_length_after_refinement(self):
        """All sequences in the refined MSA must have the same length."""
        refined_msa, _ = refine_once(self.TEST_MSA, needleman_wunsch)
        lengths = [len(s) for s in refined_msa]
        assert len(set(lengths)) == 1, f"Different lengths after refinement: {lengths}"


# 9. Profile HMM Tests

class TestProfileHMM:

    SMALL_MSA = [
        "ACGT--ACGT",
        "ACGTAAACGT",
        "ACGT--ACGG",
    ]

    def test_model_builds_without_error(self):
        """The constructor should complete without raising any exception."""
        try:
            model = ProfileHMM(self.SMALL_MSA, gap_threshold=0.5)
        except Exception as e:
            assert False, f"ProfileHMM raised: {e}"

    def test_match_columns_identified(self):
        """Column 0 (all 'A', no gaps) should be identified as a match column."""
        model = ProfileHMM(self.SMALL_MSA, gap_threshold=0.5)
        assert 0 in model.match_cols

    def test_emissions_are_valid_probabilities(self):
        """All emission values must be between 0 and 1."""
        model = ProfileHMM(self.SMALL_MSA)
        for state, probs in model.emissions.items():
            for char, prob in probs.items():
                assert 0.0 <= prob <= 1.0, f"[{state}][{char}] = {prob}"

    def test_transitions_sum_to_one(self):
        """Transition probabilities from each state must sum to approximately 1."""
        model = ProfileHMM(self.SMALL_MSA)
        for state, dests in model.transitions.items():
            total = sum(dests.values())
            assert abs(total - 1.0) < 1e-6, f"Transitions from {state} sum to {total}"

    def test_viterbi_returns_path_and_score(self):
        """Viterbi should return a (list, float) pair without crashing."""
        model    = ProfileHMM(self.SMALL_MSA)
        path, score = model.viterbi("ACGTACGT")
        assert isinstance(path,  list)
        assert isinstance(score, float)

    def test_viterbi_score_not_neg_inf(self):
        """With Laplace pseudo-counts, the Viterbi score should never be -inf."""
        import math
        model    = ProfileHMM(self.SMALL_MSA)
        _, score = model.viterbi("ACGT")
        assert not math.isinf(score), "Score is -inf — check pseudo-count smoothing."

# Helper: scoring an already-computed alignment string pair

def _score_alignment(a1, a2, match, mismatch, gap):
    """Computing the alignment score directly from two aligned strings."""
    score = 0
    for c1, c2 in zip(a1, a2):
        if c1 == '-' or c2 == '-':
            score += gap
        elif c1 == c2:
            score += match
        else:
            score += mismatch
    return score


# Test Runner

if __name__ == "__main__":
    import traceback

    test_classes = [
        TestNeedlemanWunsch,
        TestSmithWaterman,
        TestGotoh,
        TestHirschberg,
        TestBandedDP,
        TestBLASTLite,
        TestMinimizerAlign,
        TestIterativeRefinement,
        TestProfileHMM,
    ]

    total_pass = total_fail = total_error = 0

    for cls in test_classes:
        instance = cls()
        methods  = [m for m in dir(instance) if m.startswith("test_")]
        print(f"  {cls.__name__}  ({len(methods)} tests)")

        for method_name in methods:
            try:
                getattr(instance, method_name)()
                print(f"  PASS  {method_name}")
                total_pass += 1
            except AssertionError as e:
                print(f"  FAIL  {method_name}")
                print(f"         -> {e}")
                total_fail += 1
            except Exception as e:
                print(f"  ERROR {method_name}")
                print(f"         -> {e}")
                total_error += 1

    print(f"  RESULTS: {total_pass} passed, {total_fail} failed, {total_error} errors")

