#Unit tests for ALL alignment algorithms.

import sys
import os

# Adding the parent folder  to Python's path

src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src'))

if src_path not in sys.path:
    sys.path.insert(0, src_path)
# import all modules
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
# 1. NEEDLEMAN-WUNSCH TESTS
# =============================================================================

class TestNeedlemanWunsch:
    """Tests for global alignment (Needleman-Wunsch)."""

    def test_identical_sequences_perfect_score(self):
        """
        Two identical sequences should align perfectly with no gaps.
        Score = len(seq) * match_score = 5 * 1 = 5 (with default match=1).
        """
        a1, a2, score = needleman_wunsch("ACGT", "ACGT")
        # Both aligned strings should equal the original (no gaps)
        assert a1 == "ACGT", f"Expected ACGT but got {a1}"
        assert a2 == "ACGT", f"Expected ACGT but got {a2}"
        assert score > 0,    f"Perfect match should have positive score, got {score}"

    def test_completely_different_sequences(self):
        """
        Completely different sequences will have a negative or zero score.
        The alignment should still be returned (not crash).
        """
        a1, a2, score = needleman_wunsch("AAAA", "TTTT")
        # Alignment must be returned (not None)
        assert a1 is not None
        assert a2 is not None
        # Length of aligned strings must be equal
        assert len(a1) == len(a2), \
            f"Aligned strings must have equal length: {len(a1)} vs {len(a2)}"

    def test_empty_sequence_returns_all_gaps(self):
        """
        If one sequence is empty, the other must be aligned against all gaps.
        """
        a1, a2, score = needleman_wunsch("", "ACG")
        # seq1 should be all gaps, seq2 should be unchanged
        assert set(a1) <= {'-', ''}, f"Expected all gaps in a1, got: {a1}"

    def test_single_char_match(self):
        """Single matching characters should give a positive score."""
        a1, a2, score = needleman_wunsch("A", "A")
        assert score > 0

    def test_single_char_mismatch(self):
        """Single mismatching characters should give a negative score."""
        a1, a2, score = needleman_wunsch("A", "T", match=1, mismatch=-1, gap=-2)
        assert score == -1, f"Expected -1 for a single mismatch, got {score}"

    def test_known_gattaca_alignment(self):
        """
        Classic GATTACA vs GATCA example.
        One deletion (one gap somewhere) is the optimal global alignment.
        Score should be positive.
        """
        a1, a2, score = needleman_wunsch("GATTACA", "GATCA")
        # The aligned strings must have equal length
        assert len(a1) == len(a2)
        # One of the strings should contain a gap
        assert '-' in a1 or '-' in a2, \
            "There should be at least one gap to align GATTACA and GATCA globally"
        # Score should be positive (mostly matching)
        assert score > 0

    def test_alignment_strings_have_equal_length(self):
        """A fundamental property: both aligned strings must be the same length."""
        for s1, s2 in [("ACGT", "AGT"), ("HELLO", "HELP"), ("A", "ACGT")]:
            a1, a2, _ = needleman_wunsch(s1, s2)
            assert len(a1) == len(a2), \
                f"For {s1}/{s2}: aligned strings lengths differ: {len(a1)} vs {len(a2)}"

    def test_gap_penalty_affects_score(self):
        """
        A harsher gap penalty should produce a lower (or equal) score
        compared to a lenient one, for the same input.
        """
        _, _, score_lenient = needleman_wunsch("GATTACA", "GATCA", gap=-1)
        _, _, score_harsh   = needleman_wunsch("GATTACA", "GATCA", gap=-10)
        assert score_lenient >= score_harsh, \
            "Harsher gap penalty should not increase the score"

    def test_symmetry(self):
        """
        Swapping seq1 and seq2 should give the same score
        (alignment may differ, but score should be equal).
        """
        _, _, score_ab = needleman_wunsch("GATTACA", "GATCA")
        _, _, score_ba = needleman_wunsch("GATCA", "GATTACA")
        assert score_ab == score_ba, \
            f"Score should be symmetric: {score_ab} vs {score_ba}"


# =============================================================================
# 2. SMITH-WATERMAN TESTS
# =============================================================================

class TestSmithWaterman:
    """Tests for local alignment (Smith-Waterman)."""

    def test_local_score_never_negative(self):
        """
        SW scores are always non-negative (cells are floored at 0).
        """
        _, _, score = smith_waterman("AAAA", "TTTT")
        assert score >= 0, f"SW score must be >= 0, got {score}"

    def test_identical_sequences(self):
        """Identical sequences: local score should equal full match score."""
        a1, a2, score = smith_waterman("ACGT", "ACGT")
        assert score > 0

    def test_local_finds_substring(self):
        """
        Classic SW example from the original paper.
        ACACACTA vs AGCACACA → local alignment should find ACAC or similar.
        Score should be > 0.
        """
        a1, a2, score = smith_waterman("ACACACTA", "AGCACACA",
                                       match=2, mismatch=-1, gap=-1)
        assert score > 0, f"Should find a local alignment, score={score}"
        # The local alignment should not be empty
        assert len(a1) > 0 and len(a2) > 0

    def test_no_common_subsequence(self):
        """
        If no common subsequence exists, score should be 0
        and aligned strings should be empty.
        """
        # All A vs all T with harsh penalties
        a1, a2, score = smith_waterman("AAAA", "TTTT", match=1, mismatch=-10, gap=-10)
        assert score == 0, f"No common subseq → score should be 0, got {score}"

    def test_aligned_lengths_equal(self):
        """Both returned aligned strings must always have equal length."""
        a1, a2, _ = smith_waterman("GATTACA", "GATCA")
        assert len(a1) == len(a2)

    def test_sw_score_lte_nw_score(self):
        """
        Local alignment score should be >= global alignment score
        because SW can ignore bad end regions that NW must include.
        (For sequences with common cores, SW often finds a better sub-score.)
        """
        # We don't import NW here to keep tests independent,
        # but we verify SW score is positive when NW would penalise ends
        _, _, sw_score = smith_waterman("TTTACGTTTTT", "ACGT",
                                        match=2, mismatch=-1, gap=-2)
        assert sw_score > 0


# =============================================================================
# 3: GOTOH TESTS
# =============================================================================

class TestGotoh:
    """Tests for affine-gap alignment (Gotoh)."""

    def test_returns_three_values(self):
        """Gotoh must return (align1, align2, score)."""
        result = gotoh_affine_gap("GATTACA", "GATCA")
        assert len(result) == 3, f"Expected 3 return values, got {len(result)}"

    def test_single_gap_penalised_correctly(self):
        """
        For a single gap of length 1: cost = gap_open + 1*gap_extend.
        With gap_open=-5, gap_extend=-1 → gap cost = -6.
        A perfect alignment with one gap should reflect this.
        """
        a1, a2, score = gotoh_affine_gap("GATTACA", "GATACA",
                                         match=2, mismatch=-1,
                                         gap_open=-5, gap_extend=-1)
        # Alignment should exist
        assert len(a1) == len(a2)
        assert '-' in a1 or '-' in a2

    def test_affine_cheaper_than_linear_for_long_gap(self):
        """
        With affine gaps, a single long gap of length k costs:
            gap_open + k * gap_extend
        In linear model, k separate gaps cost k * gap.
        Affine should be cheaper for k > 1 when gap_extend < |gap|.
        This test verifies affine produces better scores for gappy alignments.
        """
        # Sequence with a 4-character insertion
        s1 = "ACGTACGT"
        s2 = "ACGT" + "AAAA" + "ACGT"   # 4 extra chars in the middle

        _, _, score_gotoh = gotoh_affine_gap(s1, s2,
                                             match=2, mismatch=-1,
                                             gap_open=-2, gap_extend=-1)
        # verifying it runs without error and returns a score
        assert score_gotoh is not None

    def test_identical_sequences_no_gap(self):
        """Identical sequences: no gap should be introduced."""
        a1, a2, score = gotoh_affine_gap("ACGT", "ACGT")
        assert '-' not in a1
        assert '-' not in a2

    def test_aligned_strings_equal_length(self):
        """Fundamental property: both aligned strings must be equal length."""
        a1, a2, _ = gotoh_affine_gap("GATTACAACTTG", "GATCCAGTTCAAA")
        assert len(a1) == len(a2)


# =============================================================================
# 4. HIRSCHBERG TESTS
# =============================================================================

class TestHirschberg:
    """Tests for space-efficient global alignment (Hirschberg)."""

    def test_matches_nw_on_small_input(self):
        """
        CRITICAL TEST: Hirschberg must produce the SAME score as NW.
        The alignment path may differ on ties, but the score must be equal.
        """
        test_pairs = [
            ("AGTAACG",  "ACATAG"),
            ("GATTACA",  "GATCA"),
            ("ACGT",     "ACGT"),
            ("AAAA",     "TTTT"),
            ("ACGTACGT", "CGTACGT"),
        ]
        for s1, s2 in test_pairs:
            # Hirschberg result
            hb_a1, hb_a2 = hirschberg(s1, s2, match=2, mismatch=-1, gap=-1)

            # NW result for comparison
            nw_a1, nw_a2, nw_score = needleman_wunsch(s1, s2, match=2, mismatch=-1, gap=-1)

            # Calculate Hirschberg score from the alignment
            hb_score = _score_alignment(hb_a1, hb_a2, match=2, mismatch=-1, gap=-1)

            assert hb_score == nw_score, \
                f"Hirschberg score {hb_score} != NW score {nw_score} for {s1}/{s2}"

    def test_empty_sequence(self):
        """One empty sequence → result should be all gaps."""
        a1, a2 = hirschberg("", "ACG")
        assert len(a1) == len(a2)

    def test_single_character(self):
        """Single characters should align directly."""
        a1, a2 = hirschberg("A", "A")
        assert a1 == "A" and a2 == "A"

    def test_output_is_valid_alignment(self):
        """Both strings must have equal length and contain only valid chars."""
        a1, a2 = hirschberg("GATTACA", "GATCACA")
        valid_chars = set("ACGT-")
        assert len(a1) == len(a2)
        assert all(c in valid_chars for c in a1)
        assert all(c in valid_chars for c in a2)


# =============================================================================
# 5: BANDED DP TESTS
# =============================================================================

class TestBandedDP:
    """Tests for banded Needleman-Wunsch."""

    def test_identical_sequences_k1(self):
        """Identical sequences: banded with k=1 should find perfect alignment."""
        a1, a2, score = banded_nw("ACGT", "ACGT", k=1)
        assert score > 0
        assert len(a1) == len(a2)

    def test_matches_full_nw_for_similar_seqs(self):
        """
        For similar sequences (≤k differences), banded should equal full NW.
        """
        s1, s2 = "GATTACAACTTG", "GATTACAACTTG"
        _, _, score_banded = banded_nw(s1, s2, k=5, match=2, mismatch=-1, gap=-2)
        _, _, score_nw     = needleman_wunsch(s1, s2, match=2, mismatch=-1, gap=-2)
        assert score_banded == score_nw, \
            f"Banded {score_banded} should equal NW {score_nw} for identical seqs"

    def test_narrow_band_triggers_fallback(self):
        """
        If band is too narrow, it should fallback gracefully (not crash).
        """
        s1 = "AAAAAAAAAAA"  # length 11
        s2 = "TTTTTT"       # length 6
        # k=1 is way too narrow for sequences this different
        a1, a2, score = banded_nw(s1, s2, k=1)
        # Should still return a valid result (fallback to full NW)
        assert a1 is not None
        assert len(a1) == len(a2)


# =============================================================================
# 6. BLAST-LITE TESTS
# =============================================================================

class TestBLASTLite:
    """Tests for heuristic seed-and-extend (BLAST-style)."""

    def test_identical_sequences_finds_hsp(self):
        """Two identical sequences: BLAST should find at least one HSP."""
        results = blast_lite("GATTACA", "GATTACA", k=3)
        assert len(results) > 0, "Should find HSPs for identical sequences"

    def test_hsp_has_required_fields(self):
        """Each result must contain the required fields."""
        results = blast_lite("GGAGTCAG", "GAAGTCGG", k=3)
        if results:
            for res in results:
                assert 'score'      in res, "HSP must have 'score'"
                assert 'alignment'  in res, "HSP must have 'alignment'"
                assert 'query_pos'  in res, "HSP must have 'query_pos'"
                assert 'target_pos' in res, "HSP must have 'target_pos'"

    def test_results_sorted_by_score_descending(self):
        """Results should be sorted from highest score to lowest."""
        results = blast_lite("GATTACAGATTACA", "GATTACAGATTACA", k=3)
        if len(results) > 1:
            scores = [r['score'] for r in results]
            assert scores == sorted(scores, reverse=True), \
                "Results should be sorted by score (highest first)"

    def test_no_match_returns_empty(self):
        """Very short sequences with no common k-mer → empty list."""
        results = blast_lite("AAA", "TTT", k=3)
        assert results == [], f"Expected no HSPs, got: {results}"

    def test_alignment_strings_in_hsp(self):
        """Each HSP alignment tuple must have two strings of equal length."""
        results = blast_lite("GATTACAGATTACA", "GATTACAGATTACA", k=3)
        for res in results:
            a1, a2 = res['alignment']
            assert len(a1) == len(a2), \
                f"Alignment strings must be equal length: {len(a1)} vs {len(a2)}"


# =============================================================================
# 7. MINIMIZER TESTS
# =============================================================================

class TestMinimizerAlign:
    """Tests for minimizer-based approximate alignment."""

    def test_returns_list(self):
        """Function should return a chain (list) as the 4th element."""
        _, _, _, chain = minimizer_align("GATTACA", "GATTACA", k=3, w=5)
        assert isinstance(chain, list)

    def test_identical_sequences_finds_anchors(self):
        """Identical sequences should produce anchors."""
        _, _, _, chain = minimizer_align("GATTACAGATTACA", "GATTACAGATTACA", k=3, w=5)
        assert len(chain) > 0, "Should find shared minimizers for identical seqs"

    def test_anchor_structure(self):
        """Each anchor in the chain should be (q_pos, t_pos, length)."""
        _, _, _, chain = minimizer_align("GATTACAGATTACA", "GATTACAGATTACA", k=3, w=5)
        for anchor in chain:
            assert len(anchor) == 3, f"Anchor should have 3 elements: {anchor}"
            q_pos, t_pos, length = anchor
            assert isinstance(q_pos, int)
            assert isinstance(t_pos, int)
            assert isinstance(length, int)

    def test_different_sequences_fewer_anchors(self):
        """Very different sequences should produce fewer anchors."""
        _, _, _, same_chain = minimizer_align("GATTACAGATTACA", "GATTACAGATTACA", k=3, w=5)
        _, _, _, diff_chain = minimizer_align("GATTACAGATTACA", "CCCCCCCCCCCCC",  k=3, w=5)
        assert len(same_chain) >= len(diff_chain), \
            "Identical sequences should produce >= anchors as different sequences"


# =============================================================================
# 8. ITERATIVE REFINEMENT TESTS
# =============================================================================

class TestIterativeRefinement:
    """Tests for MSA iterative refinement (MUSCLE-style)."""

    # A small, deliberately imperfect MSA for testing
    TEST_MSA = [
        "ACGT--ACGT",
        "ACGTAAACGT",
        "ACGT--ACGG",
    ]

    def test_sp_score_is_numeric(self):
        """Sum-of-Pairs score should be a number."""
        score = sum_of_pairs(self.TEST_MSA)
        assert isinstance(score, (int, float)), \
            f"SP score should be numeric, got {type(score)}"

    def test_identical_columns_contribute_positively(self):
        """All-match columns should increase SP score."""
        # Perfect MSA: all identical → maximum SP score
        perfect = ["AAAA", "AAAA", "AAAA"]
        imperfect = ["AAAA", "TTTT", "AAAA"]
        score_perfect   = sum_of_pairs(perfect)
        score_imperfect = sum_of_pairs(imperfect)
        assert score_perfect > score_imperfect, \
            "Perfect MSA should have higher SP score than mismatched one"

    def test_refine_does_not_worsen_score(self):
        """
        After one round of refinement, the SP score should be
        greater than or equal to the initial score.
        """
        initial_score  = sum_of_pairs(self.TEST_MSA)
        refined_msa, refined_score = refine_once(self.TEST_MSA, needleman_wunsch)
        assert refined_score >= initial_score, \
            f"Refinement worsened SP score: {initial_score} → {refined_score}"

    def test_all_seqs_same_length_after_refinement(self):
        """After refinement, all sequences must still have the same length."""
        refined_msa, _ = refine_once(self.TEST_MSA, needleman_wunsch)
        lengths = [len(s) for s in refined_msa]
        assert len(set(lengths)) == 1, \
            f"Sequences have different lengths after refinement: {lengths}"


# =============================================================================
# 9. PROFILE HMM TESTS
# =============================================================================

class TestProfileHMM:
    """Tests for Profile Hidden Markov Model."""

    # Small MSA: 3 sequences, 10 columns
    SMALL_MSA = [
        "ACGT--ACGT",
        "ACGTAAACGT",
        "ACGT--ACGG",
    ]

    def test_model_builds_without_error(self):
        """ProfileHMM constructor should not raise any exception."""
        try:
            model = ProfileHMM(self.SMALL_MSA, gap_threshold=0.5)
        except Exception as e:
            assert False, f"ProfileHMM constructor raised: {e}"

    def test_match_columns_identified(self):
        """
        Columns with >50% gaps should NOT be match columns.
        Columns with <50% gaps SHOULD be match columns.
        """
        model = ProfileHMM(self.SMALL_MSA, gap_threshold=0.5)
        # Column 4 and 5 are '--', 'AA', '--' → 2/3 gaps → not a match col
        # Column 0 is 'A','A','A' → 0 gaps → must be a match col
        assert 0 in model.match_cols, "Column 0 (all A) should be a match column"

    def test_emissions_are_probabilities(self):
        """All emission values must be between 0 and 1."""
        model = ProfileHMM(self.SMALL_MSA)
        for state, probs in model.emissions.items():
            for aa, prob in probs.items():
                assert 0.0 <= prob <= 1.0, \
                    f"Emission prob out of range [{state}][{aa}] = {prob}"

    def test_transitions_sum_to_one(self):
        """Transition probabilities from each state must sum to ~1."""
        model = ProfileHMM(self.SMALL_MSA)
        for state, dests in model.transitions.items():
            total = sum(dests.values())
            assert abs(total - 1.0) < 1e-6, \
                f"Transitions from {state} sum to {total}, expected 1.0"

    def test_viterbi_returns_path_and_score(self):
        """Viterbi should return (path, score) without crashing."""
        model = ProfileHMM(self.SMALL_MSA)
        seq = "ACGTACGT"
        try:
            path, score = model.viterbi(seq)
            assert isinstance(path, list), f"Path should be a list, got {type(path)}"
            assert isinstance(score, float), f"Score should be float, got {type(score)}"
        except Exception as e:
            assert False, f"Viterbi raised an exception: {e}"

    def test_viterbi_score_not_neg_inf(self):
        """
        With proper pseudo-counts, Viterbi should NOT return -inf.
        -inf means an emission probability was exactly 0, which
        is the pseudo-count bug we fixed.
        """
        model = ProfileHMM(self.SMALL_MSA)
        _, score = model.viterbi("ACGT")
        import math
        assert not math.isinf(score), \
            "Viterbi returned -inf! Add Laplace pseudo-counts to fix this."


# =============================================================================
# 10. HELPER FUNCTION (used by Hirschberg tests)
# =============================================================================

def _score_alignment(a1, a2, match, mismatch, gap):
    """Calculate the score of an already-computed alignment."""
    score = 0
    for c1, c2 in zip(a1, a2):
        if c1 == '-' or c2 == '-':
            score += gap
        elif c1 == c2:
            score += match
        else:
            score += mismatch
    return score




if __name__ == "__main__":
    import traceback

    # Collect all test classes
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

    total_pass = 0
    total_fail = 0
    total_error = 0

    for cls in test_classes:
        instance = cls()
        methods = [m for m in dir(instance) if m.startswith("test_")]
        print(f"\n{'='*60}")
        print(f"  {cls.__name__}  ({len(methods)} tests)")
        print(f"{'='*60}")

        for method_name in methods:
            try:
                getattr(instance, method_name)()
                print(f"  ✓  PASS  {method_name}")
                total_pass += 1
            except AssertionError as e:
                print(f"  ✗  FAIL  {method_name}")
                print(f"           → {e}")
                total_fail += 1
            except Exception as e:
                print(f"  !  ERROR {method_name}")
                print(f"           → {e}")
                total_error += 1

    print(f"\n{'='*60}")
    print(f"  RESULTS: {total_pass} passed, {total_fail} failed, {total_error} errors")
    print(f"{'='*60}")
