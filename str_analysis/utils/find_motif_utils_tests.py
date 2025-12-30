import math
import unittest

from str_analysis.utils.find_motif_utils import (
    compute_repeat_purity,
    find_highest_purity_motif_length,
    find_highest_purity_motif_length_for_interval,
    _compute_motif_length_quality,
    compute_most_common_motif,
    compute_motif_length_purity,
    compute_motif_length_purity_for_interval,
    compute_motif_purity_for_interval,
    compute_motif_null_quality_score_for_sequence_length,
    HAMMING_DISTANCE_METRIC,
    EDIT_DISTANCE_METRIC
)
from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif


class Tests(unittest.TestCase):

    def test_compute_repeat_purity(self):
        tests = [
            ("CACACACACA",      "CA",    0, 1,               False),
            ("CACACACACAC",     "CA",    0, 1,               False),
            ("CACACTGTGT",      "CA",    5, 0.5,             False),
            ("CA",              "CTTG",  None, float('nan'), False),
            ("CAGCAACAGCAGTT",  "CAG",   1, 11/12,           False),
            ("CAGCAACAGCAGTT",  "CAG",   3, 11/14,           True),
        ]

        for sequence, motif, expected_interruptions, expected_purity, include_partial_repeats in tests:
            purity, interruption_count = compute_repeat_purity(sequence, motif, include_partial_repeats=include_partial_repeats)
            self.assertEqual(expected_interruptions, interruption_count, msg=f"Expected interruption count {expected_interruptions} but got {interruption_count} for sequence {sequence} & motif {motif}")
            if math.isnan(expected_purity):
                self.assertTrue(math.isnan(purity), msg=f"Expected purity to be NaN but got {purity} for sequence {sequence} & motif {motif}")
            else:
                self.assertAlmostEquals(expected_purity, purity, msg=f"Expected base purity {expected_purity} but got {purity} for sequence {sequence} & motif {motif}")


    def test_compute_motif_length_quality(self):
        quality = _compute_motif_length_quality(3, {3: ("CAG", 0.5)})
        self.assertEqual(0.5 - 0.3, quality)

        quality = _compute_motif_length_quality(3, {2: ("CA", 0.1), 3: ("CAG", 0.5)})
        self.assertEqual(0.5 - 0.1, quality)


    def test_find_highest_purity_motif_length(self):
        for expected_highest_purity_motif, sequences, min_quality_score in (
            ("CAG", ["CAG"*3, "CAG"*3+"C", "CAG"*3+"CA"], 0.6),
            ("CAGA", ["CAGA"*3, "CAGA"*3+"C", "CAGA"*3+"CA"], 0.5),
            ("CAGAT", ["CAGAT"*3, "CAGAT"*3+"C", "CAGAT"*3+"CA"], 0.5),
            ("CAGATC", ["CAGATC"*30, "CAGATC"*30+"C", "CAGATC"*30+"CA", "TT"+"CAGATC"*30+"CA", "TTTT"+"CAGATC"*30+"CA"], 0.5),
        ):
            for seq in sequences:
                highest_purity_motif, highest_purity, highest_purity_motif_length_quality_score = find_highest_purity_motif_length(
                    nucleotide_sequence=seq,
                    max_motif_length=len(seq)//2,
                    distance_metric=HAMMING_DISTANCE_METRIC,
                    negligible_change_in_purity=0.025,
                    allow_partial_repeats=True,
                    verbose=True
                )
                self.assertEqual(compute_canonical_motif(highest_purity_motif), compute_canonical_motif(expected_highest_purity_motif),
                                 f"Unexpected optimal motif for sequence {seq}")
                self.assertGreater(highest_purity_motif_length_quality_score, min_quality_score,
                                   f"Unexpected optimal motif length quality score for sequence {seq}")


    def test_compute_repeat_purity_with_edit_distance(self):
        # Test with edit distance metric
        purity, distance = compute_repeat_purity("CACACACACA", "CA", distance_metric=EDIT_DISTANCE_METRIC)
        self.assertEqual(distance, 0)
        self.assertEqual(purity, 1.0)

        # Test with insertions/deletions
        purity, distance = compute_repeat_purity("CACACTGTGT", "CA", distance_metric=EDIT_DISTANCE_METRIC)
        self.assertLess(purity, 1.0)
        self.assertGreater(distance, 0)

        # Test error handling for invalid distance > sequence length
        # This should not raise an exception for valid inputs
        purity, distance = compute_repeat_purity("CAG", "CAG", distance_metric=EDIT_DISTANCE_METRIC)
        self.assertEqual(purity, 1.0)

    def test_compute_repeat_purity_error_cases(self):
        # Test with unknown distance metric
        with self.assertRaises(ValueError):
            compute_repeat_purity("CACACA", "CA", distance_metric="unknown_metric")

    def test_find_highest_purity_motif_length_edge_cases(self):
        # Test with empty sequence
        with self.assertRaises(ValueError):
            find_highest_purity_motif_length("", max_motif_length=10)

        # Test with max_motif_length < 1
        with self.assertRaises(ValueError):
            find_highest_purity_motif_length("CAGCAG", max_motif_length=0)

        # Test with min_motif_length > max_motif_length
        with self.assertRaises(ValueError):
            find_highest_purity_motif_length("CAGCAG", min_motif_length=10, max_motif_length=5)

        # Test with very long sequence (should truncate)
        long_sequence = "CAG" * 5000  # 15000bp
        motif, purity, quality = find_highest_purity_motif_length(
            long_sequence, max_motif_length=10, verbose=False
        )
        self.assertEqual(motif, "CAG")

        # Test with max_motif_length > 200 (should reduce to 200)
        motif, purity, quality = find_highest_purity_motif_length(
            "CAG" * 100, max_motif_length=300, verbose=False
        )
        self.assertEqual(motif, "CAG")

        # Test with max_motif_length > sequence length (should adjust)
        motif, purity, quality = find_highest_purity_motif_length(
            "CAGCAG", max_motif_length=100, verbose=False
        )
        self.assertIsNotNone(motif)

        # Test with allow_partial_repeats=False
        motif, purity, quality = find_highest_purity_motif_length(
            "CAGCAGC", max_motif_length=6, allow_partial_repeats=False, verbose=False
        )
        self.assertIsNotNone(motif)

        # Test with edit distance metric
        motif, purity, quality = find_highest_purity_motif_length(
            "CAGCAGCAG", max_motif_length=6, distance_metric=EDIT_DISTANCE_METRIC, verbose=False
        )
        self.assertEqual(motif, "CAG")

        # Test with unknown distance metric
        with self.assertRaises(ValueError):
            find_highest_purity_motif_length(
                "CAGCAG", max_motif_length=6, distance_metric="unknown"
            )

        # Test case where no valid motif is found (very short sequence)
        motif, purity, quality = find_highest_purity_motif_length(
            "AC", max_motif_length=10, min_motif_length=5, verbose=False
        )
        # Should return None when no valid motifs found
        self.assertTrue(motif is None or isinstance(motif, str))

    def test_compute_most_common_motif_edge_cases(self):
        # Test with sequence length < repeat_unit_length
        with self.assertRaises(ValueError):
            compute_most_common_motif("CA", 3)

        # Test with sequence length == repeat_unit_length
        result = compute_most_common_motif("CAG", 3)
        self.assertEqual(result, "CAG")

        # Test normal case
        result = compute_most_common_motif("CAGCAGCAG", 3)
        self.assertEqual(result, "CAG")

    def test_compute_motif_length_purity_edge_cases(self):
        # Test with sequence length < motif_length
        purity, motif = compute_motif_length_purity("CA", 3)
        self.assertTrue(math.isnan(purity))
        self.assertIsNone(motif)

        # Test with motif_length == sequence length
        purity, motif = compute_motif_length_purity("CAG", 3)
        self.assertEqual(purity, 1)
        self.assertEqual(motif, "CAG")

        # Test with zero-length sliced list
        purity, motif = compute_motif_length_purity("CAG", 4)
        self.assertTrue(math.isnan(purity))
        self.assertIsNone(motif)

        # Test normal case with edit distance
        purity, motif = compute_motif_length_purity("CAGCAGCAG", 3, distance_metric=EDIT_DISTANCE_METRIC)
        self.assertEqual(motif, "CAG")
        self.assertGreater(purity, 0.9)

    def test_pyfaidx_dependent_functions(self):
        # Mock pyfaidx.Fasta object for testing interval-based functions
        class MockFaidx:
            def __init__(self):
                self.one_based_attributes = False
                self.as_raw = True

        class MockChromosome:
            def __init__(self, sequence):
                self.sequence = sequence

            def __len__(self):
                return len(self.sequence)

            def __getitem__(self, key):
                if isinstance(key, slice):
                    return self.sequence[key]
                return self.sequence[key]

            def __str__(self):
                return self.sequence

        class MockFasta:
            def __init__(self, sequences):
                self.sequences = {chrom: MockChromosome(seq) for chrom, seq in sequences.items()}
                self.faidx = MockFaidx()

            def __getitem__(self, chrom):
                return self.sequences[chrom]

            def __contains__(self, chrom):
                return chrom in self.sequences

            def keys(self):
                return self.sequences.keys()

        # Test compute_motif_purity_for_interval
        reference = {"chr1": "AAAA" + "CAG"*10 + "TTTT"}
        mock_fasta = MockFasta(reference)
        purity, distance = compute_motif_purity_for_interval(
            mock_fasta, "chr1", 4, 34, "CAG", distance_metric=HAMMING_DISTANCE_METRIC
        )
        self.assertGreater(purity, 0.95)

        # Test with None reference
        purity, distance = compute_motif_purity_for_interval(
            None, "chr1", 4, 34, "CAG"
        )
        self.assertTrue(math.isnan(purity))
        self.assertIsNone(distance)

        # Test find_highest_purity_motif_length_for_interval
        motif, purity, quality = find_highest_purity_motif_length_for_interval(
            mock_fasta, "chr1", 4, 34, max_motif_length=10, verbose=False
        )
        self.assertEqual(motif, "CAG")

        # Test with None reference
        motif, purity, quality = find_highest_purity_motif_length_for_interval(
            None, "chr1", 4, 34, max_motif_length=10
        )
        self.assertIsNone(motif)
        self.assertTrue(math.isnan(purity))

        # Test with invalid pyfaidx attributes
        class BadMockFasta(MockFasta):
            def __init__(self, sequences):
                super().__init__(sequences)
                self.faidx.one_based_attributes = True

        bad_fasta = BadMockFasta(reference)
        with self.assertRaises(ValueError):
            find_highest_purity_motif_length_for_interval(
                bad_fasta, "chr1", 4, 34, max_motif_length=10
            )

        class BadMockFasta2(MockFasta):
            def __init__(self, sequences):
                super().__init__(sequences)
                self.faidx.as_raw = False

        bad_fasta2 = BadMockFasta2(reference)
        with self.assertRaises(ValueError):
            find_highest_purity_motif_length_for_interval(
                bad_fasta2, "chr1", 4, 34, max_motif_length=10
            )

        # Test compute_motif_length_purity_for_interval
        purity, motif = compute_motif_length_purity_for_interval(
            mock_fasta, "chr1", 4, 34, motif_length=3
        )
        self.assertEqual(motif, "CAG")
        self.assertGreater(purity, 0.9)

        # Test with None reference
        purity, motif = compute_motif_length_purity_for_interval(
            None, "chr1", 4, 34, motif_length=3
        )
        self.assertTrue(math.isnan(purity))
        self.assertIsNone(motif)

    def test_compute_motif_null_quality_score(self):
        # Test with HAMMING_DISTANCE_METRIC
        score = compute_motif_null_quality_score_for_sequence_length(100, distance_metric=HAMMING_DISTANCE_METRIC)
        self.assertIsInstance(score, (float, int))
        self.assertGreater(score, 0)

        # Test with EDIT_DISTANCE_METRIC
        score = compute_motif_null_quality_score_for_sequence_length(100, distance_metric=EDIT_DISTANCE_METRIC)
        self.assertIsInstance(score, (float, int))
        self.assertGreater(score, 0)

        # Test with invalid distance metric
        with self.assertRaises(ValueError):
            compute_motif_null_quality_score_for_sequence_length(100, distance_metric="unknown")

    def test_find_highest_purity_motif_additional_edge_cases(self):
        # Test sequence that triggers motif length > sequence//2 skip
        motif, purity, quality = find_highest_purity_motif_length(
            "CAGCAGCAG", max_motif_length=8, verbose=False
        )
        self.assertEqual(motif, "CAG")

        # Test sequence that reaches 100% purity and breaks early
        motif, purity, quality = find_highest_purity_motif_length(
            "CAGCAGCAGCAGCAG", max_motif_length=10, verbose=False
        )
        self.assertEqual(motif, "CAG")
        self.assertGreaterEqual(purity, 0.99)

        # Test with verbose=True to hit verbose output paths
        motif, purity, quality = find_highest_purity_motif_length(
            "CAG" * 5000, max_motif_length=300, verbose=True
        )
        self.assertEqual(motif, "CAG")


