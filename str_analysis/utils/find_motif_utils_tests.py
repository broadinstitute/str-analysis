import math
import unittest

from str_analysis.utils.find_motif_utils import compute_repeat_purity, find_optimal_motif_length, _compute_motif_length_quality, HAMMING_DISTANCE_METRIC
from utils.canonical_repeat_unit import compute_canonical_motif


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


    def test_find_optimal_motif_length(self):
        for expected_optimal_motif, sequences, min_quality_score in (
            ("CAG", ["CAG"*3, "CAG"*3+"C", "CAG"*3+"CA"], 0.6),
            ("CAGA", ["CAGA"*3, "CAGA"*3+"C", "CAGA"*3+"CA"], 0.5),
            ("CAGAT", ["CAGAT"*3, "CAGAT"*3+"C", "CAGAT"*3+"CA"], 0.5),
            ("CAGATC", ["CAGATC"*30, "CAGATC"*30+"C", "CAGATC"*30+"CA", "TT"+"CAGATC"*30+"CA", "TTTT"+"CAGATC"*30+"CA"], 0.5),
        ):
            for seq in sequences:
                optimal_motif, optimal_purity, optimal_motif_length_quality_score = find_optimal_motif_length(
                    nucleotide_sequence=seq,
                    max_motif_length=len(seq)//2,
                    distance_metric=HAMMING_DISTANCE_METRIC,
                    negligible_change_in_purity=0.025,
                    allow_partial_repeats=True,
                    verbose=True
                )
                self.assertEqual(compute_canonical_motif(optimal_motif), compute_canonical_motif(expected_optimal_motif),
                                 f"Unexpected optimal motif for sequence {seq}")
                self.assertGreater(optimal_motif_length_quality_score, min_quality_score,
                                   f"Unexpected optimal motif length quality score for sequence {seq}")


