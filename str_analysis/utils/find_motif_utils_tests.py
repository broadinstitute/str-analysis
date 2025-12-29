import math
import unittest

from str_analysis.utils.find_motif_utils import (
    compute_repeat_purity,
    find_highest_purity_motif_length,
    _compute_motif_length_quality,
    adjust_motif_and_boundaries_to_maximize_purity,
    HAMMING_DISTANCE_METRIC
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

    def test_adjust_motif_and_boundaries_to_maximize_purity(self):
        class MockChromosome:
            def __init__(self, sequence):
                self.sequence = sequence

            def __len__(self):
                return len(self.sequence)

            def __getitem__(self, key):
                if isinstance(key, slice):
                    return self.sequence[key]
                return self.sequence[key]

        class MockFasta:
            def __init__(self, sequences):
                self.sequences = {chrom: MockChromosome(seq) for chrom, seq in sequences.items()}

            def __getitem__(self, chrom):
                return self.sequences[chrom]

        # Test 1: No adjustment needed - repeat doesn't extend into flanks
        reference = {"chr1": "AAAA" + "CAG"*10 + "TTTT"}
        mock_fasta = MockFasta(reference)
        start, end, motif, adjusted = adjust_motif_and_boundaries_to_maximize_purity(
            mock_fasta, "chr1", 4, 34, "CAG"
        )
        self.assertEqual(start, 4)
        self.assertEqual(end, 34)
        self.assertEqual(motif, "CAG")
        self.assertFalse(adjusted)

        # Test 2: Extension into left flank
        reference = {"chr1": "CAG"*5 + "CAG"*10 + "TTTT"}
        mock_fasta = MockFasta(reference)
        start, end, motif, adjusted = adjust_motif_and_boundaries_to_maximize_purity(
            mock_fasta, "chr1", 15, 45, "CAG"
        )
        self.assertEqual(start, 0)  # Extended to start of chromosome
        self.assertEqual(end, 45)
        self.assertEqual(motif, "CAG")
        self.assertTrue(adjusted)

        # Test 3: Extension into right flank
        reference = {"chr1": "AAAA" + "CAG"*10 + "CAG"*5}
        mock_fasta = MockFasta(reference)
        start, end, motif, adjusted = adjust_motif_and_boundaries_to_maximize_purity(
            mock_fasta, "chr1", 4, 34, "CAG"
        )
        self.assertEqual(start, 4)
        self.assertEqual(end, 49)  # Extended to include right flank repeats
        self.assertEqual(motif, "CAG")
        self.assertTrue(adjusted)

        # Test 4: Extension into both flanks
        reference = {"chr1": "CAG"*3 + "CAG"*10 + "CAG"*3}
        mock_fasta = MockFasta(reference)
        start, end, motif, adjusted = adjust_motif_and_boundaries_to_maximize_purity(
            mock_fasta, "chr1", 9, 39, "CAG"
        )
        self.assertEqual(start, 0)  # Extended left
        self.assertEqual(end, 48)   # Extended right
        self.assertEqual(motif, "CAG")
        self.assertTrue(adjusted)

        # Test 5: Near chromosome start boundary (actual_left_flank_size < requested)
        reference = {"chr1": "CA" + "CAG"*10 + "TTTT"}
        mock_fasta = MockFasta(reference)
        start, end, motif, adjusted = adjust_motif_and_boundaries_to_maximize_purity(
            mock_fasta, "chr1", 2, 32, "CAG"
        )
        self.assertEqual(start, 2)  # Stays at position 2 (can't go to negative)
        self.assertEqual(end, 32)
        self.assertEqual(motif, "CAG")
        self.assertFalse(adjusted)

        # Test 6: Near chromosome end boundary (actual_right_flank_size < requested)
        reference = {"chr1": "AAAA" + "CAG"*10 + "TT"}
        mock_fasta = MockFasta(reference)
        start, end, motif, adjusted = adjust_motif_and_boundaries_to_maximize_purity(
            mock_fasta, "chr1", 4, 34, "CAG"
        )
        self.assertEqual(start, 4)
        self.assertEqual(end, 34)
        self.assertEqual(motif, "CAG")
        self.assertFalse(adjusted)

        # Test 7: Locus too small (less than motif length)
        reference = {"chr1": "AAAA" + "CA" + "TTTT"}
        mock_fasta = MockFasta(reference)
        start, end, motif, adjusted = adjust_motif_and_boundaries_to_maximize_purity(
            mock_fasta, "chr1", 4, 6, "CAG"
        )
        self.assertEqual(start, 4)  # Returns unchanged
        self.assertEqual(end, 6)
        self.assertEqual(motif, "CAG")
        self.assertFalse(adjusted)

        # Test 8: Motif refinement - most common motif differs from input
        # Sequence has mostly CAGCAG but one CAACAG interruption
        reference = {"chr1": "AAAA" + "CAG"*5 + "CAA" + "CAG"*4 + "TTTT"}
        mock_fasta = MockFasta(reference)
        start, end, motif, adjusted = adjust_motif_and_boundaries_to_maximize_purity(
            mock_fasta, "chr1", 4, 34, "CAG"
        )
        self.assertEqual(start, 4)
        self.assertEqual(end, 34)
        self.assertEqual(motif, "CAG")  # Should still be CAG (most common)
        # adjusted might be True or False depending on whether motif was refined

        # Test 9: Extension with motif that's already at chromosome start
        reference = {"chr1": "CAG"*20}
        mock_fasta = MockFasta(reference)
        start, end, motif, adjusted = adjust_motif_and_boundaries_to_maximize_purity(
            mock_fasta, "chr1", 30, 45, "CAG"
        )
        self.assertEqual(start, 0)  # Extended all the way to start
        self.assertTrue(adjusted)

        # Test 10: Zero-length right flank case (at chromosome end)
        reference = {"chr1": "AAAA" + "CAG"*10}
        mock_fasta = MockFasta(reference)
        start, end, motif, adjusted = adjust_motif_and_boundaries_to_maximize_purity(
            mock_fasta, "chr1", 4, 34, "CAG"
        )
        self.assertEqual(start, 4)
        self.assertEqual(end, 34)
        self.assertEqual(motif, "CAG")


