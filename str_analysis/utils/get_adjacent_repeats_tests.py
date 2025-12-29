#!/usr/bin/env python3

"""Comprehensive tests for get_adjacent_repeats.py"""

import unittest
from unittest import mock
import intervaltree

from str_analysis.utils.get_adjacent_repeats import (
    get_repeat_unit_from_fasta,
    compute_locus_id,
    get_adjacent_repeats,
    MAX_DISTANCE_BETWEEN_REPEATS,
    MAX_TOTAL_ADJACENT_REGION_SIZE,
    MAX_OVERLAP_BETWEEN_ADJACENT_REPEATS,
)


class TestGetRepeatUnitFromFasta(unittest.TestCase):
    """Test get_repeat_unit_from_fasta function."""

    def test_get_repeat_unit_from_fasta(self):
        """Test extracting repeat unit from FASTA."""
        mock_fasta = mock.Mock()
        mock_fasta.fetch.return_value = "CAGCAGCAGCAGCAG"

        repeat_unit = get_repeat_unit_from_fasta("chr1", 1001, 1015, 3, mock_fasta)

        self.assertEqual(repeat_unit, "CAG")
        mock_fasta.fetch.assert_called_once_with("chr1", 1000, 1015)  # 0-based fetch

    def test_get_repeat_unit_uppercase_conversion(self):
        """Test that lowercase sequences are converted to uppercase."""
        mock_fasta = mock.Mock()
        mock_fasta.fetch.return_value = "cagcagcag"

        repeat_unit = get_repeat_unit_from_fasta("chr2", 2001, 2009, 3, mock_fasta)

        self.assertEqual(repeat_unit, "CAG")

    def test_get_repeat_unit_different_lengths(self):
        """Test with different repeat unit lengths."""
        mock_fasta = mock.Mock()

        # 2bp repeat
        mock_fasta.fetch.return_value = "ATATAT"
        repeat_unit = get_repeat_unit_from_fasta("chr3", 3001, 3006, 2, mock_fasta)
        self.assertEqual(repeat_unit, "AT")

        # 4bp repeat
        mock_fasta.fetch.return_value = "AAAAGAAAAG"
        repeat_unit = get_repeat_unit_from_fasta("chr3", 3001, 3010, 4, mock_fasta)
        self.assertEqual(repeat_unit, "AAAA")

    def test_get_repeat_unit_mixed_case(self):
        """Test with mixed case sequences."""
        mock_fasta = mock.Mock()
        mock_fasta.fetch.return_value = "CaGcAgCaG"

        repeat_unit = get_repeat_unit_from_fasta("chr4", 4001, 4009, 3, mock_fasta)

        self.assertEqual(repeat_unit, "CAG")


class TestComputeLocusId(unittest.TestCase):
    """Test compute_locus_id function."""

    def test_compute_locus_id_basic(self):
        """Test basic locus ID computation."""
        locus_id = compute_locus_id("chr1", 1000, 1015, "CAG")
        self.assertEqual(locus_id, "chr1-1000-1015-CAG")

    def test_compute_locus_id_different_chroms(self):
        """Test locus ID with different chromosomes."""
        locus_id1 = compute_locus_id("chr1", 100, 200, "AT")
        locus_id2 = compute_locus_id("chr22", 100, 200, "AT")

        self.assertEqual(locus_id1, "chr1-100-200-AT")
        self.assertEqual(locus_id2, "chr22-100-200-AT")
        self.assertNotEqual(locus_id1, locus_id2)

    def test_compute_locus_id_different_motifs(self):
        """Test locus ID with different motifs."""
        locus_id1 = compute_locus_id("chr1", 1000, 1015, "CAG")
        locus_id2 = compute_locus_id("chr1", 1000, 1015, "CTG")

        self.assertNotEqual(locus_id1, locus_id2)

    def test_compute_locus_id_format(self):
        """Test locus ID format components."""
        locus_id = compute_locus_id("chr10", 5000, 5100, "AAAG")

        self.assertIn("chr10", locus_id)
        self.assertIn("5000", locus_id)
        self.assertIn("5100", locus_id)
        self.assertIn("AAAG", locus_id)


class TestGetAdjacentRepeats(unittest.TestCase):
    """Test get_adjacent_repeats function."""

    def setUp(self):
        """Set up mock FASTA and interval tree."""
        self.mock_fasta = mock.Mock()

        # Default behavior: return the motif repeated
        def fetch_side_effect(chrom, start, end):
            length = end - start
            if length <= 3:
                return "CAG" * ((length // 3) + 1)
            elif length <= 4:
                return "AAAG" * ((length // 4) + 1)
            else:
                return "ACGT" * ((length // 4) + 1)

        self.mock_fasta.fetch.side_effect = fetch_side_effect

    def test_get_adjacent_repeats_no_adjacent(self):
        """Test when there are no adjacent repeats."""
        # Create empty interval tree (no adjacent repeats)
        tree = intervaltree.IntervalTree()

        left, structure, right = get_adjacent_repeats(
            "chr1:1000-1015",
            "CAG",
            self.mock_fasta,
            tree
        )

        self.assertEqual(len(left), 0)
        self.assertEqual(structure, "(CAG)*")
        self.assertEqual(len(right), 0)

    def test_get_adjacent_repeats_with_left_adjacent(self):
        """Test when there's one adjacent repeat on the left."""
        # Create interval tree with one left adjacent repeat
        tree = intervaltree.IntervalTree()
        tree.add(intervaltree.Interval(990, 997, "CTG"))  # Left adjacent

        # Mock fetch to return different sequences
        def custom_fetch(chrom, start, end):
            if start < 995:
                return "CTG" * 5
            else:
                return "CAG" * 5

        self.mock_fasta.fetch.side_effect = custom_fetch

        left, structure, right = get_adjacent_repeats(
            "chr1:1000-1015",
            "CAG",
            self.mock_fasta,
            tree
        )

        self.assertEqual(len(left), 1)
        self.assertIn("CTG", left[0])
        self.assertIn("CTG", structure)
        self.assertIn("CAG", structure)
        self.assertEqual(len(right), 0)

    def test_get_adjacent_repeats_with_right_adjacent(self):
        """Test when there's one adjacent repeat on the right."""
        tree = intervaltree.IntervalTree()
        tree.add(intervaltree.Interval(1016, 1023, "CTG"))  # Right adjacent

        def custom_fetch(chrom, start, end):
            if start > 1015:
                return "CTG" * 5
            else:
                return "CAG" * 5

        self.mock_fasta.fetch.side_effect = custom_fetch

        left, structure, right = get_adjacent_repeats(
            "chr1:1000-1015",
            "CAG",
            self.mock_fasta,
            tree
        )

        self.assertEqual(len(left), 0)
        self.assertEqual(len(right), 1)
        self.assertIn("CTG", right[0])
        self.assertIn("CAG", structure)
        self.assertIn("CTG", structure)

    def test_get_adjacent_repeats_with_both_sides(self):
        """Test with adjacent repeats on both left and right."""
        tree = intervaltree.IntervalTree()
        tree.add(intervaltree.Interval(990, 997, "CTG"))   # Left adjacent
        tree.add(intervaltree.Interval(1016, 1023, "AAG"))  # Right adjacent

        def custom_fetch(chrom, start, end):
            if start < 995:
                return "CTG" * 5
            elif start > 1015:
                return "AAG" * 5
            else:
                return "CAG" * 5

        self.mock_fasta.fetch.side_effect = custom_fetch

        left, structure, right = get_adjacent_repeats(
            "chr1:1000-1015",
            "CAG",
            self.mock_fasta,
            tree
        )

        self.assertEqual(len(left), 1)
        self.assertEqual(len(right), 1)
        self.assertIn("CTG", structure)
        self.assertIn("CAG", structure)
        self.assertIn("AAG", structure)

    def test_get_adjacent_repeats_max_distance_filter(self):
        """Test that repeats beyond max distance are filtered out."""
        tree = intervaltree.IntervalTree()
        # Add a repeat far on the left (beyond max distance)
        tree.add(intervaltree.Interval(900, 907, "CTG"))  # Too far left

        left, structure, right = get_adjacent_repeats(
            "chr1:1000-1015",
            "CAG",
            self.mock_fasta,
            tree,
            max_distance_between_adjacent_repeats=10  # Only 10bp allowed
        )

        # Should not include the far left repeat
        self.assertEqual(len(left), 0)
        self.assertEqual(structure, "(CAG)*")

    def test_get_adjacent_repeats_max_adjacent_repeats_limit(self):
        """Test max_adjacent_repeats parameter limits results."""
        tree = intervaltree.IntervalTree()
        # Add multiple left adjacent repeats
        tree.add(intervaltree.Interval(980, 987, "CTG"))
        tree.add(intervaltree.Interval(990, 997, "AAG"))
        tree.add(intervaltree.Interval(970, 977, "ATG"))

        def custom_fetch(chrom, start, end):
            if 970 <= start < 980:
                return "ATG" * 5
            elif 980 <= start < 990:
                return "CTG" * 5
            elif 990 <= start < 1000:
                return "AAG" * 5
            else:
                return "CAG" * 5

        self.mock_fasta.fetch.side_effect = custom_fetch

        left, structure, right = get_adjacent_repeats(
            "chr1:1000-1015",
            "CAG",
            self.mock_fasta,
            tree,
            max_adjacent_repeats=2  # Limit to 2
        )

        # Should only get 2 adjacent repeats on the left
        self.assertLessEqual(len(left), 2)

    def test_get_adjacent_repeats_duplicate_motif_stops_search(self):
        """Test that duplicate motifs stop the search."""
        tree = intervaltree.IntervalTree()
        tree.add(intervaltree.Interval(990, 997, "CTG"))
        tree.add(intervaltree.Interval(980, 987, "CAG"))  # Same as main repeat!

        def custom_fetch(chrom, start, end):
            if 980 <= start < 990:
                return "CAG" * 5
            elif 990 <= start < 1000:
                return "CTG" * 5
            else:
                return "CAG" * 5

        self.mock_fasta.fetch.side_effect = custom_fetch

        left, structure, right = get_adjacent_repeats(
            "chr1:1000-1015",
            "CAG",
            self.mock_fasta,
            tree
        )

        # Should stop at the duplicate CAG motif
        self.assertEqual(len(left), 1)  # Only CTG, not the second CAG
        self.assertIn("CTG", structure)

    def test_get_adjacent_repeats_with_spacer_sequence(self):
        """Test that spacer sequences are included in locus structure."""
        tree = intervaltree.IntervalTree()
        tree.add(intervaltree.Interval(990, 996, "CTG"))  # Ends at 996, main starts at 1000

        def custom_fetch(chrom, start, end):
            if start == 996 and end == 1000:  # Spacer region between 996 and 1000
                return "ACGT"
            elif start < 996:
                return "CTG" * 5
            else:
                return "CAG" * 5

        self.mock_fasta.fetch.side_effect = custom_fetch

        left, structure, right = get_adjacent_repeats(
            "chr1:1000-1015",
            "CAG",
            self.mock_fasta,
            tree,
            max_distance_between_adjacent_repeats=10
        )

        # Locus structure should include the spacer
        self.assertIn("ACGT", structure)

    def test_get_adjacent_repeats_overlapping_trimmed(self):
        """Test that overlapping adjacent repeats are trimmed."""
        tree = intervaltree.IntervalTree()
        # Add a repeat that overlaps with the main repeat
        tree.add(intervaltree.Interval(995, 1002, "CTG"))

        def custom_fetch(chrom, start, end):
            if start < 1000:
                return "CTG" * 5
            else:
                return "CAG" * 5

        self.mock_fasta.fetch.side_effect = custom_fetch

        left, structure, right = get_adjacent_repeats(
            "chr1:1000-1015",
            "CAG",
            self.mock_fasta,
            tree
        )

        # Should still find the adjacent repeat but with trimmed coordinates
        self.assertEqual(len(left), 1)

    def test_get_adjacent_repeats_verbose_warning(self):
        """Test verbose warning when motif differs from reference."""
        tree = intervaltree.IntervalTree()

        # Mock fetch to return different motif than specified
        def custom_fetch(chrom, start, end):
            return "CTG" * 10  # Different from CAG

        self.mock_fasta.fetch.side_effect = custom_fetch

        with mock.patch('builtins.print') as mock_print:
            left, structure, right = get_adjacent_repeats(
                "chr1:1000-1015",
                "CAG",
                self.mock_fasta,
                tree,
                verbose=True
            )

            # Should print a warning
            mock_print.assert_called()
            warning_call = mock_print.call_args[0][0]
            self.assertIn("WARNING", warning_call)
            self.assertIn("differs from", warning_call)

    def test_get_adjacent_repeats_no_verbose_warning(self):
        """Test no warning when verbose=False."""
        tree = intervaltree.IntervalTree()

        def custom_fetch(chrom, start, end):
            return "CTG" * 10

        self.mock_fasta.fetch.side_effect = custom_fetch

        with mock.patch('builtins.print') as mock_print:
            left, structure, right = get_adjacent_repeats(
                "chr1:1000-1015",
                "CAG",
                self.mock_fasta,
                tree,
                verbose=False  # No warning expected
            )

            # Should not print anything
            mock_print.assert_not_called()


class TestConstants(unittest.TestCase):
    """Test that module constants have expected values."""

    def test_max_distance_between_repeats(self):
        """Test MAX_DISTANCE_BETWEEN_REPEATS constant."""
        self.assertEqual(MAX_DISTANCE_BETWEEN_REPEATS, 6)
        self.assertIsInstance(MAX_DISTANCE_BETWEEN_REPEATS, int)

    def test_max_total_adjacent_region_size(self):
        """Test MAX_TOTAL_ADJACENT_REGION_SIZE constant."""
        self.assertEqual(MAX_TOTAL_ADJACENT_REGION_SIZE, 1000)
        self.assertIsInstance(MAX_TOTAL_ADJACENT_REGION_SIZE, int)

    def test_max_overlap_between_adjacent_repeats(self):
        """Test MAX_OVERLAP_BETWEEN_ADJACENT_REPEATS constant."""
        self.assertEqual(MAX_OVERLAP_BETWEEN_ADJACENT_REPEATS, 3)
        self.assertIsInstance(MAX_OVERLAP_BETWEEN_ADJACENT_REPEATS, int)


if __name__ == "__main__":
    unittest.main()
