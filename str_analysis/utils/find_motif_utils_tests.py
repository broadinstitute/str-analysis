import unittest

from str_analysis.utils.find_motif_utils import compute_repeat_purity


class Tests(unittest.TestCase):

    def test_compute_repeat_purity(self):
        tests = [
            ("CACACACACA",      "CA",    0, 1,         False),
            ("CACACACACAC",     "CA",    0, 1,         False),
            ("CACACTGTGT",      "CA",    5, 0.5,       False),
            ("CA",              "CTTG",  0, 0,         False),
            ("CAGCAACAGCAGTT",  "CAG",   1, 11/12,     False),
            ("CAGCAACAGCAGTT",  "CAG",   3, 11/14,     True),
        ]

        for sequence, motif, expected_interruptions, expected_purity, include_partial_repeats in tests:
            purity, interruption_count = compute_repeat_purity(sequence, motif, include_partial_repeats=include_partial_repeats)
            self.assertEqual(expected_interruptions, interruption_count, msg=f"Expected interruption count {expected_interruptions} but got {interruption_count} for sequence {sequence} & motif {motif}")
            self.assertAlmostEquals(expected_purity, purity, msg=f"Expected base purity {expected_purity} but got {purity} for sequence {sequence} & motif {motif}")
