import unittest

from str_analysis.annotate_and_filter_str_catalog import compute_sequence_purity_stats
from intervaltree import Interval


class Tests(unittest.TestCase):

    def test_compute_sequence_purity_stats(self):
        tests = [
            ("CACACACACA",  "CA",   0, 1, 1),
            ("CACACACACAC", "CA",   0, 1, 1),
            ("CACACTGTGT",  "CA",   5, 0.5, 0.4),
            ("CA",          "CTTG", 0, 0, 0),
        ]

        for sequence, motif, i, pb, pr in tests:
            i_, pb_, pr_ = compute_sequence_purity_stats(sequence, motif)
            self.assertEqual(i, i_, msg=f"Expected interruption count {i} but got {i_} for sequence {sequence} & motif {motif}")
            self.assertAlmostEquals(pb, pb_, msg=f"Expected base purity {pb} but got {pb_} for sequence {sequence} & motif {motif}")
            self.assertAlmostEquals(pr, pr_, msg=f"Expected repeat purity {pr} but got {pr_} for sequence {sequence} & motif {motif}")
