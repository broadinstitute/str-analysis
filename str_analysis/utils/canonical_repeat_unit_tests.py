import unittest

from str_analysis.utils.canonical_repeat_unit import (
    _alphabetically_first_motif_under_shift, compute_canonical_motif,
    reverse_complement)


class Tests(unittest.TestCase):

    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("A"), "T")
        self.assertEqual(reverse_complement("C"), "G")
        self.assertEqual(reverse_complement("G"), "C")
        self.assertEqual(reverse_complement("T"), "A")
        self.assertEqual(reverse_complement("N"), "N")
        self.assertEqual(reverse_complement("CG"), "CG")
        self.assertEqual(reverse_complement("TA"), "TA")
        self.assertEqual(reverse_complement("GC"), "GC")
        self.assertEqual(reverse_complement("AT"), "AT")
        self.assertEqual(reverse_complement("CA"), "TG")
        self.assertEqual(reverse_complement("GG"), "CC")
        self.assertEqual(reverse_complement("G"*2), "C"*2)
        self.assertEqual(reverse_complement("G"*3), "C"*3)
        self.assertEqual(reverse_complement("G"*10), "C"*10)
        self.assertRaises(KeyError, lambda: reverse_complement("X"))

    def test_alphabetically_first_motif_under_shift(self):
        self.assertEqual(_alphabetically_first_motif_under_shift("C"), "C")
        self.assertEqual(_alphabetically_first_motif_under_shift("TAA"), "AAT")
        self.assertEqual(_alphabetically_first_motif_under_shift("ACA"), "AAC")

    def test_compute_canonical_motif(self):
        self.assertEqual(compute_canonical_motif("G"), "C")
        self.assertEqual(compute_canonical_motif("N"), "N")
        self.assertEqual(compute_canonical_motif("T"), "A")
        self.assertEqual(compute_canonical_motif("TGAG"), "ACTC")
        self.assertEqual(compute_canonical_motif("G"*9), "C"*9)

