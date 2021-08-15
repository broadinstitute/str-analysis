import unittest

from str_analysis.utils.canonical_repeat_unit import reverse_complement, compute_canonical_repeat_unit, \
    _minimal_unit_under_shift


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

    def test_minimal_unit_under_shift(self):
        self.assertEqual(_minimal_unit_under_shift("C"), "C")
        self.assertEqual(_minimal_unit_under_shift("TAA"), "AAT")
        self.assertEqual(_minimal_unit_under_shift("ACA"), "AAC")

    def test_compute_canonical_repeat_unit(self):
        self.assertEqual(compute_canonical_repeat_unit("G"), "C")
        self.assertEqual(compute_canonical_repeat_unit("N"), "N")
        self.assertEqual(compute_canonical_repeat_unit("T"), "A")
        self.assertEqual(compute_canonical_repeat_unit("TGAG"), "ACTC")
        self.assertEqual(compute_canonical_repeat_unit("G"*9), "C"*9)
