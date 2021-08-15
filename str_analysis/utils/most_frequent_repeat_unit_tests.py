import unittest
from str_analysis.utils.most_frequent_repeat_unit import compute_most_frequent_repeat_unit


class Test(unittest.TestCase):

    def test_compute_most_frequent_repeat_unit(self):

        self.assertEqual(
            compute_most_frequent_repeat_unit("CACAAA", repeat_unit_size=1, min_fraction_bases_covered=0.8),
            (None, 0))
        self.assertEqual(
            compute_most_frequent_repeat_unit("CACAAAAAAA", repeat_unit_size=1, min_fraction_bases_covered=0.8),
            ("A", 8))
        self.assertEqual(
            compute_most_frequent_repeat_unit("CAG", repeat_unit_size=3, min_occurrences=1, min_fraction_bases_covered=0.8),
            ("CAG", 1))
        self.assertEqual(
            compute_most_frequent_repeat_unit("TCAGCAGCAGC", repeat_unit_size=3, min_fraction_bases_covered=0.8),
            ("CAG", 3))
        self.assertEqual(
            compute_most_frequent_repeat_unit("TCAGCAGCAGC", repeat_unit_size=3, min_fraction_bases_covered=0.8),
            ("CAG", 3))
        self.assertEqual(
            compute_most_frequent_repeat_unit("GGAAGGGAAGGGAAGGGAAGGG", repeat_unit_size=5, min_fraction_bases_covered=0.8),
            ("GGAAG", 4))

        sequence = "AGGGT"+"AAAAG"*8+"AGGGG"
        self.assertEqual(
            compute_most_frequent_repeat_unit(sequence, repeat_unit_size=5, min_fraction_bases_covered=0.8),
            ("AAAAG", 8))
