import unittest

from str_analysis.utils.find_repeat_unit import find_repeat_unit


class Test(unittest.TestCase):

    def test_compute_most_frequent_repeat_unit(self):

        for repeat_unit in "A", "AC", "ACG", "ACGT", "ACGTA":
            for num_repeats in 1, 2, 3, 4, 5:
                if len(repeat_unit) == 1 and num_repeats == 1:
                    self.assertRaises(ValueError, find_repeat_unit, repeat_unit)
                    continue

                found_repeat_unit, found_num_repeats = find_repeat_unit(
                    repeat_unit * num_repeats, min_fraction_covered_by_repeat=0.95)
                self.assertEqual(found_repeat_unit, repeat_unit,
                                 f"repeat unit not detected correctly for {repeat_unit}*{num_repeats}")
                self.assertEqual(num_repeats, num_repeats,
                                 f"repeat count not detected correctly for {repeat_unit}*{num_repeats}")
