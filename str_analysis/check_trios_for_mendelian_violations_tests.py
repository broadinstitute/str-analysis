from intervaltree import Interval
import unittest

from check_trios_for_mendelian_violations import compute_mendelian_violations, \
    check_for_duplicate_keys, group_rows_by_trio, \
    intervals_overlap, \
    compute_min_distance_mendelian_ci, compute_min_distance_mendelian


class Tests(unittest.TestCase):

    def test_intervals_overlap(self):
        self.assertTrue(intervals_overlap(Interval(1,1), Interval(1,1)))
        self.assertTrue(intervals_overlap(Interval(1,1), Interval(1,2)))
        self.assertTrue(intervals_overlap(Interval(2, 3), Interval(1,2)))

        self.assertFalse(intervals_overlap(Interval(1,1), Interval(2,3)))
        self.assertFalse(intervals_overlap(Interval(1,1), Interval(2,3)))

    def test_compute_min_distance_mendelian(self):
        self.assertEqual(compute_min_distance_mendelian("23", ["11", "24"]), 1)
        self.assertEqual(compute_min_distance_mendelian("23", ["11", "10"]), 12)

    def test_compute_min_distance_mendelian_ci(self):
        self.assertEqual(compute_min_distance_mendelian_ci(Interval(20, 20), [Interval(21, 21), Interval(24, 30)]), 1)
        self.assertEqual(compute_min_distance_mendelian_ci(Interval(20, 20), [Interval(25, 50), Interval(24, 30)]), 4)
        self.assertEqual(compute_min_distance_mendelian_ci(Interval(20, 30), [Interval(25, 50), Interval(24, 30)]), 0)
        self.assertEqual(compute_min_distance_mendelian_ci(Interval(20, 30), [Interval(32, 50), Interval(43, 50)]), 2)

    def test_compute_mendelian_violations(self):
        pass

    def check_for_duplicate_keys(self):
        pass

    def group_rows_by_trio(self):
        pass
