import unittest

from check_trios_for_mendelian_violations import (
    compute_min_distance_mendelian, compute_min_distance_mendelian_ci,
    intervals_overlap, determine_transmitted_alleles)
from intervaltree import Interval


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

    def test_determine_transmitted_alleles(self):
        father_allele, mother_allele, proband_allele_from_father, proband_allele_from_mother = determine_transmitted_alleles(
            ["5", "10"], father_alleles=["3", "5"], mother_alleles=["9", "14"], is_chrX_locus=False)
        self.assertListEqual(
            [father_allele, mother_allele, proband_allele_from_father, proband_allele_from_mother],
            ["5", "9", "5", "10"])

        father_allele, mother_allele, proband_allele_from_father, proband_allele_from_mother = determine_transmitted_alleles(
            ["5", "10"], mother_alleles=["3", "5"], father_alleles=["9", "14"], is_chrX_locus=False)
        self.assertListEqual(
            [father_allele, mother_allele, proband_allele_from_father, proband_allele_from_mother],
            ["9", "5", "10", "5"])

    def test_compute_mendelian_violations(self):
        pass

    def check_for_duplicate_keys(self):
        pass

    def group_rows_by_trio(self):
        pass
