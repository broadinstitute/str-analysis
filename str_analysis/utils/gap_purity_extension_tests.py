import unittest

from str_analysis.utils.gap_purity_extension import extend_by_gap_purity, cumulative_mismatched_bases


class TestGapPurityExtension(unittest.TestCase):

    def test_pure_repeat_is_taken_whole(self):
        # 5 clean copies in the flank, nothing else
        self.assertEqual(extend_by_gap_purity("CATCATCATCATCAT", "CAT", 9, "right"), 15)

    def test_random_flank_is_not_extended(self):
        self.assertEqual(extend_by_gap_purity("GGGGTTAACCAGTGA", "CAT", 9, "right"), 0)

    def test_empty_flank_is_not_extended(self):
        self.assertEqual(extend_by_gap_purity("", "CAT", 9, "right"), 0)

    def test_gap_is_crossed_when_enough_clean_copies_follow(self):
        # one interrupted copy (CAG) then 3 clean ones: 12bp with 1 mismatch = 0.917 purity
        self.assertEqual(extend_by_gap_purity("CAGCATCATCAT", "CAT", 9, "right"), 12)

    def test_gap_is_refused_when_too_few_clean_copies_follow(self):
        # one interrupted copy then only 1 clean one: 6bp with 1 mismatch = 0.833, below the default
        self.assertEqual(extend_by_gap_purity("CAGCATTATGGG", "CAT", 9, "right"), 0)

    def test_lowering_the_purity_threshold_accepts_that_same_gap(self):
        self.assertEqual(
            extend_by_gap_purity("CAGCATTATGGG", "CAT", 9, "right", min_purity_of_new_sequence=0.83), 6)

    def test_max_repeats_in_gap_zero_requires_an_immediate_exact_copy(self):
        self.assertEqual(
            extend_by_gap_purity("CAGCATCATCAT", "CAT", 9, "right", max_repeats_in_gap=0), 0)
        self.assertEqual(
            extend_by_gap_purity("CATCATCATCAG", "CAT", 9, "right", max_repeats_in_gap=0), 9)

    def test_gap_longer_than_the_budget_stops_the_extension(self):
        # 3 interrupted copies exceeds the default budget of 2. The clean run afterwards is made long
        # enough that purity is NOT the binding constraint (99bp with 3 mismatches is 0.970), so this
        # would fail if the budget were ignored, which is what the test is here to check.
        flank = "CAG" * 3 + "CAT" * 30
        self.assertEqual(extend_by_gap_purity(flank, "CAT", 9, "right"), 0)
        self.assertEqual(extend_by_gap_purity(flank, "CAT", 9, "right", max_repeats_in_gap=3), 99)

    def test_extension_continues_past_an_accepted_chunk(self):
        # two acceptable chunks in a row, each an interrupted copy plus 3 clean ones
        self.assertEqual(extend_by_gap_purity("CAGCATCATCATCAGCATCATCAT", "CAT", 9, "right"), 24)

    def test_left_flank_is_read_outward_from_the_boundary(self):
        # the caller reverses the left flank, so the boundary-adjacent base comes first. Genomic
        # sequence ...CATCATCAT|core, reversed, should extend by all 9bp.
        self.assertEqual(extend_by_gap_purity("CATCATCAT"[::-1], "CAT", 9, "left"), 9)

    def test_phase_depends_on_the_core_length(self):
        # A core of 9bp ends on a copy boundary, so the right flank resumes at the motif's first base.
        self.assertEqual(extend_by_gap_purity("CATCATCAT", "CAT", 9, "right"), 9)
        # A core of 10bp ends one base into a copy, so the flank must resume at the motif's 2nd base.
        self.assertEqual(extend_by_gap_purity("CATCATCAT", "CAT", 10, "right"), 0)
        self.assertEqual(extend_by_gap_purity("ATCATCATC", "CAT", 10, "right"), 9)

    def test_iupac_codes_in_the_motif_match_any_base_they_represent(self):
        # N matches anything, so GCN tiles cleanly over GCA GCT GCC GCG and all 12bp are taken
        self.assertEqual(extend_by_gap_purity("GCAGCTGCCGCG", "GCN", 9, "right"), 12)
        # against a literal motif only the first copy matches, and the 3 that follow exceed the gap
        # budget, so the extension stops after 3bp
        self.assertEqual(extend_by_gap_purity("GCAGCTGCCGCG", "GCA", 9, "right"), 3)

    def test_cumulative_mismatched_bases_counts_hamming_distance(self):
        mismatches = cumulative_mismatched_bases("CATCAG", "CAT", 9, "right")
        self.assertEqual(list(mismatches), [0, 0, 0, 0, 0, 1])

    def test_side_must_be_left_or_right(self):
        with self.assertRaises(ValueError):
            extend_by_gap_purity("CATCATCAT", "CAT", 9, "middle")


if __name__ == "__main__":
    unittest.main()
