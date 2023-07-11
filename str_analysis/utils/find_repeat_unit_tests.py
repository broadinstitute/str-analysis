import unittest

from str_analysis.utils.find_repeat_unit import get_repeat_unit_regex_with_N_base, get_most_common_repeat_unit, \
    find_repeat_unit_without_allowing_interruptions, find_repeat_unit_allowing_interruptions, \
    extend_repeat_into_sequence_allowing_interruptions, extend_repeat_into_sequence_without_allowing_interruptions
from str_analysis.utils.find_repeat_unit import count_pure_repeats


def add_interruption(sequence, motif, i):
    """Insert interruption in last repeat at position i within the motif"""
    assert 0 <= i < len(motif)

    unused_bases = list(sorted(set("ACGT") - set(motif)))
    i = len(sequence) - len(motif) + i
    new_base = unused_bases[0] if unused_bases else ("C" if sequence[i] == "A" else "A")
    return sequence[:i] + new_base + sequence[i+1:]


class Test(unittest.TestCase):

    def test_add_interruption(self):
        sequence_with_interruption = add_interruption("CAG"*5, "CAG", 1)
        self.assertEqual(sequence_with_interruption, "CAG"*4 + "CTG")

    def test_get_repeat_unit_regex_with_N_base(self):
        self.assertRaises(ValueError, lambda: get_repeat_unit_regex_with_N_base("C", 0, allow_homopolymer=True))
        self.assertRaises(ValueError, lambda: get_repeat_unit_regex_with_N_base("CA", 2, allow_homopolymer=True))

        self.assertEqual(get_repeat_unit_regex_with_N_base("CA", 0, allow_homopolymer=False), "[CGT]A")
        self.assertEqual(get_repeat_unit_regex_with_N_base("CA", 1, allow_homopolymer=False), "C[AGT]")
        self.assertEqual(get_repeat_unit_regex_with_N_base("CA", 1, allow_homopolymer=True), "C[ACGT]")

        self.assertEqual(get_repeat_unit_regex_with_N_base("CAG", 0, allow_homopolymer=True), "[ACGT]AG")
        self.assertEqual(get_repeat_unit_regex_with_N_base("CAG", 1, allow_homopolymer=True), "C[ACGT]G")
        self.assertEqual(get_repeat_unit_regex_with_N_base("CAG", 2, allow_homopolymer=True), "CA[ACGT]")

        self.assertEqual(get_repeat_unit_regex_with_N_base("CACG", 3, allow_homopolymer=False), "CAC[ACGT]")
        self.assertEqual(get_repeat_unit_regex_with_N_base("CACG", 3, allow_homopolymer=True), "CAC[ACGT]")

    def test_get_most_common_repeat_unit(self):

        # Test get_most_common_repeat_unit
        self.assertEqual(get_most_common_repeat_unit("CAG"*5, repeat_unit_length=3), ("CAG", 5))
        self.assertEqual(get_most_common_repeat_unit("CAG"*5+"CAT"*5+"CAC"*5+"CAA"*5+"CAG", repeat_unit_length=3),
                         ("CAG", 6))
        self.assertEqual(get_most_common_repeat_unit("TAAAATTAAA", repeat_unit_length=1), ("A", 7))

    def test_find_repeat_unit_without_allowing_interruptions(self):
        for repeat_unit in "A", "AC", "ACG", "ACGT", "ACGTA", "ACGATA":
            for num_repeats in 1, 2, 3, 4, 5:
                for allow_partial_repeats in True, False:
                    sequence = repeat_unit * num_repeats
                    if allow_partial_repeats and len(repeat_unit) > 1 and num_repeats > 1:
                        sequence += repeat_unit[0]
                    found_repeat_unit, found_num_pure_repeats, has_partial_repeats = find_repeat_unit_without_allowing_interruptions(sequence, allow_partial_repeats=allow_partial_repeats)

                    self.assertEqual(found_repeat_unit, repeat_unit,
                                     f"repeat unit not detected correctly for {num_repeats}x{repeat_unit}")
                    self.assertEqual(found_num_pure_repeats, num_repeats,
                                     f"repeat count not detected correctly for {num_repeats}x{repeat_unit}")

    def test_find_repeat_unit_with_interruptions(self):

        # test a few examples
        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, interruption_index, has_partial_repeats = find_repeat_unit_allowing_interruptions(
            "A"*10 + "C" + "A"*20)
        self.assertEqual(found_repeat_unit, "A")
        self.assertEqual(found_num_pure_repeats, 30)
        self.assertEqual(found_num_total_repeats, 31)
        self.assertEqual(interruption_index, 0)
        self.assertFalse(has_partial_repeats)

        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, interruption_index, has_partial_repeats = find_repeat_unit_allowing_interruptions(
            "A"*10 + "CC" + "A"*20)
        self.assertEqual(found_repeat_unit, "A"*10 + "CC" + "A"*20)
        self.assertEqual(found_num_pure_repeats, 1)
        self.assertEqual(found_num_total_repeats, 1)
        self.assertEqual(interruption_index, None)
        self.assertFalse(has_partial_repeats)

        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, interruption_index, has_partial_repeats = find_repeat_unit_allowing_interruptions(
            "A"*10 + "CT" + "A"*20)
        self.assertEqual(found_repeat_unit, "A"*10 + "CT" + "A"*20)
        self.assertEqual(found_num_pure_repeats, 1)
        self.assertEqual(found_num_total_repeats, 1)
        self.assertEqual(interruption_index, None)
        self.assertFalse(has_partial_repeats)

        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, interruption_index, has_partial_repeats = find_repeat_unit_allowing_interruptions(
            "CAG"*5)
        self.assertEqual(found_repeat_unit, "CAG")
        self.assertEqual(found_num_pure_repeats, 5)
        self.assertEqual(found_num_total_repeats, 5)
        self.assertEqual(interruption_index, None)
        self.assertFalse(has_partial_repeats)

        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, interruption_index, has_partial_repeats = find_repeat_unit_allowing_interruptions(
            "CAG"*5 + "C", allow_partial_repeats=True)
        self.assertEqual(found_repeat_unit, "CAG")
        self.assertEqual(found_num_pure_repeats, 5)
        self.assertEqual(found_num_total_repeats, 5)
        self.assertEqual(interruption_index, None)
        self.assertTrue(has_partial_repeats)

        found_repeat_unit, found_num_pure_repeats, has_partial_repeats = find_repeat_unit_without_allowing_interruptions(
            "CAG"*5 + "C", allow_partial_repeats=False)
        self.assertEqual(found_repeat_unit, "CAG"*5 + "C")
        self.assertEqual(found_num_pure_repeats, 1)
        self.assertFalse(has_partial_repeats)

        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, interruption_index, has_partial_repeats = find_repeat_unit_allowing_interruptions(
            "CATCAGCAGCAGCAG")
        self.assertEqual(found_repeat_unit, "CAG")
        self.assertEqual(found_num_pure_repeats, 4)
        self.assertEqual(found_num_total_repeats, 5)
        self.assertEqual(interruption_index, 2)
        self.assertFalse(has_partial_repeats)

        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, interruption_index, has_partial_repeats = find_repeat_unit_allowing_interruptions(
            "CAGCATCAGCAGCAG")
        self.assertEqual(found_repeat_unit, "CAG")
        self.assertEqual(found_num_pure_repeats, 4)
        self.assertEqual(found_num_total_repeats, 5)
        self.assertEqual(interruption_index, 2)
        self.assertFalse(has_partial_repeats)

        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, interruption_index, has_partial_repeats = find_repeat_unit_allowing_interruptions(
            "A"*5)
        self.assertEqual(found_repeat_unit, "A")
        self.assertEqual(found_num_pure_repeats, 5)
        self.assertEqual(found_num_total_repeats, 5)
        self.assertEqual(interruption_index, None)
        self.assertFalse(has_partial_repeats)

        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, interruption_index, has_partial_repeats = find_repeat_unit_allowing_interruptions(
            "ACG"*2 + "ACA" + "ACG"*2 + "ACT"*2)
        self.assertEqual(found_repeat_unit, "ACG")
        self.assertEqual(found_num_pure_repeats, 4)
        self.assertEqual(found_num_total_repeats, 7)
        self.assertEqual(interruption_index, 2)
        self.assertFalse(has_partial_repeats)

        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, interruption_index, has_partial_repeats = find_repeat_unit_allowing_interruptions(
            "AAAAAGAAAA")
        self.assertEqual(found_repeat_unit, "A")
        self.assertEqual(found_num_pure_repeats, 9)
        self.assertEqual(found_num_total_repeats, 10)
        self.assertEqual(interruption_index, 0)
        self.assertFalse(has_partial_repeats)

        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, interruption_index, has_partial_repeats = find_repeat_unit_allowing_interruptions(
            "CTCTCTCACT")
        self.assertEqual(found_repeat_unit, "CT")
        self.assertEqual(found_num_pure_repeats, 4)
        self.assertEqual(found_num_total_repeats, 5)
        self.assertEqual(interruption_index, 1)
        self.assertFalse(has_partial_repeats)

        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, interruption_index, has_partial_repeats = find_repeat_unit_allowing_interruptions(
            "AAAAAGAAA")
        self.assertEqual(found_repeat_unit, "AAAAAGAAA")
        self.assertEqual(found_num_pure_repeats, 1)
        self.assertEqual(found_num_total_repeats, 1)
        self.assertEqual(interruption_index, None)
        self.assertFalse(has_partial_repeats)

        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, interruption_index, has_partial_repeats = find_repeat_unit_allowing_interruptions(
            "AAGAAAAAGAAAAAT")
        self.assertEqual(found_repeat_unit, "AAGAAAAAGAAAAAT")
        self.assertEqual(found_num_pure_repeats, 1)
        self.assertEqual(found_num_total_repeats, 1)
        self.assertEqual(interruption_index, None)
        self.assertFalse(has_partial_repeats)

        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, interruption_index, has_partial_repeats = find_repeat_unit_allowing_interruptions(
            "CTTCGGCAT")
        self.assertEqual(found_repeat_unit, "CTTCGGCAT")
        self.assertEqual(found_num_pure_repeats, 1)
        self.assertEqual(found_num_total_repeats, 1)
        self.assertEqual(interruption_index, None)
        self.assertFalse(has_partial_repeats)
        self.assertFalse(has_partial_repeats)

        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, interruption_index, has_partial_repeats = find_repeat_unit_allowing_interruptions(
            "CATAATGAT")
        self.assertEqual(found_repeat_unit, "CAT")
        self.assertEqual(found_num_pure_repeats, 1)
        self.assertEqual(found_num_total_repeats, 3)
        self.assertEqual(interruption_index, 0)
        self.assertFalse(has_partial_repeats)

        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, interruption_index, has_partial_repeats = find_repeat_unit_allowing_interruptions(
            "CAGGCATGCAGGCACG")
        self.assertEqual(found_repeat_unit, "CAGG")
        self.assertEqual(found_num_pure_repeats, 2)
        self.assertEqual(found_num_total_repeats, 4)
        self.assertEqual(interruption_index, 2)
        self.assertFalse(has_partial_repeats)

        # do a general test of various motif sizes and repeat counts
        for repeat_unit in "A", "AC", "ACG", "ACGT", "ACGTA", "ACGGTA", "ACGGTAC", "ACGGTACG":
            for num_repeats in range(1, 50):
                for allow_interruptions in True, False:
                    sequence = repeat_unit * num_repeats
                    if allow_interruptions:
                        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, i, has_partial_repeats = find_repeat_unit_allowing_interruptions(
                            sequence)
                    else:
                        found_repeat_unit, found_num_pure_repeats, has_partial_repeats = find_repeat_unit_without_allowing_interruptions(
                            sequence)
                        found_num_total_repeats = found_num_pure_repeats

                    self.assertEqual(found_num_pure_repeats, found_num_total_repeats, f"found_num_pure_repeats != found_num_total_repeats for {num_repeats}x{repeat_unit}: {sequence}")
                    self.assertFalse(has_partial_repeats)

                    self.assertEqual(found_repeat_unit, repeat_unit,
                                     f"repeat unit not detected correctly for {num_repeats}x{repeat_unit}")
                    self.assertEqual(found_num_total_repeats, num_repeats,
                                     f"repeat count not detected correctly for {num_repeats}x{repeat_unit}")

                    # add an interruption to the last repeat of the sequence
                    for interruption_index in range(len(repeat_unit)):
                        sequence_with_interruption = add_interruption(sequence, repeat_unit, interruption_index)

                        if allow_interruptions:
                            found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, i, has_partial_repeats = find_repeat_unit_allowing_interruptions(
                                sequence_with_interruption)
                        else:
                            found_repeat_unit, found_num_pure_repeats, has_partial_repeats = find_repeat_unit_without_allowing_interruptions(
                                sequence_with_interruption)
                            i = None
                            found_num_total_repeats = found_num_pure_repeats

                        if allow_interruptions and num_repeats > 1:
                            if len(repeat_unit) >= 3:
                                self.assertEqual(found_repeat_unit, repeat_unit, f"{num_repeats}x{repeat_unit} interrupted at {interruption_index}: {sequence_with_interruption}. found_repeat_unit: {found_repeat_unit}. Expected: {repeat_unit} ")
                                self.assertEqual(found_num_pure_repeats, found_num_total_repeats - 1, f"{num_repeats}x{repeat_unit} interrupted at {interruption_index}: {sequence_with_interruption}. found_num_pure_repeats: {found_num_pure_repeats}. Expected: {found_num_total_repeats - 1}")
                                self.assertEqual(found_num_total_repeats, num_repeats, f"{num_repeats}x{repeat_unit} interrupted at {interruption_index}: {sequence_with_interruption}. found_num_total_repeats: {found_num_total_repeats}. Expected: {num_repeats}")
                                self.assertEqual(i, interruption_index, f"{num_repeats}x{repeat_unit} interrupted at {interruption_index}: {sequence_with_interruption}. found_num_pure_repeats: {found_num_pure_repeats}. Expected: {interruption_index} ")
                        else:
                            self.assertEqual(found_repeat_unit, sequence_with_interruption, f"{num_repeats}x{repeat_unit} interrupted at {interruption_index}: {sequence_with_interruption}. found_repeat_unit: {found_repeat_unit}. Expected: {sequence_with_interruption} ")
                            self.assertEqual(found_num_pure_repeats, 1, f"{num_repeats}x{repeat_unit} interrupted at {interruption_index}: {sequence_with_interruption}. found_num_pure_repeats: {found_num_pure_repeats}. Expected: 1 ")
                            self.assertEqual(found_num_total_repeats, 1, f"{num_repeats}x{repeat_unit} interrupted at {interruption_index}: {sequence_with_interruption}. found_num_total_repeats: {found_num_total_repeats}. Expected: 1 ")
                            self.assertEqual(i, None, f"{num_repeats}x{repeat_unit} interrupted at {interruption_index}: {sequence_with_interruption}. interruption_index: {i}. Expected: None ")

    def test_extend_repeat_into_sequence_without_allowing_interruptions(self):
        num_pure_repeats = extend_repeat_into_sequence_without_allowing_interruptions(
            "CAG", "CAGCAGCATCAGTTTT")
        self.assertEqual(num_pure_repeats, 2)

        num_pure_repeats = extend_repeat_into_sequence_without_allowing_interruptions(
            "CAG", "CATCATCATCAGTTTT")
        self.assertEqual(num_pure_repeats, 0)
        num_pure_repeats = extend_repeat_into_sequence_without_allowing_interruptions(
            "A", "CATCATCATCAGTTTT")
        self.assertEqual(num_pure_repeats, 0)
        num_pure_repeats = extend_repeat_into_sequence_without_allowing_interruptions(
            "A", "ATCATCATCAGTTTT")
        self.assertEqual(num_pure_repeats, 1)

    def test_extend_repeat_into_sequence_allowing_interruptions(self):
        num_pure_repeats, num_total_repeats, repeat_unit_interruption_index = extend_repeat_into_sequence_allowing_interruptions(
            "CA", "CACAGACATCAGTTTT", repeat_unit_interruption_index=None)
        self.assertEqual(num_pure_repeats, 2)
        self.assertEqual(num_total_repeats, 2)
        self.assertEqual(repeat_unit_interruption_index, None)

        num_pure_repeats, num_total_repeats, repeat_unit_interruption_index = extend_repeat_into_sequence_allowing_interruptions(
            "CAG", "CAGCAGCATCAGTTTT", repeat_unit_interruption_index=None)
        self.assertEqual(num_pure_repeats, 3)
        self.assertEqual(num_total_repeats, 4)
        self.assertEqual(repeat_unit_interruption_index, 2)

        # extended repeat sequence must end with an exact copy of the motif
        num_pure_repeats, num_total_repeats, repeat_unit_interruption_index = extend_repeat_into_sequence_allowing_interruptions(
            "CAG", "CAGCAGCATCAT", repeat_unit_interruption_index=None)
        self.assertEqual(num_pure_repeats, 2)
        self.assertEqual(num_total_repeats, 2)
        self.assertEqual(repeat_unit_interruption_index, None)

        num_pure_repeats, num_total_repeats, repeat_unit_interruption_index = extend_repeat_into_sequence_allowing_interruptions(
            "CAG", "CAGCAGCATCAT", repeat_unit_interruption_index=2)
        self.assertEqual(num_pure_repeats, 2)
        self.assertEqual(num_total_repeats, 2)
        self.assertEqual(repeat_unit_interruption_index, 2)

        num_pure_repeats, num_total_repeats, repeat_unit_interruption_index = extend_repeat_into_sequence_allowing_interruptions(
            "CAG", "CAGCAGCATCAGTTTCAG", repeat_unit_interruption_index=None)
        self.assertEqual(num_pure_repeats, 3)
        self.assertEqual(num_total_repeats, 4)
        self.assertEqual(repeat_unit_interruption_index, 2)

        num_pure_repeats, num_total_repeats, repeat_unit_interruption_index = extend_repeat_into_sequence_allowing_interruptions(
            "CAG", "CTGCAGCATCAATTTCAG", repeat_unit_interruption_index=1)
        self.assertEqual(num_pure_repeats, 1)
        self.assertEqual(num_total_repeats, 2)
        self.assertEqual(repeat_unit_interruption_index, 1)

        num_pure_repeats, num_total_repeats, repeat_unit_interruption_index = extend_repeat_into_sequence_allowing_interruptions(
            "CACA", "CGCACGCACACATTTTT", repeat_unit_interruption_index=None)
        self.assertEqual(num_pure_repeats, 1)
        self.assertEqual(num_total_repeats, 3)
        self.assertEqual(repeat_unit_interruption_index, 1)

    def test_count_pure_repeats(self):
        # 5 test cases for count_pure_repeats
        self.assertEqual(count_pure_repeats("CAGCAGCATCAG", "CAG"), 3)
        self.assertEqual(count_pure_repeats("CAGCAGCATCAG", "CAGCAG"), 1)
        self.assertEqual(count_pure_repeats("CAGCAGCATCAG", "CAGCAGCATCAG"), 1)
        self.assertEqual(count_pure_repeats("CAGCAGCATCAG", "CAT"), 1)
        self.assertEqual(count_pure_repeats("CAGCAGCATCAG", "AGC"), 0)

        # 5 test cases for count_pure_repeats with a 4bp repeat unit
        self.assertEqual(count_pure_repeats("CAGACATACAGA", "CAGA"), 2)
        self.assertEqual(count_pure_repeats("CAGACATACAGA", "GACA"), 0)

