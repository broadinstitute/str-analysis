import unittest

from str_analysis.utils.find_repeat_unit import find_repeat_unit
from utils.find_repeat_unit import get_repeat_unit_regex_with_N_base, extend_repeat_into_sequence


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
        self.assertEqual(get_repeat_unit_regex_with_N_base("CAG", 0, allow_homopolymer=True), "[ACGT]AG")
        self.assertEqual(get_repeat_unit_regex_with_N_base("CAG", 1, allow_homopolymer=True), "C[ACGT]G")
        self.assertEqual(get_repeat_unit_regex_with_N_base("CAG", 2, allow_homopolymer=True), "CA[ACGT]")

    def test_find_repeat_unit(self):

        for repeat_unit in "A", "AC", "ACG", "ACGT", "ACGTA", "ACGATA":
            for num_repeats in 1, 2, 3, 4, 5:
                for allow_interruptions in True, False:
                    sequence = repeat_unit * num_repeats
                    found_repeat_unit, found_num_pure_repeats, found_num_repeats, i = find_repeat_unit(
                        sequence, allow_interruptions=allow_interruptions)
                    if not allow_interruptions:
                        self.assertEqual(found_num_pure_repeats, found_num_repeats)
                        self.assertIsNone(i)

                    self.assertEqual(found_repeat_unit, repeat_unit,
                                     f"repeat unit not detected correctly for {num_repeats}x{repeat_unit}")
                    self.assertEqual(found_num_repeats, num_repeats,
                                     f"repeat count not detected correctly for {num_repeats}x{repeat_unit}")


    def test_find_repeat_unit(self):
        # test a few examples
        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, interruption_index = find_repeat_unit(
            "CAG"*5, allow_interruptions=True)
        self.assertEqual(found_repeat_unit, "CAG")
        self.assertEqual(found_num_pure_repeats, 5)
        self.assertEqual(found_num_total_repeats, 5)
        self.assertEqual(interruption_index, None)

        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, interruption_index = find_repeat_unit(
            "CAG"*5 + "C", allow_interruptions=True, allow_partial_repeats=True)
        self.assertEqual(found_repeat_unit, "CAG")
        self.assertEqual(found_num_pure_repeats, 5)
        self.assertEqual(found_num_total_repeats, 5)
        self.assertEqual(interruption_index, None)

        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, interruption_index = find_repeat_unit(
            "CAG"*5 + "C", allow_interruptions=True, allow_partial_repeats=False)
        self.assertEqual(found_repeat_unit, "CAG"*5 + "C")
        self.assertEqual(found_num_pure_repeats, 1)
        self.assertEqual(found_num_total_repeats, 1)
        self.assertEqual(interruption_index, None)

        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, interruption_index = find_repeat_unit(
            "CATCAGCAGCAGCAG", allow_interruptions=True)
        self.assertEqual(found_repeat_unit, "CAT")
        self.assertEqual(found_num_pure_repeats, 1)
        self.assertEqual(found_num_total_repeats, 5)
        self.assertEqual(interruption_index, 2)

        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, interruption_index = find_repeat_unit(
            "CAGCATCAGCAGCAG", allow_interruptions=True)
        self.assertEqual(found_repeat_unit, "CAG")
        self.assertEqual(found_num_pure_repeats, 4)
        self.assertEqual(found_num_total_repeats, 5)
        self.assertEqual(interruption_index, 2)

        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, interruption_index = find_repeat_unit(
            "A"*5, allow_interruptions=True)
        self.assertEqual(found_repeat_unit, "A")
        self.assertEqual(found_num_pure_repeats, 5)
        self.assertEqual(found_num_total_repeats, 5)
        self.assertEqual(interruption_index, None)

        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, interruption_index = find_repeat_unit(
            "ACG"*2 + "ACA" + "ACG"*2 + "ACT"*2, allow_interruptions=True)
        self.assertEqual(found_repeat_unit, "ACG")
        self.assertEqual(found_num_pure_repeats, 4)
        self.assertEqual(found_num_total_repeats, 7)
        self.assertEqual(interruption_index, 2)

        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, interruption_index = find_repeat_unit(
            "AAAAAGAAAA", allow_interruptions=True)
        self.assertEqual(found_repeat_unit, "AAAAAGAAAA")
        self.assertEqual(found_num_pure_repeats, 1)
        self.assertEqual(found_num_total_repeats, 1)
        self.assertEqual(interruption_index, None)

        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, interruption_index = find_repeat_unit(
            "CTCTCTCACT", allow_interruptions=True)
        self.assertEqual(found_repeat_unit, "CTCTCTCACT")
        self.assertEqual(found_num_pure_repeats, 1)
        self.assertEqual(found_num_total_repeats, 1)
        self.assertEqual(interruption_index, None)

        # do a general test of various motif sizes and repeat counts
        for repeat_unit in "A", "AC", "ACG", "ACGT", "ACGTA", "ACGGTA", "ACGGTAC", "ACGGTACG":
            for num_repeats in 1, 2, 3, 4, 5:
                for allow_interruptions in True, False:
                    sequence = repeat_unit * num_repeats
                    found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, i = find_repeat_unit(
                        sequence, allow_interruptions=allow_interruptions)

                    self.assertEqual(found_num_pure_repeats, found_num_total_repeats)
                    self.assertIsNone(i)

                    self.assertEqual(found_repeat_unit, repeat_unit,
                                     f"repeat unit not detected correctly for {num_repeats}x{repeat_unit}")
                    self.assertEqual(found_num_total_repeats, num_repeats,
                                     f"repeat count not detected correctly for {num_repeats}x{repeat_unit}")

                    # add an interruption to the last repeat of the sequence
                    for interruption_index in range(len(repeat_unit)):
                        sequence_with_interruption = add_interruption(sequence, repeat_unit, interruption_index)

                        found_repeat_unit, found_num_pure_repeats, found_num_total_repeats, i = find_repeat_unit(
                            sequence_with_interruption, allow_interruptions=allow_interruptions)

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


    def test_extend_repeat_into_sequence(self):
        num_pure_repeats, num_total_repeats, repeat_unit_interruption_index = extend_repeat_into_sequence(
            "CAG", "CAGCAGCATCAG", allow_interruptions=False, repeat_unit_interruption_index=None)
        self.assertEqual(num_pure_repeats, 2)
        self.assertEqual(num_total_repeats, 2)
        self.assertEqual(repeat_unit_interruption_index, None)

        num_pure_repeats, num_total_repeats, repeat_unit_interruption_index = extend_repeat_into_sequence(
            "CAG", "CAGCAGCATCAG", allow_interruptions=True, repeat_unit_interruption_index=None)
        self.assertEqual(num_pure_repeats, 2)
        self.assertEqual(num_total_repeats, 4)
        self.assertEqual(repeat_unit_interruption_index, 2)

        # extended repeat sequence must end with an exact copy of the motif
        num_pure_repeats, num_total_repeats, repeat_unit_interruption_index = extend_repeat_into_sequence(
            "CAG", "CAGCAGCATCAT", allow_interruptions=True, repeat_unit_interruption_index=None)
        self.assertEqual(num_pure_repeats, 2)
        self.assertEqual(num_total_repeats, 2)
        self.assertEqual(repeat_unit_interruption_index, None)

        num_pure_repeats, num_total_repeats, repeat_unit_interruption_index = extend_repeat_into_sequence(
            "CAG", "CAGCAGCATCAT", allow_interruptions=True, repeat_unit_interruption_index=2)
        self.assertEqual(num_pure_repeats, 2)
        self.assertEqual(num_total_repeats, 2)
        self.assertEqual(repeat_unit_interruption_index, 2)

        num_pure_repeats, num_total_repeats, repeat_unit_interruption_index = extend_repeat_into_sequence(
            "CAG", "CAGCAGCATCAGTTTCAG", allow_interruptions=True, repeat_unit_interruption_index=None)
        self.assertEqual(num_pure_repeats, 2)
        self.assertEqual(num_total_repeats, 4)
        self.assertEqual(repeat_unit_interruption_index, 2)

        num_pure_repeats, num_total_repeats, repeat_unit_interruption_index = extend_repeat_into_sequence(
            "CAG", "CTGCAGCATCAATTTCAG", allow_interruptions=True, repeat_unit_interruption_index=1)
        self.assertEqual(num_pure_repeats, 0)
        self.assertEqual(num_total_repeats, 2)
        self.assertEqual(repeat_unit_interruption_index, 1)