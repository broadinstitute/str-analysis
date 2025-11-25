import unittest

import collections
import random

from str_analysis.utils.misc_utils import reverse_complement
from parse_motif_composition_from_short_reads import LocusParser, _parse_motif_ids_from_processed_sequence


class LocusParserTest(unittest.TestCase):

    def test_parses_known_motifs_in_forward_orientation(self):
        reference_motif_frequency_dict = {"AC": 5, "GT": 3}
        locus_parser = LocusParser(reference_motif_frequency_dict)

        parsed_sequence = locus_parser.convert_nucleotide_seq_to_motif_seq(
            "ACGTAC",
            record_reference_motif_counts=True,
        )

        self.assertEqual(parsed_sequence, "|0||1||0|")
        self.assertEqual(locus_parser._observed_motif_frequency_dict["AC"], 2)
        self.assertEqual(locus_parser._observed_motif_frequency_dict["GT"], 1)


    def test_prefers_reverse_complement_when_it_matches(self):
        reference_motif_frequency_dict = {"GTA": 3, "TCG": 2}
        locus_parser = LocusParser(reference_motif_frequency_dict)

        parsed_sequence = locus_parser.convert_nucleotide_seq_to_motif_seq(
            reverse_complement("ACGTAGTAGTAAC"),
            record_reference_motif_counts=True,
        )

        self.assertEqual(parsed_sequence, "AC|0||0||0|AC")
        self.assertEqual(locus_parser._observed_motif_frequency_dict["GTA"], 3)


    def test_records_novel_motifs_matching_known_lengths(self):
        reference_motif_frequency_dict = {
            "AAAAG": 10,
            "AAGGG": 5,
            "AAAGG": 3,
            "ACAGG": 3,
        }

        novel_motif = "ACCGG"

        seq = ""
        expected_parsed_seq = ""
        for i, (key, count) in enumerate(reference_motif_frequency_dict.items()):
            seq += key * count
            expected_parsed_seq += f"|{i}|"*count
        seq += novel_motif
        expected_parsed_seq += novel_motif
        seq *= 2
        expected_parsed_seq *= 2
        #print("seq length:", len(seq), collections.Counter(seq))

        locus_parser = LocusParser(reference_motif_frequency_dict)
        parsed_sequence = locus_parser.convert_nucleotide_seq_to_motif_seq(
            #reverse_complement(seq),
            seq,
            record_reference_motif_counts=True,
            record_novel_motif_counts=True,
        )

        self.assertEqual(parsed_sequence, expected_parsed_seq)

    def test_probability_of_result_by_chance(self):
        reference_motif_frequency_dict = {
            "AAAAG": 10,
            "AAGGG": 5,
            "AAAGG": 3,
            "ACAGG": 3,
        }

        novel_motif = "ACCGG"

        sequence = ""
        for key, count in reference_motif_frequency_dict.items():
            sequence += key * count
        sequence += novel_motif
        sequence = sequence * 2

        random_sequence = list(sequence)
        random.shuffle(random_sequence)
        random_sequence = "".join(random_sequence)

        for is_random, seq in [(False, sequence), (True, random_sequence)]:
            locus_parser = LocusParser(reference_motif_frequency_dict)
            parsed_sequence = locus_parser.convert_nucleotide_seq_to_motif_seq(
                seq,
                chance_occurrence_threshold=None,
                check_reverse_complement=True,
                record_reference_motif_counts=False,
                record_novel_motif_counts=False,
                record_motif_pair_counts=False,
            )
            motif_ids = _parse_motif_ids_from_processed_sequence(parsed_sequence)
            probability = locus_parser._probability_of_result_by_chance(seq, motif_ids)
            parsed_sequence2 = locus_parser.convert_nucleotide_seq_to_motif_seq(
                seq,
                record_reference_motif_counts=True,
                record_novel_motif_counts=True,
                record_motif_pair_counts=True,
            )

            #print(probability)
            #print(parsed_sequence)
            #print(parsed_sequence2)
            #print("-----")

            if is_random:
                self.assertGreater(probability, 0.001)
                self.assertIsNone(parsed_sequence2)
            else:
                self.assertLess(probability, 0.001)
                self.assertIsNotNone(parsed_sequence2)
                self.assertDictEqual(locus_parser.get_observed_motif_frequency_dict(), {
                    motif: count * 2 for motif, count in reference_motif_frequency_dict.items()
                })
                self.assertDictEqual(locus_parser.get_novel_motif_frequency_dict(), {
                    novel_motif: 2
                })
                self.assertDictEqual(
                    {
                        ("AAAAG", "AAAAG"): (10 - 1) * 2,
                        ("AAGGG", "AAGGG"): (5 - 1) * 2,
                        ("AAAGG", "AAAGG"): (3 - 1) * 2,
                        ("ACAGG", "ACAGG"): (3 - 1) * 2,

                        ("AAAAG", "AAGGG"): 2,
                        ("AAGGG", "AAAGG"): 2,
                        ("AAAGG", "ACAGG"): 2,
                        #("ACAGG", novel_motif): 2,  # novel motifs are not counted in pairs
                    },
                    locus_parser.get_observed_motif_pair_frequency_dict(),
                )
