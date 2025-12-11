import unittest

import collections
import random
random.seed(0)

from str_analysis.utils.misc_utils import reverse_complement
from str_analysis.parse_motif_composition import LocusParser, _parse_motif_ids_from_processed_sequence


class LocusParserTest(unittest.TestCase):

    def test_parses_known_motifs_in_forward_orientation(self):
        reference_motif_frequency_dict = {"AC": 5, "GT": 3}
        locus_parser = LocusParser(reference_motif_frequency_dict)

        parsed_sequence = locus_parser.convert_nucleotide_seq_to_motif_seq(
            "ACGTAC",
            record_reference_motif_counts=True,
        )

        self.assertEqual(parsed_sequence, "|0||1||0|")
        self.assertEqual(locus_parser._observed_motif_id_frequency_dict[0], 2)
        self.assertEqual(locus_parser._observed_motif_id_frequency_dict[1], 1)


    def test_prefers_reverse_complement_when_it_matches(self):
        reference_motif_frequency_dict = {"GTA": 3, "TCG": 2}
        locus_parser = LocusParser(reference_motif_frequency_dict)

        parsed_sequence = locus_parser.convert_nucleotide_seq_to_motif_seq(
            reverse_complement("ACGTAGTAGTAAC"),
            record_reference_motif_counts=True,
        )

        self.assertEqual(parsed_sequence, "AC|0||0||0|AC")
        self.assertEqual(locus_parser._observed_motif_id_frequency_dict[0], 3)


    def test_records_novel_motifs_matching_known_lengths(self):
        reference_motif_frequency_dict = {
            "AAAAG": 10,
            "AAGGG": 5,
            "AAAGG": 3,
            "ACAGG": 3,
        }
        novel_motif = "ACCGG"

        locus_parser = LocusParser(reference_motif_frequency_dict)

        seq = ""
        expected_parsed_seq = ""
        motif_to_motif_id = locus_parser.get_reference_motif_to_motif_id_dict()
        for motif, count in reference_motif_frequency_dict.items():
            seq += motif * count
            i = motif_to_motif_id[motif]
            expected_parsed_seq += f"|{i}|"*count
        seq += novel_motif
        expected_parsed_seq += novel_motif
        seq *= 2
        expected_parsed_seq *= 2
        #print("seq length:", len(seq), collections.Counter(seq))

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
        motif_sequence = []
        for motif, count in reference_motif_frequency_dict.items():
            sequence += motif * count
            motif_sequence += [motif] * count
        sequence += novel_motif
        motif_sequence.append(novel_motif)
        sequence = sequence * 2
        motif_sequence = motif_sequence * 2

        expected_pair_counts = collections.defaultdict(int)
        expected_triplet_counts = collections.defaultdict(int)
        for i in range(len(motif_sequence)):
            for count_dict, tuple_size in [(expected_pair_counts, 2), (expected_triplet_counts, 3)]:
                if i - tuple_size + 1 >= 0:
                    t = tuple([motif_sequence[x] for x in range(i - tuple_size + 1, i + 1)])
                    if all(m in reference_motif_frequency_dict for m in t):
                        count_dict[t] += 1  # since pairs and tripplets will only be counted once in the second call to locus_parser.convert_nucleotide_sequence



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
                record_motif_triplet_counts=False,
            )
            motif_ids = _parse_motif_ids_from_processed_sequence(parsed_sequence)
            probability = locus_parser._probability_of_result_by_chance(seq, motif_ids)
            parsed_sequence2 = locus_parser.convert_nucleotide_seq_to_motif_seq(
                seq,
                record_reference_motif_counts=True,
                record_novel_motif_counts=True,
                record_motif_pair_counts=True,
                record_motif_triplet_counts=True,
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
                self.assertDictEqual(expected_pair_counts, locus_parser.get_observed_motif_pair_frequency_dict())
                print("get_observed_motif_pair_frequency_dict:", len(locus_parser.get_observed_motif_pair_frequency_dict()))

                self.assertDictEqual(expected_triplet_counts, locus_parser.get_observed_motif_triplet_frequency_dict())
                print("get_observed_motif_triplet_frequency_dict:", len(locus_parser.get_observed_motif_triplet_frequency_dict()))
                print(locus_parser.get_observed_motif_frequency_dict())
                print("motifs:", locus_parser.get_observed_motif_frequency_dict(as_string=True))
                print("motif pairs:", locus_parser.get_observed_motif_pair_frequency_dict(as_string=True))
                print("triplets:", locus_parser.get_observed_motif_triplet_frequency_dict(as_string=True))
