import unittest

import collections
import random
random.seed(0)

from str_analysis.utils.misc_utils import reverse_complement
from str_analysis.parse_motif_composition import (
    LocusParser,
    _parse_motif_ids_from_processed_sequence,
    _is_nucleotide_sequence,
    parse_motif_composition_from_alignment_file,
)


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


class TestIsNucleotideSequence(unittest.TestCase):
    """Tests for _is_nucleotide_sequence helper function."""

    def test_sequences_with_nucleotides(self):
        """Test that sequences containing nucleotides return True."""
        self.assertTrue(_is_nucleotide_sequence("ACGT"))
        self.assertTrue(_is_nucleotide_sequence("A"))
        self.assertTrue(_is_nucleotide_sequence("C"))
        self.assertTrue(_is_nucleotide_sequence("G"))
        self.assertTrue(_is_nucleotide_sequence("T"))
        self.assertTrue(_is_nucleotide_sequence("N"))
        self.assertTrue(_is_nucleotide_sequence("ACGTN"))
        self.assertTrue(_is_nucleotide_sequence("123A456"))

    def test_sequences_without_nucleotides(self):
        """Test that sequences without nucleotides return False."""
        self.assertFalse(_is_nucleotide_sequence("123"))
        self.assertFalse(_is_nucleotide_sequence(""))
        self.assertFalse(_is_nucleotide_sequence("XYZ"))
        self.assertFalse(_is_nucleotide_sequence("||"))


class TestParseMotifIdsFromProcessedSequence(unittest.TestCase):
    """Tests for _parse_motif_ids_from_processed_sequence helper function."""

    def test_parse_single_motif_id(self):
        """Test parsing a single motif ID."""
        result = _parse_motif_ids_from_processed_sequence("|0|")
        self.assertEqual(result, [0])

    def test_parse_multiple_motif_ids(self):
        """Test parsing multiple motif IDs."""
        result = _parse_motif_ids_from_processed_sequence("|0||1||0|")
        self.assertEqual(result, [0, 1, 0])

    def test_parse_with_nucleotides(self):
        """Test parsing sequence with both motif IDs and nucleotides."""
        result = _parse_motif_ids_from_processed_sequence("AC|0||1|GT")
        self.assertEqual(result, [0, 1])

    def test_parse_empty_sequence(self):
        """Test parsing empty sequence."""
        result = _parse_motif_ids_from_processed_sequence("")
        self.assertEqual(result, [])

    def test_parse_nucleotides_only(self):
        """Test parsing sequence with only nucleotides."""
        result = _parse_motif_ids_from_processed_sequence("ACGTACGT")
        self.assertEqual(result, [])


class TestLocusParserInit(unittest.TestCase):
    """Tests for LocusParser.__init__."""

    def test_init_with_valid_dict(self):
        """Test initialization with valid motif frequency dictionary."""
        motif_dict = {"CAG": 10, "CCG": 5, "GGC": 3}
        parser = LocusParser(motif_dict)

        self.assertEqual(parser._reference_motif_frequency_dict, motif_dict)
        self.assertIsNotNone(parser._sorted_reference_motif_list)
        self.assertIsNotNone(parser._reference_motif_id_to_motif)

    def test_init_empty_dict_raises_error(self):
        """Test that initialization with empty dict raises ValueError."""
        with self.assertRaises(ValueError) as context:
            LocusParser({})
        self.assertIn("reference_motif_frequency_dict not specified", str(context.exception))

    def test_init_none_raises_error(self):
        """Test that initialization with None raises ValueError."""
        with self.assertRaises(ValueError) as context:
            LocusParser(None)
        self.assertIn("reference_motif_frequency_dict not specified", str(context.exception))

    def test_motif_sorting_by_frequency(self):
        """Test that motifs are sorted by frequency (descending)."""
        motif_dict = {"AAA": 1, "BBB": 10, "CCC": 5}
        parser = LocusParser(motif_dict)

        # BBB should be motif_id 0 (highest frequency)
        # CCC should be motif_id 1 (middle frequency)
        # AAA should be motif_id 2 (lowest frequency)
        self.assertEqual(parser._reference_motif_id_to_motif[0], "BBB")
        self.assertEqual(parser._reference_motif_id_to_motif[1], "CCC")
        self.assertEqual(parser._reference_motif_id_to_motif[2], "AAA")

    def test_internal_dicts_initialized(self):
        """Test that internal dictionaries are properly initialized."""
        parser = LocusParser({"CAG": 10})

        self.assertIsInstance(parser._observed_motif_id_frequency_dict, collections.defaultdict)
        self.assertIsInstance(parser._observed_motif_id_pair_frequency_dict, collections.defaultdict)
        self.assertIsInstance(parser._observed_motif_id_triplet_frequency_dict, collections.defaultdict)
        self.assertIsInstance(parser._novel_motif_frequency_dict, collections.defaultdict)


class TestConvertNucleotideSeqEdgeCases(unittest.TestCase):
    """Tests for edge cases in convert_nucleotide_seq_to_motif_seq."""

    def setUp(self):
        """Set up test fixtures."""
        self.motif_dict = {"CAG": 10, "CCG": 5}
        self.parser = LocusParser(self.motif_dict)

    def test_no_motifs_returns_none(self):
        """Test that sequence with no motifs returns None."""
        sequence = "ATGATGATG"
        result = self.parser.convert_nucleotide_seq_to_motif_seq(
            sequence,
            check_reverse_complement=False,
            chance_occurrence_threshold=None
        )
        self.assertIsNone(result)

    def test_lowercase_sequence(self):
        """Test that lowercase sequences are handled."""
        sequence = "cagcagcag"
        result = self.parser.convert_nucleotide_seq_to_motif_seq(
            sequence,
            check_reverse_complement=False,
            chance_occurrence_threshold=None
        )
        self.assertIsNotNone(result)
        self.assertIn("|0|", result)

    def test_invalid_sequence_raises_error(self):
        """Test that invalid nucleotide sequence raises ValueError."""
        sequence = "CAGXYZ"
        with self.assertRaises(ValueError) as context:
            self.parser.convert_nucleotide_seq_to_motif_seq(
                sequence,
                check_reverse_complement=False
            )
        self.assertIn("not a valid nucleotide sequence", str(context.exception))

    def test_sequence_with_n_bases(self):
        """Test that sequences with N bases are handled."""
        sequence = "CAGNAGCAG"
        result = self.parser.convert_nucleotide_seq_to_motif_seq(
            sequence,
            check_reverse_complement=False,
            chance_occurrence_threshold=None
        )
        # Should still process the sequence
        self.assertIn("|0|", result)

    def test_min_matches_threshold(self):
        """Test min_matches_threshold filtering."""
        sequence = "CAGCAG"  # Only 2 matches

        # With threshold of 3, should return None
        result = self.parser.convert_nucleotide_seq_to_motif_seq(
            sequence,
            check_reverse_complement=False,
            chance_occurrence_threshold=None,
            min_matches_threshold=3
        )
        self.assertIsNone(result)

        # With threshold of 2, should return result
        result = self.parser.convert_nucleotide_seq_to_motif_seq(
            sequence,
            check_reverse_complement=False,
            chance_occurrence_threshold=None,
            min_matches_threshold=2
        )
        self.assertIsNotNone(result)

    def test_mixed_motifs_with_flanking(self):
        """Test sequence with multiple different motifs and flanking sequence."""
        sequence = "ATCAGCCGCAGGT"
        result = self.parser.convert_nucleotide_seq_to_motif_seq(
            sequence,
            check_reverse_complement=False,
            chance_occurrence_threshold=None
        )

        self.assertIn("AT", result)
        self.assertIn("|0|", result)  # CAG
        self.assertIn("|1|", result)  # CCG
        self.assertIn("GT", result)


class TestGetterMethods(unittest.TestCase):
    """Tests for LocusParser getter methods."""

    def setUp(self):
        """Set up test fixtures."""
        self.motif_dict = {"CAG": 10, "CCG": 5}
        self.parser = LocusParser(self.motif_dict)

        # Process a sequence to populate internal dicts
        self.parser.convert_nucleotide_seq_to_motif_seq(
            "CAGCAGCCGAAACAGCCG",
            check_reverse_complement=False,
            chance_occurrence_threshold=None,
            record_reference_motif_counts=True,
            record_novel_motif_counts=True,
            record_motif_pair_counts=True,
            record_motif_triplet_counts=True
        )

    def test_get_reference_motif_id_to_motif_dict(self):
        """Test get_reference_motif_id_to_motif_dict returns correct mapping."""
        result = self.parser.get_reference_motif_id_to_motif_dict()

        self.assertIsInstance(result, dict)
        self.assertIn(0, result)
        self.assertIn(1, result)
        self.assertEqual(result[0], "CAG")
        self.assertEqual(result[1], "CCG")

    def test_get_reference_motif_to_motif_id_dict(self):
        """Test get_reference_motif_to_motif_id_dict returns correct mapping."""
        result = self.parser.get_reference_motif_to_motif_id_dict()

        self.assertIsInstance(result, dict)
        self.assertIn("CAG", result)
        self.assertIn("CCG", result)
        self.assertEqual(result["CAG"], 0)
        self.assertEqual(result["CCG"], 1)

    def test_get_observed_motif_frequency_dict(self):
        """Test get_observed_motif_frequency_dict returns motif counts."""
        result = self.parser.get_observed_motif_frequency_dict()

        self.assertIsInstance(result, dict)
        self.assertIn("CAG", result)
        self.assertIn("CCG", result)
        self.assertEqual(result["CAG"], 3)
        self.assertEqual(result["CCG"], 2)

    def test_get_observed_motif_frequency_dict_as_string(self):
        """Test get_observed_motif_frequency_dict with as_string=True."""
        result = self.parser.get_observed_motif_frequency_dict(as_string=True)

        self.assertIsInstance(result, str)
        self.assertIn("[CAG]:", result)
        self.assertIn("[CCG]:", result)

    def test_get_observed_motif_pair_frequency_dict(self):
        """Test get_observed_motif_pair_frequency_dict returns pair counts."""
        result = self.parser.get_observed_motif_pair_frequency_dict()

        self.assertIsInstance(result, dict)
        # Should have some pairs
        self.assertGreater(len(result), 0)

        # All keys should be tuples of motifs
        for key in result.keys():
            self.assertIsInstance(key, tuple)
            self.assertEqual(len(key), 2)

    def test_get_observed_motif_pair_frequency_dict_as_string(self):
        """Test get_observed_motif_pair_frequency_dict with as_string=True."""
        result = self.parser.get_observed_motif_pair_frequency_dict(as_string=True)

        self.assertIsInstance(result, str)
        # Format should be [motif1][motif2]:count
        if result:  # Only check if there are pairs
            self.assertIn("]:", result)

    def test_get_observed_motif_triplet_frequency_dict(self):
        """Test get_observed_motif_triplet_frequency_dict returns triplet counts."""
        result = self.parser.get_observed_motif_triplet_frequency_dict()

        self.assertIsInstance(result, dict)

        # All keys should be tuples of motifs
        for key in result.keys():
            self.assertIsInstance(key, tuple)
            self.assertEqual(len(key), 3)

    def test_get_observed_motif_triplet_frequency_dict_as_string(self):
        """Test get_observed_motif_triplet_frequency_dict with as_string=True."""
        result = self.parser.get_observed_motif_triplet_frequency_dict(as_string=True)

        self.assertIsInstance(result, str)
        # Format should be [motif1][motif2][motif3]:count
        if result:  # Only check if there are triplets
            self.assertIn("]:", result)

    def test_get_novel_motif_frequency_dict(self):
        """Test get_novel_motif_frequency_dict returns novel motifs."""
        result = self.parser.get_novel_motif_frequency_dict()

        self.assertIsInstance(result, dict)
        self.assertIn("AAA", result)
        self.assertEqual(result["AAA"], 1)

    def test_get_novel_motif_frequency_dict_as_string(self):
        """Test get_novel_motif_frequency_dict with as_string=True."""
        result = self.parser.get_novel_motif_frequency_dict(as_string=True)

        self.assertIsInstance(result, str)
        self.assertIn("[AAA]:", result)

    def test_novel_motif_not_in_reference_raises_error(self):
        """Test that novel motifs in reference dict raises Exception."""
        # Manually add a reference motif to novel dict (shouldn't happen normally)
        self.parser._novel_motif_frequency_dict["CAG"] = 1

        with self.assertRaises(Exception) as context:
            self.parser.get_novel_motif_frequency_dict()
        self.assertIn("is not novel", str(context.exception))


class TestParseMotifCompositionFromAlignmentFile(unittest.TestCase):
    """Tests for parse_motif_composition_from_alignment_file function."""

    def test_nucleotide_sequence_input(self):
        """Test processing a nucleotide sequence (not a file)."""
        import tempfile
        import os

        motif_dict = {"CAG": 10, "CCG": 5}
        sequence = "CAGCAGCCG"

        with tempfile.TemporaryDirectory() as tmpdir:
            output_prefix = os.path.join(tmpdir, "test_output")

            result = parse_motif_composition_from_alignment_file(
                input_sequence_or_path=sequence,
                motif_frequency_dict=motif_dict,
                output_prefix=output_prefix,
                output_format="JSON",
                check_reverse_complement=False,
            )

            self.assertIsInstance(result, dict)
            self.assertEqual(result["input"], sequence)
            self.assertIsNone(result["input_file_size"])
            self.assertIsNone(result["total_mapped_reads"])
            self.assertIn("motif_frequency", result)
            self.assertIn("[CAG]:", result["motif_frequency"])

    def test_invalid_nucleotide_sequence_raises_error(self):
        """Test that invalid nucleotide sequence raises ValueError."""
        import tempfile
        import os

        motif_dict = {"CAG": 10}
        sequence = "CAGXYZ"

        with tempfile.TemporaryDirectory() as tmpdir:
            output_prefix = os.path.join(tmpdir, "test_output")

            with self.assertRaises(ValueError) as context:
                parse_motif_composition_from_alignment_file(
                    input_sequence_or_path=sequence,
                    motif_frequency_dict=motif_dict,
                    output_prefix=output_prefix,
                    output_format="JSON",
                )
            self.assertIn("Invalid characters detected", str(context.exception))

    def test_no_output_prefix_raises_error(self):
        """Test that missing output_prefix raises ValueError."""
        motif_dict = {"CAG": 10}
        sequence = "CAGCAGCAG"

        with self.assertRaises(ValueError) as context:
            parse_motif_composition_from_alignment_file(
                input_sequence_or_path=sequence,
                motif_frequency_dict=motif_dict,
                output_prefix=None,
                output_format="JSON",
            )
        self.assertIn("output_prefix not specified", str(context.exception))

    def test_json_output_format(self):
        """Test that JSON output format creates .json.gz file."""
        import tempfile
        import os

        motif_dict = {"CAG": 10}
        sequence = "CAGCAGCAG"

        with tempfile.TemporaryDirectory() as tmpdir:
            output_prefix = os.path.join(tmpdir, "test_output")

            result = parse_motif_composition_from_alignment_file(
                input_sequence_or_path=sequence,
                motif_frequency_dict=motif_dict,
                output_prefix=output_prefix,
                output_format="JSON",
                check_reverse_complement=False,
            )

            output_file = output_prefix + ".json.gz"
            self.assertTrue(os.path.exists(output_file))

    def test_tsv_output_format(self):
        """Test that TSV output format creates .tsv.gz file."""
        import tempfile
        import os

        motif_dict = {"CAG": 10}
        sequence = "CAGCAGCAG"

        with tempfile.TemporaryDirectory() as tmpdir:
            output_prefix = os.path.join(tmpdir, "test_output")

            result = parse_motif_composition_from_alignment_file(
                input_sequence_or_path=sequence,
                motif_frequency_dict=motif_dict,
                output_prefix=output_prefix,
                output_format="TSV",
                check_reverse_complement=False,
            )

            output_file = output_prefix + ".tsv.gz"
            self.assertTrue(os.path.exists(output_file))

    def test_unsupported_output_format_raises_error(self):
        """Test that unsupported output format raises ValueError."""
        import tempfile
        import os

        motif_dict = {"CAG": 10}
        sequence = "CAGCAGCAG"

        with tempfile.TemporaryDirectory() as tmpdir:
            output_prefix = os.path.join(tmpdir, "test_output")

            with self.assertRaises(ValueError) as context:
                parse_motif_composition_from_alignment_file(
                    input_sequence_or_path=sequence,
                    motif_frequency_dict=motif_dict,
                    output_prefix=output_prefix,
                    output_format="XML",
                    check_reverse_complement=False,
                )
            self.assertIn("is not supported", str(context.exception))

    def test_output_contains_expected_fields(self):
        """Test that output dictionary contains all expected fields."""
        import tempfile
        import os

        motif_dict = {"CAG": 10}
        sequence = "CAGCAGCAG"

        with tempfile.TemporaryDirectory() as tmpdir:
            output_prefix = os.path.join(tmpdir, "test_output")

            result = parse_motif_composition_from_alignment_file(
                input_sequence_or_path=sequence,
                motif_frequency_dict=motif_dict,
                output_prefix=output_prefix,
                output_format="JSON",
                check_reverse_complement=False,
            )

            expected_fields = [
                "input",
                "input_file_size",
                "total_mapped_reads",
                "motif_frequency",
                "motif_pair_frequency",
                "motif_triplet_frequency",
                "novel_motif_frequency",
                "min_mapq",
                "include_low_quality_alignments",
            ]

            for field in expected_fields:
                self.assertIn(field, result)

    def test_verbose_output(self):
        """Test verbose output mode."""
        import tempfile
        import os
        import sys
        from io import StringIO

        motif_dict = {"CAG": 10}
        sequence = "CAGCAGCAG"

        with tempfile.TemporaryDirectory() as tmpdir:
            output_prefix = os.path.join(tmpdir, "test_output")

            # Capture stdout
            old_stdout = sys.stdout
            sys.stdout = captured_output = StringIO()

            try:
                result = parse_motif_composition_from_alignment_file(
                    input_sequence_or_path=sequence,
                    motif_frequency_dict=motif_dict,
                    output_prefix=output_prefix,
                    output_format="JSON",
                    check_reverse_complement=False,
                    verbose=True,
                )
            finally:
                sys.stdout = old_stdout

            output = captured_output.getvalue()

            # Check that verbose output was produced
            self.assertIn("Parsed motif sequence", output)
            self.assertIn("Observed motif frequencies", output)
