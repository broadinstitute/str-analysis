#!/usr/bin/env python3

"""Test script for the filter_vcf_to_catalog_tandem_repeats.py script."""

import argparse
import collections
import os
import pkgutil
import shutil
import tempfile
import unittest
from unittest import mock

import pyfaidx

from str_analysis.filter_vcf_to_tandem_repeats import Allele, TandemRepeatAllele, ReferenceTandemRepeat, \
    GenotypedTandemRepeat, \
    DETECTION_MODE_PURE_REPEATS, DETECTION_MODE_ALLOW_INTERRUPTIONS, DETECTION_MODE_TRF, \
    FILTER_ALLELE_INDEL_WITHOUT_REPEATS, FILTER_TR_ALLELE_NOT_ENOUGH_REPEATS, \
    FILTER_TR_ALLELE_DOESNT_SPAN_ENOUGH_BASE_PAIRS, FILTER_TR_ALLELE_REPEAT_UNIT_TOO_SHORT, \
    FILTER_TR_ALLELE_REPEAT_UNIT_TOO_LONG, FILTER_ALLELE_WITH_N_BASES, \
    FILTER_TR_ALLELE_NOT_ENOUGH_REPEATS_IN_REFERENCE, FILTER_TR_ALLELE_TOO_MANY_REPEATS, \
    FILTER_TR_ALLELE_SPANS_TOO_MANY_BASE_PAIRS, FILTER_TR_ALLELE_PURITY_IS_TOO_LOW, \
    TRF_MAX_REPEATS_IN_REFERENCE_THRESHOLD, TRF_MAX_SPAN_IN_REFERENCE_THRESHOLD, \
    merge_overlapping_tandem_repeat_loci, detect_perfect_and_almost_perfect_tandem_repeats, \
    detect_tandem_repeats_using_trf, check_if_allele_is_tandem_repeat, \
    check_if_tandem_repeat_allele_failed_filters, compute_repeat_unit_id, are_repeat_units_similar, \
    need_to_reprocess_allele_with_extended_flanking_sequence, \
    detect_chromosome_naming_convention, normalize_chromosome_name, \
    open_vcf_for_genotyping, get_overlapping_vcf_variants, \
    convert_variants_to_haplotype_sequence, extract_haplotype_sequences_from_vcf, \
    compute_repeat_counts_from_sequence, genotype_single_locus


class TestAllele(unittest.TestCase):
    """Test the Allele class."""

    def setUp(self):
        """Set up the test case."""
        # create file-like object around a binary string
        fasta_data = pkgutil.get_data("str_analysis", "data/tests/chr22_11Mb.fa.gz")
        # write it to a temporary named file
        with tempfile.NamedTemporaryFile(suffix=".fa.gz", delete=False) as fasta_file:
            self._temp_fasta_path = fasta_file.name
            fasta_file.write(fasta_data)
            fasta_file.flush()

        self._fasta_obj = pyfaidx.Fasta(self._temp_fasta_path, one_based_attributes=False, as_raw=True)
        self._fasta_obj_with_one_based_attributes = pyfaidx.Fasta(self._temp_fasta_path, one_based_attributes=True, as_raw=True)

        self._poly_a_insertion1 = Allele("chr22", 10513201, "T", "TAAAAA", self._fasta_obj)
        self._poly_a_insertion2 = Allele("chr22", 10513201, "T", "TAAAA", self._fasta_obj)
        self._poly_a_deletion1 = Allele("chr22", 10513201, "TAAAAA", "T", self._fasta_obj)
        self._poly_a_deletion2 = Allele("chr22", 10513201, "TAAAA", "T", self._fasta_obj)

        self._AAGA_insertion = Allele("chr22", 10515040, "T", "TAAGA", self._fasta_obj)
        self._AAGA_deletion = Allele("chr22", 10515040, "TAAGA", "T", self._fasta_obj)

        self._AAGA_tandem_repeat_expansion = TandemRepeatAllele(self._AAGA_insertion, "AAGA", False, 0, 4, 37, DETECTION_MODE_PURE_REPEATS)
        self._AAGA_tandem_repeat_contraction = TandemRepeatAllele(self._AAGA_deletion, "AAGA", False, 0, 4, 33, DETECTION_MODE_PURE_REPEATS)

        self._trf_working_dir = tempfile.mkdtemp()
        self._default_args = argparse.Namespace(
            min_repeat_unit_length=1,
            max_repeat_unit_length=1000,
            show_progress_bar=False,
            min_repeats=3,
            min_tandem_repeat_length=9,
            debug=False,
            trf_working_dir=self._trf_working_dir,
            input_vcf_prefix="test",
            trf_executable_path="trf",
            trf_threads=2,
            verbose=False,
            allow_multiple_trf_results_per_locus=False,
            dont_allow_interruptions=False,
            dont_run_trf=False,
            min_indel_size_to_run_trf=7,
            trf_min_repeats_in_reference=2,
            trf_min_purity=0.2,
        )

    def test_getters(self):
        """Test the Allele class initialization."""
        for insertion_allele in [self._poly_a_insertion1, self._poly_a_insertion2]:
            self.assertEqual(insertion_allele.chrom, "chr22")
            self.assertEqual(insertion_allele.pos, 10513201)
            self.assertEqual(insertion_allele.ref, "T")
            self.assertTrue(insertion_allele.alt.startswith("TAAAA"))
            self.assertEqual(insertion_allele.ins_or_del, "INS")

            left_flanking_sequence = insertion_allele.get_left_flanking_sequence()
            self.assertEqual(left_flanking_sequence[-6:], "TAGATT")

            right_flanking_sequence = insertion_allele.get_right_flanking_sequence()
            self.assertEqual(right_flanking_sequence[:15], "A"*9 + "CTGTGG")
            self.assertEqual(insertion_allele.get_left_flank_end(), 10513201)
            self.assertEqual(insertion_allele.get_right_flank_start_0based(), 10513201)
            self.assertEqual(len(insertion_allele.get_left_flanking_sequence()), 300)
            self.assertEqual(len(insertion_allele.get_right_flanking_sequence()), 300)

        for deletion_allele in [self._poly_a_deletion1, self._poly_a_deletion2]:
            self.assertEqual(deletion_allele.chrom, "chr22")
            self.assertEqual(deletion_allele.pos, 10513201)
            self.assertTrue(deletion_allele.ref.startswith("TAAAA"))
            self.assertEqual(deletion_allele.alt, "T")
            self.assertEqual(deletion_allele.ins_or_del, "DEL")
            self.assertTrue(deletion_allele.variant_bases.startswith("AAAA"))

            left_flanking_sequence = deletion_allele.get_left_flanking_sequence()
            self.assertEqual(left_flanking_sequence[-6:], "TAGATT")

            right_flanking_sequence = deletion_allele.get_right_flanking_sequence()
            self.assertTrue(right_flanking_sequence.startswith("A"*(9 - len(deletion_allele.variant_bases)) + "CTGTGG"))

            self.assertEqual(deletion_allele.get_left_flank_end(), 10513201)
            self.assertEqual(deletion_allele.get_right_flank_start_0based(), deletion_allele.get_left_flank_end() + len(deletion_allele.variant_bases))
            self.assertEqual(len(deletion_allele.get_left_flanking_sequence()), 300)
            self.assertEqual(len(deletion_allele.get_right_flanking_sequence()), 300)


        self.assertRaises(ValueError, Allele, "chr1", 100, "A", "T", self._fasta_obj)  # SNV
        self.assertRaises(ValueError, Allele, "chr1", 100, "AT", "TAG", self._fasta_obj)  # Complex MNP
        self.assertRaises(ValueError, Allele, "chr1", 100, "A", "AT", self._fasta_obj_with_one_based_attributes)  # incorrectly initialized fasta

    def test_getters2(self):
        """Test the Allele class getters."""
        self.assertEqual(self._AAGA_insertion.get_left_flank_end(), 10515040)
        self.assertTrue(self._AAGA_insertion.get_left_flanking_sequence().endswith("AAAT"))

        self.assertEqual(self._AAGA_insertion.get_right_flank_start_0based(), 10515040)
        self.assertTrue(self._AAGA_insertion.get_right_flanking_sequence().startswith("AAGAAAGA"))

    def test_increment_flanking_sequence_size(self):
        """Test the Allele class increment_flanking_sequence_size method."""
        for expected_length in [1000, 3000, 10000, 30000, 100000]:
            for allele in [self._poly_a_insertion1, self._poly_a_insertion2, self._poly_a_deletion1, self._poly_a_deletion2]:
                allele.increase_flanking_sequence_size()
                left_flank = allele.get_left_flanking_sequence()
                right_flank = allele.get_right_flanking_sequence()

                # Flanking sequences may be shorter than expected if N's are encountered
                self.assertLessEqual(len(left_flank), expected_length)
                self.assertLessEqual(len(right_flank), expected_length)

                # Verify no N's in flanking sequences
                self.assertNotIn("N", left_flank)
                self.assertNotIn("N", right_flank)

                self.assertEqual(allele.get_right_flank_end(), allele.get_right_flank_start_0based() + len(right_flank))
                self.assertEqual(allele.get_left_flank_end(), allele.get_left_flank_start_0based() + len(left_flank))

    def test_allele_variant_id(self):
        """Test variant ID generation."""
        allele = self._poly_a_insertion1
        variant_id = allele.variant_id
        self.assertIn("chr22", variant_id)
        self.assertIn("10513201", variant_id)

        shortened_id = allele.shortened_variant_id
        self.assertIsNotNone(shortened_id)
        self.assertIn("chr22", shortened_id)

    def test_allele_lazy_loading_flanking_sequences(self):
        """Test that flanking sequences are loaded lazily."""
        allele = Allele("chr22", 10513201, "T", "TAAAAA", self._fasta_obj)
        # Access left flanking sequence
        left_seq = allele.get_left_flanking_sequence()
        self.assertIsNotNone(left_seq)
        # Access again - should use cached value
        left_seq2 = allele.get_left_flanking_sequence()
        self.assertEqual(left_seq, left_seq2)

    def test_allele_info_field_dict(self):
        """Test info field dict handling."""
        info_dict = {"AC": 1, "AF": 0.5}
        allele = Allele("chr22", 10513201, "T", "TAAAAA", self._fasta_obj, info_field_dict=info_dict)
        self.assertEqual(allele.info_field_dict, info_dict)

    def test_allele_previously_increased_flanking_sequence_size(self):
        """Test tracking of flanking sequence size increases."""
        allele = Allele("chr22", 10513201, "T", "TAAAAA", self._fasta_obj)
        self.assertFalse(allele.previously_increased_flanking_sequence_size)
        allele.increase_flanking_sequence_size()
        self.assertTrue(allele.previously_increased_flanking_sequence_size)
        # increase_flanking_sequence_size() increases both left and right, so count is 2
        self.assertEqual(allele.number_of_times_flanking_sequence_size_was_increased, 2)

    def test_flanking_sequences_exclude_ns_left(self):
        """Test that N's in left flanking sequence are excluded."""
        # Create mock fasta object with N's in the left flank
        mock_fasta = mock.MagicMock()
        mock_fasta.faidx.one_based_attributes = False
        mock_chrom = mock.MagicMock()

        # Sequence with N's: ...NNNACGT[variant at 1000]...
        # When retrieving left flank, it should stop at the N
        left_sequence = "ACGT" + "N" * 10 + "ACGT"
        right_sequence = "G" * 300

        def mock_getitem(slice_obj):
            if slice_obj.start < 1000:
                # Left flanking region
                start_offset = max(0, 1000 - slice_obj.stop)
                return left_sequence[start_offset:]
            else:
                # Right flanking region
                return right_sequence[:slice_obj.stop - slice_obj.start]

        mock_chrom.__getitem__.side_effect = mock_getitem
        mock_chrom.__len__.return_value = 100000
        mock_fasta.__getitem__.return_value = mock_chrom

        allele = Allele("chr1", 1000, "A", "AAAAA", mock_fasta)
        left_flank = allele.get_left_flanking_sequence()

        # Left flanking sequence should not contain any N's
        self.assertNotIn("N", left_flank)
        # Should have been truncated to just "ACGT" (after the N's)
        self.assertTrue(left_flank.endswith("ACGT"))

    def test_flanking_sequences_exclude_ns_right(self):
        """Test that N's in right flanking sequence are excluded."""
        # Create mock fasta object with N's in the right flank
        mock_fasta = mock.MagicMock()
        mock_fasta.faidx.one_based_attributes = False
        mock_chrom = mock.MagicMock()

        # Sequence with N's: ...[variant at 1000]ACGTNNN...
        # When retrieving right flank, it should stop at the N
        left_sequence = "G" * 300
        right_sequence = "ACGT" + "N" * 10 + "ACGT"

        def mock_getitem(slice_obj):
            if slice_obj.start < 1000:
                # Left flanking region
                return left_sequence[-(slice_obj.stop - slice_obj.start):]
            else:
                # Right flanking region
                return right_sequence[:slice_obj.stop - slice_obj.start]

        mock_chrom.__getitem__.side_effect = mock_getitem
        mock_chrom.__len__.return_value = 100000
        mock_fasta.__getitem__.return_value = mock_chrom

        allele = Allele("chr1", 1000, "A", "AAAAA", mock_fasta)
        right_flank = allele.get_right_flanking_sequence()

        # Right flanking sequence should not contain any N's
        self.assertNotIn("N", right_flank)
        # Should have been truncated to just "ACGT" (before the N's)
        self.assertEqual(right_flank, "ACGT")

    def test_flanking_sequences_no_ns(self):
        """Test that flanking sequences without N's are unchanged."""
        # Create mock fasta object without N's
        mock_fasta = mock.MagicMock()
        mock_fasta.faidx.one_based_attributes = False
        mock_chrom = mock.MagicMock()

        left_sequence = "A" * 1000
        right_sequence = "G" * 1000

        def mock_getitem(slice_obj):
            if slice_obj.start < 1000:
                # Left flanking region
                return left_sequence[-(slice_obj.stop - slice_obj.start):]
            else:
                # Right flanking region
                return right_sequence[:slice_obj.stop - slice_obj.start]

        mock_chrom.__getitem__.side_effect = mock_getitem
        mock_chrom.__len__.return_value = 100000
        mock_fasta.__getitem__.return_value = mock_chrom

        allele = Allele("chr1", 1000, "A", "AAAAA", mock_fasta)
        left_flank = allele.get_left_flanking_sequence()
        right_flank = allele.get_right_flanking_sequence()

        # Both flanking sequences should be 300bp (default size)
        self.assertEqual(len(left_flank), 300)
        self.assertEqual(len(right_flank), 300)
        # Should not contain any N's
        self.assertNotIn("N", left_flank)
        self.assertNotIn("N", right_flank)

    def test_tandem_repeat_allele_getters(self):
        """Test the TandemRepeatAllele class getters"""
        for tandem_repeat_allele in [self._AAGA_tandem_repeat_expansion, self._AAGA_tandem_repeat_contraction]:
            self.assertEqual(tandem_repeat_allele.start_0based, 10515040)
            self.assertEqual(tandem_repeat_allele.end_1based, 10515077)
            self.assertEqual(tandem_repeat_allele.detection_mode, DETECTION_MODE_PURE_REPEATS)
            self.assertEqual(tandem_repeat_allele.repeat_unit, "AAGA")
            self.assertEqual(tandem_repeat_allele.num_repeats_in_variant, 1)
            self.assertEqual(tandem_repeat_allele.num_repeat_bases_in_variant, 4)
            self.assertEqual(tandem_repeat_allele.num_repeat_bases_in_left_flank, 0)

            self.assertEqual(tandem_repeat_allele.num_repeats_ref, 9)
            if tandem_repeat_allele.ins_or_del == "INS":
                self.assertEqual(tandem_repeat_allele.end_1based, tandem_repeat_allele.start_0based + tandem_repeat_allele.num_repeat_bases_in_left_flank + tandem_repeat_allele.num_repeat_bases_in_right_flank)
                self.assertEqual(tandem_repeat_allele.num_repeats_alt, 10)
                self.assertEqual(tandem_repeat_allele.num_repeat_bases_in_right_flank, 37)
                self.assertEqual(tandem_repeat_allele.num_repeats_in_right_flank, 9)
                self.assertEqual(tandem_repeat_allele.ref_allele_repeat_sequence, "AAGA"*9 + "A")
                self.assertEqual(tandem_repeat_allele.alt_allele_repeat_sequence, "AAGA"*10 + "A")
            else:
                self.assertEqual(tandem_repeat_allele.end_1based, tandem_repeat_allele.start_0based + tandem_repeat_allele.num_repeat_bases_in_left_flank + tandem_repeat_allele.num_repeat_bases_in_variant + tandem_repeat_allele.num_repeat_bases_in_right_flank)
                self.assertEqual(tandem_repeat_allele.num_repeats_alt, 8)
                self.assertEqual(tandem_repeat_allele.num_repeat_bases_in_right_flank, 33)
                self.assertEqual(tandem_repeat_allele.num_repeats_in_right_flank, 8)
                self.assertEqual(tandem_repeat_allele.ref_allele_repeat_sequence, "AAGA"*9 + "A")
                self.assertEqual(tandem_repeat_allele.alt_allele_repeat_sequence, "AAGA"*8 + "A")

            self.assertEqual(tandem_repeat_allele.num_repeats_in_left_flank, 0)

            self.assertFalse(tandem_repeat_allele.do_repeats_cover_entire_left_flanking_sequence())
            self.assertFalse(tandem_repeat_allele.do_repeats_cover_entire_right_flanking_sequence())

    def tearDown(self):
        """Tear down the test case."""
        self._fasta_obj.close()
        self._fasta_obj_with_one_based_attributes.close()
        if os.path.exists(self._temp_fasta_path):
            os.unlink(self._temp_fasta_path)
        fai_path = self._temp_fasta_path + ".fai"
        if os.path.exists(fai_path):
            os.unlink(fai_path)
        if os.path.exists(self._trf_working_dir):
            shutil.rmtree(self._trf_working_dir)


class TestTandemRepeatAllele(unittest.TestCase):
    """Test the TandemRepeatAllele class comprehensive functionality."""

    def setUp(self):
        """Set up test case."""
        fasta_data = pkgutil.get_data("str_analysis", "data/tests/chr22_11Mb.fa.gz")
        with tempfile.NamedTemporaryFile(suffix=".fa.gz", delete=False) as fasta_file:
            self._temp_fasta_path = fasta_file.name
            fasta_file.write(fasta_data)
            fasta_file.flush()

        self._fasta_obj = pyfaidx.Fasta(self._temp_fasta_path, one_based_attributes=False, as_raw=True)
        self._poly_a_insertion = Allele("chr22", 10513201, "T", "TAAAAA", self._fasta_obj)

    def test_tandem_repeat_allele_with_purity_adjustment(self):
        """Test TandemRepeatAllele with purity adjustment enabled."""
        allele = Allele("chr22", 10515040, "T", "TAAGA", self._fasta_obj)
        tr_allele = TandemRepeatAllele(allele, "AAGA", True, 0, 4, 37, DETECTION_MODE_PURE_REPEATS)
        self.assertEqual(tr_allele.repeat_unit, "AAGA")
        self.assertGreater(tr_allele.repeat_purity, 0.9)

    def test_tandem_repeat_allele_without_purity_adjustment(self):
        """Test TandemRepeatAllele without purity adjustment."""
        allele = Allele("chr22", 10515040, "T", "TAAGA", self._fasta_obj)
        tr_allele = TandemRepeatAllele(allele, "AAGA", False, 0, 4, 37, DETECTION_MODE_PURE_REPEATS)
        self.assertEqual(tr_allele.repeat_unit, "AAGA")

    def test_repeat_purity_calculation(self):
        """Test repeat purity computation."""
        allele = Allele("chr22", 10515040, "T", "TAAGA", self._fasta_obj)
        tr_allele = TandemRepeatAllele(allele, "AAGA", False, 0, 4, 37, DETECTION_MODE_PURE_REPEATS)
        purity = tr_allele.repeat_purity
        self.assertGreater(purity, 0.0)
        self.assertLessEqual(purity, 1.0)

    def test_is_pure_repeat(self):
        """Test is_pure_repeat property."""
        allele = Allele("chr22", 10515040, "T", "TAAGA", self._fasta_obj)
        tr_allele = TandemRepeatAllele(allele, "AAGA", False, 0, 4, 37, DETECTION_MODE_PURE_REPEATS)
        # This should be a pure repeat based on the test data
        is_pure = tr_allele.is_pure_repeat
        self.assertIsInstance(is_pure, bool)

    def test_canonical_repeat_unit(self):
        """Test canonical repeat unit computation."""
        allele = Allele("chr22", 10515040, "T", "TAAGA", self._fasta_obj)
        tr_allele = TandemRepeatAllele(allele, "AAGA", False, 0, 4, 37, DETECTION_MODE_PURE_REPEATS)
        canonical = tr_allele.canonical_repeat_unit
        self.assertIsNotNone(canonical)
        self.assertIsInstance(canonical, str)

    def test_locus_id_generation(self):
        """Test locus ID format."""
        allele = Allele("chr22", 10515040, "T", "TAAGA", self._fasta_obj)
        tr_allele = TandemRepeatAllele(allele, "AAGA", False, 0, 4, 37, DETECTION_MODE_PURE_REPEATS)
        locus_id = tr_allele.locus_id
        self.assertIn("chr22", locus_id)
        self.assertIn("AAGA", locus_id)

    def test_summary_string(self):
        """Test summary string generation."""
        allele = Allele("chr22", 10515040, "T", "TAAGA", self._fasta_obj)
        tr_allele = TandemRepeatAllele(allele, "AAGA", False, 0, 4, 37, DETECTION_MODE_PURE_REPEATS)
        summary = tr_allele.summary_string
        self.assertIsNotNone(summary)
        self.assertIsInstance(summary, str)

    def test_variant_and_flanks_repeat_sequence(self):
        """Test variant_and_flanks_repeat_sequence property."""
        allele = Allele("chr22", 10515040, "T", "TAAGA", self._fasta_obj)
        tr_allele = TandemRepeatAllele(allele, "AAGA", False, 0, 4, 37, DETECTION_MODE_PURE_REPEATS)
        seq = tr_allele.variant_and_flanks_repeat_sequence
        self.assertIsInstance(seq, str)
        self.assertGreater(len(seq), 0)

    def test_flank_coverage_checks(self):
        """Test methods checking if repeats cover entire flanking sequences."""
        allele = Allele("chr22", 10515040, "T", "TAAGA", self._fasta_obj)
        tr_allele = TandemRepeatAllele(allele, "AAGA", False, 0, 4, 37, DETECTION_MODE_PURE_REPEATS)

        covers_left = tr_allele.do_repeats_cover_entire_left_flanking_sequence()
        covers_right = tr_allele.do_repeats_cover_entire_right_flanking_sequence()
        covers_both = tr_allele.do_repeats_cover_entire_flanking_sequence()

        self.assertIsInstance(covers_left, bool)
        self.assertIsInstance(covers_right, bool)
        self.assertIsInstance(covers_both, bool)

    def tearDown(self):
        """Tear down test case."""
        self._fasta_obj.close()
        if os.path.exists(self._temp_fasta_path):
            os.unlink(self._temp_fasta_path)
        fai_path = self._temp_fasta_path + ".fai"
        if os.path.exists(fai_path):
            os.unlink(fai_path)


class TestReferenceTandemRepeat(unittest.TestCase):
    """Test the ReferenceTandemRepeat class."""

    def test_constructor_valid(self):
        """Test constructor with valid inputs."""
        ref_tr = ReferenceTandemRepeat("chr1", 1000, 1100, "CAG", DETECTION_MODE_PURE_REPEATS)
        self.assertEqual(ref_tr.chrom, "chr1")
        self.assertEqual(ref_tr.start_0based, 1000)
        self.assertEqual(ref_tr.end_1based, 1100)
        self.assertEqual(ref_tr.repeat_unit, "CAG")
        self.assertEqual(ref_tr.detection_mode, DETECTION_MODE_PURE_REPEATS)

    def test_constructor_invalid_coordinates(self):
        """Test constructor rejects start > end."""
        with self.assertRaises(ValueError):
            ReferenceTandemRepeat("chr1", 1100, 1000, "CAG")

    def test_repeat_unit_length(self):
        """Test repeat_unit_length property."""
        ref_tr = ReferenceTandemRepeat("chr1", 1000, 1100, "CAG")
        self.assertEqual(ref_tr.repeat_unit_length, 3)

    def test_ref_interval_size(self):
        """Test ref_interval_size calculation."""
        ref_tr = ReferenceTandemRepeat("chr1", 1000, 1100, "CAG")
        self.assertEqual(ref_tr.ref_interval_size, 100)

    def test_num_repeats_ref(self):
        """Test num_repeats_ref calculation."""
        ref_tr = ReferenceTandemRepeat("chr1", 1000, 1099, "CAG")  # 99 bp = 33 repeats
        self.assertEqual(ref_tr.num_repeats_ref, 33)

    def test_locus_id_format(self):
        """Test locus_id format."""
        ref_tr = ReferenceTandemRepeat("chr1", 1000, 1100, "CAG")
        locus_id = ref_tr.locus_id
        self.assertEqual(locus_id, "chr1-1000-1100-CAG")

    def test_canonical_repeat_unit(self):
        """Test canonical_repeat_unit computation."""
        ref_tr = ReferenceTandemRepeat("chr1", 1000, 1100, "CAG")
        canonical = ref_tr.canonical_repeat_unit
        self.assertIsNotNone(canonical)
        self.assertIsInstance(canonical, str)

    def test_summary_string_short_motif(self):
        """Test summary string for short motif."""
        ref_tr = ReferenceTandemRepeat("chr1", 1000, 1099, "CAG", DETECTION_MODE_PURE_REPEATS)
        summary = ref_tr.summary_string
        self.assertIn("CAG", summary)
        self.assertIn("pure", summary)

    def test_summary_string_long_motif(self):
        """Test summary string truncates long motifs."""
        long_motif = "A" * 50
        ref_tr = ReferenceTandemRepeat("chr1", 1000, 1100, long_motif, DETECTION_MODE_TRF)
        summary = ref_tr.summary_string
        self.assertIn("...", summary)  # Should be truncated

    def test_str_and_repr(self):
        """Test __str__ and __repr__ methods."""
        ref_tr = ReferenceTandemRepeat("chr1", 1000, 1100, "CAG")
        str_rep = str(ref_tr)
        repr_rep = repr(ref_tr)
        self.assertEqual(str_rep, "chr1-1000-1100-CAG")
        self.assertEqual(repr_rep, "chr1-1000-1100-CAG")


class TestDetectionFunctions(unittest.TestCase):
    """Test tandem repeat detection functions."""

    def setUp(self):
        """Set up test case."""
        fasta_data = pkgutil.get_data("str_analysis", "data/tests/chr22_11Mb.fa.gz")
        with tempfile.NamedTemporaryFile(suffix=".fa.gz", delete=False) as fasta_file:
            self._temp_fasta_path = fasta_file.name
            fasta_file.write(fasta_data)
            fasta_file.flush()

        self._fasta_obj = pyfaidx.Fasta(self._temp_fasta_path, one_based_attributes=False, as_raw=True)

        self._args = argparse.Namespace(
            min_repeat_unit_length=1,
            max_repeat_unit_length=1000,
            min_repeats=3,
            min_tandem_repeat_length=9,
            debug=False,
            trf_min_repeats_in_reference=2,
            trf_min_purity=0.2,
        )

    def test_check_if_allele_is_tandem_repeat_pure_mode(self):
        """Test check_if_allele_is_tandem_repeat with DETECTION_MODE_PURE_REPEATS."""
        allele = Allele("chr22", 10515040, "T", "TAAGA", self._fasta_obj)
        tr_allele, filter_reason = check_if_allele_is_tandem_repeat(allele, self._args, DETECTION_MODE_PURE_REPEATS)
        # Should detect this as a tandem repeat
        self.assertTrue(tr_allele is not None or filter_reason is not None)

    def test_check_if_allele_is_tandem_repeat_interrupted_mode(self):
        """Test check_if_allele_is_tandem_repeat with DETECTION_MODE_ALLOW_INTERRUPTIONS."""
        allele = Allele("chr22", 10515040, "T", "TAAGA", self._fasta_obj)
        tr_allele, filter_reason = check_if_allele_is_tandem_repeat(allele, self._args, DETECTION_MODE_ALLOW_INTERRUPTIONS)
        self.assertTrue(tr_allele is not None or filter_reason is not None)

    def test_check_if_allele_is_tandem_repeat_invalid_mode(self):
        """Test check_if_allele_is_tandem_repeat rejects invalid mode."""
        allele = Allele("chr22", 10515040, "T", "TAAGA", self._fasta_obj)
        with self.assertRaises(ValueError):
            check_if_allele_is_tandem_repeat(allele, self._args, "invalid_mode")

    def test_check_if_allele_is_tandem_repeat_return_tuple(self):
        """Test return tuple structure from check_if_allele_is_tandem_repeat."""
        allele = Allele("chr22", 10515040, "T", "TAAGA", self._fasta_obj)
        result = check_if_allele_is_tandem_repeat(allele, self._args, DETECTION_MODE_PURE_REPEATS)
        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 2)

    def test_detect_perfect_and_almost_perfect_tandem_repeats(self):
        """Test detect_perfect_and_almost_perfect_tandem_repeats function."""
        alleles = [Allele("chr22", 10515040, "T", "TAAGA", self._fasta_obj)]
        counters = collections.defaultdict(int)

        args = argparse.Namespace(
            min_repeat_unit_length=1,
            max_repeat_unit_length=1000,
            min_repeats=3,
            min_tandem_repeat_length=9,
            debug=False,
            show_progress_bar=False,
            verbose=False,
            dont_allow_interruptions=False,
            dont_run_trf=True,
            min_indel_size_to_run_trf=7,
        )

        tr_alleles, trf_queue = detect_perfect_and_almost_perfect_tandem_repeats(alleles, counters, args)

        self.assertIsInstance(tr_alleles, list)
        self.assertIsInstance(trf_queue, list)

    def tearDown(self):
        """Tear down test case."""
        self._fasta_obj.close()
        if os.path.exists(self._temp_fasta_path):
            os.unlink(self._temp_fasta_path)
        fai_path = self._temp_fasta_path + ".fai"
        if os.path.exists(fai_path):
            os.unlink(fai_path)


class TestFilterFunctions(unittest.TestCase):
    """Test filter functions for tandem repeat alleles."""

    def setUp(self):
        """Set up test case."""
        fasta_data = pkgutil.get_data("str_analysis", "data/tests/chr22_11Mb.fa.gz")
        with tempfile.NamedTemporaryFile(suffix=".fa.gz", delete=False) as fasta_file:
            self._temp_fasta_path = fasta_file.name
            fasta_file.write(fasta_data)
            fasta_file.flush()

        self._fasta_obj = pyfaidx.Fasta(self._temp_fasta_path, one_based_attributes=False, as_raw=True)

        self._args = argparse.Namespace(
            min_repeat_unit_length=2,
            max_repeat_unit_length=100,
            min_repeats=3,
            min_tandem_repeat_length=9,
            trf_min_repeats_in_reference=2,
            trf_min_purity=0.5,
        )

    def test_filter_passes_all(self):
        """Test allele that passes all filters."""
        allele = Allele("chr22", 10515040, "T", "TAAGA", self._fasta_obj)
        tr_allele = TandemRepeatAllele(allele, "AAGA", False, 0, 4, 37, DETECTION_MODE_PURE_REPEATS)
        result = check_if_tandem_repeat_allele_failed_filters(self._args, tr_allele, detected_by_trf=False)
        self.assertIsNone(result)  # None means it passed

    def test_filter_not_enough_repeats(self):
        """Test filter for minimum repeats."""
        allele = Allele("chr22", 10515040, "T", "TAAGA", self._fasta_obj)
        tr_allele = TandemRepeatAllele(allele, "AAGA", False, 0, 4, 4, DETECTION_MODE_PURE_REPEATS)  # Only 2 total repeats
        result = check_if_tandem_repeat_allele_failed_filters(self._args, tr_allele, detected_by_trf=False)
        self.assertIsNotNone(result)
        self.assertIn("contains <", result)

    def test_filter_repeat_unit_too_short(self):
        """Test filter for minimum repeat unit length."""
        allele = Allele("chr22", 10515040, "T", "TA", self._fasta_obj)
        tr_allele = TandemRepeatAllele(allele, "A", False, 0, 1, 10, DETECTION_MODE_PURE_REPEATS)
        result = check_if_tandem_repeat_allele_failed_filters(self._args, tr_allele, detected_by_trf=False)
        self.assertIsNotNone(result)
        self.assertIn("repeat unit <", result)

    def test_filter_repeat_unit_too_long(self):
        """Test filter for maximum repeat unit length."""
        long_motif = "A" * 150
        allele = Allele("chr22", 10515040, "T", "T" + long_motif, self._fasta_obj)
        # Note: The TandemRepeatAllele may normalize the repeat unit, so we use a longer motif
        # and expect it to be filtered for being too long
        tr_allele = TandemRepeatAllele(allele, long_motif, False, 0, len(long_motif), len(long_motif), DETECTION_MODE_PURE_REPEATS)
        result = check_if_tandem_repeat_allele_failed_filters(self._args, tr_allele, detected_by_trf=False)
        self.assertIsNotNone(result)
        # Should fail either for being too long or not enough repeats (since it's a single "A")
        self.assertTrue("repeat unit >" in result or "contains <" in result or "INDEL without repeats" in result)

    def test_filter_trf_specific_purity(self):
        """Test TRF-specific purity filter."""
        allele = Allele("chr22", 10515040, "T", "TAAGA", self._fasta_obj)
        tr_allele = TandemRepeatAllele(allele, "AAGA", False, 0, 4, 37, DETECTION_MODE_TRF)

        # Mock low purity
        with mock.patch.object(TandemRepeatAllele, 'repeat_purity', new_callable=mock.PropertyMock) as mock_purity:
            mock_purity.return_value = 0.1  # Below threshold of 0.5
            result = check_if_tandem_repeat_allele_failed_filters(self._args, tr_allele, detected_by_trf=True)
            self.assertIsNotNone(result)
            self.assertIn("purity", result)

    def tearDown(self):
        """Tear down test case."""
        self._fasta_obj.close()
        if os.path.exists(self._temp_fasta_path):
            os.unlink(self._temp_fasta_path)
        fai_path = self._temp_fasta_path + ".fai"
        if os.path.exists(fai_path):
            os.unlink(fai_path)


class TestMergeFunctions(unittest.TestCase):
    """Test merge and utility functions."""

    def setUp(self):
        """Set up test case."""
        fasta_data = pkgutil.get_data("str_analysis", "data/tests/chr22_11Mb.fa.gz")
        with tempfile.NamedTemporaryFile(suffix=".fa.gz", delete=False) as fasta_file:
            self._temp_fasta_path = fasta_file.name
            fasta_file.write(fasta_data)
            fasta_file.flush()

        self._fasta_obj = pyfaidx.Fasta(self._temp_fasta_path, one_based_attributes=False, as_raw=True)

    def test_compute_repeat_unit_id_short_motif(self):
        """Test compute_repeat_unit_id for short motifs."""
        result = compute_repeat_unit_id("CAG")
        self.assertEqual(result, "CAG")

    def test_compute_repeat_unit_id_long_motif(self):
        """Test compute_repeat_unit_id for long motifs (>6bp)."""
        long_motif = "ACGTACGT"  # 8bp
        result = compute_repeat_unit_id(long_motif)
        self.assertEqual(result, 8)

    def test_compute_repeat_unit_id_boundary(self):
        """Test compute_repeat_unit_id at 6bp boundary."""
        motif = "ACGTAC"  # Exactly 6bp
        result = compute_repeat_unit_id(motif)
        self.assertEqual(result, "ACGTAC")

    def test_are_repeat_units_similar_same_short(self):
        """Test are_repeat_units_similar for same short motifs."""
        result = are_repeat_units_similar("CAG", "CAG")
        self.assertTrue(result)

    def test_are_repeat_units_similar_different_short(self):
        """Test are_repeat_units_similar for different short motifs."""
        result = are_repeat_units_similar("CAG", "CTG")
        self.assertFalse(result)

    def test_are_repeat_units_similar_same_length_long(self):
        """Test are_repeat_units_similar for same length long motifs."""
        motif1 = "ACGTACGT"  # 8bp
        motif2 = "TGCATGCA"  # 8bp
        result = are_repeat_units_similar(motif1, motif2)
        self.assertTrue(result)  # Same length, so similar

    def test_are_repeat_units_similar_different_length_long(self):
        """Test are_repeat_units_similar for different length long motifs."""
        motif1 = "ACGTACGT"  # 8bp
        motif2 = "ACGTACGTA"  # 9bp
        result = are_repeat_units_similar(motif1, motif2)
        self.assertFalse(result)

    def test_merge_overlapping_tandem_repeat_loci_non_overlapping(self):
        """Test merge with non-overlapping loci."""
        allele1 = Allele("chr22", 10515040, "T", "TAAGA", self._fasta_obj)
        allele2 = Allele("chr22", 10516040, "T", "TAAGA", self._fasta_obj)

        tr1 = TandemRepeatAllele(allele1, "AAGA", False, 0, 4, 37, DETECTION_MODE_PURE_REPEATS)
        tr2 = TandemRepeatAllele(allele2, "AAGA", False, 0, 4, 37, DETECTION_MODE_PURE_REPEATS)

        result = merge_overlapping_tandem_repeat_loci([tr1, tr2], self._fasta_obj, verbose=False)

        self.assertIsInstance(result, list)
        # Non-overlapping should result in separate loci
        self.assertGreaterEqual(len(result), 1)

    def test_merge_overlapping_tandem_repeat_loci_single_allele(self):
        """Test merge with single allele (should pass through)."""
        allele = Allele("chr22", 10515040, "T", "TAAGA", self._fasta_obj)
        tr = TandemRepeatAllele(allele, "AAGA", False, 0, 4, 37, DETECTION_MODE_PURE_REPEATS)

        result = merge_overlapping_tandem_repeat_loci([tr], self._fasta_obj, verbose=False)

        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        # When there's only one allele, it returns the original TandemRepeatAllele, not a ReferenceTandemRepeat
        self.assertIsInstance(result[0], TandemRepeatAllele)

    def test_need_to_reprocess_allele_no_coverage(self):
        """Test need_to_reprocess when repeats don't cover flanks."""
        allele = Allele("chr22", 10515040, "T", "TAAGA", self._fasta_obj)
        tr_allele = TandemRepeatAllele(allele, "AAGA", False, 0, 4, 37, DETECTION_MODE_PURE_REPEATS)

        result = need_to_reprocess_allele_with_extended_flanking_sequence(tr_allele)

        self.assertFalse(result)  # Repeats don't cover entire flank

    def tearDown(self):
        """Tear down test case."""
        self._fasta_obj.close()
        if os.path.exists(self._temp_fasta_path):
            os.unlink(self._temp_fasta_path)
        fai_path = self._temp_fasta_path + ".fai"
        if os.path.exists(fai_path):
            os.unlink(fai_path)


class TestTRFIntegration(unittest.TestCase):
    """Test TRF integration with mocking."""

    def setUp(self):
        """Set up test case."""
        fasta_data = pkgutil.get_data("str_analysis", "data/tests/chr22_11Mb.fa.gz")
        with tempfile.NamedTemporaryFile(suffix=".fa.gz", delete=False) as fasta_file:
            self._temp_fasta_path = fasta_file.name
            fasta_file.write(fasta_data)
            fasta_file.flush()

        self._fasta_obj = pyfaidx.Fasta(self._temp_fasta_path, one_based_attributes=False, as_raw=True)

        self._trf_working_dir = tempfile.mkdtemp()
        self._args = argparse.Namespace(
            min_repeat_unit_length=1,
            max_repeat_unit_length=1000,
            min_repeats=3,
            min_tandem_repeat_length=9,
            debug=False,
            trf_working_dir=self._trf_working_dir,
            input_vcf_prefix="test",
            trf_executable_path="trf",
            trf_threads=2,
            verbose=False,
            show_progress_bar=False,
            allow_multiple_trf_results_per_locus=False,
            dont_allow_interruptions=False,
            dont_run_trf=False,
            min_indel_size_to_run_trf=7,
            trf_min_repeats_in_reference=2,
            trf_min_purity=0.2,
        )

    @mock.patch('str_analysis.filter_vcf_to_tandem_repeats.TRFRunner')
    def test_detect_tandem_repeats_using_trf_with_mock(self, mock_trf_runner_class):
        """Test detect_tandem_repeats_using_trf with mocked TRFRunner."""
        # Mock TRFRunner instance
        mock_trf_instance = mock.Mock()
        mock_trf_runner_class.return_value = mock_trf_instance

        # Mock the run_trf_on_fasta_file method
        mock_trf_instance.run_trf_on_fasta_file.return_value = None

        # Mock parse_html_results to return empty list (no TRF results)
        mock_trf_instance.parse_html_results.return_value = []

        alleles = [Allele("chr22", 10515040, "T", "TAAGAAAGA", self._fasta_obj)]
        counters = collections.defaultdict(int)

        try:
            result = detect_tandem_repeats_using_trf(alleles, counters, self._args)
            self.assertIsInstance(result, list)
        except (FileNotFoundError, OSError) as e:
            # TRF executable not available - skip the test
            # OSError covers cases where the binary exists but can't execute
            self.skipTest(f"TRF not available: {e}")

    def tearDown(self):
        """Tear down test case."""
        self._fasta_obj.close()
        if os.path.exists(self._temp_fasta_path):
            os.unlink(self._temp_fasta_path)
        fai_path = self._temp_fasta_path + ".fai"
        if os.path.exists(fai_path):
            os.unlink(fai_path)
        if os.path.exists(self._trf_working_dir):
            shutil.rmtree(self._trf_working_dir)


class TestChromosomeNamingFunctions(unittest.TestCase):
    """Test chromosome naming detection and normalization functions."""

    def test_normalize_chromosome_name_add_chr_prefix(self):
        """Test adding chr prefix when target convention is 'chr'."""
        self.assertEqual(normalize_chromosome_name("1", "chr"), "chr1")
        self.assertEqual(normalize_chromosome_name("X", "chr"), "chrX")
        self.assertEqual(normalize_chromosome_name("22", "chr"), "chr22")

    def test_normalize_chromosome_name_already_has_chr(self):
        """Test chromosome name that already has chr prefix."""
        self.assertEqual(normalize_chromosome_name("chr1", "chr"), "chr1")
        self.assertEqual(normalize_chromosome_name("chrX", "chr"), "chrX")

    def test_normalize_chromosome_name_remove_chr_prefix(self):
        """Test removing chr prefix when target convention is 'no_chr'."""
        self.assertEqual(normalize_chromosome_name("chr1", "no_chr"), "1")
        self.assertEqual(normalize_chromosome_name("chrX", "no_chr"), "X")
        self.assertEqual(normalize_chromosome_name("chr22", "no_chr"), "22")

    def test_normalize_chromosome_name_no_chr_already(self):
        """Test chromosome name that already doesn't have chr prefix."""
        self.assertEqual(normalize_chromosome_name("1", "no_chr"), "1")
        self.assertEqual(normalize_chromosome_name("X", "no_chr"), "X")

    def test_normalize_chromosome_name_unknown_convention(self):
        """Test normalization with unknown convention returns as-is."""
        self.assertEqual(normalize_chromosome_name("chr1", None), "chr1")
        self.assertEqual(normalize_chromosome_name("1", None), "1")


class TestVCFFunctions(unittest.TestCase):
    """Test VCF-related functions for genotyping."""

    def setUp(self):
        """Set up test fixtures with a small VCF file."""
        # Create a temporary VCF file for testing
        self._temp_vcf_path = None
        self._temp_vcf_gz_path = None

    def _create_temp_vcf(self, vcf_content, compress=False):
        """Helper to create temporary VCF file."""
        import gzip
        suffix = ".vcf.gz" if compress else ".vcf"
        with tempfile.NamedTemporaryFile(suffix=suffix, delete=False, mode='wb' if compress else 'w') as f:
            if compress:
                with gzip.open(f.name, 'wt') as gz:
                    gz.write(vcf_content)
            else:
                f.write(vcf_content)
            return f.name

    def test_open_vcf_for_genotyping_valid_single_sample(self):
        """Test opening a valid single-sample VCF."""
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t1000\t.\tA\tAT\t.\tPASS\t.\tGT\t0/1
"""
        vcf_path = self._create_temp_vcf(vcf_content)
        try:
            vcf_file, sample_name, chrom_convention = open_vcf_for_genotyping(vcf_path)
            self.assertEqual(sample_name, "SAMPLE1")
            self.assertEqual(chrom_convention, "chr")
            vcf_file.close()
        finally:
            os.unlink(vcf_path)

    def test_open_vcf_for_genotyping_no_chr_prefix(self):
        """Test detecting chromosome naming convention without chr prefix."""
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=1,length=1000000>
##contig=<ID=2,length=1000000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
1\t1000\t.\tA\tAT\t.\tPASS\t.\tGT\t0/1
"""
        vcf_path = self._create_temp_vcf(vcf_content)
        try:
            vcf_file, sample_name, chrom_convention = open_vcf_for_genotyping(vcf_path)
            self.assertEqual(chrom_convention, "no_chr")
            vcf_file.close()
        finally:
            os.unlink(vcf_path)

    def test_open_vcf_for_genotyping_multi_sample_error(self):
        """Test that multi-sample VCF raises error."""
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2
chr1\t1000\t.\tA\tAT\t.\tPASS\t.\tGT\t0/1\t1/1
"""
        vcf_path = self._create_temp_vcf(vcf_content)
        try:
            with self.assertRaises(ValueError) as context:
                open_vcf_for_genotyping(vcf_path)
            self.assertIn("single-sample", str(context.exception))
        finally:
            os.unlink(vcf_path)

    def test_open_vcf_for_genotyping_file_not_found(self):
        """Test that missing VCF raises FileNotFoundError."""
        with self.assertRaises(FileNotFoundError):
            open_vcf_for_genotyping("/nonexistent/path/to/file.vcf")

    def test_get_overlapping_vcf_variants_basic(self):
        """Test fetching overlapping variants from a VCF."""
        import pysam
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t1000\t.\tA\tAT\t.\tPASS\t.\tGT\t0/1
chr1\t1050\t.\tG\tC\t.\tPASS\t.\tGT\t1/1
chr1\t2000\t.\tC\tCA\t.\tPASS\t.\tGT\t0/1
"""
        vcf_path = self._create_temp_vcf(vcf_content)
        try:
            # Need to compress and index for fetch to work
            import subprocess
            vcf_gz_path = vcf_path + ".gz"
            with open(vcf_gz_path, 'wb') as gz_out:
                subprocess.run(["bgzip", "-c", vcf_path], stdout=gz_out, check=True)
            subprocess.run(["tabix", "-p", "vcf", vcf_gz_path], check=True)

            vcf_file = pysam.VariantFile(vcf_gz_path)

            # Fetch variants overlapping region 990-1060 (should get 2 variants)
            variants, has_multiallelic = get_overlapping_vcf_variants(
                vcf_file, "chr1", 990, 1060, vcf_chrom_convention="chr"
            )

            self.assertEqual(len(variants), 2)
            self.assertFalse(has_multiallelic)
            self.assertEqual(variants[0].pos, 1000)
            self.assertEqual(variants[1].pos, 1050)

            vcf_file.close()
        except FileNotFoundError:
            # bgzip/tabix not available - skip test
            self.skipTest("bgzip/tabix not available")
        finally:
            os.unlink(vcf_path)
            if os.path.exists(vcf_path + ".gz"):
                os.unlink(vcf_path + ".gz")
            if os.path.exists(vcf_path + ".gz.tbi"):
                os.unlink(vcf_path + ".gz.tbi")

    def test_get_overlapping_vcf_variants_multiallelic(self):
        """Test detection of multi-allelic variants."""
        import pysam
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t1000\t.\tA\tAT,ATT,ATTT\t.\tPASS\t.\tGT\t1/2
"""
        vcf_path = self._create_temp_vcf(vcf_content)
        try:
            # Need to compress and index for fetch to work
            import subprocess
            vcf_gz_path = vcf_path + ".gz"
            with open(vcf_gz_path, 'wb') as gz_out:
                subprocess.run(["bgzip", "-c", vcf_path], stdout=gz_out, check=True)
            subprocess.run(["tabix", "-p", "vcf", vcf_gz_path], check=True)

            vcf_file = pysam.VariantFile(vcf_gz_path)

            variants, has_multiallelic = get_overlapping_vcf_variants(
                vcf_file, "chr1", 990, 1010, vcf_chrom_convention="chr"
            )

            self.assertEqual(len(variants), 1)
            self.assertTrue(has_multiallelic)  # Should detect multi-allelic

            vcf_file.close()
        except FileNotFoundError:
            self.skipTest("bgzip/tabix not available")
        finally:
            os.unlink(vcf_path)
            if os.path.exists(vcf_path + ".gz"):
                os.unlink(vcf_path + ".gz")
            if os.path.exists(vcf_path + ".gz.tbi"):
                os.unlink(vcf_path + ".gz.tbi")

    def test_get_overlapping_vcf_variants_no_overlap(self):
        """Test fetching from region with no overlapping variants."""
        import pysam
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t1000\t.\tA\tAT\t.\tPASS\t.\tGT\t0/1
"""
        vcf_path = self._create_temp_vcf(vcf_content)
        try:
            import subprocess
            vcf_gz_path = vcf_path + ".gz"
            with open(vcf_gz_path, 'wb') as gz_out:
                subprocess.run(["bgzip", "-c", vcf_path], stdout=gz_out, check=True)
            subprocess.run(["tabix", "-p", "vcf", vcf_gz_path], check=True)

            vcf_file = pysam.VariantFile(vcf_gz_path)

            # Fetch from region that has no variants
            variants, has_multiallelic = get_overlapping_vcf_variants(
                vcf_file, "chr1", 2000, 3000, vcf_chrom_convention="chr"
            )

            self.assertEqual(len(variants), 0)
            self.assertFalse(has_multiallelic)

            vcf_file.close()
        except FileNotFoundError:
            self.skipTest("bgzip/tabix not available")
        finally:
            os.unlink(vcf_path)
            if os.path.exists(vcf_path + ".gz"):
                os.unlink(vcf_path + ".gz")
            if os.path.exists(vcf_path + ".gz.tbi"):
                os.unlink(vcf_path + ".gz.tbi")

    def test_get_overlapping_vcf_variants_chromosome_not_in_vcf(self):
        """Test fetching from chromosome not in VCF returns empty list."""
        import pysam
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t1000\t.\tA\tAT\t.\tPASS\t.\tGT\t0/1
"""
        vcf_path = self._create_temp_vcf(vcf_content)
        try:
            import subprocess
            vcf_gz_path = vcf_path + ".gz"
            with open(vcf_gz_path, 'wb') as gz_out:
                subprocess.run(["bgzip", "-c", vcf_path], stdout=gz_out, check=True)
            subprocess.run(["tabix", "-p", "vcf", vcf_gz_path], check=True)

            vcf_file = pysam.VariantFile(vcf_gz_path)

            # Fetch from chromosome not in VCF
            variants, has_multiallelic = get_overlapping_vcf_variants(
                vcf_file, "chr2", 1000, 2000, vcf_chrom_convention="chr"
            )

            self.assertEqual(len(variants), 0)
            self.assertFalse(has_multiallelic)

            vcf_file.close()
        except FileNotFoundError:
            self.skipTest("bgzip/tabix not available")
        finally:
            os.unlink(vcf_path)
            if os.path.exists(vcf_path + ".gz"):
                os.unlink(vcf_path + ".gz")
            if os.path.exists(vcf_path + ".gz.tbi"):
                os.unlink(vcf_path + ".gz.tbi")

    def test_get_overlapping_vcf_variants_chromosome_normalization(self):
        """Test that chromosome names are normalized correctly."""
        import pysam
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=1,length=1000000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
1\t1000\t.\tA\tAT\t.\tPASS\t.\tGT\t0/1
"""
        vcf_path = self._create_temp_vcf(vcf_content)
        try:
            import subprocess
            vcf_gz_path = vcf_path + ".gz"
            with open(vcf_gz_path, 'wb') as gz_out:
                subprocess.run(["bgzip", "-c", vcf_path], stdout=gz_out, check=True)
            subprocess.run(["tabix", "-p", "vcf", vcf_gz_path], check=True)

            vcf_file = pysam.VariantFile(vcf_gz_path)

            # Query using chr1 but VCF has "1" - should normalize correctly
            variants, has_multiallelic = get_overlapping_vcf_variants(
                vcf_file, "chr1", 990, 1010, vcf_chrom_convention="no_chr"
            )

            self.assertEqual(len(variants), 1)
            self.assertEqual(variants[0].pos, 1000)

            vcf_file.close()
        except FileNotFoundError:
            self.skipTest("bgzip/tabix not available")
        finally:
            os.unlink(vcf_path)
            if os.path.exists(vcf_path + ".gz"):
                os.unlink(vcf_path + ".gz")
            if os.path.exists(vcf_path + ".gz.tbi"):
                os.unlink(vcf_path + ".gz.tbi")


class TestConvertVariantsToHaplotypeSequence(unittest.TestCase):
    """Test the convert_variants_to_haplotype_sequence function."""

    def test_single_snp(self):
        """Test applying a single SNP to a reference sequence."""
        # Reference: ACGTACGT at position 100
        # Position 100=A, 101=C, 102=G, 103=T, 104=A, 105=C, 106=G, 107=T
        # SNP at position 103: T -> A
        # Expected: ACGAACGT
        result = convert_variants_to_haplotype_sequence(
            pos_1based=100,
            reference_sequence="ACGTACGT",
            variants=[(103, "T", "A")]
        )
        self.assertEqual(result, "ACGAACGT")

    def test_single_insertion(self):
        """Test applying a single insertion to a reference sequence."""
        # Reference: ACGTACGT at position 100
        # Position 100=A, 101=C, 102=G, 103=T, 104=A, 105=C, 106=G, 107=T
        # Insertion at position 102: G -> GAA (insert AA after G)
        # Expected: AC + GAA + TACGT = ACGAATACGT
        result = convert_variants_to_haplotype_sequence(
            pos_1based=100,
            reference_sequence="ACGTACGT",
            variants=[(102, "G", "GAA")]
        )
        self.assertEqual(result, "ACGAATACGT")

    def test_single_deletion(self):
        """Test applying a single deletion to a reference sequence."""
        # Reference: ACGTACGT at position 100
        # Position 100=A, 101=C, 102=G, 103=T, 104=A, 105=C, 106=G, 107=T
        # Deletion at position 102: GT -> G (delete T)
        # Expected: AC + G + ACGT = ACGACGT
        result = convert_variants_to_haplotype_sequence(
            pos_1based=100,
            reference_sequence="ACGTACGT",
            variants=[(102, "GT", "G")]
        )
        self.assertEqual(result, "ACGACGT")

    def test_multiple_variants(self):
        """Test applying multiple variants in order."""
        # Reference: ACGTACGT at position 100
        # Position 100=A, 101=C, 102=G, 103=T, 104=A, 105=C, 106=G, 107=T
        # SNP at position 101: C -> T
        # Insertion at position 105: C -> CAA
        # Expected: A + T + GTA + CAA + GT = ATGTACAAGT
        result = convert_variants_to_haplotype_sequence(
            pos_1based=100,
            reference_sequence="ACGTACGT",
            variants=[(101, "C", "T"), (105, "C", "CAA")]
        )
        self.assertEqual(result, "ATGTACAAGT")

    def test_shared_suffix_trimming(self):
        """Test that shared suffixes are properly trimmed."""
        # Reference: CAGCAGCAG at position 100
        # Variant at position 100: CAGCAG -> CAG (6bp deletion represented with shared suffix)
        # After suffix trimming: CAG -> empty, so we're deleting CAG
        # Expected: CAG + CAG = CAGCAG (one repeat deleted)
        result = convert_variants_to_haplotype_sequence(
            pos_1based=100,
            reference_sequence="CAGCAGCAG",
            variants=[(100, "CAGCAG", "CAG")]
        )
        # After suffix trimming: CAGCAG -> CAG means we remove the suffix CAG from both
        # ref becomes "CAG", alt becomes "" (empty)
        # So we replace first 3bp (CAG) with nothing, leaving CAGCAG
        self.assertEqual(result, "CAGCAG")

    def test_empty_variant_list(self):
        """Test with no variants - should return reference unchanged."""
        result = convert_variants_to_haplotype_sequence(
            pos_1based=100,
            reference_sequence="ACGTACGT",
            variants=[]
        )
        self.assertEqual(result, "ACGTACGT")

    def test_variant_at_start(self):
        """Test variant at the very start of the sequence."""
        result = convert_variants_to_haplotype_sequence(
            pos_1based=100,
            reference_sequence="ACGTACGT",
            variants=[(100, "A", "T")]
        )
        self.assertEqual(result, "TCGTACGT")

    def test_variant_at_end(self):
        """Test variant at the very end of the sequence."""
        result = convert_variants_to_haplotype_sequence(
            pos_1based=100,
            reference_sequence="ACGTACGT",
            variants=[(107, "T", "A")]
        )
        self.assertEqual(result, "ACGTACGA")

    def test_ref_mismatch_raises_error(self):
        """Test that mismatched reference raises ValueError."""
        # Reference: ACGTACGT at position 100
        # Position 103 has T, not G - so this should raise an error
        with self.assertRaises(ValueError) as context:
            convert_variants_to_haplotype_sequence(
                pos_1based=100,
                reference_sequence="ACGTACGT",
                variants=[(103, "G", "A")]  # Position 103 has T, not G
            )
        self.assertIn("does not match", str(context.exception))

    def test_out_of_order_variants_raises_error(self):
        """Test that out-of-order variants raise ValueError."""
        with self.assertRaises(ValueError) as context:
            convert_variants_to_haplotype_sequence(
                pos_1based=100,
                reference_sequence="ACGTACGT",
                variants=[(105, "C", "T"), (103, "G", "A")]  # Out of order
            )
        self.assertIn("before or overlaps", str(context.exception))

    def test_variant_beyond_sequence_raises_error(self):
        """Test that variant beyond reference sequence raises ValueError."""
        with self.assertRaises(ValueError) as context:
            convert_variants_to_haplotype_sequence(
                pos_1based=100,
                reference_sequence="ACGT",  # Only 4bp, ends at position 103
                variants=[(110, "A", "T")]  # Beyond the end
            )
        self.assertIn("beyond the end", str(context.exception))

    def test_invalid_ref_allele_raises_error(self):
        """Test that invalid DNA bases in ref raise ValueError."""
        with self.assertRaises(ValueError) as context:
            convert_variants_to_haplotype_sequence(
                pos_1based=100,
                reference_sequence="ACGTACGT",
                variants=[(100, "X", "A")]  # X is not a valid DNA base
            )
        self.assertIn("Invalid ref allele", str(context.exception))

    def test_invalid_alt_allele_raises_error(self):
        """Test that invalid DNA bases in alt raise ValueError."""
        with self.assertRaises(ValueError) as context:
            convert_variants_to_haplotype_sequence(
                pos_1based=100,
                reference_sequence="ACGTACGT",
                variants=[(100, "A", "Z")]  # Z is not a valid DNA base
            )
        self.assertIn("Invalid alt allele", str(context.exception))

    def test_case_insensitivity(self):
        """Test that lowercase sequences are handled correctly."""
        # Position 103 has 't' (lowercase), variant has 't' -> 'a'
        result = convert_variants_to_haplotype_sequence(
            pos_1based=100,
            reference_sequence="acgtacgt",
            variants=[(103, "t", "a")]
        )
        self.assertEqual(result, "ACGAACGT")

    def test_complex_example_from_docstring(self):
        """Test the example from the docstring."""
        # From docstring: reference at pos 100 is "ACAGCAG", variants are:
        # (103, "G", "A"), (106, "G", "AT")
        # Expected: "ACAACAAT"
        result = convert_variants_to_haplotype_sequence(
            pos_1based=100,
            reference_sequence="ACAGCAG",
            variants=[(103, "G", "A"), (106, "G", "AT")]
        )
        self.assertEqual(result, "ACAACAAT")


class TestExtractHaplotypeSequencesFromVcf(unittest.TestCase):
    """Test the extract_haplotype_sequences_from_vcf function."""

    def setUp(self):
        """Set up test fixtures with a mock fasta object."""
        # Default: return a CAG repeat sequence chr1:100-120 = "CAGCAGCAGCAGCAGCAGCA"
        # This is 20bp spanning 6.67 CAG repeats
        self._reference_seq = "CAGCAGCAGCAGCAGCAGCA"

    def _create_mock_fasta(self, ref_seq=None):
        """Create a mock pyfaidx Fasta object.

        Args:
            ref_seq (str): Reference sequence to return. If None, uses default.

        Returns:
            mock.MagicMock: A mock fasta object
        """
        if ref_seq is None:
            ref_seq = self._reference_seq

        class MockChrom:
            def __init__(self, seq):
                self.seq = seq

            def __getitem__(self, s):
                return self.seq[s.start:s.stop]

        class MockFasta:
            def __init__(self, seq):
                self.chrom = MockChrom(seq)
                self.raise_keyerror = False

            def __getitem__(self, chrom):
                if self.raise_keyerror:
                    raise KeyError(f"Chromosome {chrom} not found")
                return self.chrom

        return MockFasta(ref_seq)

    def _create_mock_variant(self, pos, ref, alt, gt, phased=True):
        """Helper to create a mock pysam variant record.

        Args:
            pos (int): 1-based position
            ref (str): Reference allele
            alt (str): Alternate allele
            gt (tuple): Genotype tuple, e.g., (0, 1) for 0/1
            phased (bool): Whether the genotype is phased

        Returns:
            mock.MagicMock: A mock variant record
        """
        variant = mock.MagicMock()
        variant.pos = pos
        variant.ref = ref
        variant.alleles = (ref, alt)

        # Mock samples[0] for single-sample VCF
        sample = mock.MagicMock()
        sample.get = mock.MagicMock(return_value=gt)
        sample.phased = phased
        variant.samples = [sample]

        return variant

    def test_no_variants_returns_reference(self):
        """Test that no variants returns reference sequence for both haplotypes."""
        mock_fasta = self._create_mock_fasta()

        result = extract_haplotype_sequences_from_vcf(
            chrom="chr1",
            start_0based=0,
            end=9,  # First 9 bp: "CAGCAGCAG"
            fasta_obj=mock_fasta,
            vcf_variants=[]
        )

        self.assertEqual(result[0], "CAGCAGCAG")
        self.assertEqual(result[1], "CAGCAGCAG")

    def test_single_het_variant_phased(self):
        """Test a single heterozygous phased variant (0|1)."""
        # Reference: CAGCAGCAGCAGCAGCAGCA (20bp starting at position 0)
        # Position (1-based): 1 2 3 4 5 6 7 8 9...
        # Character:          C A G C A G C A G...
        # Variant at position 4 (1-based): C -> T on haplotype 1
        # Haplotype 0: CAGCAGCAG...  (reference)
        # Haplotype 1: CAGTAG... (with SNP at position 4)
        mock_fasta = self._create_mock_fasta()

        variant = self._create_mock_variant(
            pos=4,  # 1-based position 4 (0-based index 3) = "C"
            ref="C",
            alt="T",
            gt=(0, 1),  # Het: ref on hap0, alt on hap1
            phased=True
        )

        result = extract_haplotype_sequences_from_vcf(
            chrom="chr1",
            start_0based=0,
            end=12,  # First 12 bp
            fasta_obj=mock_fasta,
            vcf_variants=[variant]
        )

        # Haplotype 0 should be reference: CAGCAGCAGCAG
        self.assertEqual(result[0], "CAGCAGCAGCAG")
        # Haplotype 1 should have SNP at position 4 (idx 3): CAGTAGCAGCAG
        self.assertEqual(result[1], "CAGTAGCAGCAG")

    def test_single_hom_alt_variant(self):
        """Test a homozygous alternate variant (1|1)."""
        mock_fasta = self._create_mock_fasta()

        variant = self._create_mock_variant(
            pos=4,  # 1-based position 4 (idx 3) = "C"
            ref="C",
            alt="T",
            gt=(1, 1),  # Hom alt
            phased=True
        )

        result = extract_haplotype_sequences_from_vcf(
            chrom="chr1",
            start_0based=0,
            end=12,
            fasta_obj=mock_fasta,
            vcf_variants=[variant]
        )

        # Both haplotypes should have the variant: CAGTAGCAGCAG
        self.assertEqual(result[0], "CAGTAGCAGCAG")
        self.assertEqual(result[1], "CAGTAGCAGCAG")

    def test_single_hom_ref_variant(self):
        """Test a homozygous reference variant (0|0) returns reference."""
        mock_fasta = self._create_mock_fasta()

        variant = self._create_mock_variant(
            pos=4,  # 1-based position 4 (idx 3) = "C"
            ref="C",
            alt="T",
            gt=(0, 0),  # Hom ref
            phased=True
        )

        result = extract_haplotype_sequences_from_vcf(
            chrom="chr1",
            start_0based=0,
            end=12,
            fasta_obj=mock_fasta,
            vcf_variants=[variant]
        )

        # Both should be reference
        self.assertEqual(result[0], "CAGCAGCAGCAG")
        self.assertEqual(result[1], "CAGCAGCAGCAG")

    def test_insertion_variant(self):
        """Test an insertion variant."""
        mock_fasta = self._create_mock_fasta()

        # Insert CAG at position 4 (after "CAGC")
        variant = self._create_mock_variant(
            pos=4,  # 1-based position
            ref="C",
            alt="CCAG",  # Insert CAG after C
            gt=(0, 1),
            phased=True
        )

        result = extract_haplotype_sequences_from_vcf(
            chrom="chr1",
            start_0based=0,
            end=12,
            fasta_obj=mock_fasta,
            vcf_variants=[variant]
        )

        # Haplotype 0: reference
        self.assertEqual(result[0], "CAGCAGCAGCAG")
        # Haplotype 1: CAGCCAGAGCAGCAG (insert CAG after position 3)
        # Reference[0:12] = "CAGCAGCAGCAG"
        # Position 4 (1-based) = index 3 = "C"
        # After variant: "CAG" + "CCAG" + "AGCAGCAG" = "CAGCCAGAGCAGCAG"
        self.assertEqual(result[1], "CAGCCAGAGCAGCAG")

    def test_deletion_variant(self):
        """Test a deletion variant."""
        mock_fasta = self._create_mock_fasta()

        # Delete CAG at positions 4-7 (1-based 4-6)
        variant = self._create_mock_variant(
            pos=4,  # 1-based
            ref="CAGC",  # Delete AGC
            alt="C",
            gt=(0, 1),
            phased=True
        )

        result = extract_haplotype_sequences_from_vcf(
            chrom="chr1",
            start_0based=0,
            end=12,
            fasta_obj=mock_fasta,
            vcf_variants=[variant]
        )

        # Haplotype 0: reference "CAGCAGCAGCAG"
        self.assertEqual(result[0], "CAGCAGCAGCAG")
        # Haplotype 1: "CAG" + "C" + "AGCAGCAG" but wait...
        # Reference: C A G C A G C A G C A G
        # Index:     0 1 2 3 4 5 6 7 8 9 10 11
        # Position:  1 2 3 4 5 6 7 8 9 10 11 12
        # Variant at pos 4 (idx 3): CAGC -> C, deletes AGC
        # Result: C A G + C + A G C A G = "CAGCAGCAG" (9 bp)
        self.assertEqual(result[1], "CAGCAGCAG")

    def test_multiple_unphased_variants_returns_missing(self):
        """Test that multiple unphased variants returns missing genotype (None, None).

        Per Design Decision #3: If multiple variants overlap AND any GT is unphased,
        return (None, None) to indicate missing genotype.
        """
        mock_fasta = self._create_mock_fasta()

        variant1 = self._create_mock_variant(
            pos=4, ref="C", alt="T", gt=(0, 1), phased=False  # Unphased!
        )
        variant2 = self._create_mock_variant(
            pos=7, ref="A", alt="G", gt=(1, 0), phased=True
        )

        result = extract_haplotype_sequences_from_vcf(
            chrom="chr1",
            start_0based=0,
            end=12,
            fasta_obj=mock_fasta,
            vcf_variants=[variant1, variant2]
        )

        # Should return (None, None) due to unphased multiple variants
        self.assertIsNone(result[0])
        self.assertIsNone(result[1])

    def test_single_unphased_variant_works(self):
        """Test that a single unphased variant still works (no phasing ambiguity)."""
        mock_fasta = self._create_mock_fasta()

        # Single variant - phasing doesn't matter
        # Position 4 (1-based, idx 3) = "C"
        variant = self._create_mock_variant(
            pos=4, ref="C", alt="T", gt=(0, 1), phased=False
        )

        result = extract_haplotype_sequences_from_vcf(
            chrom="chr1",
            start_0based=0,
            end=12,
            fasta_obj=mock_fasta,
            vcf_variants=[variant]
        )

        # Should work - single variant doesn't need phasing
        self.assertIsNotNone(result[0])
        self.assertIsNotNone(result[1])

    def test_missing_gt_returns_none_for_haplotype(self):
        """Test that missing genotype (.) returns None for that haplotype."""
        mock_fasta = self._create_mock_fasta()

        # Variant with missing genotype on haplotype 1
        # Position 4 (1-based, idx 3) = "C"
        variant = self._create_mock_variant(
            pos=4, ref="C", alt="T", gt=(0, None), phased=True
        )

        result = extract_haplotype_sequences_from_vcf(
            chrom="chr1",
            start_0based=0,
            end=12,
            fasta_obj=mock_fasta,
            vcf_variants=[variant]
        )

        # Haplotype 0 should be reference (gt=0)
        self.assertEqual(result[0], "CAGCAGCAGCAG")
        # Haplotype 1 should be None (missing genotype)
        self.assertIsNone(result[1])

    def test_star_allele_skipped(self):
        """Test that star alleles (*) are properly skipped."""
        mock_fasta = self._create_mock_fasta()

        # Create a variant with star allele
        # Position 4 (1-based, idx 3) = "C"
        variant = mock.MagicMock()
        variant.pos = 4
        variant.ref = "C"
        variant.alleles = ("C", "*")  # Star allele

        sample = mock.MagicMock()
        sample.get = mock.MagicMock(return_value=(0, 1))
        sample.phased = True
        variant.samples = [sample]

        result = extract_haplotype_sequences_from_vcf(
            chrom="chr1",
            start_0based=0,
            end=12,
            fasta_obj=mock_fasta,
            vcf_variants=[variant]
        )

        # Both should be reference since star allele is skipped
        self.assertEqual(result[0], "CAGCAGCAGCAG")
        self.assertEqual(result[1], "CAGCAGCAGCAG")

    def test_multiple_phased_variants(self):
        """Test multiple phased variants on same haplotype."""
        mock_fasta = self._create_mock_fasta()

        # Two variants on haplotype 1
        # Reference: C A G C A G C A G C A G
        # Position:  1 2 3 4 5 6 7 8 9 10 11 12
        # Index:     0 1 2 3 4 5 6 7 8 9 10 11
        # Position 2 (idx 1) = "A", Position 5 (idx 4) = "A"
        variant1 = self._create_mock_variant(
            pos=2, ref="A", alt="T", gt=(0, 1), phased=True
        )
        variant2 = self._create_mock_variant(
            pos=5, ref="A", alt="G", gt=(0, 1), phased=True
        )

        result = extract_haplotype_sequences_from_vcf(
            chrom="chr1",
            start_0based=0,
            end=12,
            fasta_obj=mock_fasta,
            vcf_variants=[variant1, variant2]
        )

        # Haplotype 0: reference
        self.assertEqual(result[0], "CAGCAGCAGCAG")
        # Haplotype 1: two SNPs at positions 2 and 5
        # Position 2 (idx 1): A -> T, Position 5 (idx 4): A -> G
        # Reference: C A G C A G C A G C A G
        # With SNPs: C T G C G G C A G C A G
        self.assertEqual(result[1], "CTGCGGCAGCAG")

    def test_variant_spanning_beyond_locus_is_trimmed(self):
        """Test that variants spanning beyond locus boundaries are properly handled.

        Per Design Decision #4: Fetch reference with padding to cover all variants,
        build full haplotype, trim to locus boundaries.
        """
        # Create a longer reference sequence for this test
        # Use DNA bases: AT prefix, CAG repeats, then GC suffix
        extended_ref = "ATCAGCAGCAGCAGCAGCAGCAGC"
        mock_fasta = self._create_mock_fasta(extended_ref)

        # Variant at position 1 (1-based) that spans before our locus start at index 2
        # Locus is index 2-14 (positions 3-15, 1-based)
        variant = self._create_mock_variant(
            pos=1,  # Variant starts before locus
            ref="AT",  # Spans 2 bases (indices 0-1)
            alt="TTT",  # Replace with 3 bases
            gt=(0, 1),
            phased=True
        )

        result = extract_haplotype_sequences_from_vcf(
            chrom="chr1",
            start_0based=2,  # Locus starts at index 2
            end=14,  # Locus ends at index 14
            fasta_obj=mock_fasta,
            vcf_variants=[variant]
        )

        # Haplotype 0: reference from index 2-14 = "CAGCAGCAGCAG"
        self.assertEqual(result[0], "CAGCAGCAGCAG")

        # Haplotype 1: After applying variant "AT" -> "TTT" at position 1,
        # the full sequence becomes "TTTCAGCAGCAGCAGCAGCAGCAGC"
        # Original locus was indices 2-14 (12 bp)
        # After the +1bp insertion at the start, the content that was at index 2
        # is now at index 3 in the modified sequence.
        # But we trim back to the original locus boundaries (indices 2-14 in output space)
        # The trimming logic adjusts: output_offset_at_locus_start += length_change (+1)
        # So we extract from index 3 to 15 in the modified sequence
        # Modified sequence: T T T C A G C A G C A G C A G ...
        # Indices:           0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
        # Extracting [3:15] = "TCAGCAGCAGCA" (wait, that's not right either)
        # Actually index 3 is 'C', so [3:15] = "CAGCAGCAGCAG"
        # Hmm, let me reconsider...
        # Full modified: "TTTCAGCAGCAGCAGCAGCAGCAGC" (25 chars)
        # We want 12bp starting where the original index 2 content now lives
        # Original index 2 had 'C' (first char of CAG repeat)
        # After +1bp insertion, that 'C' is now at modified index 3
        # So we want modified[3:15] = "CAGCAGCAGCAG" (unchanged!)
        self.assertEqual(result[1], "CAGCAGCAGCAG")

    def test_chromosome_not_found_returns_none(self):
        """Test that chromosome not in fasta returns (None, None)."""
        mock_fasta = self._create_mock_fasta()
        mock_fasta.raise_keyerror = True

        result = extract_haplotype_sequences_from_vcf(
            chrom="chrUNKNOWN",
            start_0based=0,
            end=12,
            fasta_obj=mock_fasta,
            vcf_variants=[]
        )

        self.assertIsNone(result[0])
        self.assertIsNone(result[1])


class TestGenotypeSingleLocus(unittest.TestCase):
    """Test the genotype_single_locus function."""

    def setUp(self):
        """Set up the test case."""
        # Default: return a CAG repeat sequence chr1:0-12 = "CAGCAGCAGCAG"
        # This is 12bp spanning 4 CAG repeats
        self._reference_seq = "CAGCAGCAGCAG"

    def _create_mock_fasta(self, ref_seq=None):
        """Create a mock pyfaidx Fasta object."""
        if ref_seq is None:
            ref_seq = self._reference_seq

        class MockChrom:
            def __init__(self, seq):
                self.seq = seq

            def __getitem__(self, s):
                return self.seq[s.start:s.stop]

        class MockFasta:
            def __init__(self, seq):
                self.chrom = MockChrom(seq)
                self.raise_keyerror = False

            def __getitem__(self, chrom):
                if self.raise_keyerror:
                    raise KeyError(f"Chromosome {chrom} not found")
                return self.chrom

        return MockFasta(ref_seq)

    def _create_mock_variant(self, pos, ref, alt, gt, phased=True):
        """Helper to create a mock pysam variant record."""
        variant = mock.MagicMock()
        variant.pos = pos
        variant.ref = ref
        variant.alleles = (ref, alt)

        sample = mock.MagicMock()
        sample.get = mock.MagicMock(return_value=gt)
        sample.phased = phased
        variant.samples = [sample]

        return variant

    def _create_mock_vcf_file(self, variants):
        """Helper to create a mock VCF file that returns given variants on fetch."""
        vcf_file = mock.MagicMock()
        vcf_file.fetch = mock.MagicMock(return_value=variants)
        return vcf_file

    def test_no_overlapping_variants_returns_hom_ref(self):
        """Test that no overlapping variants returns HOM with reference sequence."""
        tr_locus = ReferenceTandemRepeat(
            chrom="chr1",
            start_0based=0,
            end_1based=12,
            repeat_unit="CAG"
        )
        mock_fasta = self._create_mock_fasta()
        mock_vcf = self._create_mock_vcf_file([])

        result = genotype_single_locus(tr_locus, mock_vcf, mock_fasta)

        self.assertEqual(result.zygosity, "HOM")
        self.assertEqual(result.num_repeats_allele1, 4)
        self.assertEqual(result.num_repeats_allele2, 4)
        self.assertEqual(result.allele1_sequence, "CAGCAGCAGCAG")
        self.assertEqual(result.allele2_sequence, "CAGCAGCAGCAG")
        self.assertEqual(result.num_overlapping_variants, 0)

    def test_het_expansion_returns_het(self):
        """Test heterozygous expansion (one allele larger than reference)."""
        # Reference: CAGCAGCAGCAG (4 repeats)
        # Variant at pos 1: CAG -> CAGCAG (insert CAG)
        # GT: 0|1 means haplotype 0 = ref, haplotype 1 = alt
        variant = self._create_mock_variant(
            pos=1, ref="C", alt="CCAG", gt=(0, 1), phased=True
        )

        tr_locus = ReferenceTandemRepeat(
            chrom="chr1",
            start_0based=0,
            end_1based=12,
            repeat_unit="CAG"
        )
        mock_fasta = self._create_mock_fasta()
        mock_vcf = self._create_mock_vcf_file([variant])

        result = genotype_single_locus(tr_locus, mock_vcf, mock_fasta)

        self.assertEqual(result.zygosity, "HET")
        # Allele 1 (haplotype 0): reference = 4 repeats
        self.assertEqual(result.num_repeats_allele1, 4)
        # Allele 2 (haplotype 1): expansion = 5 repeats (15bp // 3bp)
        self.assertEqual(result.num_repeats_allele2, 5)
        self.assertEqual(result.num_repeats_short_allele, 4)
        self.assertEqual(result.num_repeats_long_allele, 5)
        self.assertEqual(result.num_overlapping_variants, 1)

    def test_hom_alt_returns_hom(self):
        """Test homozygous alternate (both alleles same non-ref)."""
        # Reference: CAGCAGCAGCAG (4 repeats)
        # Variant: insert CAG on both alleles
        variant = self._create_mock_variant(
            pos=1, ref="C", alt="CCAG", gt=(1, 1), phased=True
        )

        tr_locus = ReferenceTandemRepeat(
            chrom="chr1",
            start_0based=0,
            end_1based=12,
            repeat_unit="CAG"
        )
        mock_fasta = self._create_mock_fasta()
        mock_vcf = self._create_mock_vcf_file([variant])

        result = genotype_single_locus(tr_locus, mock_vcf, mock_fasta)

        self.assertEqual(result.zygosity, "HOM")
        self.assertEqual(result.num_repeats_allele1, 5)
        self.assertEqual(result.num_repeats_allele2, 5)

    def test_multiallelic_variant_returns_missing(self):
        """Test that multi-allelic variants cause missing genotype."""
        # Multi-allelic variant has >2 alleles
        variant = mock.MagicMock()
        variant.pos = 1
        variant.ref = "C"
        variant.alleles = ("C", "T", "G")  # 3 alleles = multi-allelic

        sample = mock.MagicMock()
        sample.get = mock.MagicMock(return_value=(0, 1))
        sample.phased = True
        variant.samples = [sample]

        tr_locus = ReferenceTandemRepeat(
            chrom="chr1",
            start_0based=0,
            end_1based=12,
            repeat_unit="CAG"
        )
        mock_fasta = self._create_mock_fasta()
        mock_vcf = self._create_mock_vcf_file([variant])

        result = genotype_single_locus(tr_locus, mock_vcf, mock_fasta)

        # Should be missing genotype due to multi-allelic
        self.assertIsNone(result.zygosity)
        self.assertIsNone(result.num_repeats_allele1)
        self.assertIsNone(result.num_repeats_allele2)
        self.assertIsNone(result.allele1_sequence)
        self.assertIsNone(result.allele2_sequence)
        # Variants should still be recorded
        self.assertEqual(result.num_overlapping_variants, 1)

    def test_multiple_unphased_variants_returns_missing(self):
        """Test that multiple unphased variants cause missing genotype."""
        # Two variants with unphased genotypes
        variant1 = self._create_mock_variant(
            pos=1, ref="C", alt="T", gt=(0, 1), phased=False
        )
        variant2 = self._create_mock_variant(
            pos=4, ref="C", alt="T", gt=(0, 1), phased=False
        )

        tr_locus = ReferenceTandemRepeat(
            chrom="chr1",
            start_0based=0,
            end_1based=12,
            repeat_unit="CAG"
        )
        mock_fasta = self._create_mock_fasta()
        mock_vcf = self._create_mock_vcf_file([variant1, variant2])

        result = genotype_single_locus(tr_locus, mock_vcf, mock_fasta)

        # Should be missing genotype due to unphased multiple variants
        self.assertIsNone(result.zygosity)
        self.assertIsNone(result.allele1_sequence)
        self.assertIsNone(result.allele2_sequence)

    def test_single_unphased_variant_works(self):
        """Test that a single unphased variant still works."""
        # Single variant with unphased genotype should work (no phasing ambiguity)
        variant = self._create_mock_variant(
            pos=1, ref="C", alt="CCAG", gt=(0, 1), phased=False
        )

        tr_locus = ReferenceTandemRepeat(
            chrom="chr1",
            start_0based=0,
            end_1based=12,
            repeat_unit="CAG"
        )
        mock_fasta = self._create_mock_fasta()
        mock_vcf = self._create_mock_vcf_file([variant])

        result = genotype_single_locus(tr_locus, mock_vcf, mock_fasta)

        # Should work (single variant has no phasing ambiguity)
        self.assertEqual(result.zygosity, "HET")
        self.assertIsNotNone(result.allele1_sequence)
        self.assertIsNotNone(result.allele2_sequence)

    def test_hemi_with_missing_haplotype(self):
        """Test HEMI zygosity when one haplotype is missing."""
        # Create variant with one missing haplotype (./1 genotype represented as (None, 1))
        variant = self._create_mock_variant(
            pos=1, ref="C", alt="CCAG", gt=(None, 1), phased=True
        )

        tr_locus = ReferenceTandemRepeat(
            chrom="chr1",
            start_0based=0,
            end_1based=12,
            repeat_unit="CAG"
        )
        mock_fasta = self._create_mock_fasta()
        mock_vcf = self._create_mock_vcf_file([variant])

        result = genotype_single_locus(tr_locus, mock_vcf, mock_fasta)

        # Should be HEMI (one allele missing)
        self.assertEqual(result.zygosity, "HEMI")
        self.assertIsNone(result.allele1_sequence)
        self.assertIsNotNone(result.allele2_sequence)

    def test_locus_properties_preserved(self):
        """Test that locus properties are correctly preserved in result."""
        tr_locus = ReferenceTandemRepeat(
            chrom="chr1",
            start_0based=0,
            end_1based=12,
            repeat_unit="CAG"
        )
        mock_fasta = self._create_mock_fasta()
        mock_vcf = self._create_mock_vcf_file([])

        result = genotype_single_locus(tr_locus, mock_vcf, mock_fasta)

        self.assertEqual(result.chrom, "chr1")
        self.assertEqual(result.start_0based, 0)
        self.assertEqual(result.end, 12)
        self.assertEqual(result.motif, "CAG")
        self.assertEqual(result.motif_size, 3)
        self.assertEqual(result.locus, "chr1:0-12")
        self.assertEqual(result.locus_id, "chr1-0-12-CAG")
        self.assertEqual(result.num_repeats_in_reference, 4)


class TestGenotypingPipeline(unittest.TestCase):
    """Integration tests for the genotyping pipeline.

    These tests create actual temporary VCF and BED files to test the full
    genotyping workflow including file I/O, variant parsing, and genotype computation.
    """

    def setUp(self):
        """Set up test fixtures."""
        self._temp_files = []

    def tearDown(self):
        """Clean up temporary files."""
        for f in self._temp_files:
            if os.path.exists(f):
                os.unlink(f)
            # Clean up associated index files
            for ext in [".tbi", ".csi", ".fai"]:
                idx_file = f + ext
                if os.path.exists(idx_file):
                    os.unlink(idx_file)

    def _create_temp_file(self, content, suffix, compress=False):
        """Create a temporary file with given content.

        Args:
            content (str): File content
            suffix (str): File suffix (e.g., '.vcf', '.bed')
            compress (bool): If True, compress with gzip

        Returns:
            str: Path to the temporary file
        """
        import gzip
        if compress:
            suffix = suffix + ".gz"

        with tempfile.NamedTemporaryFile(suffix=suffix, delete=False, mode='wb' if compress else 'w') as f:
            if compress:
                with gzip.open(f.name, 'wt') as gz:
                    gz.write(content)
            else:
                f.write(content)
            self._temp_files.append(f.name)
            return f.name

    def _create_test_vcf_and_index(self, vcf_content):
        """Create a bgzipped and indexed VCF file.

        Args:
            vcf_content (str): VCF content string

        Returns:
            str: Path to the bgzipped VCF file, or None if bgzip/tabix unavailable
        """
        import subprocess
        vcf_path = self._create_temp_file(vcf_content, ".vcf")
        vcf_gz_path = vcf_path + ".gz"
        self._temp_files.append(vcf_gz_path)
        self._temp_files.append(vcf_gz_path + ".tbi")

        try:
            with open(vcf_gz_path, 'wb') as gz_out:
                subprocess.run(["bgzip", "-c", vcf_path], stdout=gz_out, check=True)
            subprocess.run(["tabix", "-p", "vcf", vcf_gz_path], check=True)
            return vcf_gz_path
        except (FileNotFoundError, subprocess.CalledProcessError):
            return None

    def _create_test_fasta(self, seq_dict):
        """Create a temporary FASTA file and index.

        Args:
            seq_dict (dict): Dictionary mapping chromosome names to sequences

        Returns:
            str: Path to the FASTA file, or None if pyfaidx fails
        """
        fasta_content = ""
        for chrom, seq in seq_dict.items():
            fasta_content += f">{chrom}\n{seq}\n"

        fasta_path = self._create_temp_file(fasta_content, ".fa")

        try:
            import pyfaidx
            # Create index
            pyfaidx.Fasta(fasta_path)
            self._temp_files.append(fasta_path + ".fai")
            return fasta_path
        except (ImportError, IOError, OSError):
            # ImportError - pyfaidx not installed
            # IOError/OSError - file system errors during indexing
            return None

    def test_genotype_simple_expansion_hom(self):
        """Test genotyping a homozygous expansion.

        Create a locus with a CAG repeat where both alleles have an expansion
        compared to the reference.
        """
        import pysam
        import pyfaidx

        # Reference: 4 CAG repeats = CAGCAGCAGCAG
        fasta_path = self._create_test_fasta({"chr1": "AAACAGCAGCAGCAGAAA"})
        if fasta_path is None:
            self.skipTest("pyfaidx unavailable")

        # VCF with homozygous expansion: insert CAG at position 4 (1-based)
        # This adds one CAG repeat on both alleles (1|1)
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=18>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t4\t.\tC\tCCAG\t.\tPASS\t.\tGT\t1|1
"""
        vcf_gz_path = self._create_test_vcf_and_index(vcf_content)
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        # Create TR locus: positions 4-15 (0-based: 3-15) = 12bp = 4 CAG repeats
        tr_locus = ReferenceTandemRepeat(
            chrom="chr1",
            start_0based=3,
            end_1based=15,  # 12bp
            repeat_unit="CAG"
        )

        fasta_obj = pyfaidx.Fasta(fasta_path, one_based_attributes=False, as_raw=True)
        vcf_file = pysam.VariantFile(vcf_gz_path)

        try:
            result = genotype_single_locus(tr_locus, vcf_file, fasta_obj, vcf_chrom_convention="chr")

            # Both alleles should have 5 repeats (reference 4 + insertion 1)
            self.assertEqual(result.zygosity, "HOM")
            self.assertEqual(result.num_repeats_allele1, 5)
            self.assertEqual(result.num_repeats_allele2, 5)
            self.assertEqual(result.num_overlapping_variants, 1)
        finally:
            vcf_file.close()
            fasta_obj.close()

    def test_genotype_simple_contraction_het(self):
        """Test genotyping a heterozygous contraction.

        Create a locus where one allele has a contraction (fewer repeats).
        """
        import pysam
        import pyfaidx

        # Reference: 4 CAG repeats = CAGCAGCAGCAG
        fasta_path = self._create_test_fasta({"chr1": "AAACAGCAGCAGCAGAAA"})
        if fasta_path is None:
            self.skipTest("pyfaidx unavailable")

        # VCF with het deletion: delete one CAG on haplotype 1 (0|1)
        # Position 4 (1-based), CAGC -> C (delete AGC which is part of repeat)
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=18>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t4\t.\tCAGC\t.\tC\t.\tPASS\t.\tGT\t0|1
"""
        # Note: This VCF has a deletion represented as CAGC -> C
        # Actually let's use a simpler representation
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=18>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t4\t.\tCAGC\tC\t.\tPASS\t.\tGT\t0|1
"""
        vcf_gz_path = self._create_test_vcf_and_index(vcf_content)
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        tr_locus = ReferenceTandemRepeat(
            chrom="chr1",
            start_0based=3,
            end_1based=15,
            repeat_unit="CAG"
        )

        fasta_obj = pyfaidx.Fasta(fasta_path, one_based_attributes=False, as_raw=True)
        vcf_file = pysam.VariantFile(vcf_gz_path)

        try:
            result = genotype_single_locus(tr_locus, vcf_file, fasta_obj, vcf_chrom_convention="chr")

            # Haplotype 0: reference = 4 repeats
            # Haplotype 1: deletion = 3 repeats (9bp)
            self.assertEqual(result.zygosity, "HET")
            self.assertEqual(result.num_repeats_allele1, 4)
            self.assertEqual(result.num_repeats_allele2, 3)
            self.assertEqual(result.num_repeats_short_allele, 3)
            self.assertEqual(result.num_repeats_long_allele, 4)
        finally:
            vcf_file.close()
            fasta_obj.close()

    def test_genotype_no_overlapping_variants(self):
        """Test genotyping when no VCF variants overlap the locus.

        Should return reference genotype (HOM with reference repeat count).
        """
        import pysam
        import pyfaidx

        fasta_path = self._create_test_fasta({"chr1": "AAACAGCAGCAGCAGAAA"})
        if fasta_path is None:
            self.skipTest("pyfaidx unavailable")

        # VCF with variant at a different position (doesn't overlap locus)
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=18>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t1\t.\tA\tT\t.\tPASS\t.\tGT\t0|1
"""
        vcf_gz_path = self._create_test_vcf_and_index(vcf_content)
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        tr_locus = ReferenceTandemRepeat(
            chrom="chr1",
            start_0based=3,
            end_1based=15,
            repeat_unit="CAG"
        )

        fasta_obj = pyfaidx.Fasta(fasta_path, one_based_attributes=False, as_raw=True)
        vcf_file = pysam.VariantFile(vcf_gz_path)

        try:
            result = genotype_single_locus(tr_locus, vcf_file, fasta_obj, vcf_chrom_convention="chr")

            self.assertEqual(result.zygosity, "HOM")
            self.assertEqual(result.num_repeats_allele1, 4)
            self.assertEqual(result.num_repeats_allele2, 4)
            self.assertEqual(result.num_overlapping_variants, 0)
            self.assertEqual(result.allele1_sequence, "CAGCAGCAGCAG")
            self.assertEqual(result.allele2_sequence, "CAGCAGCAGCAG")
        finally:
            vcf_file.close()
            fasta_obj.close()

    def test_genotype_complex_locus_multiple_variants(self):
        """Test genotyping a locus with multiple overlapping phased variants."""
        import pysam
        import pyfaidx

        # Reference with 5 CAG repeats
        fasta_path = self._create_test_fasta({"chr1": "AAACAGCAGCAGCAGCAGAAA"})
        if fasta_path is None:
            self.skipTest("pyfaidx unavailable")

        # Two phased variants on different haplotypes
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=21>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t4\t.\tC\tCCAG\t.\tPASS\t.\tGT\t1|0
chr1\t10\t.\tC\tCCAG\t.\tPASS\t.\tGT\t0|1
"""
        vcf_gz_path = self._create_test_vcf_and_index(vcf_content)
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        # Locus: positions 4-18 (0-based: 3-18) = 15bp = 5 CAG repeats
        tr_locus = ReferenceTandemRepeat(
            chrom="chr1",
            start_0based=3,
            end_1based=18,
            repeat_unit="CAG"
        )

        fasta_obj = pyfaidx.Fasta(fasta_path, one_based_attributes=False, as_raw=True)
        vcf_file = pysam.VariantFile(vcf_gz_path)

        try:
            result = genotype_single_locus(tr_locus, vcf_file, fasta_obj, vcf_chrom_convention="chr")

            # Both haplotypes get one insertion, so both have 6 repeats
            self.assertEqual(result.zygosity, "HOM")
            self.assertEqual(result.num_repeats_allele1, 6)
            self.assertEqual(result.num_repeats_allele2, 6)
            self.assertEqual(result.num_overlapping_variants, 2)
        finally:
            vcf_file.close()
            fasta_obj.close()

    def test_genotype_missing_genotype_in_vcf(self):
        """Test handling of missing genotype (./.) in VCF."""
        import pysam
        import pyfaidx

        fasta_path = self._create_test_fasta({"chr1": "AAACAGCAGCAGCAGAAA"})
        if fasta_path is None:
            self.skipTest("pyfaidx unavailable")

        # VCF with missing genotype (./.)
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=18>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t4\t.\tC\tCCAG\t.\tPASS\t.\tGT\t./.
"""
        vcf_gz_path = self._create_test_vcf_and_index(vcf_content)
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        tr_locus = ReferenceTandemRepeat(
            chrom="chr1",
            start_0based=3,
            end_1based=15,
            repeat_unit="CAG"
        )

        fasta_obj = pyfaidx.Fasta(fasta_path, one_based_attributes=False, as_raw=True)
        vcf_file = pysam.VariantFile(vcf_gz_path)

        try:
            result = genotype_single_locus(tr_locus, vcf_file, fasta_obj, vcf_chrom_convention="chr")

            # Both haplotypes should be None (missing)
            self.assertIsNone(result.allele1_sequence)
            self.assertIsNone(result.allele2_sequence)
            self.assertIsNone(result.zygosity)
        finally:
            vcf_file.close()
            fasta_obj.close()

    def test_genotype_unphased_multiple_variants_returns_missing(self):
        """Test that multiple unphased variants cause missing genotype.

        Per Design Decision #3: If multiple variants overlap AND any GT is unphased,
        return missing genotype.
        """
        import pysam
        import pyfaidx

        fasta_path = self._create_test_fasta({"chr1": "AAACAGCAGCAGCAGAAA"})
        if fasta_path is None:
            self.skipTest("pyfaidx unavailable")

        # Two variants with unphased genotypes (using /)
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=18>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t4\t.\tC\tT\t.\tPASS\t.\tGT\t0/1
chr1\t7\t.\tC\tT\t.\tPASS\t.\tGT\t0/1
"""
        vcf_gz_path = self._create_test_vcf_and_index(vcf_content)
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        tr_locus = ReferenceTandemRepeat(
            chrom="chr1",
            start_0based=3,
            end_1based=15,
            repeat_unit="CAG"
        )

        fasta_obj = pyfaidx.Fasta(fasta_path, one_based_attributes=False, as_raw=True)
        vcf_file = pysam.VariantFile(vcf_gz_path)

        try:
            result = genotype_single_locus(tr_locus, vcf_file, fasta_obj, vcf_chrom_convention="chr")

            # Should be missing due to unphased multiple variants
            self.assertIsNone(result.zygosity)
            self.assertIsNone(result.allele1_sequence)
            self.assertIsNone(result.allele2_sequence)
            # But variants should still be recorded
            self.assertEqual(result.num_overlapping_variants, 2)
        finally:
            vcf_file.close()
            fasta_obj.close()

    def test_genotype_multiallelic_variant_returns_missing(self):
        """Test that multi-allelic variants (>2 alleles) cause missing genotype.

        Per Design Decision #15: Mark entire locus as missing if any overlapping
        variant has >2 alleles.
        """
        import pysam
        import pyfaidx

        fasta_path = self._create_test_fasta({"chr1": "AAACAGCAGCAGCAGAAA"})
        if fasta_path is None:
            self.skipTest("pyfaidx unavailable")

        # Multi-allelic variant: REF with two ALT alleles
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=18>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t4\t.\tC\tCCAG,CCAGCAG\t.\tPASS\t.\tGT\t1/2
"""
        vcf_gz_path = self._create_test_vcf_and_index(vcf_content)
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        tr_locus = ReferenceTandemRepeat(
            chrom="chr1",
            start_0based=3,
            end_1based=15,
            repeat_unit="CAG"
        )

        fasta_obj = pyfaidx.Fasta(fasta_path, one_based_attributes=False, as_raw=True)
        vcf_file = pysam.VariantFile(vcf_gz_path)

        try:
            result = genotype_single_locus(tr_locus, vcf_file, fasta_obj, vcf_chrom_convention="chr")

            # Should be missing due to multi-allelic variant
            self.assertIsNone(result.zygosity)
            self.assertIsNone(result.num_repeats_allele1)
            self.assertIsNone(result.num_repeats_allele2)
            # Variant should still be recorded
            self.assertEqual(result.num_overlapping_variants, 1)
        finally:
            vcf_file.close()
            fasta_obj.close()

    def test_genotype_variant_spanning_beyond_locus(self):
        """Test that variants spanning beyond locus boundaries are handled correctly.

        Per Design Decision #4: Fetch reference with padding to cover all variants,
        build full haplotype, trim to locus boundaries.
        """
        import pysam
        import pyfaidx

        # Extended reference to allow variant before locus
        fasta_path = self._create_test_fasta({"chr1": "GGGGCAGCAGCAGCAGTTTT"})
        if fasta_path is None:
            self.skipTest("pyfaidx unavailable")

        # Variant at position 3 (1-based) that spans before our locus (which starts at 0-based 4)
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=20>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t3\t.\tGGC\tGC\t.\tPASS\t.\tGT\t0|1
"""
        vcf_gz_path = self._create_test_vcf_and_index(vcf_content)
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        # Locus: positions 5-16 (0-based: 4-16) = 12bp = 4 CAG repeats
        tr_locus = ReferenceTandemRepeat(
            chrom="chr1",
            start_0based=4,
            end_1based=16,
            repeat_unit="CAG"
        )

        fasta_obj = pyfaidx.Fasta(fasta_path, one_based_attributes=False, as_raw=True)
        vcf_file = pysam.VariantFile(vcf_gz_path)

        try:
            result = genotype_single_locus(tr_locus, vcf_file, fasta_obj, vcf_chrom_convention="chr")

            # Haplotype 0 should be reference
            self.assertEqual(result.allele1_sequence, "CAGCAGCAGCAG")
            # Haplotype 1: The variant is before the locus, so after trimming
            # the locus content should still be the same
            self.assertEqual(result.zygosity, "HOM")
        finally:
            vcf_file.close()
            fasta_obj.close()

    def test_chromosome_naming_normalization_no_chr(self):
        """Test that chromosome names are normalized correctly (chr1 vs 1).

        Per Design Decision #14: Auto-detect convention from each file and
        normalize internally.
        """
        import pysam
        import pyfaidx

        # VCF uses "1" without chr prefix
        fasta_path = self._create_test_fasta({"1": "AAACAGCAGCAGCAGAAA"})
        if fasta_path is None:
            self.skipTest("pyfaidx unavailable")

        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=1,length=18>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
1\t4\t.\tC\tCCAG\t.\tPASS\t.\tGT\t0|1
"""
        vcf_gz_path = self._create_test_vcf_and_index(vcf_content)
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        # Locus uses "chr1" with prefix - but should work via normalization
        tr_locus = ReferenceTandemRepeat(
            chrom="chr1",  # Uses chr prefix
            start_0based=3,
            end_1based=15,
            repeat_unit="CAG"
        )

        fasta_obj = pyfaidx.Fasta(fasta_path, one_based_attributes=False, as_raw=True)
        vcf_file = pysam.VariantFile(vcf_gz_path)

        try:
            result = genotype_single_locus(tr_locus, vcf_file, fasta_obj, vcf_chrom_convention="no_chr")

            # Should successfully genotype despite naming mismatch
            self.assertEqual(result.zygosity, "HET")
            self.assertEqual(result.num_repeats_allele1, 4)
            self.assertEqual(result.num_repeats_allele2, 5)
        finally:
            vcf_file.close()
            fasta_obj.close()

    def test_genotype_partial_haplotype_returns_hemi(self):
        """Test HEMI zygosity when one haplotype is missing (e.g., chrX in males)."""
        import pysam
        import pyfaidx

        fasta_path = self._create_test_fasta({"chr1": "AAACAGCAGCAGCAGAAA"})
        if fasta_path is None:
            self.skipTest("pyfaidx unavailable")

        # VCF with partial genotype (one allele missing): .|1
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=18>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t4\t.\tC\tCCAG\t.\tPASS\t.\tGT\t.|1
"""
        vcf_gz_path = self._create_test_vcf_and_index(vcf_content)
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        tr_locus = ReferenceTandemRepeat(
            chrom="chr1",
            start_0based=3,
            end_1based=15,
            repeat_unit="CAG"
        )

        fasta_obj = pyfaidx.Fasta(fasta_path, one_based_attributes=False, as_raw=True)
        vcf_file = pysam.VariantFile(vcf_gz_path)

        try:
            result = genotype_single_locus(tr_locus, vcf_file, fasta_obj, vcf_chrom_convention="chr")

            # Should be HEMI with one allele present
            self.assertEqual(result.zygosity, "HEMI")
            # Allele 1 (haplotype 0) should be None
            self.assertIsNone(result.allele1_sequence)
            # Allele 2 (haplotype 1) should have the expansion
            self.assertIsNotNone(result.allele2_sequence)
            self.assertEqual(result.num_repeats_allele2, 5)
        finally:
            vcf_file.close()
            fasta_obj.close()

    def test_genotype_pure_vs_impure_repeats(self):
        """Test that repeat purity is correctly computed for genotyped alleles."""
        import pysam
        import pyfaidx

        fasta_path = self._create_test_fasta({"chr1": "AAACAGCAGCAGCAGAAA"})
        if fasta_path is None:
            self.skipTest("pyfaidx unavailable")

        # Variant introduces an interruption: C -> T at position 7 (within repeat)
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=18>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t7\t.\tC\tT\t.\tPASS\t.\tGT\t0|1
"""
        vcf_gz_path = self._create_test_vcf_and_index(vcf_content)
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        tr_locus = ReferenceTandemRepeat(
            chrom="chr1",
            start_0based=3,
            end_1based=15,
            repeat_unit="CAG"
        )

        fasta_obj = pyfaidx.Fasta(fasta_path, one_based_attributes=False, as_raw=True)
        vcf_file = pysam.VariantFile(vcf_gz_path)

        try:
            result = genotype_single_locus(tr_locus, vcf_file, fasta_obj, vcf_chrom_convention="chr")

            # Allele 1 (reference) should be pure
            self.assertIsNotNone(result.allele1_purity)
            self.assertGreater(result.allele1_purity, 0.99)  # Pure repeat

            # Allele 2 (with SNP) should be impure
            self.assertIsNotNone(result.allele2_purity)
            self.assertLess(result.allele2_purity, 1.0)  # Impure due to SNP

            # Overall is_pure_repeat should be False (one allele is impure)
            self.assertFalse(result.is_pure_repeat)
        finally:
            vcf_file.close()
            fasta_obj.close()

    def test_genotype_to_tsv_dict_method(self):
        """Test that the to_tsv_dict method produces correct output."""
        import pysam
        import pyfaidx

        fasta_path = self._create_test_fasta({"chr1": "AAACAGCAGCAGCAGAAA"})
        if fasta_path is None:
            self.skipTest("pyfaidx unavailable")

        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=18>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t4\t.\tC\tCCAG\t.\tPASS\t.\tGT\t0|1
"""
        vcf_gz_path = self._create_test_vcf_and_index(vcf_content)
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        tr_locus = ReferenceTandemRepeat(
            chrom="chr1",
            start_0based=3,
            end_1based=15,
            repeat_unit="CAG"
        )

        fasta_obj = pyfaidx.Fasta(fasta_path, one_based_attributes=False, as_raw=True)
        vcf_file = pysam.VariantFile(vcf_gz_path)

        try:
            result = genotype_single_locus(tr_locus, vcf_file, fasta_obj, vcf_chrom_convention="chr")
            tsv_dict = result.to_tsv_dict()

            # Check required columns are present
            self.assertEqual(tsv_dict["Chrom"], "chr1")
            self.assertEqual(tsv_dict["Start0Based"], 3)
            self.assertEqual(tsv_dict["End"], 15)
            self.assertEqual(tsv_dict["Motif"], "CAG")
            self.assertEqual(tsv_dict["MotifSize"], 3)
            self.assertEqual(tsv_dict["Zygosity"], "HET")
            self.assertEqual(tsv_dict["NumRepeatsShortAllele"], 4)
            self.assertEqual(tsv_dict["NumRepeatsLongAllele"], 5)
            self.assertEqual(tsv_dict["NumOverlappingVariants"], 1)
        finally:
            vcf_file.close()
            fasta_obj.close()


class TestRepeatCounting(unittest.TestCase):
    """Test the compute_repeat_counts_from_sequence function."""

    def test_count_repeats_pure_str(self):
        """Test counting repeats in a pure CAG repeat sequence."""
        # Pure CAG repeat: CAGCAGCAGCAG = 4 complete repeats
        result = compute_repeat_counts_from_sequence("CAGCAGCAGCAG", "CAG")

        self.assertEqual(result["num_repeats"], 4)
        self.assertEqual(result["repeat_size_bp"], 12)
        self.assertAlmostEqual(result["purity"], 1.0, places=2)
        self.assertTrue(result["is_pure"])

    def test_count_repeats_impure_str(self):
        """Test counting repeats in an STR with interruptions."""
        # CAG repeat with one interruption: CAGCAACAGCAG
        # Has 4 motif-sized chunks but one is CAA instead of CAG
        result = compute_repeat_counts_from_sequence("CAGCAACAGCAG", "CAG")

        self.assertEqual(result["num_repeats"], 4)
        self.assertEqual(result["repeat_size_bp"], 12)
        # Purity should be < 1.0 due to interruption (1 base mismatch out of 12)
        self.assertLess(result["purity"], 1.0)
        self.assertGreater(result["purity"], 0.8)  # Still mostly pure
        self.assertFalse(result["is_pure"])

    def test_count_repeats_vntr(self):
        """Test counting repeats with a longer VNTR motif."""
        # AAGGG repeat (5bp motif): AAGGGAAGGGAAGGG = 3 complete repeats
        result = compute_repeat_counts_from_sequence("AAGGGAAGGGAAGGG", "AAGGG")

        self.assertEqual(result["num_repeats"], 3)
        self.assertEqual(result["repeat_size_bp"], 15)
        self.assertAlmostEqual(result["purity"], 1.0, places=2)
        self.assertTrue(result["is_pure"])

    def test_count_repeats_partial_repeat(self):
        """Test counting when sequence is not divisible by motif length."""
        # 10bp sequence with 3bp motif = 3 full repeats (10 // 3 = 3)
        # CAGCAGCAGC = 3 full CAG + partial C
        result = compute_repeat_counts_from_sequence("CAGCAGCAGC", "CAG")

        self.assertEqual(result["num_repeats"], 3)  # Integer division: 10 // 3 = 3
        self.assertEqual(result["repeat_size_bp"], 10)
        # Purity should still be reasonable since most bases match
        self.assertIsNotNone(result["purity"])
        self.assertGreater(result["purity"], 0.9)

    def test_count_repeats_empty_sequence(self):
        """Test counting repeats with an empty sequence."""
        result = compute_repeat_counts_from_sequence("", "CAG")

        self.assertEqual(result["num_repeats"], 0)
        self.assertEqual(result["repeat_size_bp"], 0)
        # Purity is undefined for empty sequence
        self.assertIsNone(result["purity"])
        self.assertIsNone(result["is_pure"])

    def test_count_repeats_none_sequence(self):
        """Test handling of None sequence (missing genotype)."""
        result = compute_repeat_counts_from_sequence(None, "CAG")

        self.assertIsNone(result["num_repeats"])
        self.assertIsNone(result["repeat_size_bp"])
        self.assertIsNone(result["purity"])
        self.assertIsNone(result["is_pure"])

    def test_count_repeats_single_motif(self):
        """Test counting when sequence is exactly one motif."""
        result = compute_repeat_counts_from_sequence("CAG", "CAG")

        self.assertEqual(result["num_repeats"], 1)
        self.assertEqual(result["repeat_size_bp"], 3)
        self.assertAlmostEqual(result["purity"], 1.0, places=2)
        self.assertTrue(result["is_pure"])

    def test_count_repeats_sequence_shorter_than_motif(self):
        """Test when sequence is shorter than the motif."""
        # Sequence "CA" is shorter than motif "CAG"
        result = compute_repeat_counts_from_sequence("CA", "CAG")

        self.assertEqual(result["num_repeats"], 0)  # 2 // 3 = 0
        self.assertEqual(result["repeat_size_bp"], 2)
        # Purity may be None or undefined when sequence < motif
        # The compute_repeat_purity function may return NaN

    def test_count_repeats_homopolymer(self):
        """Test counting repeats in a homopolymer (1bp motif)."""
        # Poly-A: AAAAAAAA = 8 A's
        result = compute_repeat_counts_from_sequence("AAAAAAAA", "A")

        self.assertEqual(result["num_repeats"], 8)
        self.assertEqual(result["repeat_size_bp"], 8)
        self.assertAlmostEqual(result["purity"], 1.0, places=2)
        self.assertTrue(result["is_pure"])

    def test_count_repeats_dinucleotide(self):
        """Test counting repeats in a dinucleotide repeat."""
        # AT repeat: ATATAT = 3 complete AT repeats
        result = compute_repeat_counts_from_sequence("ATATAT", "AT")

        self.assertEqual(result["num_repeats"], 3)
        self.assertEqual(result["repeat_size_bp"], 6)
        self.assertAlmostEqual(result["purity"], 1.0, places=2)
        self.assertTrue(result["is_pure"])

    def test_count_repeats_severe_interruption(self):
        """Test counting repeats with severe interruptions."""
        # Multiple interruptions: CAGTAGCAGAAT
        # Expected pattern: CAGCAGCAGCAG
        # Mismatches at positions: 3 (T instead of C), 10 (A instead of C), 11 (T instead of G)
        result = compute_repeat_counts_from_sequence("CAGTAGCAGAAT", "CAG")

        self.assertEqual(result["num_repeats"], 4)  # 12 // 3 = 4
        self.assertEqual(result["repeat_size_bp"], 12)
        # Purity should be notably lower with multiple mismatches
        self.assertLess(result["purity"], 0.9)
        self.assertFalse(result["is_pure"])

    def test_count_repeats_case_insensitivity(self):
        """Test that repeat counting handles lowercase sequences."""
        result = compute_repeat_counts_from_sequence("cagcagcagcag", "CAG")

        self.assertEqual(result["num_repeats"], 4)
        self.assertEqual(result["repeat_size_bp"], 12)
        # Purity calculation should handle case differences
        self.assertIsNotNone(result["purity"])


class TestEndToEndGenotyping(unittest.TestCase):
    """End-to-end tests for the genotype subcommand using realistic test data.

    These tests create a small set of realistic TR loci modeled after known pathogenic
    repeat loci, along with corresponding VCF variants and reference sequences, to test
    the full genotyping workflow.

    Test data design:
    - Uses simplified versions of known TR loci (CAG repeats similar to Huntingtin, etc.)
    - Tests various genotype scenarios: expansions, contractions, reference genotypes
    - Verifies TSV and JSON output format and content

    Note: These tests require bgzip and tabix to be available in the PATH.
    """

    def setUp(self):
        """Set up test fixtures with realistic TR loci and variants."""
        self._temp_files = []
        self._temp_dir = None

    def tearDown(self):
        """Clean up temporary files and directories."""
        for f in self._temp_files:
            if os.path.exists(f):
                os.unlink(f)
            # Clean up associated index files
            for ext in [".tbi", ".csi", ".fai"]:
                idx_file = f + ext
                if os.path.exists(idx_file):
                    os.unlink(idx_file)
        if self._temp_dir and os.path.exists(self._temp_dir):
            shutil.rmtree(self._temp_dir)

    def _create_temp_file(self, content, suffix, compress=False):
        """Create a temporary file with given content."""
        import gzip
        if compress:
            suffix = suffix + ".gz"

        with tempfile.NamedTemporaryFile(suffix=suffix, delete=False, mode='wb' if compress else 'w') as f:
            if compress:
                with gzip.open(f.name, 'wt') as gz:
                    gz.write(content)
            else:
                f.write(content)
            self._temp_files.append(f.name)
            return f.name

    def _create_test_vcf_and_index(self, vcf_content):
        """Create a bgzipped and indexed VCF file."""
        import subprocess
        vcf_path = self._create_temp_file(vcf_content, ".vcf")
        vcf_gz_path = vcf_path + ".gz"
        self._temp_files.append(vcf_gz_path)
        self._temp_files.append(vcf_gz_path + ".tbi")

        try:
            with open(vcf_gz_path, 'wb') as gz_out:
                subprocess.run(["bgzip", "-c", vcf_path], stdout=gz_out, check=True)
            subprocess.run(["tabix", "-p", "vcf", vcf_gz_path], check=True)
            return vcf_gz_path
        except (FileNotFoundError, subprocess.CalledProcessError):
            return None

    def _create_test_fasta(self, seq_dict):
        """Create a temporary FASTA file and index."""
        fasta_content = ""
        for chrom, seq in seq_dict.items():
            fasta_content += f">{chrom}\n{seq}\n"

        fasta_path = self._create_temp_file(fasta_content, ".fa")

        try:
            import pyfaidx
            pyfaidx.Fasta(fasta_path)
            self._temp_files.append(fasta_path + ".fai")
            return fasta_path
        except (ImportError, IOError, OSError):
            # ImportError - pyfaidx not installed
            # IOError/OSError - file system errors during indexing
            return None

    def _create_test_bed_and_index(self, bed_content):
        """Create a bgzipped and indexed BED file."""
        import subprocess
        bed_path = self._create_temp_file(bed_content, ".bed")
        bed_gz_path = bed_path + ".gz"
        self._temp_files.append(bed_gz_path)
        self._temp_files.append(bed_gz_path + ".tbi")

        try:
            with open(bed_gz_path, 'wb') as gz_out:
                subprocess.run(["bgzip", "-c", bed_path], stdout=gz_out, check=True)
            subprocess.run(["tabix", "-p", "bed", bed_gz_path], check=True)
            return bed_gz_path
        except (FileNotFoundError, subprocess.CalledProcessError):
            return None

    def test_end_to_end_genotype_multiple_loci(self):
        """End-to-end test with multiple TR loci simulating realistic data.

        Test data description:
        - chr1:100-124: CAG repeat (8 repeats in reference) - heterozygous expansion
        - chr1:200-218: AT repeat (9 repeats in reference) - homozygous reference
        - chr1:300-327: CTG repeat (9 repeats in reference) - heterozygous contraction
        - chr1:400-430: AAGGG repeat (6 repeats in reference) - homozygous expansion
        - chr1:500-518: GCC repeat (6 repeats in reference) - no overlapping variants
        - chr1:600-636: CAG repeat (12 repeats in reference) - complex multi-variant

        This tests:
        - Various motif sizes (2bp, 3bp, 5bp)
        - Various genotype scenarios (HOM, HET, expansions, contractions)
        - Correct repeat count computation
        - TSV output format and column values
        """
        import pyfaidx
        import gzip

        # Create reference with embedded tandem repeats
        # Positions are 0-based for BED coordinates
        reference_seq = (
            # Positions 0-99: flanking region before first locus
            "A" * 100 +
            # Position 100-123: CAG repeat (8 CAG = 24bp)
            "CAG" * 8 +
            # Positions 124-199: flanking between loci
            "T" * 76 +
            # Position 200-217: AT repeat (9 AT = 18bp)
            "AT" * 9 +
            # Positions 218-299: flanking
            "G" * 82 +
            # Position 300-326: CTG repeat (9 CTG = 27bp)
            "CTG" * 9 +
            # Positions 327-399: flanking
            "C" * 73 +
            # Position 400-429: AAGGG repeat (6 AAGGG = 30bp)
            "AAGGG" * 6 +
            # Positions 430-499: flanking
            "A" * 70 +
            # Position 500-517: GCC repeat (6 GCC = 18bp)
            "GCC" * 6 +
            # Positions 518-599: flanking
            "T" * 82 +
            # Position 600-635: CAG repeat (12 CAG = 36bp)
            "CAG" * 12 +
            # Positions 636-699: trailing flanking
            "G" * 64
        )

        fasta_path = self._create_test_fasta({"chr1": reference_seq})
        if fasta_path is None:
            self.skipTest("pyfaidx unavailable")

        # Create catalog BED file with 6 TR loci
        # Format: chrom start end name (name contains motif)
        bed_content = """chr1\t100\t124\tCAG
chr1\t200\t218\tAT
chr1\t300\t327\tCTG
chr1\t400\t430\tAAGGG
chr1\t500\t518\tGCC
chr1\t600\t636\tCAG
"""
        bed_gz_path = self._create_test_bed_and_index(bed_content)
        if bed_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        # Create VCF with variants at specific loci
        # Note: VCF uses 1-based positions. REF alleles must match the reference exactly.
        #
        # Variant 1: Heterozygous expansion at locus 1 (CAG at 100-123)
        #   - Position 101 (0-based 100) is 'C', insert CAGCAG after C -> C to CCAGCAG (adds 2 CAG)
        # Variant 2: No variant at locus 2 (AT) - should get reference genotype
        # Variant 3: Heterozygous contraction at locus 3 (CTG at 300-326)
        #   - Position 301 (0-based 300) is 'CTGC', delete CTG -> CTGC to C (removes 1 CTG)
        # Variant 4: Homozygous expansion at locus 4 (AAGGG at 400-429)
        #   - Position 401 (0-based 400) is 'A', insert AGGGA after A -> A to AAGGGA (adds 1 AAGGG)
        # Variant 5: No variant at locus 5 (GCC) - should get reference genotype
        # Variant 6: Two variants at locus 6 (CAG at 600-635) - one on each haplotype
        #   - Position 601 (0-based 600) is 'C', insert AGCAGCAG -> C to CAGCAGCAG (adds 3 CAG on hap0)
        #   - Position 619 (0-based 618) is 'C', insert AG -> C to CAG (adds 1 CAG on hap1)
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=700>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t101\t.\tC\tCCAGCAG\t.\tPASS\t.\tGT\t0|1
chr1\t301\t.\tCTGC\tC\t.\tPASS\t.\tGT\t1|0
chr1\t401\t.\tA\tAAGGGA\t.\tPASS\t.\tGT\t1|1
chr1\t601\t.\tC\tCAGCAGCAG\t.\tPASS\t.\tGT\t1|0
chr1\t619\t.\tC\tCAG\t.\tPASS\t.\tGT\t0|1
"""
        vcf_gz_path = self._create_test_vcf_and_index(vcf_content)
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        # Create temporary output directory
        self._temp_dir = tempfile.mkdtemp()
        output_prefix = os.path.join(self._temp_dir, "test_output")

        # Create args object
        args = argparse.Namespace(
            reference_fasta_path=fasta_path,
            catalog_bed=bed_gz_path,
            input_vcf_path=vcf_gz_path,
            input_vcf_prefix="test",
            output_prefix=output_prefix,
            interval=None,
            verbose=False,
            show_progress_bar=False,
            write_vcf=False,
            write_json=False,
            add_motif_counts=None,
            trf_executable_path=None,
        )

        # Import and run the subcommand
        from str_analysis.filter_vcf_to_tandem_repeats import do_genotype_subcommand

        # Run the genotype subcommand
        do_genotype_subcommand(args)

        # Verify output TSV was created
        tsv_path = f"{output_prefix}.tandem_repeat_genotypes.tsv.gz"
        self._temp_files.append(tsv_path)
        self.assertTrue(os.path.exists(tsv_path), f"Expected output TSV not found: {tsv_path}")

        # Read and parse the output TSV
        with gzip.open(tsv_path, 'rt') as f:
            lines = f.readlines()

        self.assertGreater(len(lines), 1, "TSV should have header + data rows")

        # Parse header
        header = lines[0].strip().split('\t')
        expected_columns = [
            "Chrom", "Start0Based", "End", "Locus", "LocusId", "Motif", "CanonicalMotif",
            "MotifSize", "NumRepeatsInReference", "NumRepeatsShortAllele", "NumRepeatsLongAllele",
            "RepeatSizeShortAlleleBp", "RepeatSizeLongAlleleBp", "Zygosity", "IsPureRepeat",
            "RepeatPurity", "Allele1Sequence", "Allele2Sequence", "NumOverlappingVariants",
            "VariantPositions"
        ]
        for col in expected_columns:
            self.assertIn(col, header, f"Missing column: {col}")

        # Parse data rows
        data_rows = []
        for line in lines[1:]:
            values = line.strip().split('\t')
            row_dict = dict(zip(header, values))
            data_rows.append(row_dict)

        self.assertEqual(len(data_rows), 6, "Should have 6 genotyped loci")

        # Verify specific loci
        # Sort by Start0Based to ensure consistent order
        data_rows.sort(key=lambda x: int(x["Start0Based"]))

        # Locus 1 (CAG at 100-124): Heterozygous expansion
        locus1 = data_rows[0]
        self.assertEqual(locus1["Chrom"], "chr1")
        self.assertEqual(int(locus1["Start0Based"]), 100)
        self.assertEqual(int(locus1["End"]), 124)
        self.assertEqual(locus1["Motif"], "CAG")
        self.assertEqual(int(locus1["NumRepeatsInReference"]), 8)
        self.assertEqual(locus1["Zygosity"], "HET")
        # One allele is reference (8), one has expansion (+2 = 10)
        self.assertEqual(int(locus1["NumRepeatsShortAllele"]), 8)
        self.assertEqual(int(locus1["NumRepeatsLongAllele"]), 10)
        self.assertEqual(int(locus1["NumOverlappingVariants"]), 1)

        # Locus 2 (AT at 200-218): Reference genotype (no variants)
        locus2 = data_rows[1]
        self.assertEqual(int(locus2["Start0Based"]), 200)
        self.assertEqual(locus2["Motif"], "AT")
        self.assertEqual(locus2["Zygosity"], "HOM")
        self.assertEqual(int(locus2["NumRepeatsInReference"]), 9)
        self.assertEqual(int(locus2["NumRepeatsShortAllele"]), 9)
        self.assertEqual(int(locus2["NumRepeatsLongAllele"]), 9)
        self.assertEqual(int(locus2["NumOverlappingVariants"]), 0)

        # Locus 3 (CTG at 300-327): Heterozygous contraction
        locus3 = data_rows[2]
        self.assertEqual(int(locus3["Start0Based"]), 300)
        self.assertEqual(locus3["Motif"], "CTG")
        self.assertEqual(locus3["Zygosity"], "HET")
        self.assertEqual(int(locus3["NumRepeatsInReference"]), 9)
        # One allele has contraction (-1 = 8), one is reference (9)
        self.assertEqual(int(locus3["NumRepeatsShortAllele"]), 8)
        self.assertEqual(int(locus3["NumRepeatsLongAllele"]), 9)

        # Locus 4 (AAGGG at 400-430): Homozygous expansion
        locus4 = data_rows[3]
        self.assertEqual(int(locus4["Start0Based"]), 400)
        self.assertEqual(locus4["Motif"], "AAGGG")
        self.assertEqual(int(locus4["MotifSize"]), 5)
        self.assertEqual(locus4["Zygosity"], "HOM")
        self.assertEqual(int(locus4["NumRepeatsInReference"]), 6)
        # Both alleles have +1 expansion = 7
        self.assertEqual(int(locus4["NumRepeatsShortAllele"]), 7)
        self.assertEqual(int(locus4["NumRepeatsLongAllele"]), 7)

        # Locus 5 (GCC at 500-518): Reference genotype (no variants)
        locus5 = data_rows[4]
        self.assertEqual(int(locus5["Start0Based"]), 500)
        self.assertEqual(locus5["Motif"], "GCC")
        self.assertEqual(locus5["Zygosity"], "HOM")
        self.assertEqual(int(locus5["NumOverlappingVariants"]), 0)

        # Locus 6 (CAG at 600-636): Complex with multiple variants
        locus6 = data_rows[5]
        self.assertEqual(int(locus6["Start0Based"]), 600)
        self.assertEqual(locus6["Motif"], "CAG")
        self.assertEqual(int(locus6["NumRepeatsInReference"]), 12)
        self.assertEqual(int(locus6["NumOverlappingVariants"]), 2)
        # Haplotype 0 gets +3 CAG (from variant at 601): 12 + 3*3/3 = 12 + 8bp/3bp = ~14 repeats
        # Haplotype 1 gets +1 CAG (from variant at 619): 12 + 2bp/3bp = ~12 repeats
        # The repeat counting is len(seq) // len(motif), so:
        # - Original: 36bp / 3 = 12 repeats
        # - Haplotype 0: (36 + 8)bp / 3 = 14 repeats (44bp)
        # - Haplotype 1: (36 + 2)bp / 3 = 12 repeats (38bp, since 38//3 = 12)
        self.assertEqual(locus6["Zygosity"], "HET")
        self.assertIn(int(locus6["NumRepeatsShortAllele"]), [12, 13])
        self.assertIn(int(locus6["NumRepeatsLongAllele"]), [14, 15])

    def test_end_to_end_genotype_with_json_output(self):
        """Test end-to-end genotyping with JSON output enabled.

        Verifies that the JSON output contains all expected fields and
        matches the TSV output.
        """
        import pyfaidx
        import gzip
        import json

        # Simple reference with one CAG repeat locus
        reference_seq = "A" * 50 + "CAG" * 10 + "T" * 50  # 10 CAG repeats at position 50-79

        fasta_path = self._create_test_fasta({"chr1": reference_seq})
        if fasta_path is None:
            self.skipTest("pyfaidx unavailable")

        bed_content = "chr1\t50\t80\tCAG\n"
        bed_gz_path = self._create_test_bed_and_index(bed_content)
        if bed_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        # VCF with heterozygous expansion
        # Position 51 (0-based 50) is 'C', inserting AGCAG -> C to CAGCAG (adds 2 CAG)
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=130>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t51\t.\tC\tCCAGCAG\t.\tPASS\t.\tGT\t0|1
"""
        vcf_gz_path = self._create_test_vcf_and_index(vcf_content)
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        self._temp_dir = tempfile.mkdtemp()
        output_prefix = os.path.join(self._temp_dir, "test_json_output")

        args = argparse.Namespace(
            reference_fasta_path=fasta_path,
            catalog_bed=bed_gz_path,
            input_vcf_path=vcf_gz_path,
            input_vcf_prefix="test",
            output_prefix=output_prefix,
            interval=None,
            verbose=False,
            show_progress_bar=False,
            write_vcf=False,
            write_json=True,  # Enable JSON output
            add_motif_counts=None,
            trf_executable_path=None,
        )

        from str_analysis.filter_vcf_to_tandem_repeats import do_genotype_subcommand
        do_genotype_subcommand(args)

        # Verify JSON output was created
        json_path = f"{output_prefix}.tandem_repeat_genotypes.json.gz"
        self._temp_files.append(json_path)
        self.assertTrue(os.path.exists(json_path), f"Expected JSON output not found: {json_path}")

        # Parse JSON
        with gzip.open(json_path, 'rt') as f:
            json_data = json.load(f)

        self.assertIsInstance(json_data, list)
        self.assertEqual(len(json_data), 1)

        locus = json_data[0]
        self.assertEqual(locus["Chrom"], "chr1")
        self.assertEqual(locus["Start0Based"], 50)
        self.assertEqual(locus["End"], 80)
        self.assertEqual(locus["Motif"], "CAG")
        self.assertEqual(locus["Zygosity"], "HET")
        self.assertEqual(locus["NumRepeatsInReference"], 10)
        self.assertEqual(locus["NumRepeatsShortAllele"], 10)
        self.assertEqual(locus["NumRepeatsLongAllele"], 12)

        # Verify allele sequences are present
        self.assertIn("Allele1Sequence", locus)
        self.assertIn("Allele2Sequence", locus)
        # One allele should be reference (30bp = 10 CAG), one expanded (36bp = 12 CAG)
        allele_lengths = sorted([len(locus["Allele1Sequence"]), len(locus["Allele2Sequence"])])
        self.assertEqual(allele_lengths, [30, 36])

    def test_end_to_end_genotype_with_basic_motif_counts(self):
        """Test JSON output with basic motif counting enabled.

        Verifies that motif counts are correctly computed when
        --add-motif-counts basic is specified.
        """
        import pyfaidx
        import gzip
        import json

        # Reference with pure CAG repeat
        reference_seq = "A" * 30 + "CAG" * 6 + "T" * 30  # 6 CAG repeats

        fasta_path = self._create_test_fasta({"chr1": reference_seq})
        if fasta_path is None:
            self.skipTest("pyfaidx unavailable")

        bed_content = "chr1\t30\t48\tCAG\n"
        bed_gz_path = self._create_test_bed_and_index(bed_content)
        if bed_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        # VCF with no variants - both alleles are pure reference
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=78>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t1\t.\tA\tT\t.\tPASS\t.\tGT\t0|0
"""
        vcf_gz_path = self._create_test_vcf_and_index(vcf_content)
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        self._temp_dir = tempfile.mkdtemp()
        output_prefix = os.path.join(self._temp_dir, "test_motif_counts")

        args = argparse.Namespace(
            reference_fasta_path=fasta_path,
            catalog_bed=bed_gz_path,
            input_vcf_path=vcf_gz_path,
            input_vcf_prefix="test",
            output_prefix=output_prefix,
            interval=None,
            verbose=False,
            show_progress_bar=False,
            write_vcf=False,
            write_json=True,
            add_motif_counts="basic",  # Enable basic motif counting
            trf_executable_path=None,
        )

        from str_analysis.filter_vcf_to_tandem_repeats import do_genotype_subcommand
        do_genotype_subcommand(args)

        json_path = f"{output_prefix}.tandem_repeat_genotypes.json.gz"
        self._temp_files.append(json_path)

        with gzip.open(json_path, 'rt') as f:
            json_data = json.load(f)

        locus = json_data[0]

        # Verify motif counts are present
        self.assertIn("Allele1MotifCounts", locus)
        self.assertIn("Allele2MotifCounts", locus)

        # Both alleles should be pure CAG repeats (6 each)
        self.assertEqual(locus["Allele1MotifCounts"].get("CAG", 0), 6)
        self.assertEqual(locus["Allele2MotifCounts"].get("CAG", 0), 6)

    def test_end_to_end_genotype_missing_genotype_output(self):
        """Test that missing genotypes are correctly represented in output.

        When a locus has multiple unphased variants or multiallelic variants,
        the genotype should be missing and output as empty/null values.
        """
        import pyfaidx
        import gzip

        reference_seq = "A" * 30 + "CAG" * 6 + "T" * 30

        fasta_path = self._create_test_fasta({"chr1": reference_seq})
        if fasta_path is None:
            self.skipTest("pyfaidx unavailable")

        bed_content = "chr1\t30\t48\tCAG\n"
        bed_gz_path = self._create_test_bed_and_index(bed_content)
        if bed_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        # VCF with multiple unphased variants overlapping locus
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=78>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t31\t.\tA\tT\t.\tPASS\t.\tGT\t0/1
chr1\t34\t.\tC\tT\t.\tPASS\t.\tGT\t0/1
"""
        vcf_gz_path = self._create_test_vcf_and_index(vcf_content)
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        self._temp_dir = tempfile.mkdtemp()
        output_prefix = os.path.join(self._temp_dir, "test_missing")

        args = argparse.Namespace(
            reference_fasta_path=fasta_path,
            catalog_bed=bed_gz_path,
            input_vcf_path=vcf_gz_path,
            input_vcf_prefix="test",
            output_prefix=output_prefix,
            interval=None,
            verbose=False,
            show_progress_bar=False,
            write_vcf=False,
            write_json=False,
            add_motif_counts=None,
            trf_executable_path=None,
        )

        from str_analysis.filter_vcf_to_tandem_repeats import do_genotype_subcommand
        do_genotype_subcommand(args)

        tsv_path = f"{output_prefix}.tandem_repeat_genotypes.tsv.gz"
        self._temp_files.append(tsv_path)

        with gzip.open(tsv_path, 'rt') as f:
            lines = f.readlines()

        header = lines[0].strip().split('\t')
        values = lines[1].strip().split('\t')
        row = dict(zip(header, values))

        # Zygosity should be empty for missing genotype
        self.assertEqual(row["Zygosity"], "")
        # Allele sequences should be empty
        self.assertEqual(row["Allele1Sequence"], "")
        self.assertEqual(row["Allele2Sequence"], "")
        # But NumOverlappingVariants should still be counted
        self.assertEqual(int(row["NumOverlappingVariants"]), 2)


if __name__ == "__main__":
    unittest.main()
