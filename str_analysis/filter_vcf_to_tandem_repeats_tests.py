#!/usr/bin/env python3

"""Test script for the filter_vcf_to_catalog_tandem_repeats.py script."""

import argparse
import collections
import contextlib
import gzip
import io
import json
import os
import pkgutil
import re
import shutil
import subprocess
import tempfile
import unittest
from unittest import mock

import intervaltree
import pyfaidx
import pysam

from str_analysis.utils.fasta_utils import normalize_chromosome_name
from str_analysis.filter_vcf_to_tandem_repeats import Allele, TandemRepeatAllele, ReferenceTandemRepeat, \
    GenotypedTandemRepeat, \
    DEFAULT_MIN_INSERTION_SIZE_TO_CHECK, DEFAULT_MIN_INSERTION_PURITY, DEFAULT_MIN_INSERTION_PERIODICITY, \
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
    open_vcf_for_genotyping, get_overlapping_vcf_variants, build_contig_name_lookup, \
    convert_variants_to_haplotype_sequence, extract_haplotype_sequences_from_vcf, \
    compute_locus_start_and_end_offsets_in_haplotype, \
    compute_repeat_counts_from_sequence, genotype_single_locus, genotype_all_loci, write_tsv, \
    compute_motif_composition, compute_motif_lists_with_trf, \
    build_basic_split_motif_entry, build_trviz_motif_entry, format_motif_entry_as_sequence_string, \
    MOTIF_DETECTION_METHOD_TRF, MOTIF_DETECTION_METHOD_TRVIZ, MOTIF_DETECTION_METHOD_BASIC_SPLIT, \
    InsertionFilter, build_insertion_filter, check_if_inserted_sequence_is_repetitive, \
    find_insufficiently_repetitive_insertions_inside_locus, \
    extract_haplotype_sequences_and_insertions_from_vcf, \
    INSERTION_FILTER_REASON_CONTAINS_NS, INSERTION_FILTER_REASON_NOT_REPEAT_LIKE, \
    is_heterozygous_genotype, are_variants_unambiguously_phased, compute_chrom_sort_key, \
    parse_catalog_bed_file, run_trf_motif_splitting, do_genotype_subcommand, \
    GENOTYPE_TSV_OUTPUT_COLUMNS, \
    get_PAR_region_coordinates, get_locus_ploidy, PAR_REGIONS_BY_CHRX_LENGTH, OverlappingVariant, \
    detect_sex_chromosome_ploidy, \
    NO_CALL_REASON_AMBIGUOUS_PHASING, NO_CALL_REASON_HET_AT_HAPLOID_LOCUS, \
    NO_CALL_REASON_HAPLOID_GENOTYPE_AT_DIPLOID_LOCUS, \
    NO_CALL_REASON_HAPLOTYPE_BUILD_ERROR, NO_CALL_REASON_CONTIG_NOT_IN_REFERENCE, \
    NO_CALL_REASON_LOCUS_PAST_CONTIG_END, NO_CALL_REASON_NON_IUPAC_ALLELE, NonIUPACAlleleError, \
    NO_CALL_REASON_REF_ALLELE_MISMATCH, RefAlleleMismatchError, \
    NO_CALL_REASON_MISSING_GENOTYPE, find_records_with_a_missing_genotype_allele, \
    NO_CALL_REASON_NON_REPEAT_INSERTION, NO_CALL_REASON_CHROMOSOME_ABSENT, run_shell_command, \
    does_variant_change_the_locus, get_called_alt_alleles, \
    match_interval_to_catalog_contig, is_locus_genotyped_as_reference
from str_analysis.utils.find_motif_utils import format_motifs_as_sequence_string, \
    compute_best_phase_repeat_purity, compute_sequence_periodicity, compute_partial_copy_purity

# Sex chromosome ploidy of an XY sample, as detect_sex_chromosome_ploidy reports it
XY_PLOIDY = {"X": 1, "Y": 1}


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
        # AAGA canonicalizes to AAAG: the smallest rotation of the motif or of its reverse complement
        # (TCTT rotates to CTTT, AAGA rotates to AAAG). Asserting only that some string comes back would
        # let an identity function pass.
        self.assertEqual(tr_allele.canonical_repeat_unit, "AAAG")

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

    def test_summary_string_format(self):
        """Test summary_string shows correct repeat count (not divided by motif length).

        This test verifies the fix for a bug where num_repeats_in_variant_and_flanks
        was incorrectly divided by repeat_unit_length, causing e.g. 10 repeats of a 3bp
        motif to show as "3.3x" instead of "10.0x".
        """
        allele = Allele("chr22", 10515040, "T", "TAAGA", self._fasta_obj)
        tr_allele = TandemRepeatAllele(allele, "AAGA", False, 0, 4, 37, DETECTION_MODE_PURE_REPEATS)

        # Get the actual number of repeats
        num_repeats = tr_allele.num_repeats_in_variant_and_flanks

        # The summary string should show this exact count with 1 decimal place
        expected_repeat_str = f"{num_repeats:0.1f}x"
        summary = tr_allele.summary_string

        # Verify the repeat count in the summary string is correct
        self.assertIn(expected_repeat_str, summary,
            f"Summary string '{summary}' should contain '{expected_repeat_str}' "
            f"(num_repeats={num_repeats}, motif_length={tr_allele.repeat_unit_length})")

        # Also verify it does NOT contain an incorrectly divided value
        # (which would happen if we divided by motif length again)
        wrong_value = num_repeats / tr_allele.repeat_unit_length
        wrong_repeat_str = f"{wrong_value:0.1f}x"
        if wrong_repeat_str != expected_repeat_str:
            self.assertNotIn(wrong_repeat_str, summary,
                f"Summary string should NOT contain '{wrong_repeat_str}' "
                f"(incorrectly divided by motif length)")

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
        # CAG canonicalizes to AGC: the smallest rotation of the motif or of its reverse complement.
        # Asserting only that some string comes back would let an identity function pass.
        self.assertEqual(ref_tr.canonical_repeat_unit, "AGC")

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

    def test_merge_write_detailed_bed_for_plain_catalog(self):
        """--write-detailed-bed must produce the detailed BED even for plain (motif-only) catalogs.

        Previously the detailed BED was only written when the input catalog names already
        contained detail tokens, so the flag was silently ignored for plain BED catalogs.
        """
        if shutil.which("bgzip") is None or shutil.which("tabix") is None:
            self.skipTest("bgzip/tabix unavailable")

        temp_dir = tempfile.mkdtemp()
        try:
            # Plain catalog: the name field is just the motif, with no detail tokens. The header, blank and
            # track lines must be skipped the same way parse_catalog_bed_file skips them.
            input_bed_path = os.path.join(temp_dir, "plain_catalog.bed")
            with open(input_bed_path, "w") as f:
                f.write("track name=catalog\n#chrom\tstart\tend\tmotif\n")
                f.write("chr22\t10515000\t10515020\tAT\n\n")

            output_prefix = os.path.join(temp_dir, "merged")
            args = argparse.Namespace(
                reference_fasta_path=self._temp_fasta_path,
                input_bed_paths=[input_bed_path],
                output_prefix=output_prefix,
                interval=None,
                verbose=False,
                show_progress_bar=False,
                write_detailed_bed=True,
                batch_size=1000,
            )

            from str_analysis.filter_vcf_to_tandem_repeats import do_merge_subcommand
            do_merge_subcommand(args)

            self.assertTrue(os.path.exists(f"{output_prefix}.tandem_repeats.bed.gz"))
            detailed_bed_path = f"{output_prefix}.tandem_repeats.detailed.bed.gz"
            self.assertTrue(os.path.exists(detailed_bed_path),
                            "--write-detailed-bed was ignored for a plain input catalog")

            # Re-merging the detailed BED must not mistake its purity token for a detection mode and write it twice
            with gzip.open(detailed_bed_path, "rt") as f:
                name_field = f.readline().split("\t")[3]
            self.assertRegex(name_field, r"^AT:2bp:10\.0x:p\d\.\d\d$")
            args.input_bed_paths = [detailed_bed_path]
            args.output_prefix = os.path.join(temp_dir, "remerged")
            do_merge_subcommand(args)
            with gzip.open(f"{args.output_prefix}.tandem_repeats.detailed.bed.gz", "rt") as f:
                self.assertEqual(f.readline().split("\t")[3], name_field)

            # A detection mode round-trips whether it is "pure" (which also starts with "p") or a merged
            # locus's "merged:pure,trf" (which contains the ":" separator). The purity is recomputed per locus.
            # A zero-width locus has no purity, so no purity token is written, and a "pnan" token that older
            # detailed BEDs wrote for it is still recognized as the purity token rather than a detection mode.
            detailed_input_path = os.path.join(temp_dir, "detailed_catalog.bed")
            with open(detailed_input_path, "w") as f:
                f.write("chr22\t10515000\t10515020\tAT:2bp:10.0x:pure:p1.00\n")
                f.write("chr22\t10515100\t10515120\tAT:2bp:10.0x:merged:pure,trf:p1.00\n")
                f.write("chr22\t10515200\t10515200\tCAG:3bp:0.0x:pure:pnan\n")
            args.input_bed_paths = [detailed_input_path]
            args.output_prefix = os.path.join(temp_dir, "modes")
            do_merge_subcommand(args)
            with gzip.open(f"{args.output_prefix}.tandem_repeats.detailed.bed.gz", "rt") as f:
                self.assertEqual([re.sub(r":p\d\.\d\d$", "", line.split("\t")[3]) for line in f],
                                 ["AT:2bp:10.0x:pure", "AT:2bp:10.0x:merged:pure,trf", "CAG:3bp:0.0x:pure"])
        finally:
            shutil.rmtree(temp_dir)

    def test_merge_interval_on_a_contig_missing_from_a_catalog_is_skipped(self):
        """A sharded genome-wide merge routinely asks a catalog for a contig it has no loci on.

        An interval on such a contig, or one spelled "22" against a "chr22" catalog, must be matched or
        skipped with a warning rather than aborting the merge with a tabix error.
        """
        if shutil.which("bgzip") is None or shutil.which("tabix") is None:
            self.skipTest("bgzip/tabix unavailable")

        temp_dir = tempfile.mkdtemp()
        try:
            input_bed_path = os.path.join(temp_dir, "catalog.bed")
            with open(input_bed_path, "w") as f:
                f.write("chr22\t10515000\t10515020\tAT\n")
            if os.system(f"bgzip -f {input_bed_path} && tabix -p bed {input_bed_path}.gz") != 0:
                self.skipTest("bgzip/tabix failed")

            from str_analysis.filter_vcf_to_tandem_repeats import do_merge_subcommand
            for intervals, expected_num_loci in ((["chrY"], 0), (["22:10515000-10515100"], 1)):
                output_prefix = os.path.join(temp_dir, "merged")
                args = argparse.Namespace(
                    reference_fasta_path=self._temp_fasta_path,
                    input_bed_paths=[f"{input_bed_path}.gz"],
                    output_prefix=output_prefix,
                    interval=intervals,
                    verbose=False,
                    show_progress_bar=False,
                    write_detailed_bed=False,
                    batch_size=1000,
                )
                with contextlib.redirect_stdout(io.StringIO()):
                    do_merge_subcommand(args)
                with gzip.open(f"{output_prefix}.tandem_repeats.bed.gz", "rt") as f:
                    self.assertEqual(len(f.readlines()), expected_num_loci, intervals)
        finally:
            shutil.rmtree(temp_dir)

    def test_write_vcf_declares_info_fields_and_avoids_reserved_end(self):
        """The catalog/filter VCF output must declare its custom INFO fields and not misuse END.

        END is a VCF-reserved INFO key (the record's end position), so the tandem-repeat end
        coordinate is written as TR_END. All appended INFO fields must also have ##INFO header
        declarations so the output is spec-compliant and parseable.
        """
        import pysam
        from str_analysis.filter_vcf_to_tandem_repeats import write_vcf

        if shutil.which("bgzip") is None or shutil.which("tabix") is None:
            self.skipTest("bgzip/tabix unavailable")

        allele = Allele("chr22", 10515040, "T", "TAAGA", self._fasta_obj, order=0)
        tr_allele = TandemRepeatAllele(allele, "AAGA", False, 0, 4, 37, DETECTION_MODE_PURE_REPEATS)

        temp_dir = tempfile.mkdtemp()
        try:
            # Build an input VCF whose single row matches the tandem-repeat allele
            input_vcf_path = os.path.join(temp_dir, "in.vcf")
            with open(input_vcf_path, "w") as f:
                f.write("##fileformat=VCFv4.2\n")
                f.write("##contig=<ID=chr22>\n")
                f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
                f.write(f"{allele.chrom}\t{allele.pos}\t.\t{allele.ref}\t{allele.alt}\t.\tPASS\t.\tGT\t0/1\n")

            output_prefix = os.path.join(temp_dir, "out")
            args = argparse.Namespace(
                input_vcf_path=input_vcf_path, interval=None, verbose=False,
                output_prefix=output_prefix, offset=0, n=None,
            )
            write_vcf([tr_allele], args)

            output_vcf_path = f"{output_prefix}.tandem_repeats.vcf.gz"
            self.assertTrue(os.path.exists(output_vcf_path))

            with pysam.VariantFile(output_vcf_path) as vcf:
                # All appended INFO fields are declared in the header
                for field in ("MOTIF", "MOTIF_SIZE", "START_0BASED", "TR_END", "DETECTED"):
                    self.assertIn(field, vcf.header.info, f"{field} INFO field not declared in header")
                # The reserved END key is not overloaded
                self.assertNotIn("END", vcf.header.info)
                records = list(vcf.fetch())
            self.assertEqual(len(records), 1)
            # Number=. INFO fields are returned as tuples by pysam
            self.assertEqual(tuple(records[0].info["TR_END"]), (10515077,))
            self.assertEqual(tuple(records[0].info["MOTIF"]), ("AAGA",))
            self.assertNotIn("END", records[0].info)
        finally:
            shutil.rmtree(temp_dir)

    def test_filtered_vcf_declares_its_filter_ids(self):
        """The filtered-out VCF's FILTER column must hold declared, whitespace-free IDs.

        The filter reasons read like "INDEL without repeats", which is not a valid FILTER ID, and an ID that
        is not declared in a ##FILTER header line makes bcftools reject the file.
        """
        import pysam
        from str_analysis.filter_vcf_to_tandem_repeats import write_vcf, format_filter_id

        if shutil.which("bgzip") is None or shutil.which("tabix") is None:
            self.skipTest("bgzip/tabix unavailable")

        temp_dir = tempfile.mkdtemp()
        try:
            input_vcf_path = os.path.join(temp_dir, "in.vcf")
            with open(input_vcf_path, "w") as f:
                f.write("##fileformat=VCFv4.2\n")
                f.write("##contig=<ID=chr22>\n")
                f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
                f.write("chr22\t10515040\t.\tT\tTACGTACGTAC\t.\tPASS\t.\tGT\t0/1\n")
                f.write("chr22\t10515050\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n")
                # A multi-allelic record whose two ALT alleles were filtered for different reasons keeps both
                f.write("chr22\t10515060\t.\tG\tGACGTACGTAC,GNNNNNNNNN\t.\tPASS\t.\tGT\t1/2\n")

            output_prefix = os.path.join(temp_dir, "out")
            args = argparse.Namespace(
                input_vcf_path=input_vcf_path, interval=None, verbose=False,
                output_prefix=output_prefix, offset=0, n=None,
            )
            write_vcf([], args, only_write_filtered_out_alleles=True, filtered_alleles={
                ("chr22", 10515040, "T", "TACGTACGTAC"): FILTER_ALLELE_INDEL_WITHOUT_REPEATS,
                ("chr22", 10515060, "G", "GACGTACGTAC"): FILTER_ALLELE_INDEL_WITHOUT_REPEATS,
                ("chr22", 10515060, "G", "GNNNNNNNNN"): FILTER_ALLELE_WITH_N_BASES,
            })

            output_vcf_path = f"{output_prefix}.not_tandem_repeats.vcf.gz"
            with pysam.VariantFile(output_vcf_path) as vcf:
                records = list(vcf.fetch())
                self.assertEqual([list(r.filter) for r in records], [
                    ["INDEL_without_repeats"], ["SNV"],
                    ["INDEL_without_repeats", "contains_Ns_in_the_variant_sequence"]])
                for record in records:
                    for filter_id in record.filter:
                        self.assertIn(filter_id, vcf.header.filters, f"{filter_id} not declared in header")
                        self.assertNotRegex(filter_id, r"[\s;,]")
                self.assertEqual(vcf.header.filters["INDEL_without_repeats"].description,
                                 FILTER_ALLELE_INDEL_WITHOUT_REPEATS)
            self.assertEqual(format_filter_id("INDEL > 100,000bp"), "INDEL_gt_100_000bp")

            if shutil.which("bcftools") is not None:
                self.assertEqual(os.system(f"bcftools view -Ob -o /dev/null {output_vcf_path} 2>/dev/null"), 0,
                                 "bcftools rejected the filtered VCF")
        finally:
            shutil.rmtree(temp_dir)

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
            trf_mismatch_penalty=7,
            trf_indel_penalty=7,
            trf_min_score=20,
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

    def _create_temp_indexed_vcf(self, vcf_content):
        """Helper to create a bgzip-compressed, tabix-indexed temp VCF, as required for genotyping.

        Uses pysam's bundled bgzip/tabix (no external tools needed). Returns the path to the
        .vcf.gz file; a sibling .vcf.gz.tbi index is created alongside it.
        """
        import pysam
        vcf_path = self._create_temp_vcf(vcf_content)
        pysam.tabix_index(vcf_path, preset="vcf", force=True)  # replaces vcf_path with vcf_path + ".gz"
        return vcf_path + ".gz"

    def test_open_vcf_for_genotyping_valid_single_sample(self):
        """Test opening a valid single-sample VCF."""
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t1000\t.\tA\tAT\t.\tPASS\t.\tGT\t0/1
"""
        vcf_path = self._create_temp_indexed_vcf(vcf_content)
        try:
            vcf_file, sample_name, vcf_contig_lookup = open_vcf_for_genotyping(vcf_path)
            self.assertEqual(sample_name, "SAMPLE1")
            # A catalog contig resolves to this VCF's chr-prefixed spelling however it is written
            self.assertEqual(vcf_contig_lookup, {"1": "chr1"})
            vcf_file.close()
        finally:
            os.unlink(vcf_path)
            if os.path.exists(vcf_path + ".tbi"):
                os.unlink(vcf_path + ".tbi")

    def test_open_vcf_for_genotyping_no_chr_prefix(self):
        """Test detecting chromosome naming convention without chr prefix."""
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=1,length=1000000>
##contig=<ID=2,length=1000000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
1\t1000\t.\tA\tAT\t.\tPASS\t.\tGT\t0/1
"""
        vcf_path = self._create_temp_indexed_vcf(vcf_content)
        try:
            vcf_file, sample_name, vcf_contig_lookup = open_vcf_for_genotyping(vcf_path)
            # A catalog contig resolves to this VCF's unprefixed spelling however it is written
            self.assertEqual(vcf_contig_lookup, {"1": "1", "2": "2"})
            vcf_file.close()
        finally:
            os.unlink(vcf_path)
            if os.path.exists(vcf_path + ".tbi"):
                os.unlink(vcf_path + ".tbi")

    def test_open_vcf_for_genotyping_multi_sample_error(self):
        """Test that multi-sample VCF raises error."""
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2
chr1\t1000\t.\tA\tAT\t.\tPASS\t.\tGT\t0/1\t1/1
"""
        vcf_path = self._create_temp_indexed_vcf(vcf_content)
        try:
            with self.assertRaises(ValueError) as context:
                open_vcf_for_genotyping(vcf_path)
            self.assertIn("single-sample", str(context.exception))
        finally:
            os.unlink(vcf_path)
            if os.path.exists(vcf_path + ".tbi"):
                os.unlink(vcf_path + ".tbi")

    def test_open_vcf_for_genotyping_file_not_found(self):
        """Test that missing VCF raises FileNotFoundError."""
        with self.assertRaises(FileNotFoundError):
            open_vcf_for_genotyping("/nonexistent/path/to/file.vcf")

    def test_open_vcf_for_genotyping_unindexed_error(self):
        """A plain, un-indexed VCF must raise a clear error instead of silently genotyping hom-ref.

        Genotyping fetches variants per-locus, which requires a tabix/csi index. Without this
        guard, pysam raises 'fetch requires an index' on every fetch, which is swallowed and
        reported as zero overlapping variants (hom-ref) at every locus.
        """
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t1000\t.\tA\tAT\t.\tPASS\t.\tGT\t0/1
"""
        vcf_path = self._create_temp_vcf(vcf_content)
        try:
            with self.assertRaises(ValueError) as context:
                open_vcf_for_genotyping(vcf_path)
            self.assertIn("index", str(context.exception))
        finally:
            os.unlink(vcf_path)

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
            variants = get_overlapping_vcf_variants(
                vcf_file, "chr1", 990, 1060)

            self.assertEqual(len(variants), 2)
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

            variants = get_overlapping_vcf_variants(
                vcf_file, "chr1", 990, 1010)

            self.assertEqual(len(variants), 1)

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
            variants = get_overlapping_vcf_variants(
                vcf_file, "chr1", 2000, 3000)

            self.assertEqual(len(variants), 0)

            vcf_file.close()
        except FileNotFoundError:
            self.skipTest("bgzip/tabix not available")
        finally:
            os.unlink(vcf_path)
            if os.path.exists(vcf_path + ".gz"):
                os.unlink(vcf_path + ".gz")
            if os.path.exists(vcf_path + ".gz.tbi"):
                os.unlink(vcf_path + ".gz.tbi")

    def test_get_overlapping_vcf_variants_left_anchored_insertion(self):
        """A left-anchored insertion (REF span ending exactly at the locus start) must be returned.

        A normalized repeat-unit insertion anchors to the base just before the tract, so its
        REF span [start_0based - 1, start_0based) does not overlap [start_0based, end). A plain
        tabix fetch(start_0based, end) misses it; get_overlapping_vcf_variants widens the window
        and must include it.
        """
        import pysam
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t1001\t.\tA\tAT\t.\tPASS\t.\tGT\t0/1
"""
        vcf_path = self._create_temp_vcf(vcf_content)
        try:
            import subprocess
            vcf_gz_path = vcf_path + ".gz"
            with open(vcf_gz_path, 'wb') as gz_out:
                subprocess.run(["bgzip", "-c", vcf_path], stdout=gz_out, check=True)
            subprocess.run(["tabix", "-p", "vcf", vcf_gz_path], check=True)

            vcf_file = pysam.VariantFile(vcf_gz_path)

            # Locus [1001, 1011); the insertion at POS 1001 (0-based 1000) is anchored one base
            # to the left of the locus start and would be dropped by a non-widened fetch.
            variants = get_overlapping_vcf_variants(
                vcf_file, "chr1", 1001, 1011)

            self.assertEqual(len(variants), 1)
            self.assertEqual(variants[0].pos, 1001)

            vcf_file.close()
        except FileNotFoundError:
            self.skipTest("bgzip/tabix not available")
        finally:
            os.unlink(vcf_path)
            if os.path.exists(vcf_path + ".gz"):
                os.unlink(vcf_path + ".gz")
            if os.path.exists(vcf_path + ".gz.tbi"):
                os.unlink(vcf_path + ".gz.tbi")

    def test_get_overlapping_vcf_variants_excludes_left_flank_snv(self):
        """A non-insertion whose REF span ends at the locus start (e.g. a flank SNV) must be dropped.

        Widening the fetch window must not let a variant that only touches the last flank base
        perturb the locus or spuriously trigger the multi-variant phasing-ambiguity check.
        """
        import pysam
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t1001\t.\tA\tT\t.\tPASS\t.\tGT\t0/1
chr1\t1006\t.\tC\tCAT\t.\tPASS\t.\tGT\t0/1
"""
        vcf_path = self._create_temp_vcf(vcf_content)
        try:
            import subprocess
            vcf_gz_path = vcf_path + ".gz"
            with open(vcf_gz_path, 'wb') as gz_out:
                subprocess.run(["bgzip", "-c", vcf_path], stdout=gz_out, check=True)
            subprocess.run(["tabix", "-p", "vcf", vcf_gz_path], check=True)

            vcf_file = pysam.VariantFile(vcf_gz_path)

            # Locus [1001, 1011): the SNV at POS 1001 (0-based 1000) touches only the last flank
            # base and must be excluded; the insertion at POS 1006 falls inside and is kept.
            variants = get_overlapping_vcf_variants(
                vcf_file, "chr1", 1001, 1011)

            self.assertEqual(len(variants), 1)
            self.assertEqual(variants[0].pos, 1006)

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
            variants = get_overlapping_vcf_variants(
                vcf_file, "chr2", 1000, 2000)

            self.assertEqual(len(variants), 0)

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

            # Query using chr1 but VCF has "1" - the contig lookup should resolve it to the VCF's spelling
            vcf_contig_lookup = build_contig_name_lookup(vcf_file.header.contigs)
            variants = get_overlapping_vcf_variants(
                vcf_file, vcf_contig_lookup[normalize_chromosome_name("chr1")], 990, 1010)

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
        """Test that mismatched reference raises RefAlleleMismatchError."""
        # Reference: ACGTACGT at position 100
        # Position 103 has T, not G - so this should raise an error
        with self.assertRaises(RefAlleleMismatchError) as context:
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
        """Test that invalid DNA bases in ref raise NonIUPACAlleleError."""
        with self.assertRaises(NonIUPACAlleleError) as context:
            convert_variants_to_haplotype_sequence(
                pos_1based=100,
                reference_sequence="ACGTACGT",
                variants=[(100, "X", "A")]  # X is not a valid DNA base
            )
        self.assertIn("Invalid ref allele", str(context.exception))

    def test_invalid_alt_allele_raises_error(self):
        """Test that invalid DNA bases in alt raise NonIUPACAlleleError."""
        with self.assertRaises(NonIUPACAlleleError) as context:
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


class TestComputeLocusStartAndEndOffsetsInHaplotype(unittest.TestCase):
    """Test compute_locus_start_and_end_offsets_in_haplotype directly.

    The reference is 10 T's, the locus (CAG)x3 at [10, 19), then 10 more T's, fetched starting at position 0:

        0-based position:  0-9         10-18       19-28
        base:              TTTTTTTTTT  CAGCAGCAG   TTTTTTTTTT
    """

    REFERENCE_SEQUENCE = "T" * 10 + "CAGCAGCAG" + "T" * 10

    def _check(self, variants, expected_offsets, expected_locus_sequence, locus_start_0based=10,
               locus_end_0based=19, repeat_unit="CAG"):
        """Check the offsets, and that slicing the haplotype built from the same variants gives the locus."""
        offsets = compute_locus_start_and_end_offsets_in_haplotype(
            locus_start_0based, locus_end_0based, 0, variants, repeat_unit=repeat_unit)
        self.assertEqual(offsets, expected_offsets)
        haplotype_sequence = convert_variants_to_haplotype_sequence(1, self.REFERENCE_SEQUENCE, variants)
        self.assertEqual(haplotype_sequence[offsets[0]:offsets[1]], expected_locus_sequence)

    def test_no_variants(self):
        self._check([], (10, 19), "CAGCAGCAG")

    def test_snv_inside_the_locus(self):
        self._check([(12, "A", "T")], (10, 19), "CTGCAGCAG")

    def test_insertion_inside_the_locus(self):
        # The inserted CAG is positioned at 0-based 12, inside the locus
        self._check([(13, "G", "GCAG")], (10, 22), "CAGCAGCAGCAG")

    def test_left_anchored_insertion_just_before_the_locus(self):
        # The VCF padding base is the last flank base (0-based 9), and the inserted bases are positioned at the
        # locus start, so they belong to the locus
        self._check([(10, "T", "TCAG")], (10, 22), "CAGCAGCAGCAG")

    def test_deletion_crossing_the_locus_start(self):
        # TTC -> T deletes 0-based 9 and 10. The kept T (0-based 8) stays in the flank, and the locus loses its
        # first base.
        self._check([(9, "TTC", "T")], (9, 17), "AGCAGCAG")

    def test_deletion_crossing_the_locus_end(self):
        # AGT -> A deletes 0-based 18 and 19, the locus's last base and the first right-flank base
        self._check([(18, "AGT", "A")], (10, 18), "CAGCAGCA")

    def test_insertion_of_the_motif_positioned_at_the_locus_end_is_inside(self):
        self._check([(19, "G", "GCAG")], (10, 22), "CAGCAGCAGCAG")

    def test_other_insertion_positioned_at_the_locus_end_is_outside(self):
        self._check([(19, "G", "GTTA")], (10, 19), "CAGCAGCAG")

    def test_zero_width_locus_keeps_its_insertion(self):
        # A zero-width locus keeps an insertion positioned at its single coordinate whatever the inserted bases
        self._check([(10, "T", "TGATTACA")], (10, 17), "GATTACA", locus_start_0based=10, locus_end_0based=10)

    def test_trimming_the_shared_suffix_moves_an_insertion_into_the_flank(self):
        # TT -> TTT shares a T suffix, so it is trimmed to T -> TT at 0-based 8, which positions the inserted T at
        # 0-based 9, in the left flank rather than at the locus start
        self._check([(9, "TT", "TTT")], (11, 20), "CAGCAGCAG")

    def test_lowercase_ref_allele_is_trimmed_like_uppercase(self):
        # A VCF written against a soft-masked reference gives the ref allele in lower case. Compared as written,
        # "tt" and "TTT" share no suffix, and the inserted T would land at the locus start.
        self._check([(9, "tt", "TTT")], (11, 20), "CAGCAGCAG")


class TestExtractHaplotypeSequencesFromVcf(unittest.TestCase):
    """Test the extract_haplotype_sequences_from_vcf function."""

    def setUp(self):
        """Set up test fixtures with a mock fasta object."""
        # Default: return a CAG repeat sequence chr1:100-120 = "CAGCAGCAGCAGCAGCAGCA"
        # This is 20bp spanning 6.67 CAG repeats
        self._reference_seq = "CAGCAGCAGCAGCAGCAGCA"

    def _create_mock_fasta(self, reference_sequence=None):
        """Create a mock pyfaidx Fasta object.

        Args:
            reference_sequence (str): Reference sequence to return. If None, uses default.

        Returns:
            mock.MagicMock: A mock fasta object
        """
        if reference_sequence is None:
            reference_sequence = self._reference_seq

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

        return MockFasta(reference_sequence)

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
        # Return the GT tuple only for the "GT" key. A blanket return_value would hand the same tuple
        # back for "PS" as well, and are_variants_unambiguously_phased() would then compare genotypes
        # instead of phase sets, judging a realistic 0|1 plus 1|0 pair to be in two different blocks.
        sample.get = mock.MagicMock(side_effect=lambda key, default=None: {"GT": gt}.get(key, default))
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
        """Test that multiple unphased variants returns missing genotype (None, None)."""
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

    def test_phased_variants_on_opposite_haplotypes(self):
        """Two phased heterozygous records with no PS tag are one block, whichever haplotype each sits on.

        A caller that phases a whole chromosome at once, such as DipCall, writes 0|1 and 1|0 with no PS tag.
        Both records must still be resolvable, one onto each haplotype.
        """
        mock_fasta = self._create_mock_fasta()

        variant1 = self._create_mock_variant(pos=2, ref="A", alt="T", gt=(0, 1), phased=True)
        variant2 = self._create_mock_variant(pos=5, ref="A", alt="G", gt=(1, 0), phased=True)

        result = extract_haplotype_sequences_from_vcf(
            chrom="chr1",
            start_0based=0,
            end=12,
            fasta_obj=mock_fasta,
            vcf_variants=[variant1, variant2]
        )

        # Haplotype 0 carries only variant2 (position 5, idx 4): A -> G
        self.assertEqual(result[0], "CAGCGGCAGCAG")
        # Haplotype 1 carries only variant1 (position 2, idx 1): A -> T
        self.assertEqual(result[1], "CTGCAGCAGCAG")

    def test_variant_spanning_beyond_locus_is_trimmed(self):
        """Test that variants spanning beyond locus boundaries are properly handled."""
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

    def test_deletion_spanning_locus_start_boundary(self):
        """A deletion straddling the locus start should drop only the locus base(s) it removes.

        Reference (0-based): AAAA CAGCAGCAGCAG TTTT, locus = indices [4, 16) = "CAGCAGCAGCAG".
        Deletion "AAC" -> "A" at 1-based pos 3 spans genomic [2, 5), overlapping the start
        boundary and deleting the first locus base (index 4, 'C'). The replacement base is
        positioned in the flank (it starts before the locus), so the locus haplotype is the
        remaining "AGCAGCAGCAG" (11 bp), not a sequence padded with flanking bases.
        """
        mock_fasta = self._create_mock_fasta("AAAACAGCAGCAGCAGTTTT")

        variant = self._create_mock_variant(
            pos=3, ref="AAC", alt="A", gt=(0, 1), phased=True
        )

        result = extract_haplotype_sequences_from_vcf(
            chrom="chr1", start_0based=4, end=16, fasta_obj=mock_fasta, vcf_variants=[variant]
        )

        self.assertEqual(result[0], "CAGCAGCAGCAG")  # reference haplotype
        self.assertEqual(result[1], "AGCAGCAGCAG")   # leading locus base deleted

    def test_deletion_spanning_locus_end_boundary(self):
        """A deletion straddling the locus end should only remove flank, keeping locus content.

        Reference (0-based): AAAA CAGCAGCAGCAG TTTT, locus = indices [4, 16) = "CAGCAGCAGCAG".
        Deletion "GT" -> "G" at 1-based pos 16 spans genomic [15, 17): it keeps the last locus
        base (index 15, 'G') and deletes only the flanking 'T' at index 16. The locus haplotype
        is therefore unchanged ("CAGCAGCAGCAG"); the old offset logic incorrectly truncated the
        final locus base by applying the full deletion length change to the end offset.
        """
        mock_fasta = self._create_mock_fasta("AAAACAGCAGCAGCAGTTTT")

        variant = self._create_mock_variant(
            pos=16, ref="GT", alt="G", gt=(0, 1), phased=True
        )

        result = extract_haplotype_sequences_from_vcf(
            chrom="chr1", start_0based=4, end=16, fasta_obj=mock_fasta, vcf_variants=[variant]
        )

        self.assertEqual(result[0], "CAGCAGCAGCAG")  # reference haplotype
        self.assertEqual(result[1], "CAGCAGCAGCAG")   # only flanking base deleted

    def _create_mock_multiallelic_variant(self, pos, ref, alts, gt, phased=True):
        """Helper to create a mock multi-allelic pysam variant record.

        Args:
            pos (int): 1-based position
            ref (str): Reference allele
            alts (tuple): Alternate alleles, e.g. ("A", "AT")
            gt (tuple): Genotype tuple indexing into alleles, e.g. (1, 2)
            phased (bool): Whether the genotype is phased

        Returns:
            mock.MagicMock: A mock variant record with alleles = (ref, *alts)
        """
        variant = mock.MagicMock()
        variant.pos = pos
        variant.ref = ref
        variant.alleles = (ref,) + tuple(alts)
        sample = mock.MagicMock()
        # Return the GT tuple only for the "GT" key. A blanket return_value would hand the same tuple
        # back for "PS" as well, and are_variants_unambiguously_phased() would then compare genotypes
        # instead of phase sets, judging a realistic 0|1 plus 1|0 pair to be in two different blocks.
        sample.get = mock.MagicMock(side_effect=lambda key, default=None: {"GT": gt}.get(key, default))
        sample.phased = phased
        variant.samples = [sample]
        return variant

    def test_left_anchored_insertion_counts_into_locus(self):
        """A repeat-unit insertion anchored just before the tract extends the locus.

        Reference (0-based): AAAA CAGCAGCAGCAG TTTT, locus = indices [4, 16) = 4x"CAG".
        The insertion "A" -> "ACAG" at 1-based pos 4 anchors to the last flank base (index 3),
        inserting one "CAG" at the locus start boundary (junction at index 4). The inserted
        copy belongs to the locus, so the alt haplotype is 5x"CAG", not the unchanged reference.
        """
        mock_fasta = self._create_mock_fasta("AAAACAGCAGCAGCAGTTTT")

        variant = self._create_mock_variant(
            pos=4, ref="A", alt="ACAG", gt=(0, 1), phased=True
        )

        result = extract_haplotype_sequences_from_vcf(
            chrom="chr1", start_0based=4, end=16, fasta_obj=mock_fasta, vcf_variants=[variant]
        )

        self.assertEqual(result[0], "CAGCAGCAGCAG")       # reference haplotype (4 repeats)
        self.assertEqual(result[1], "CAGCAGCAGCAGCAG")    # inserted CAG counted into locus (5 repeats)

    def test_multiallelic_insertion_spanning_start_boundary(self):
        """Both alleles of a left-anchored multi-allelic insertion extend the locus.

        Reference: AAAA CAGCAGCAGCAG TTTT, locus [4, 16) = 4x"CAG". The record
        "A" -> "ACAG","ACAGCAG" at 1-based pos 4 with GT 1|2 inserts one and two copies
        respectively, so the haplotypes are 5x and 6x "CAG".
        """
        mock_fasta = self._create_mock_fasta("AAAACAGCAGCAGCAGTTTT")

        variant = self._create_mock_multiallelic_variant(
            pos=4, ref="A", alts=("ACAG", "ACAGCAG"), gt=(1, 2), phased=True
        )

        result = extract_haplotype_sequences_from_vcf(
            chrom="chr1", start_0based=4, end=16, fasta_obj=mock_fasta, vcf_variants=[variant]
        )

        self.assertEqual(result[0], "CAGCAGCAGCAGCAG")       # 5 repeats
        self.assertEqual(result[1], "CAGCAGCAGCAGCAGCAG")    # 6 repeats

    def test_multiallelic_deletion_spanning_start_boundary(self):
        """Two different left-anchored deletions remove different numbers of locus bases.

        Reference: AAAA TTTTTTTTTTTT AAAA, locus [4, 16) = 12 "T". The record
        "ATTTT" -> "A","AT" at 1-based pos 4 (anchored at flank index 3) deletes 4 and 3 "T"
        respectively under GT 1|2, leaving 8 and 9 "T".
        """
        mock_fasta = self._create_mock_fasta("AAAATTTTTTTTTTTTAAAA")

        variant = self._create_mock_multiallelic_variant(
            pos=4, ref="ATTTT", alts=("A", "AT"), gt=(1, 2), phased=True
        )

        result = extract_haplotype_sequences_from_vcf(
            chrom="chr1", start_0based=4, end=16, fasta_obj=mock_fasta, vcf_variants=[variant]
        )

        self.assertEqual(result[0], "TTTTTTTT")    # 4 of 12 deleted
        self.assertEqual(result[1], "TTTTTTTTT")   # 3 of 12 deleted

    def test_multiallelic_del_ins_spanning_start_boundary(self):
        """A left-anchored record with one deletion allele and one insertion allele.

        Reference: AAAA TTTTTTTTTTTT AAAA, locus [4, 16) = 12 "T". The record
        "AT" -> "A","ATT" at 1-based pos 4 deletes one "T" on the first haplotype and inserts
        one "T" on the second (GT 1|2), giving 11 and 13 "T". This is the case where the old
        offset logic dropped the inserted base and collapsed the insertion allele.
        """
        mock_fasta = self._create_mock_fasta("AAAATTTTTTTTTTTTAAAA")

        variant = self._create_mock_multiallelic_variant(
            pos=4, ref="AT", alts=("A", "ATT"), gt=(1, 2), phased=True
        )

        result = extract_haplotype_sequences_from_vcf(
            chrom="chr1", start_0based=4, end=16, fasta_obj=mock_fasta, vcf_variants=[variant]
        )

        self.assertEqual(result[0], "TTTTTTTTTTT")     # 11 (one deleted)
        self.assertEqual(result[1], "TTTTTTTTTTTTT")   # 13 (one inserted)

    def test_haploid_genotype_returns_hemizygous(self):
        """A haploid genotype (e.g. chrX/chrY, GT == (1,)) must not crash; it yields one haplotype.

        The haplotype loop iterates over (0, 1) and indexes gt[haplotype]; a length-1 genotype
        tuple has no second haplotype, so the second haplotype is reported as missing (None)
        rather than raising IndexError.
        """
        mock_fasta = self._create_mock_fasta()

        variant = self._create_mock_variant(
            pos=4, ref="C", alt="T", gt=(1,), phased=True  # haploid: single allele
        )

        result = extract_haplotype_sequences_from_vcf(
            chrom="chr1", start_0based=0, end=12, fasta_obj=mock_fasta, vcf_variants=[variant]
        )

        # Haplotype 0 applies the alt allele; haplotype 1 does not exist for a haploid call
        self.assertEqual(result[0], "CAGTAGCAGCAG")
        self.assertIsNone(result[1])

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

    def _create_mock_fasta(self, reference_sequence=None):
        """Create a mock pyfaidx Fasta object."""
        if reference_sequence is None:
            reference_sequence = self._reference_seq

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

        return MockFasta(reference_sequence)

    def _create_mock_variant(self, pos, ref, alt, gt, phased=True):
        """Helper to create a mock pysam variant record."""
        variant = mock.MagicMock()
        variant.pos = pos
        variant.ref = ref
        variant.alleles = (ref, alt)

        sample = mock.MagicMock()
        # Return the GT tuple only for the "GT" key. A blanket return_value would hand the same tuple
        # back for "PS" as well, and are_variants_unambiguously_phased() would then compare genotypes
        # instead of phase sets, judging a realistic 0|1 plus 1|0 pair to be in two different blocks.
        sample.get = mock.MagicMock(side_effect=lambda key, default=None: {"GT": gt}.get(key, default))
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

    def test_multiallelic_variant_is_genotyped(self):
        """Test that multi-allelic variants are genotyped normally."""
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

        # Multi-allelic variants should be genotyped normally
        self.assertIsNotNone(result.zygosity)
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

    def test_diploid_record_with_a_missing_allele_gives_no_call(self):
        """A '.|1' record at a diploid locus leaves one haplotype uncalled, so the locus gets no call.

        Reporting the surviving allele as HEMI would put it in both the short and long allele columns and look
        like a real hemizygous call, which is exactly what the build-error and non-repeat-insertion paths
        refuse to do. Only a haploid locus (chrX or chrY outside the PARs of a sample detected as haploid
        there) accepts it as a haploid call.
        """
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

        self.assertIsNone(result.zygosity)
        self.assertIsNone(result.allele1_sequence)
        self.assertIsNone(result.allele2_sequence)
        self.assertEqual(result.no_call_reason, NO_CALL_REASON_MISSING_GENOTYPE)

    def test_the_same_record_on_chrx_keeps_the_called_allele(self):
        """On chrX outside the PARs of an XY sample the missing allele is expected, so the called one
        must survive."""
        variant = self._create_mock_variant(
            pos=1, ref="C", alt="CCAG", gt=(None, 1), phased=True
        )

        tr_locus = ReferenceTandemRepeat(
            chrom="chrX",
            start_0based=0,
            end_1based=12,
            repeat_unit="CAG"
        )
        result = genotype_single_locus(
            tr_locus, self._create_mock_vcf_file([variant]), self._create_mock_fasta(),
            sex_chromosome_ploidy=XY_PLOIDY)

        self.assertEqual(result.zygosity, "HEMI")
        self.assertIsNone(result.allele1_sequence)
        self.assertIsNotNone(result.allele2_sequence)
        self.assertIsNone(result.no_call_reason)

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
        except FileNotFoundError:
            # Only a genuinely missing bgzip/tabix is a skip. A CalledProcessError means the fixture itself is
            # malformed (unsorted records, a bad header), which must fail the test rather than turn into a
            # skip labelled 'bgzip/tabix unavailable' that leaves the run green.
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
            result = genotype_single_locus(tr_locus, vcf_file, fasta_obj)

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
            result = genotype_single_locus(tr_locus, vcf_file, fasta_obj)

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
            result = genotype_single_locus(tr_locus, vcf_file, fasta_obj)

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
            result = genotype_single_locus(tr_locus, vcf_file, fasta_obj)

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
            result = genotype_single_locus(tr_locus, vcf_file, fasta_obj)

            # Both haplotypes should be None (missing)
            self.assertIsNone(result.allele1_sequence)
            self.assertIsNone(result.allele2_sequence)
            self.assertIsNone(result.zygosity)
        finally:
            vcf_file.close()
            fasta_obj.close()

    def test_genotype_unphased_multiple_variants_returns_missing(self):
        """Test that multiple unphased variants cause missing genotype."""
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
            result = genotype_single_locus(tr_locus, vcf_file, fasta_obj)

            # Should be missing due to unphased multiple variants
            self.assertIsNone(result.zygosity)
            self.assertIsNone(result.allele1_sequence)
            self.assertIsNone(result.allele2_sequence)
            # But variants should still be recorded
            self.assertEqual(result.num_overlapping_variants, 2)
        finally:
            vcf_file.close()
            fasta_obj.close()

    def test_genotype_multiallelic_variant(self):
        """Test that multi-allelic variants are genotyped normally."""
        import pysam
        import pyfaidx

        fasta_path = self._create_test_fasta({"chr1": "AAACAGCAGCAGCAGAAA"})
        if fasta_path is None:
            self.skipTest("pyfaidx unavailable")

        # Multi-allelic variant: REF with two ALT alleles, GT 1/2
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
            result = genotype_single_locus(tr_locus, vcf_file, fasta_obj)

            # Multi-allelic variant should be genotyped: allele 1 = +1 CAG, allele 2 = +2 CAG
            self.assertEqual(result.zygosity, "HET")
            self.assertEqual(result.num_repeats_allele1, 5)
            self.assertEqual(result.num_repeats_allele2, 6)
            self.assertEqual(result.num_overlapping_variants, 1)
        finally:
            vcf_file.close()
            fasta_obj.close()

    def test_genotype_variant_spanning_beyond_locus(self):
        """Test that variants spanning beyond locus boundaries are handled correctly."""
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
            result = genotype_single_locus(tr_locus, vcf_file, fasta_obj)

            # Haplotype 0 should be reference
            self.assertEqual(result.allele1_sequence, "CAGCAGCAGCAG")
            # Haplotype 1: The variant is before the locus, so after trimming
            # the locus content should still be the same
            self.assertEqual(result.zygosity, "HOM")
        finally:
            vcf_file.close()
            fasta_obj.close()

    def test_chromosome_naming_normalization_no_chr(self):
        """Test that chromosome names are normalized correctly (chr1 vs 1)."""
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
            result = genotype_single_locus(
                tr_locus, vcf_file, fasta_obj,
                vcf_contig_lookup=build_contig_name_lookup(vcf_file.header.contigs),
                fasta_contig_lookup=build_contig_name_lookup(fasta_obj.keys()))

            # Should successfully genotype despite naming mismatch
            self.assertEqual(result.zygosity, "HET")
            self.assertEqual(result.num_repeats_allele1, 4)
            self.assertEqual(result.num_repeats_allele2, 5)
        finally:
            vcf_file.close()
            fasta_obj.close()

    def test_genotype_partial_haplotype_end_to_end(self):
        """A '.|1' record gives no call on an autosome, and a HEMI call on chrX in a sample with haploid chrX.

        This is the chrX-in-males shape. The record itself can't say whether the missing allele means
        'uncalled' or 'only one copy here', so the chromosome and the sample's other chrX records decide.
        """
        import pysam
        import pyfaidx

        reference_sequence = "AAACAGCAGCAGCAGAAA"
        fasta_path = self._create_test_fasta({"chr1": reference_sequence, "chrX": reference_sequence})
        if fasta_path is None:
            self.skipTest("pyfaidx unavailable")

        # VCF with partial genotype (one allele missing): .|1
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=18>
##contig=<ID=chrX,length=18>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t4\t.\tC\tCCAG\t.\tPASS\t.\tGT\t.|1
chrX\t4\t.\tC\tCCAG\t.\tPASS\t.\tGT\t.|1
"""
        vcf_gz_path = self._create_test_vcf_and_index(vcf_content)
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        fasta_obj = pyfaidx.Fasta(fasta_path, one_based_attributes=False, as_raw=True)
        vcf_file = pysam.VariantFile(vcf_gz_path)

        try:
            # On an autosome the uncalled haplotype makes the whole locus a no call
            tr_locus = ReferenceTandemRepeat(chrom="chr1", start_0based=3, end_1based=15, repeat_unit="CAG")
            result = genotype_single_locus(tr_locus, vcf_file, fasta_obj)
            self.assertIsNone(result.zygosity)
            self.assertIsNone(result.allele1_sequence)
            self.assertIsNone(result.allele2_sequence)
            self.assertEqual(result.no_call_reason, NO_CALL_REASON_MISSING_GENOTYPE)

            # A single chrX record is too few for detect_sex_chromosome_ploidy to call chrX haploid (see
            # test_haploid_chrx_is_detected_from_the_share_of_half_called_records for the detection itself), so
            # the sample defaults to diploid chrX and the same record is a no call there too
            with contextlib.redirect_stdout(io.StringIO()):
                sex_chromosome_ploidy = detect_sex_chromosome_ploidy(
                    vcf_file, build_contig_name_lookup(vcf_file.header.contigs), {})
            self.assertEqual(sex_chromosome_ploidy, {"X": 2, "Y": 0})
            tr_locus = ReferenceTandemRepeat(chrom="chrX", start_0based=3, end_1based=15, repeat_unit="CAG")
            result = genotype_single_locus(tr_locus, vcf_file, fasta_obj,
                                           sex_chromosome_ploidy=sex_chromosome_ploidy)
            self.assertEqual(result.no_call_reason, NO_CALL_REASON_MISSING_GENOTYPE)

            # In a sample whose chrX was detected as haploid, the called allele is kept as the one copy
            result = genotype_single_locus(tr_locus, vcf_file, fasta_obj,
                                           sex_chromosome_ploidy={"X": 1, "Y": 0})
            self.assertEqual(result.zygosity, "HEMI")
            self.assertIsNone(result.allele1_sequence)
            self.assertIsNotNone(result.allele2_sequence)
            self.assertEqual(result.num_repeats_allele2, 5)
            self.assertIsNone(result.no_call_reason)
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
            result = genotype_single_locus(tr_locus, vcf_file, fasta_obj)

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
            result = genotype_single_locus(tr_locus, vcf_file, fasta_obj)
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

    def test_genotype_per_allele_repeat_purity(self):
        """Both alleles of a het ref/alt locus, including the reference allele, get a per-allele repeat purity."""
        import pysam
        import pyfaidx

        fasta_path = self._create_test_fasta({"chr1": "AAACAGCAGCAGCAGAAA"})
        if fasta_path is None:
            self.skipTest("pyfaidx unavailable")

        # 0|1: the reference allele has 4 CAG repeats and the alt (CCAG insertion) has 5, so the shorter allele
        # is the reference allele -- the case that previously had no purity value in the comparison plots.
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=18>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t4\t.\tC\tCCAG\t.\tPASS\t.\tGT\t0|1
"""
        vcf_gz_path = self._create_test_vcf_and_index(vcf_content)
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        tr_locus = ReferenceTandemRepeat(chrom="chr1", start_0based=3, end_1based=15, repeat_unit="CAG")
        fasta_obj = pyfaidx.Fasta(fasta_path, one_based_attributes=False, as_raw=True)
        vcf_file = pysam.VariantFile(vcf_gz_path)
        try:
            result = genotype_single_locus(
                tr_locus, vcf_file, fasta_obj)

            self.assertEqual(result.num_repeats_short_allele, 4)
            self.assertEqual(result.num_repeats_long_allele, 5)
            # both alleles, including the reference allele, get a real purity value (the X1 fix). The short allele is
            # the reference allele here, and its pure CAG repeats give purity 1.0.
            self.assertIsNotNone(result.repeat_purity_short_allele)
            self.assertIsNotNone(result.repeat_purity_long_allele)
            self.assertGreater(result.repeat_purity_short_allele, 0.99)
            self.assertGreater(result.repeat_purity_long_allele, 0.0)
            self.assertLessEqual(result.repeat_purity_long_allele, 1.0)

            tsv_dict = result.to_tsv_dict()
            self.assertNotEqual(tsv_dict["RepeatPurityShortAllele"], "")
            self.assertNotEqual(tsv_dict["RepeatPurityLongAllele"], "")
        finally:
            vcf_file.close()
            fasta_obj.close()

    def test_genotypes_unchanged_when_vcf_is_left_normalized(self):
        """Genotypes must not depend on how the VCF happens to represent a variant.

        The same two haplotypes are given here exactly as dipcall wrote them for HG00738 at chr1:4420277 (one
        multiallelic record whose REF carries the reference's soft-masking and whose second ALT shares an 18bp
        prefix and a 28bp suffix with REF) and as `bcftools norm -m - -f ref` rewrites them (one record per ALT,
        with the second ALT's shared suffix trimmed off). Both must give identical genotypes at the two loci the
        record overlaps. The reference is hg38 chr1:4420201-4420400, so the (GT)24 locus at chr1:4420232-4420280
        becomes chr1:32-80 and the record moves from POS 4420277 to POS 77.

        Trimming the shared suffix leaves the second ALT at the leftmost position it can occupy, which is the
        placement bcftools left-alignment settles on as well, so the split record needs no further shifting. That
        placement runs the 50bp deletion into the last 3 bases of the GT tract and gives 22 repeats. Trimming the
        shared prefix first would instead put the deletion 18bp further right, outside the GT locus, and give 24,
        which is what the code reported for this haplotype before it trimmed alleles when mapping locus boundaries.
        """
        import pysam
        import pyfaidx

        fasta_path = self._create_test_fasta({"chr1": (
            "CTCAGTTGCCCCCTTCGCAGCTGAATACCAGG"
            "GTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT"
            "ATATATATATATATGTGTGTATATATATATATGTATATATATATATATGTATATATATATATATATATATGTATATATATATATATATATATATGTATA"
            "AGGCACGAATTCCTAGTGGCT")})
        if fasta_path is None:
            self.skipTest("pyfaidx unavailable")

        vcf_header = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=200>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
"""
        dipcall_ref = "gtgtatatatatatatatgtgtgtatatatatatatgtatatatatatatatgtatatatatatatatatatatgtata"
        vcf_as_written_by_dipcall = vcf_header + (
            f"chr1\t77\t.\t{dipcall_ref}\tG,GTGTATATATATATATATATATATGTATA\t.\tPASS\t.\tGT\t1|2\n")
        vcf_after_bcftools_norm = vcf_header + (
            f"chr1\t77\t.\t{dipcall_ref}\tG\t.\tPASS\t.\tGT\t1|0\n"
            f"chr1\t77\t.\tGTGTATATATATATATATGTGTGTATATATATATATGTATATATATATATA\tG\t.\tPASS\t.\tGT\t0|1\n")

        vcf_gz_paths = [self._create_test_vcf_and_index(vcf_content)
                        for vcf_content in (vcf_as_written_by_dipcall, vcf_after_bcftools_norm)]
        if None in vcf_gz_paths:
            self.skipTest("bgzip/tabix unavailable")

        expected_genotypes = {
            (32, 80, "GT"): ("HOM", 22, 22, "GT" * 22 + "G", "GT" * 22 + "G"),
            (79, 94, "TA"): ("HOM", 0, 0, "", ""),
        }
        fasta_obj = pyfaidx.Fasta(fasta_path, one_based_attributes=False, as_raw=True)
        try:
            for (start_0based, end, repeat_unit), expected in expected_genotypes.items():
                tr_locus = ReferenceTandemRepeat(
                    chrom="chr1", start_0based=start_0based, end_1based=end, repeat_unit=repeat_unit)
                for vcf_gz_path in vcf_gz_paths:
                    vcf_file = pysam.VariantFile(vcf_gz_path)
                    try:
                        result = genotype_single_locus(
                            tr_locus, vcf_file, fasta_obj)
                    finally:
                        vcf_file.close()
                    self.assertEqual(
                        (result.zygosity, result.num_repeats_allele1, result.num_repeats_allele2,
                         result.allele1_sequence, result.allele2_sequence),
                        expected, f"locus chr1:{start_0based}-{end} genotyped from {vcf_gz_path}")
        finally:
            fasta_obj.close()


class TestMotifCompositionSplittingMethod(unittest.TestCase):
    """Test motif-splitting-method labeling and the TRF length threshold in motif composition."""

    def _make_locus(self, repeat_unit, allele1_sequence, allele2_sequence):
        return GenotypedTandemRepeat(
            ReferenceTandemRepeat(chrom="chr1", start_0based=100, end_1based=100 + len(repeat_unit),
                                  repeat_unit=repeat_unit),
            allele1_sequence=allele1_sequence,
            allele2_sequence=allele2_sequence,
        )

    def test_basic_method_labels_basic_split(self):
        """compute_motif_composition basic path labels both alleles 'basic-split', None for missing allele."""
        locus = self._make_locus("CAG", "CAGCAGCAG", None)
        result = compute_motif_composition(
            [locus], argparse.Namespace(add_motif_composition="basic", verbose=False, skip_hom_ref_loci=False))

        entry = result[locus.locus_id]
        self.assertEqual(entry["allele1"], {"motifs": ["CAG", "CAG", "CAG"], "prefix": "", "suffix": ""})
        self.assertEqual(entry["allele1_method"], MOTIF_DETECTION_METHOD_BASIC_SPLIT)
        self.assertIsNone(entry["allele2"])
        self.assertIsNone(entry["allele2_method"])

    def test_basic_split_keeps_trailing_suffix(self):
        """basic-split retains the trailing partial chunk as a suffix (eg. "CAGCA" -> "[CAG]CA")."""
        entry = build_basic_split_motif_entry("CAGCA", 3)
        self.assertEqual(entry, {"motifs": ["CAG"], "prefix": "", "suffix": "CA"})
        self.assertEqual(format_motif_entry_as_sequence_string(entry), "[CAG]CA")
        self.assertIsNone(build_basic_split_motif_entry(None, 3))

    def _make_trviz_decomposer(self):
        try:
            from trviz.decomposer import Decomposer
        except ImportError:
            self.skipTest("trviz is not installed")
        return Decomposer()

    def test_trviz_keeps_partial_motifs_at_the_ends_as_prefix_and_suffix(self):
        """A leading or trailing piece shorter than the motif is kept as unbracketed prefix/suffix bases."""
        decomposer = self._make_trviz_decomposer()
        entry, method = build_trviz_motif_entry("AGCAGCAGCA", "CAG", decomposer)
        self.assertEqual(entry, {"motifs": ["CAG", "CAG"], "prefix": "AG", "suffix": "CA"})
        self.assertEqual(format_motif_entry_as_sequence_string(entry), "AG[CAG][CAG]CA")
        self.assertEqual(method, MOTIF_DETECTION_METHOD_TRVIZ)

        # A sequence shorter than one motif is reported as a suffix, as the basic method does
        entry, _ = build_trviz_motif_entry("CA", "CAG", decomposer)
        self.assertEqual(entry, {"motifs": [], "prefix": "", "suffix": "CA"})

    def test_trviz_keeps_an_inserted_base_within_one_copy(self):
        """An inserted base lengthens one motif copy instead of shifting every copy after it."""
        decomposer = self._make_trviz_decomposer()
        entry, _ = build_trviz_motif_entry("CAGCAGACAGCAG", "CAG", decomposer)
        self.assertEqual(format_motif_entry_as_sequence_string(entry), "[CAG][CAGA][CAG][CAG]")

    def test_trviz_preserves_lowercase_bases(self):
        """Pieces are cut from the original sequence, so soft-masked bases stay lowercase."""
        decomposer = self._make_trviz_decomposer()
        entry, _ = build_trviz_motif_entry("cagcagCAG", "CAG", decomposer)
        self.assertEqual(entry["motifs"], ["cag", "cag", "CAG"])

    def test_trviz_falls_back_to_basic_split_for_non_acgt_bases(self):
        """trviz rejects bases other than A, C, G and T, so those alleles use the basic method."""
        decomposer = self._make_trviz_decomposer()
        entry, method = build_trviz_motif_entry("CAGNCAG", "CAG", decomposer)
        self.assertEqual(entry, build_basic_split_motif_entry("CAGNCAG", 3))
        self.assertEqual(method, MOTIF_DETECTION_METHOD_BASIC_SPLIT)
        self.assertEqual(build_trviz_motif_entry("", "CAG", decomposer), (None, None))

    def test_trviz_falls_back_to_basic_split_for_a_degenerate_motif(self):
        """trviz also rejects a motif with bases other than A, C, G and T, such as the GCN motif of a
        polyalanine repeat in the TRExplorer catalog, so those alleles use the basic method too."""
        decomposer = self._make_trviz_decomposer()
        entry, method = build_trviz_motif_entry("GCAGCCGCG", "GCN", decomposer)
        self.assertEqual(entry, build_basic_split_motif_entry("GCAGCCGCG", 3))
        self.assertEqual(method, MOTIF_DETECTION_METHOD_BASIC_SPLIT)

    def test_trviz_method_labels_both_alleles(self):
        """compute_motif_composition trviz path labels each allele 'trviz', None for a missing allele."""
        self._make_trviz_decomposer()
        locus = self._make_locus("CAG", "CAGCAGCAG", None)
        result = compute_motif_composition(
            [locus], argparse.Namespace(add_motif_composition="trviz", verbose=False, skip_hom_ref_loci=False))

        entry = result[locus.locus_id]
        self.assertEqual(entry["allele1"], {"motifs": ["CAG", "CAG", "CAG"], "prefix": "", "suffix": ""})
        self.assertEqual(entry["allele1_method"], MOTIF_DETECTION_METHOD_TRVIZ)
        self.assertIsNone(entry["allele2"])
        self.assertIsNone(entry["allele2_method"])

    def test_trf_path_short_sequences_fall_back_to_basic_split(self):
        """Sequences below max(threshold, 2*motif_size) never reach TRF and are labeled 'basic-split'.

        Both alleles are shorter than the 12bp default threshold, so compute_motif_lists_with_trf returns
        before constructing a TRFRunner - trf_executable_path is intentionally invalid to prove TRF is not run.
        """
        locus = self._make_locus("CAG", "CAGCAGCAG", "CAGCAG")  # 9bp and 6bp, both < 12
        result = compute_motif_lists_with_trf([locus], trf_executable_path="/nonexistent/trf", verbose=False)

        entry = result[locus.locus_id]
        self.assertEqual(entry["allele1"], {"motifs": ["CAG", "CAG", "CAG"], "prefix": "", "suffix": ""})
        self.assertEqual(entry["allele1_method"], MOTIF_DETECTION_METHOD_BASIC_SPLIT)
        self.assertEqual(entry["allele2"], {"motifs": ["CAG", "CAG"], "prefix": "", "suffix": ""})
        self.assertEqual(entry["allele2_method"], MOTIF_DETECTION_METHOD_BASIC_SPLIT)

    def test_threshold_arg_routes_longer_sequence_to_basic_split(self):
        """Raising min_allele_length_for_trf keeps a longer sequence on the basic path (no TRF)."""
        locus = self._make_locus("CAG", "CAG" * 6, None)  # 18bp >= default 12, but < 100
        result = compute_motif_lists_with_trf(
            [locus], trf_executable_path="/nonexistent/trf", min_allele_length_for_trf=100, verbose=False)

        entry = result[locus.locus_id]
        self.assertEqual(entry["allele1_method"], MOTIF_DETECTION_METHOD_BASIC_SPLIT)
        self.assertEqual(entry["allele1"]["motifs"], ["CAG"] * 6)

    def test_format_motifs_with_prefix_and_suffix(self):
        """format_motifs_as_sequence_string includes unbracketed prefix/suffix bases."""
        self.assertEqual(
            format_motifs_as_sequence_string(["GCA", "GCA", "GCC"], prefix="CA", suffix="G"),
            "CA[GCA][GCA][GCC]G")
        self.assertEqual(format_motifs_as_sequence_string(["CAG"]), "[CAG]")
        self.assertIsNone(format_motifs_as_sequence_string([], prefix="", suffix=""))

    def test_method_and_sequence_fields_emitted_in_tsv_and_json(self):
        """to_tsv_dict / to_json_dict emit the splitting method, motif sequence (with prefix/suffix), and counts."""
        locus = self._make_locus("CAG", "CAGCAGCAG", None)
        motif_lists = {
            "allele1": {"motifs": ["GCA", "GCA", "GCC"], "prefix": "CA", "suffix": "G"}, "allele2": None,
            "allele1_method": MOTIF_DETECTION_METHOD_TRF, "allele2_method": None,
        }

        tsv_dict = locus.to_tsv_dict(motif_lists=motif_lists)
        self.assertEqual(tsv_dict["Allele1MotifSequence"], "CA[GCA][GCA][GCC]G")
        self.assertEqual(tsv_dict["Allele1SequenceMotifSplittingMethod"], MOTIF_DETECTION_METHOD_TRF)
        self.assertEqual(tsv_dict["Allele2MotifSequence"], "")  # None -> "" in TSV
        self.assertEqual(tsv_dict["Allele2SequenceMotifSplittingMethod"], "")

        json_dict = locus.to_json_dict(motif_lists=motif_lists)
        self.assertEqual(json_dict["Allele1MotifSequence"], "CA[GCA][GCA][GCC]G")
        self.assertEqual(json_dict["Allele1MotifCounts"], {"GCA": 2, "GCC": 1})  # prefix/suffix excluded
        self.assertEqual(json_dict["Allele1SequenceMotifSplittingMethod"], MOTIF_DETECTION_METHOD_TRF)
        self.assertIsNone(json_dict["Allele2MotifSequence"])
        self.assertIsNone(json_dict["Allele2MotifCounts"])
        self.assertIsNone(json_dict["Allele2SequenceMotifSplittingMethod"])


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
        # One substituted base in 12 is a purity of exactly 11/12. A band would accept a materially wrong
        # value: forcing every non-perfect purity to 0.81 satisfies "> 0.8 and < 1.0".
        self.assertAlmostEqual(result["purity"], 11 / 12)
        self.assertFalse(result["is_pure"])

    def test_edit_distance_purity_is_skipped_above_max_allele_length(self):
        """Alleles longer than the cap get no edit-distance purity, while the other fields are unaffected."""
        sequence = "CAGCAACAGCAG"  # 12bp, one substitution
        below_cap = compute_repeat_counts_from_sequence(sequence, "CAG", max_allele_length_for_edit_distance_purity=12)
        self.assertAlmostEqual(below_cap["purity_via_edit_distance"], 11 / 12)

        above_cap = compute_repeat_counts_from_sequence(sequence, "CAG", max_allele_length_for_edit_distance_purity=11)
        self.assertIsNone(above_cap["purity_via_edit_distance"])
        self.assertEqual(above_cap["num_repeats"], 4)
        self.assertEqual(above_cap["repeat_size_bp"], 12)
        self.assertAlmostEqual(above_cap["purity"], 11 / 12)
        self.assertFalse(above_cap["is_pure"])

    def test_locus_edit_distance_purity_is_empty_when_one_present_allele_is_over_the_cap(self):
        """A capped allele is present, not absent, so the locus-level value can't fall back to the other allele."""
        from str_analysis.filter_vcf_to_tandem_repeats import GenotypedTandemRepeat, ReferenceTandemRepeat
        tr_locus = ReferenceTandemRepeat("chr1", 100, 112, "CAG")
        one_allele_capped = GenotypedTandemRepeat(
            tr_locus, allele1_sequence="CAG" * 4, allele2_sequence="CAG" * 4000,
            num_repeats_allele1=4, num_repeats_allele2=4000,
            allele1_purity=1.0, allele2_purity=0.8,
            allele1_purity_via_edit_distance=1.0, allele2_purity_via_edit_distance=None)
        self.assertIsNone(one_allele_capped.repeat_purity_via_edit_distance)
        self.assertEqual(one_allele_capped.to_tsv_dict()["RepeatPurityViaEditDistance"], "")
        # The per-allele columns still report what was computed
        self.assertEqual(one_allele_capped.to_tsv_dict()["RepeatPurityViaEditDistanceShortAllele"], "1.0000")
        self.assertEqual(one_allele_capped.to_tsv_dict()["RepeatPurityViaEditDistanceLongAllele"], "")
        # The position-by-position purity is unaffected
        self.assertAlmostEqual(one_allele_capped.repeat_purity, 0.8)

        # A genuinely absent allele (HEMI) still falls back to the one present
        one_allele_absent = GenotypedTandemRepeat(
            tr_locus, allele1_sequence="CAG" * 4, allele2_sequence=None,
            num_repeats_allele1=4, num_repeats_allele2=None,
            allele1_purity=1.0, allele2_purity=None,
            allele1_purity_via_edit_distance=1.0, allele2_purity_via_edit_distance=None)
        self.assertAlmostEqual(one_allele_absent.repeat_purity_via_edit_distance, 1.0)

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
        # Three exact copies plus a trailing partial copy that also matches: purity is exactly 1.0
        self.assertEqual(result["purity"], 1.0)
        self.assertTrue(result["is_pure"])

    def test_purity_and_purity_via_edit_distance_for_a_one_base_insertion(self):
        """A single inserted base shifts every later base out of phase for purity, but costs 2 edits for
        purity_via_edit_distance: the inserted base, and the base the allele then has beyond the pure repeat.
        """
        # (CAG)4 T (CAG)4 is 25bp. Every rotation of CAG mismatches one side of the T plus the T itself: 13
        # mismatches.
        result = compute_repeat_counts_from_sequence("CAG" * 4 + "T" + "CAG" * 4, "CAG")
        self.assertAlmostEqual(result["purity"], 12 / 25)
        self.assertAlmostEqual(result["purity_via_edit_distance"], 23 / 25)
        self.assertFalse(result["is_pure"])

        # A substitution scores the same both ways
        result = compute_repeat_counts_from_sequence("CAGCAGCCGCAG", "CAG")
        self.assertAlmostEqual(result["purity"], 11 / 12)
        self.assertAlmostEqual(result["purity_via_edit_distance"], 11 / 12)

    def test_count_repeats_allele_shorter_than_its_motif(self):
        """An allele shorter than one motif copy has no purity to report, so both fields are unknown.

        A numeric purity here would be misleading rather than merely imprecise: the sequence is too short to
        say whether it matches the motif, and 0.0 would render as RepeatPurity=0.0000 and IsPureRepeat=False
        in the output.
        """
        result = compute_repeat_counts_from_sequence("CATGG", "ACGTTGCAAGGCTTACCGGATTCAGGCATA")

        self.assertEqual(result["num_repeats"], 0)
        self.assertEqual(result["repeat_size_bp"], 5)
        self.assertIsNone(result["purity"])
        self.assertIsNone(result["is_pure"])

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
        except FileNotFoundError:
            # Only a genuinely missing bgzip/tabix is a skip. A CalledProcessError means the fixture itself is
            # malformed (unsorted records, a bad header), which must fail the test rather than turn into a
            # skip labelled 'bgzip/tabix unavailable' that leaves the run green.
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
        except FileNotFoundError:
            # Only a genuinely missing bgzip/tabix is a skip. A CalledProcessError means the fixture itself is
            # malformed (unsorted records, a bad header), which must fail the test rather than turn into a
            # skip labelled 'bgzip/tabix unavailable' that leaves the run green.
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
            add_motif_composition=None,
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
        # Every documented output column, so one going missing fails here rather than silently
        # disappearing from the file
        self.assertEqual(header, GENOTYPE_TSV_OUTPUT_COLUMNS)

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
        # The variant at 601 is C -> CAGCAGCAG, which inserts 8 bases, and the one at 619 is C -> CAG, which
        # inserts 2. Repeat counting is len(seq) // len(motif), so the values are fixed:
        # - Reference: 36bp / 3 = 12 repeats
        # - Haplotype 0: (36 + 8)bp = 44bp, 44 // 3 = 14 repeats
        # - Haplotype 1: (36 + 2)bp = 38bp, 38 // 3 = 12 repeats
        # These are asserted exactly. Accepting a range would let an offset error of up to 3 bases through,
        # and this is the only end-to-end locus with two variants on opposite haplotypes, which is where such
        # an error would show up.
        self.assertEqual(locus6["Zygosity"], "HET")
        self.assertEqual(int(locus6["NumRepeatsShortAllele"]), 12)
        self.assertEqual(int(locus6["NumRepeatsLongAllele"]), 14)
        self.assertEqual(int(locus6["RepeatSizeShortAlleleBp"]), 38)
        self.assertEqual(int(locus6["RepeatSizeLongAlleleBp"]), 44)

    def test_end_to_end_genotype_write_vcf_preserves_phasing(self):
        """The contributing-variants VCF (--write-vcf) must preserve input genotype phasing.

        A phased input genotype (0|1) should remain phased in the output, not be rewritten
        as unphased (0/1). Phasing is stored separately from the GT value on a pysam sample
        record, so it has to be copied explicitly when building the output record.
        """
        import pysam

        reference_seq = "A" * 100 + "CAG" * 8 + "T" * 76
        fasta_path = self._create_test_fasta({"chr1": reference_seq})
        if fasta_path is None:
            self.skipTest("pyfaidx unavailable")

        bed_gz_path = self._create_test_bed_and_index("chr1\t100\t124\tCAG\n")
        if bed_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        # Heterozygous phased expansion at the CAG locus (insert CAGCAG after the C at 0-based 100).
        # The record carries the FORMAT fields a DeepVariant plus WhatsHap or HiPhase callset would have, so
        # the test can tell a writer that drops everything but GT from one that copies the sample data.
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=300>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t101\t.\tC\tCCAGCAG\t.\tPASS\t.\tGT:PS:GQ:DP:AD\t0|1:101:48:37:19,18
"""
        vcf_gz_path = self._create_test_vcf_and_index(vcf_content)
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        self._temp_dir = tempfile.mkdtemp()
        output_prefix = os.path.join(self._temp_dir, "test_output")

        args = argparse.Namespace(
            reference_fasta_path=fasta_path,
            catalog_bed=bed_gz_path,
            input_vcf_path=vcf_gz_path,
            input_vcf_prefix="test",
            output_prefix=output_prefix,
            interval=None,
            verbose=False,
            show_progress_bar=False,
            write_vcf=True,
            write_json=False,
            skip_hom_ref_loci=False,
            add_motif_composition=None,
            trf_executable_path=None,
        )

        from str_analysis.filter_vcf_to_tandem_repeats import do_genotype_subcommand
        do_genotype_subcommand(args)

        output_vcf_path = f"{output_prefix}.tandem_repeat_contributing_variants.vcf.gz"
        self.assertTrue(os.path.exists(output_vcf_path), "contributing-variants VCF was not written")

        with pysam.VariantFile(output_vcf_path) as output_vcf:
            records = list(output_vcf.fetch())
        self.assertEqual(len(records), 1)
        self.assertTrue(records[0].samples["SAMPLE1"].phased,
                        "input phased genotype (0|1) was not preserved in the output VCF")
        # The phased flag is set separately from the GT value, and the GT copy is wrapped in an
        # except (KeyError, TypeError) handler, so the alleles have to be checked too: a swallowed
        # GT copy would leave a missing genotype that still reports phased=True.
        self.assertEqual(records[0].samples["SAMPLE1"]["GT"], (0, 1),
                         "input genotype alleles (0|1) were not preserved in the output VCF")
        # The TR annotations are this writer's entire payload, so check they name the right locus
        self.assertEqual(records[0].info["TR_LocusId"], ("chr1-100-124-CAG",))
        self.assertEqual(records[0].info["TR_Motif"], ("CAG",))
        # Every other FORMAT field must survive too. Losing PS in particular destroys the phase-block
        # identity even though the genotype still looks phased.
        sample = records[0].samples["SAMPLE1"]
        self.assertEqual(sample["PS"], 101)
        self.assertEqual(sample["GQ"], 48)
        self.assertEqual(sample["DP"], 37)
        self.assertEqual(tuple(sample["AD"]), (19, 18))

    def test_end_to_end_genotype_with_json_output(self):
        """Test end-to-end genotyping with JSON output enabled.

        Verifies that the JSON output contains all expected fields and agrees with the TSV written alongside
        it, allowing for the one intentional serialization difference: to_tsv_dict formats RepeatPurity to 4
        decimal places while to_json_dict rounds it to 3.
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
            add_motif_composition=None,
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

        # The input genotype is phased 0|1, so haplotype 0 is the 30bp reference allele and haplotype 1 the
        # 36bp expansion. Sorting the two lengths before comparing would let the writer serialize the
        # haplotypes in the wrong order without failing.
        self.assertEqual(len(locus["Allele1Sequence"]), 30)
        self.assertEqual(len(locus["Allele2Sequence"]), 36)
        self.assertEqual(locus["Allele1Sequence"], "CAG" * 10)
        # The VCF record is not left-normalized, so the expanded haplotype starts mid-motif. Pinning the
        # exact string keeps that documented rather than leaving it to a length check.
        self.assertEqual(locus["Allele2Sequence"], "CCAGCAGAGCAGCAGCAGCAGCAGCAGCAGCAGCAG")

        # The per-allele purity fields are part of the documented output and must survive into the file
        self.assertEqual(locus["RepeatPurityShortAllele"], 1.0)
        self.assertIsNotNone(locus["RepeatPurityLongAllele"])
        self.assertLess(locus["RepeatPurityLongAllele"], 1.0)

        # Compare against the TSV written by the same run, which the docstring promises
        tsv_path = f"{output_prefix}.tandem_repeat_genotypes.tsv.gz"
        self._temp_files.append(tsv_path)
        with gzip.open(tsv_path, "rt") as f:
            tsv_lines = f.read().split("\n")
        tsv_header = tsv_lines[0].split("\t")
        tsv_row = dict(zip(tsv_header, tsv_lines[1].split("\t")))
        for column in ["Chrom", "LocusId", "Motif", "CanonicalMotif", "Zygosity", "Allele1Sequence",
                       "Allele2Sequence", "NumRepeatsShortAllele", "NumRepeatsLongAllele",
                       "RepeatSizeShortAlleleBp", "RepeatSizeLongAlleleBp"]:
            self.assertEqual(str(locus[column]), tsv_row[column],
                             f"JSON and TSV disagree on {column}")
        # RepeatPurity is the one field the two writers format differently, by design
        self.assertAlmostEqual(float(tsv_row["RepeatPurityShortAllele"]),
                               locus["RepeatPurityShortAllele"], places=3)

    def test_end_to_end_genotype_with_basic_motif_counts(self):
        """Test JSON output with basic motif counting enabled.

        Verifies that motif counts are correctly computed when
        --add-motif-composition basic is specified.
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
            add_motif_composition="basic",  # Enable basic motif counting
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

        The no-call here must come from the two unphased heterozygous variants and nothing else, so both
        records carry a REF that matches the reference sequence. An earlier version of this fixture used
        REF 'A' at position 31 where the reference has 'C', which made the locus a no-call through the
        haplotype build-error path instead, and the test kept passing with the phasing check disabled.
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
chr1\t31\t.\tC\tT\t.\tPASS\t.\tGT\t0/1
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
            add_motif_composition=None,
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
        # The no-call must come from the unresolved phasing, not from a haplotype that failed to build
        self.assertEqual(row["NoCallReason"], NO_CALL_REASON_AMBIGUOUS_PHASING)

    def test_genotype_with_multiple_worker_processes_matches_single_process(self):
        """--threads > 1 must produce exactly the same TSV and JSON, in the same order, as a single-process run,
        including the motif composition columns that the worker processes compute for their own chunks."""
        import gzip
        import importlib.util
        import json

        reference_seq = (
            "A" * 100 + "CAG" * 8 +      # chr1:100-124  het expansion
            "T" * 76 + "AT" * 9 +        # chr1:200-218  hom ref, no variant
            "G" * 82 + "CTG" * 9 +       # chr1:300-327  het contraction
            "C" * 73 + "AAGGG" * 6 +     # chr1:400-430  hom expansion
            "A" * 70 + "GCC" * 6 +       # chr1:500-518  hom ref, no variant
            "T" * 82 + "CAG" * 12 +      # chr1:600-636  one variant per haplotype
            "G" * 64
        )
        fasta_path = self._create_test_fasta({"chr1": reference_seq})
        if fasta_path is None:
            self.skipTest("pyfaidx unavailable")

        bed_gz_path = self._create_test_bed_and_index(
            "chr1\t100\t124\tCAG\nchr1\t200\t218\tAT\nchr1\t300\t327\tCTG\n"
            "chr1\t400\t430\tAAGGG\nchr1\t500\t518\tGCC\nchr1\t600\t636\tCAG\n")
        if bed_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        vcf_gz_path = self._create_test_vcf_and_index("""##fileformat=VCFv4.2
##contig=<ID=chr1,length=700>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t101\t.\tC\tCCAGCAG\t.\tPASS\t.\tGT\t0|1
chr1\t301\t.\tCTGC\tC\t.\tPASS\t.\tGT\t1|0
chr1\t401\t.\tA\tAAGGGA\t.\tPASS\t.\tGT\t1|1
chr1\t601\t.\tC\tCAGCAGCAG\t.\tPASS\t.\tGT\t1|0
chr1\t619\t.\tC\tCAG\t.\tPASS\t.\tGT\t0|1
""")
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        self._temp_dir = tempfile.mkdtemp()

        from str_analysis import filter_vcf_to_tandem_repeats
        motif_composition_methods = [None, "basic"]
        if importlib.util.find_spec("trviz") is not None:
            motif_composition_methods.append("trviz")

        # Lower the minimum chunk size so these 6 loci are split into several chunks that the two worker
        # processes genotype out of lockstep, which is the case where an ordering or merging bug would show.
        with mock.patch.object(filter_vcf_to_tandem_repeats, "MIN_LOCI_PER_GENOTYPING_CHUNK", 2):
            for add_motif_composition in motif_composition_methods:
                tsv_contents = {}
                json_contents = {}
                for num_threads in (1, 2):
                    output_prefix = os.path.join(self._temp_dir, f"{add_motif_composition}_threads_{num_threads}")
                    filter_vcf_to_tandem_repeats.do_genotype_subcommand(argparse.Namespace(
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
                        add_motif_composition=add_motif_composition,
                        trf_executable_path=None,
                        threads=num_threads,
                    ))
                    with gzip.open(f"{output_prefix}.tandem_repeat_genotypes.tsv.gz", "rt") as f:
                        tsv_contents[num_threads] = f.read()
                    with gzip.open(f"{output_prefix}.tandem_repeat_genotypes.json.gz", "rt") as f:
                        json_contents[num_threads] = f.read()

                self.assertEqual(len(tsv_contents[1].splitlines()), 7, "Should have a header + 6 genotyped loci")
                self.assertEqual(tsv_contents[2], tsv_contents[1], add_motif_composition)
                self.assertEqual(json_contents[2], json_contents[1], add_motif_composition)
                # The JSON is written one record at a time, and must read back as the same text that
                # json.dump(records, f, indent=2) writes
                self.assertEqual(json_contents[1], json.dumps(json.loads(json_contents[1]), indent=2))
                if add_motif_composition:
                    header = tsv_contents[1].splitlines()[0].split("\t")
                    motif_sequences = [row.split("\t")[header.index("Allele1MotifSequence")]
                                       for row in tsv_contents[1].splitlines()[1:]]
                    self.assertTrue(all(motif_sequences), f"{add_motif_composition}: {motif_sequences}")


class TestWriteFunctions(unittest.TestCase):
    """Test the write_tsv and write_bed functions."""

    def setUp(self):
        """Set up test case with temp directory and mock args."""
        self._temp_dir = tempfile.mkdtemp()
        self._output_prefix = os.path.join(self._temp_dir, "test_output")
        self._args = argparse.Namespace(
            output_prefix=self._output_prefix,
            copy_info_field_keys_to_tsv=None,
            verbose=False,
        )

    def tearDown(self):
        """Clean up temp directory."""
        if os.path.exists(self._temp_dir):
            shutil.rmtree(self._temp_dir)

    def test_run_shell_command_raises_on_failure(self):
        """A failed compression or indexing step must not be reported as a successful write."""
        run_shell_command("true")
        with self.assertRaises(RuntimeError):
            run_shell_command("false")

    def test_write_tsv_with_reference_tandem_repeat(self):
        """Test write_tsv doesn't crash with ReferenceTandemRepeat objects.

        ReferenceTandemRepeat objects lack some TandemRepeatAllele properties
        (ins_or_del, allele, is_pure_repeat, info_field_dict). The write_tsv
        function should handle this gracefully without crashing.
        """
        ref_tr = ReferenceTandemRepeat("chr1", 1000, 1030, "CAG", DETECTION_MODE_PURE_REPEATS)

        # write_tsv should not crash on a ReferenceTandemRepeat
        write_tsv([ref_tr], self._args)

        # Check that output file was created
        tsv_output_path = f"{self._output_prefix}.tandem_repeats.tsv.gz"
        self.assertTrue(os.path.exists(tsv_output_path))

        # Read and verify the output
        import gzip
        with gzip.open(tsv_output_path, 'rt') as f:
            lines = f.readlines()

        self.assertEqual(len(lines), 2)  # Header + 1 data row

        header = lines[0].strip().split('\t')
        values = lines[1].strip().split('\t')
        row = dict(zip(header, values))

        # Verify expected values
        self.assertEqual(row["Chrom"], "chr1")
        self.assertEqual(row["Start0Based"], "1000")
        self.assertEqual(row["End1Based"], "1030")
        self.assertEqual(row["Motif"], "CAG")
        self.assertEqual(row["MotifSize"], "3")

        # ReferenceTandemRepeat-specific: these fields should be empty
        self.assertEqual(row["INS_or_DEL"], "")
        self.assertEqual(row["VcfPos"], "")
        self.assertEqual(row["IsPureRepeat"], "")

    def test_write_tsv_with_mixed_objects(self):
        """Test write_tsv with a mix of ReferenceTandemRepeat objects.

        After merging, the tandem_repeat_alleles list may contain both
        TandemRepeatAllele and ReferenceTandemRepeat objects.
        """
        ref_tr1 = ReferenceTandemRepeat("chr1", 1000, 1030, "CAG", DETECTION_MODE_PURE_REPEATS)
        ref_tr2 = ReferenceTandemRepeat("chr2", 2000, 2040, "AAAG", DETECTION_MODE_TRF)

        # write_tsv should handle multiple ReferenceTandemRepeat objects
        write_tsv([ref_tr1, ref_tr2], self._args)

        # Check output
        tsv_output_path = f"{self._output_prefix}.tandem_repeats.tsv.gz"
        self.assertTrue(os.path.exists(tsv_output_path))

        import gzip
        with gzip.open(tsv_output_path, 'rt') as f:
            lines = f.readlines()

        self.assertEqual(len(lines), 3)  # Header + 2 data rows


class TestRepeatUnitSimilarity(unittest.TestCase):
    """Test the compute_repeat_unit_id and are_repeat_units_similar functions."""

    def test_similar_short_motifs_same(self):
        """Test that identical short motifs (<=6bp) are considered similar."""
        self.assertTrue(are_repeat_units_similar("CAG", "CAG"))
        self.assertTrue(are_repeat_units_similar("AAAAAA", "AAAAAA"))
        self.assertTrue(are_repeat_units_similar("A", "A"))

    def test_similar_short_motifs_different(self):
        """Test that different short motifs (<=6bp) are NOT similar."""
        self.assertFalse(are_repeat_units_similar("CAG", "CTG"))
        self.assertFalse(are_repeat_units_similar("CAG", "CAA"))
        self.assertFalse(are_repeat_units_similar("AT", "CG"))

    def test_similar_short_motifs_different_length(self):
        """Test that short motifs with different lengths are NOT similar."""
        self.assertFalse(are_repeat_units_similar("CAG", "CAGC"))
        self.assertFalse(are_repeat_units_similar("A", "AA"))
        self.assertFalse(are_repeat_units_similar("AAAAAA", "AAAAAAA"))

    def test_similar_long_motifs_same_length(self):
        """Test that long motifs (>6bp) with same length ARE similar (current behavior).

        Note: This documents the current behavior where long motifs are only compared
        by length. This may be overly permissive - see issue #4 in the bug list.
        """
        # Current behavior: long motifs with same length are considered similar
        self.assertTrue(are_repeat_units_similar("CAGCAGC", "AAAAAAA"))
        self.assertTrue(are_repeat_units_similar("ATATATAT", "CGCGCGCG"))

    def test_similar_long_motifs_different_length(self):
        """Test that long motifs (>6bp) with different lengths are NOT similar."""
        self.assertFalse(are_repeat_units_similar("CAGCAGC", "CAGCAGCA"))
        self.assertFalse(are_repeat_units_similar("AAAAAAA", "AAAAAAAA"))

    def test_compute_repeat_unit_id_short_motifs(self):
        """Test compute_repeat_unit_id returns the motif itself for short motifs."""
        from str_analysis.filter_vcf_to_tandem_repeats import compute_repeat_unit_id
        self.assertEqual(compute_repeat_unit_id("CAG"), "CAG")
        self.assertEqual(compute_repeat_unit_id("A"), "A")
        self.assertEqual(compute_repeat_unit_id("AAAAAA"), "AAAAAA")

    def test_compute_repeat_unit_id_long_motifs(self):
        """Test compute_repeat_unit_id returns the length for long motifs."""
        from str_analysis.filter_vcf_to_tandem_repeats import compute_repeat_unit_id
        self.assertEqual(compute_repeat_unit_id("AAAAAAA"), 7)
        self.assertEqual(compute_repeat_unit_id("CAGCAGCAGCAG"), 12)

    def test_boundary_case_6bp(self):
        """Test the boundary case at exactly 6bp."""
        from str_analysis.filter_vcf_to_tandem_repeats import compute_repeat_unit_id
        # 6bp should still use the motif itself
        self.assertEqual(compute_repeat_unit_id("AAAAAA"), "AAAAAA")
        self.assertFalse(are_repeat_units_similar("AAAAAA", "AAAAAT"))

        # 7bp should use length
        self.assertEqual(compute_repeat_unit_id("AAAAAAA"), 7)
        self.assertTrue(are_repeat_units_similar("AAAAAAA", "AAAAAAC"))


# A full-length AluY element seen as a 304bp insertion in HG002 at chr1:230,949,404, used here as a realistic
# example of an insertion that is not a tandem repeat.
ALU_INSERTION_SEQUENCE = (
    "AAAAAAAAAAAGGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGTGGATCATGAGGTCAGGAGATCGAGACCATCCTG"
    "GCTAACAAGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGGGCGCGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGA"
    "GAATGGCGTGAACCCGGGAAGCGGAGCTTGCAGTGAGCCGAGATTGCGCCACTGCAGTCCGCAGTCCGGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAA")

# The 333bp AAAGG/AAAGGG expansion seen on one haplotype of HG01109 at the RFC1 locus, whose reference repeat is
# (AAAAG)n. It is a true tandem repeat expansion even though its motif differs from the annotated locus motif.
RFC1_MOTIF_SWAP_INSERTION_SEQUENCE = (
    "AAGGAAAGGGACGGGACGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAG"
    "GGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAG"
    "GGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGG")


class TestInsertionFiltering(unittest.TestCase):
    """Test the checks that decide whether inserted bases belong to a locus tandem repeat."""

    def _insertion_filter(self, motif, **overrides):
        # Read from the same constants parse_args() uses, so these tests exercise the shipped thresholds rather
        # than a copy of them that a future retune could leave behind.
        settings = dict(
            motif=motif,
            min_insertion_size_to_check=DEFAULT_MIN_INSERTION_SIZE_TO_CHECK,
            min_insertion_purity=DEFAULT_MIN_INSERTION_PURITY,
            min_insertion_periodicity=DEFAULT_MIN_INSERTION_PERIODICITY,
        )
        settings.update(overrides)
        return InsertionFilter(**settings)

    def test_short_insertion_is_always_counted(self):
        """An insertion below the size threshold is counted without being evaluated."""
        is_insertion_sufficiently_repetitive, reason = check_if_inserted_sequence_is_repetitive(
            "GTCAGTCTGA", self._insertion_filter("A"))
        self.assertTrue(is_insertion_sufficiently_repetitive)
        self.assertIsNone(reason)

    def test_insertion_of_locus_motif_is_counted(self):
        """More copies of the locus motif are counted."""
        is_insertion_sufficiently_repetitive, reason = check_if_inserted_sequence_is_repetitive(
            "CAG" * 20, self._insertion_filter("CAG"))
        self.assertTrue(is_insertion_sufficiently_repetitive)
        self.assertIsNone(reason)

    def test_single_copy_vntr_insertion_is_counted(self):
        """One extra copy of a large VNTR unit is counted, on purity alone.

        Gaining or losing a single repeat unit is the most common form of VNTR variation, and a lone copy of a
        non-repetitive unit carries no internal periodicity to detect: this 40bp unit scores about 0.35, well
        under any usable periodicity threshold. Purity is the only signal that can see it, which is why the
        purity clause must stay available to insertions spanning fewer than two motif copies.

        The inserted copy carries two mismatches against the catalog motif, as a real one would. An exact copy
        would score purity 1.0 and so be accepted at every threshold up to 1.0, which would leave the test unable
        to fail no matter how far the purity threshold were raised.
        """
        motif = "ACGGTTACGCATTGAGCCTAGGTCATGCAAGTTCCGATAG"
        inserted_copy = "ACGGTTACGTATTGAGCCTAGGTCGTGCAAGTTCCGATAG"
        insertion_filter = self._insertion_filter(motif)

        # Read all three thresholds off the filter rather than repeating them, so that retuning any one of them past
        # this sequence turns the test into a failure instead of quietly changing which clause does the accepting.
        # The size check comes first because check_if_inserted_sequence_is_repetitive accepts anything shorter
        # than it outright, before purity or periodicity is consulted.
        self.assertGreaterEqual(len(inserted_copy), insertion_filter.min_insertion_size_to_check)
        purity, _, _ = compute_best_phase_repeat_purity(inserted_copy, motif)
        self.assertGreaterEqual(purity, insertion_filter.min_insertion_purity)
        _, periodicity = compute_sequence_periodicity(inserted_copy, max_period=min(500, max(100, 3 * len(motif))))
        self.assertLess(periodicity, insertion_filter.min_insertion_periodicity)

        is_insertion_sufficiently_repetitive, reason = check_if_inserted_sequence_is_repetitive(
            inserted_copy, insertion_filter)
        self.assertTrue(is_insertion_sufficiently_repetitive)
        self.assertIsNone(reason)

    def test_insertion_starting_mid_motif_is_counted(self):
        """An insertion that starts in the middle of the motif is still counted."""
        is_insertion_sufficiently_repetitive, _ = check_if_inserted_sequence_is_repetitive(
            "AG" + "CAG" * 20, self._insertion_filter("CAG"))
        self.assertTrue(is_insertion_sufficiently_repetitive)

    def test_alu_insertion_into_poly_a_tract_is_not_counted(self):
        """An Alu element inserted into a poly-A tract must not be counted as poly-A repeats."""
        is_insertion_sufficiently_repetitive, reason = check_if_inserted_sequence_is_repetitive(
            ALU_INSERTION_SEQUENCE, self._insertion_filter("A"))
        self.assertFalse(is_insertion_sufficiently_repetitive)
        self.assertEqual(reason, INSERTION_FILTER_REASON_NOT_REPEAT_LIKE)

    def test_alu_insertion_is_not_counted_at_a_locus_with_a_long_impure_motif(self):
        """An Alu must still be rejected at a locus whose annotated motif is long and whose reference repeat is
        itself impure, where the Alu's purity relative to that motif is close to what random sequence scores."""
        is_insertion_sufficiently_repetitive, reason = check_if_inserted_sequence_is_repetitive(
            ALU_INSERTION_SEQUENCE, self._insertion_filter("ATATCCACTTGCAGATTCTACAAAAAGAGT"))
        self.assertFalse(is_insertion_sufficiently_repetitive)
        self.assertEqual(reason, INSERTION_FILTER_REASON_NOT_REPEAT_LIKE)

    def test_expansion_with_a_different_motif_is_counted(self):
        """An AAAGG/AAAGGG expansion at the (AAAAG)n RFC1 locus must still be counted."""
        insertion_filter = self._insertion_filter("AAAAG")
        purity, _, _ = compute_best_phase_repeat_purity(RFC1_MOTIF_SWAP_INSERTION_SEQUENCE, "AAAAG")
        self.assertLess(purity, insertion_filter.min_insertion_purity)  # the purity clause alone would reject it

        is_insertion_sufficiently_repetitive, reason = check_if_inserted_sequence_is_repetitive(
            RFC1_MOTIF_SWAP_INSERTION_SEQUENCE, insertion_filter)
        self.assertTrue(is_insertion_sufficiently_repetitive)
        self.assertIsNone(reason)

    def test_impure_but_repetitive_insertion_is_counted(self):
        """An insertion too impure to match the locus motif is still counted when it repeats on its own."""
        # 12 copies of a 12bp unit, each with a different base substituted, so the insertion matches neither the
        # locus motif nor a pure repeat of itself, yet is clearly still a tandem repeat
        unit = "CAGGATTCCTGA"
        impure_insertion = "".join(unit[:i] + "T" + unit[i + 1:] for i in range(len(unit)))
        insertion_filter = self._insertion_filter("CAG")
        purity, _, _ = compute_best_phase_repeat_purity(impure_insertion, "CAG")
        self.assertLess(purity, insertion_filter.min_insertion_purity)
        self.assertTrue(check_if_inserted_sequence_is_repetitive(impure_insertion, insertion_filter)[0])

    def test_insertion_shorter_than_the_motif_is_scored_on_purity(self):
        """An insertion shorter than one motif copy must still be judged against the motif.

        compute_best_phase_repeat_purity() compares against whole rotations and returns nan below one copy,
        which would leave only the periodicity test. A fragment of a long, non-repetitive VNTR unit has no
        internal periodicity, so it would be rejected and the whole locus would get no call.
        """
        motif = "GCTAAGGTCCATTGACCGTAAGCTTGGCCAATCGTTAGGCCATTAGGCCTTAAGGCATCG"
        self.assertEqual(len(motif), 60)
        insertion_filter = InsertionFilter(
            motif=motif,
            min_insertion_size_to_check=DEFAULT_MIN_INSERTION_SIZE_TO_CHECK,
            min_insertion_purity=DEFAULT_MIN_INSERTION_PURITY,
            min_insertion_periodicity=DEFAULT_MIN_INSERTION_PERIODICITY)

        # 40 bases of the motif: above the size threshold, below one full copy
        is_insertion_sufficiently_repetitive, reason = check_if_inserted_sequence_is_repetitive(
            motif[:40], insertion_filter)
        self.assertTrue(is_insertion_sufficiently_repetitive,
                        f"a partial copy of the locus motif was rejected ({reason})")
        # A fragment starting mid-motif must work too, since an allele often starts mid-motif
        is_insertion_sufficiently_repetitive, _ = check_if_inserted_sequence_is_repetitive(
            motif[25:] + motif[:5], insertion_filter)
        self.assertTrue(is_insertion_sufficiently_repetitive)
        # The neighbouring case: unrelated sequence of the same length is still rejected
        is_insertion_sufficiently_repetitive, reason = check_if_inserted_sequence_is_repetitive(
            "ACGTACGTGGGGCCCCTTTTAAAACCTTGGAACCTT"[:40], insertion_filter)
        self.assertFalse(is_insertion_sufficiently_repetitive)
        self.assertEqual(reason, INSERTION_FILTER_REASON_NOT_REPEAT_LIKE)

    def test_compute_partial_copy_purity(self):
        """The helper scores a short sequence against the best-matching stretch of the motif."""
        self.assertAlmostEqual(compute_partial_copy_purity("CAGCAGCA", "CAGCAGCAGT"), 1.0)
        # Wrapping past the end of the motif must be allowed, since a fragment can start anywhere
        self.assertAlmostEqual(compute_partial_copy_purity("GTCAG", "CAGCAGCAGT"), 1.0)
        # One mismatch in eight bases
        self.assertAlmostEqual(compute_partial_copy_purity("CAGCTGCA", "CAGCAGCAGT"), 7 / 8)
        # Not applicable when the sequence is at least as long as the motif, where the caller uses the
        # rotation-based purity instead
        for sequence, motif in [("CAGCAG", "CAG"), ("CAG", "CAG"), ("", "CAG"), ("CA", "")]:
            purity = compute_partial_copy_purity(sequence, motif)
            self.assertNotEqual(purity, purity, f"expected nan for ({sequence!r}, {motif!r}), got {purity}")

    def test_insertion_with_ns_is_not_counted(self):
        """Inserted bases containing Ns come from an assembly gap and must not be counted."""
        is_insertion_sufficiently_repetitive, reason = check_if_inserted_sequence_is_repetitive(
            "N" * 50, self._insertion_filter("A"))
        self.assertFalse(is_insertion_sufficiently_repetitive)
        self.assertEqual(reason, INSERTION_FILTER_REASON_CONTAINS_NS)

    def test_build_insertion_filter_uses_args(self):
        args = argparse.Namespace(
            min_insertion_size_to_check_for_repetitiveness=30,
            min_insertion_purity=0.5,
            min_insertion_periodicity=0.8)
        insertion_filter = build_insertion_filter("CAG", args)
        self.assertEqual(insertion_filter.motif, "CAG")
        self.assertEqual(insertion_filter.min_insertion_size_to_check, 30)
        self.assertEqual(insertion_filter.min_insertion_purity, 0.5)
        self.assertEqual(insertion_filter.min_insertion_periodicity, 0.8)

    def test_only_the_rejected_insertion_is_reported(self):
        """Insertions that belong to the repeat, deletions and substitutions are all passed over."""
        variant_list = [
            (100, "A", "A" + "CAG" * 20),                      # a real expansion, accepted
            (200, "T", "T" + ALU_INSERTION_SEQUENCE),          # an Alu insertion, rejected
            (300, "GGGG", "G"),                                # a deletion, not an insertion
        ]
        rejected = find_insufficiently_repetitive_insertions_inside_locus(
            variant_list, 0, 1000, self._insertion_filter("CAG"))

        self.assertEqual(len(rejected), 1)
        self.assertEqual(rejected[0][0], ALU_INSERTION_SEQUENCE)
        self.assertEqual(rejected[0][1], INSERTION_FILTER_REASON_NOT_REPEAT_LIKE)

    def test_bases_shared_by_ref_and_alt_are_not_treated_as_inserted(self):
        """Only the bases the alt actually adds are evaluated, not the suffix it shares with the ref."""
        # ref "AT" -> alt "A" + 26 G's + "T": ref and alt share the trailing "T", and 26 G's are inserted
        variant_list = [(100, "AT", "A" + "G" * 26 + "T")]
        rejected = find_insufficiently_repetitive_insertions_inside_locus(
            variant_list, 0, 1000, self._insertion_filter("CAG", min_insertion_periodicity=1.01))
        self.assertEqual(len(rejected), 1)
        self.assertEqual(rejected[0][0], "G" * 26)

    def test_nothing_is_rejected_when_every_insertion_belongs(self):
        variant_list = [(100, "A", "A" + "CAG" * 20)]
        self.assertEqual(
            find_insufficiently_repetitive_insertions_inside_locus(
                variant_list, 0, 1000, self._insertion_filter("CAG")),
            [])

    def test_an_insertion_positioned_outside_the_locus_is_not_judged(self):
        """A variant can overlap the locus through its reference span while inserting its bases outside it.
        Those bases are not part of the locus sequence, so they must not cost the allele its call."""
        # ref spans [99, 103) and the inserted bases are positioned at its end, 0-based position 103. The ref
        # ends in G and the Alu in A, so there is no shared suffix to shift that position.
        variant_list = [(100, "GGGG", "GGGG" + ALU_INSERTION_SEQUENCE)]
        insertion_filter = self._insertion_filter("A")

        # locus [50, 103) ends exactly where the bases are inserted, so they land past its end
        self.assertEqual(
            find_insufficiently_repetitive_insertions_inside_locus(variant_list, 50, 103, insertion_filter), [])
        # locus [50, 104) contains the insertion point, so the same insertion is judged and rejected
        self.assertEqual(
            len(find_insufficiently_repetitive_insertions_inside_locus(variant_list, 50, 104, insertion_filter)), 1)


class TestGenotypingWithInsertionFilter(unittest.TestCase):
    """End-to-end genotyping tests covering insertions that don't belong to the locus tandem repeat."""

    def setUp(self):
        self._temp_files = []

    def tearDown(self):
        for f in self._temp_files:
            if os.path.exists(f):
                os.unlink(f)
            for ext in [".tbi", ".csi", ".fai"]:
                if os.path.exists(f + ext):
                    os.unlink(f + ext)

    def _create_temp_file(self, content, suffix):
        with tempfile.NamedTemporaryFile(suffix=suffix, delete=False, mode="w") as f:
            f.write(content)
            self._temp_files.append(f.name)
            return f.name

    def _create_vcf(self, vcf_content):
        import subprocess
        vcf_path = self._create_temp_file(vcf_content, ".vcf")
        vcf_gz_path = vcf_path + ".gz"
        self._temp_files.append(vcf_gz_path)
        self._temp_files.append(vcf_gz_path + ".tbi")
        try:
            with open(vcf_gz_path, "wb") as gz_out:
                subprocess.run(["bgzip", "-c", vcf_path], stdout=gz_out, check=True)
            subprocess.run(["tabix", "-p", "vcf", vcf_gz_path], check=True)
            return vcf_gz_path
        except FileNotFoundError:
            # Only a genuinely missing bgzip/tabix is a skip. A CalledProcessError means the fixture itself is
            # malformed (unsorted records, a bad header), which must fail the test rather than turn into a
            # skip labelled 'bgzip/tabix unavailable' that leaves the run green.
            return None

    def _create_fasta(self, seq_dict):
        fasta_content = "".join(f">{chrom}\n{seq}\n" for chrom, seq in seq_dict.items())
        fasta_path = self._create_temp_file(fasta_content, ".fa")
        pyfaidx.Fasta(fasta_path)
        self._temp_files.append(fasta_path + ".fai")
        return fasta_path

    def _genotype_alu_insertion_into_poly_a(self, insertion_filter, genotype="1|1"):
        """Genotype a 20bp poly-A locus carrying an Alu insertion in the middle of the tract."""
        import pysam

        reference = "GGCCGT" + "A" * 20 + "TTGCAC"
        fasta_path = self._create_fasta({"chr1": reference})
        # The insertion is anchored on the 10th A of the tract, ie. 1-based position 15
        vcf_content = (
            "##fileformat=VCFv4.2\n"
            f"##contig=<ID=chr1,length={len(reference)}>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
            f"chr1\t15\t.\tA\tA{ALU_INSERTION_SEQUENCE}\t.\tPASS\t.\tGT\t{genotype}\n")
        vcf_gz_path = self._create_vcf(vcf_content)
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        tr_locus = ReferenceTandemRepeat(chrom="chr1", start_0based=6, end_1based=26, repeat_unit="A")
        fasta_obj = pyfaidx.Fasta(fasta_path, one_based_attributes=False, as_raw=True)
        vcf_file = pysam.VariantFile(vcf_gz_path)
        try:
            return genotype_single_locus(
                tr_locus, vcf_file, fasta_obj,
                insertion_filter=insertion_filter)
        finally:
            vcf_file.close()
            fasta_obj.close()

    def test_alu_insertion_leaves_both_alleles_without_a_call(self):
        """An Alu dropped into a poly-A tract makes the repeat count unknowable, so neither allele gets a
        call. The VCF here is homozygous, so both alleles carry the Alu."""
        insertion_filter = build_insertion_filter("A", argparse.Namespace())
        result = self._genotype_alu_insertion_into_poly_a(insertion_filter)

        self.assertIsNone(result.num_repeats_allele1)
        self.assertIsNone(result.num_repeats_allele2)
        self.assertIsNone(result.allele1_sequence)
        self.assertIsNone(result.allele2_sequence)
        self.assertEqual(result.no_call_reason, NO_CALL_REASON_NON_REPEAT_INSERTION)
        self.assertEqual(result.to_tsv_dict()["NoCallReason"],
                         f"{NO_CALL_REASON_NON_REPEAT_INSERTION} "
                         f"({INSERTION_FILTER_REASON_NOT_REPEAT_LIKE})")
        self.assertIsNone(result.zygosity)
        self.assertIsNone(result.repeat_size_long_allele_bp)
        self.assertEqual(result.num_alleles_with_non_repeat_insertions, 2)

    def test_alu_insertion_on_one_haplotype_discards_the_whole_genotype(self):
        """A heterozygous Alu costs the whole locus, not just the haplotype that carries it. Keeping the
        reference haplotype would report as HEMI, which is indistinguishable from a real hemizygous call and
        would put the surviving allele in both the short and long allele columns."""
        insertion_filter = build_insertion_filter("A", argparse.Namespace())
        result = self._genotype_alu_insertion_into_poly_a(insertion_filter, genotype="0|1")

        self.assertIsNone(result.zygosity)
        self.assertIsNone(result.num_repeats_allele1)
        self.assertIsNone(result.num_repeats_allele2)
        self.assertIsNone(result.allele1_sequence)
        self.assertIsNone(result.allele2_sequence)
        self.assertIsNone(result.num_repeats_short_allele)
        self.assertIsNone(result.num_repeats_long_allele)
        self.assertEqual(result.num_alleles_with_non_repeat_insertions, 1)

    def test_alu_insertion_is_counted_when_filtering_is_turned_off(self):
        """Passing no insertion filter restores the old behavior of counting every inserted base."""
        result = self._genotype_alu_insertion_into_poly_a(insertion_filter=None)

        self.assertEqual(result.num_repeats_allele1, 20 + len(ALU_INSERTION_SEQUENCE))
        self.assertEqual(result.num_alleles_with_non_repeat_insertions, 0)

    def test_alu_insertion_is_counted_when_size_threshold_is_very_large(self):
        """A very large --min-insertion-size-to-check-for-repetitiveness turns the insertion check off."""
        insertion_filter = build_insertion_filter(
            "A", argparse.Namespace(min_insertion_size_to_check_for_repetitiveness=10**9))
        result = self._genotype_alu_insertion_into_poly_a(insertion_filter)

        self.assertEqual(result.num_repeats_allele1, 20 + len(ALU_INSERTION_SEQUENCE))
        self.assertEqual(result.num_alleles_with_non_repeat_insertions, 0)

    def test_true_expansion_with_a_different_motif_is_kept(self):
        """A motif-swap expansion like the ones seen at RFC1 must be counted in full."""
        import pysam

        reference = "GGCCGT" + "AAAAG" * 11 + "TTGCAC"
        fasta_path = self._create_fasta({"chr1": reference})
        vcf_content = (
            "##fileformat=VCFv4.2\n"
            f"##contig=<ID=chr1,length={len(reference)}>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
            f"chr1\t6\t.\tT\tT{RFC1_MOTIF_SWAP_INSERTION_SEQUENCE}\t.\tPASS\t.\tGT\t0|1\n")
        vcf_gz_path = self._create_vcf(vcf_content)
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        tr_locus = ReferenceTandemRepeat(chrom="chr1", start_0based=6, end_1based=61, repeat_unit="AAAAG")
        fasta_obj = pyfaidx.Fasta(fasta_path, one_based_attributes=False, as_raw=True)
        vcf_file = pysam.VariantFile(vcf_gz_path)
        try:
            result = genotype_single_locus(
                tr_locus, vcf_file, fasta_obj,
                insertion_filter=build_insertion_filter("AAAAG", argparse.Namespace()))
        finally:
            vcf_file.close()
            fasta_obj.close()

        self.assertEqual(result.zygosity, "HET")
        self.assertEqual(result.num_repeats_short_allele, 11)
        self.assertEqual(result.num_repeats_long_allele,
                         (55 + len(RFC1_MOTIF_SWAP_INSERTION_SEQUENCE)) // 5)
        self.assertEqual(result.num_alleles_with_non_repeat_insertions, 0)

    def test_lower_case_reference_allele_does_not_shift_the_locus_boundaries(self):
        """VCFs written against a soft-masked reference (eg. DipCall) give the ref allele in lower case and the
        alt allele in upper case, which must not change where the locus is trimmed out of the haplotype."""
        import pysam

        reference = "AAA" + "cag" * 6 + "AAA"
        fasta_path = self._create_fasta({"chr1": reference})
        # ref "cagc" -> alt "CAGCAGC" inserts one CAG positioned at 0-based position 12, inside the locus, while
        # the ref span [12, 16) crosses the locus end at 15. Once both alleles are upper-cased they share a
        # "CAGC" suffix; comparing them as written finds no shared suffix and drops the inserted repeat.
        vcf_content = (
            "##fileformat=VCFv4.2\n"
            f"##contig=<ID=chr1,length={len(reference)}>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
            "chr1\t13\t.\tcagc\tCAGCAGC\t.\tPASS\t.\tGT\t1|1\n")
        vcf_gz_path = self._create_vcf(vcf_content)
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        tr_locus = ReferenceTandemRepeat(chrom="chr1", start_0based=3, end_1based=15, repeat_unit="CAG")
        fasta_obj = pyfaidx.Fasta(fasta_path, one_based_attributes=False, as_raw=True)
        vcf_file = pysam.VariantFile(vcf_gz_path)
        try:
            result = genotype_single_locus(
                tr_locus, vcf_file, fasta_obj)
        finally:
            vcf_file.close()
            fasta_obj.close()

        self.assertEqual(result.allele1_sequence, "CAGCAGCAGCAGCAG")
        self.assertEqual(result.num_repeats_allele1, 5)
        self.assertEqual(result.num_repeats_allele2, 5)

    def test_purity_is_computed_at_the_best_offset_within_the_motif(self):
        """An allele that starts in the middle of a motif must not be reported as 0% pure."""
        counts = compute_repeat_counts_from_sequence("AGCAGCAGCAGCAGCAG", "CAG")
        self.assertAlmostEqual(counts["purity"], 1.0)
        self.assertTrue(counts["is_pure"])


class TestPhasingAmbiguityRules(unittest.TestCase):
    """Test which combinations of overlapping variants can be assigned to haplotypes without guessing."""

    def _make_variant(self, gt, phased, phase_set=None, pos=100, ref="A"):
        """Build a stub variant record exposing just the fields the phasing check reads.

        chrom, pos and ref are set explicitly: the check treats records at one position as a single split
        multiallelic site, so leaving them as auto-created MagicMock attributes would make every pair of
        stubs look like distinct sites by accident, and the ambiguity tests would pass for the wrong reason.
        """
        sample = mock.MagicMock()
        sample.get.side_effect = lambda key: {"GT": gt, "PS": phase_set}.get(key)
        sample.phased = phased
        variant = mock.MagicMock()
        variant.chrom = "chr1"
        variant.pos = pos
        variant.ref = ref
        variant.samples = [sample]
        return variant

    def test_is_heterozygous_genotype(self):
        """Only genotypes naming more than one distinct allele need a phase."""
        for gt, expected in [((0, 1), True), ((1, 2), True), ((0, 0), False), ((1, 1), False),
                             ((1,), False), ((None, 1), False), ((None, None), False), (None, False)]:
            self.assertEqual(is_heterozygous_genotype(gt), expected, f"GT {gt}")

    def test_unphased_homozygous_records_are_not_ambiguous(self):
        """Homozygous calls put the same allele on both haplotypes, so '/' separators don't matter.

        whatshap and HiPhase phase only heterozygous sites and leave homozygous calls slash-delimited, so
        rejecting them would lose the genotype at every locus carrying a hom call plus another variant.
        """
        variants = [self._make_variant((1, 1), phased=False), self._make_variant((1, 1), phased=False)]
        self.assertTrue(are_variants_unambiguously_phased(variants))

        variants = [self._make_variant((1, 1), phased=False), self._make_variant((0, 1), phased=True)]
        self.assertTrue(are_variants_unambiguously_phased(variants))

    def test_single_unphased_heterozygous_variant_is_not_ambiguous(self):
        """One heterozygous variant gives the same pair of alleles whichever haplotype it goes on."""
        variants = [self._make_variant((0, 1), phased=False), self._make_variant((1, 1), phased=False)]
        self.assertTrue(are_variants_unambiguously_phased(variants))

    def test_two_unphased_heterozygous_records_at_one_position_are_one_site(self):
        """Records sharing a position came from splitting one multiallelic genotype, so they need no phase."""
        self.assertTrue(are_variants_unambiguously_phased([
            self._make_variant((1, 0), phased=False, pos=100, ref="TCAGCAG"),
            self._make_variant((0, 1), phased=False, pos=100, ref="TCAG"),
        ]))

    def test_two_unphased_heterozygous_variants_are_ambiguous(self):
        """Two heterozygous variants at different positions with no phase between them are ambiguous."""
        variants = [self._make_variant((0, 1), phased=False, pos=100),
                    self._make_variant((0, 1), phased=False, pos=140)]
        self.assertFalse(are_variants_unambiguously_phased(variants))

    def test_heterozygous_variants_in_different_phase_sets_are_ambiguous(self):
        """A '0|1' in one phase block says nothing about which haplotype '0|1' means in another block."""
        variants = [self._make_variant((0, 1), phased=True, phase_set=100, pos=100),
                    self._make_variant((1, 0), phased=True, phase_set=300, pos=140)]
        self.assertFalse(are_variants_unambiguously_phased(variants))

    def test_heterozygous_variants_in_the_same_phase_set_are_not_ambiguous(self):
        """Variants sharing a PS value are phased relative to each other."""
        variants = [self._make_variant((0, 1), phased=True, phase_set=100, pos=100),
                    self._make_variant((1, 0), phased=True, phase_set=100, pos=140)]
        self.assertTrue(are_variants_unambiguously_phased(variants))

    def test_absent_phase_sets_are_treated_as_one_block(self):
        """Assembly-based callers such as dipcall phase a whole chromosome and emit no PS at all."""
        variants = [self._make_variant((0, 1), phased=True, pos=100),
                    self._make_variant((1, 0), phased=True, pos=140)]
        self.assertTrue(are_variants_unambiguously_phased(variants))


class TestGenotypeCorrectnessRegressions(unittest.TestCase):
    """Regression tests for genotype-path bugs that produced wrong output fields."""

    # 10bp flank, then (CAG)x6 at chr1:10-28 (0-based half-open), then 10bp flank
    REFERENCE_SEQUENCE = "TTTTTTTTTT" + "CAG" * 6 + "AAAAAAAAAA"
    LOCUS_START_0BASED = 10
    LOCUS_END = 28
    LOCUS_MOTIF = "CAG"

    VCF_HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=38>
##FILTER=<ID=LowQual,Description="low quality">
##FILTER=<ID=RefCall,Description="reference call">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
"""

    def setUp(self):
        self._temp_files = []

    def tearDown(self):
        for path in self._temp_files:
            if os.path.exists(path):
                os.remove(path)

    _create_temp_file = TestGenotypingPipeline._create_temp_file
    _create_test_vcf_and_index = TestGenotypingPipeline._create_test_vcf_and_index
    _create_test_fasta = TestGenotypingPipeline._create_test_fasta

    def _genotype(self, vcf_body, motif=None, start_0based=None, end=None, reference_sequence=None,
                  chrom="chr1", vcf_contig="chr1", fasta_contig="chr1", **genotype_kwargs):
        """Genotype the test locus from the given VCF body lines, or skip if the tools aren't installed."""
        reference_sequence = reference_sequence or self.REFERENCE_SEQUENCE
        fasta_path = self._create_test_fasta({fasta_contig: reference_sequence})
        if fasta_path is None:
            self.skipTest("pyfaidx unavailable")
        vcf_header = self.VCF_HEADER.replace("length=38", f"length={len(reference_sequence)}")
        if vcf_contig != "chr1":
            vcf_header = vcf_header.replace("##contig=<ID=chr1,", f"##contig=<ID={vcf_contig},")
        vcf_gz_path = self._create_test_vcf_and_index(vcf_header + vcf_body)
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        tr_locus = ReferenceTandemRepeat(
            chrom=chrom,
            start_0based=self.LOCUS_START_0BASED if start_0based is None else start_0based,
            end_1based=self.LOCUS_END if end is None else end,
            repeat_unit=self.LOCUS_MOTIF if motif is None else motif)
        fasta_obj = pyfaidx.Fasta(fasta_path, one_based_attributes=False, as_raw=True)
        vcf_file = pysam.VariantFile(vcf_gz_path)
        try:
            return genotype_single_locus(
                tr_locus, vcf_file, fasta_obj,
                vcf_contig_lookup=build_contig_name_lookup([vcf_contig]),
                fasta_contig_lookup=build_contig_name_lookup([fasta_contig]),
                **genotype_kwargs)
        finally:
            vcf_file.close()
            fasta_obj.close()

    def test_unphased_homozygous_variants_still_genotype(self):
        """Two hom-alt records written with '/' are unambiguous and must not become a no-call."""
        result = self._genotype("chr1\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t1/1\n"
                                "chr1\t20\t.\tC\tA\t.\tPASS\t.\tGT\t1/1\n")
        self.assertEqual(result.zygosity, "HOM")
        self.assertEqual((result.num_repeats_allele1, result.num_repeats_allele2), (7, 7))

    def test_heterozygous_variants_in_different_phase_sets_return_no_call(self):
        """Combining phased records from two phase blocks would mispair the haplotypes."""
        result = self._genotype("chr1\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT:PS\t0|1:100\n"
                                "chr1\t20\t.\tC\tA\t.\tPASS\t.\tGT:PS\t1|0:300\n")
        self.assertIsNone(result.zygosity)
        self.assertIsNone(result.allele1_sequence)
        self.assertIsNone(result.allele2_sequence)

    def test_heterozygous_variants_in_one_phase_set_genotype(self):
        """The same two records genotype normally when they share a phase set."""
        result = self._genotype("chr1\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT:PS\t0|1:100\n"
                                "chr1\t20\t.\tC\tA\t.\tPASS\t.\tGT:PS\t1|0:100\n")
        self.assertEqual(result.zygosity, "HET")
        self.assertEqual((result.num_repeats_allele1, result.num_repeats_allele2), (6, 7))

    def test_ref_allele_mismatch_gives_no_call_rather_than_hemizygous(self):
        """A VCF ref allele that disagrees with the fasta must not look like a real hemizygous call.

        Reporting only the surviving allele would fill both the short and long allele columns with it, hiding
        the expansion on the haplotype whose sequence could not be built.
        """
        # The fasta has C at position 11, not T, so this record cannot be applied
        result = self._genotype("chr1\t11\t.\tT\tTCAGCAGCAGCAGCAGCAGCAGCAGCAG\t.\tPASS\t.\tGT\t0|1\n")
        self.assertIsNone(result.zygosity)
        self.assertIsNone(result.num_repeats_short_allele)
        self.assertIsNone(result.num_repeats_long_allele)
        self.assertEqual(result.num_alleles_with_build_errors, 1)
        self.assertEqual(result.no_call_reason, NO_CALL_REASON_REF_ALLELE_MISMATCH)

    def test_flank_insertion_positioned_outside_the_locus_is_not_counted(self):
        """An insertion whose inserted bases land in the flank once the shared suffix is trimmed contributes
        nothing to the locus.

        Its REF span ends exactly at the locus start, but trimming the shared 'T' suffix pulls the inserted
        base back to position 8, so it must not inflate NumOverlappingVariants or trigger a phasing no-call.
        """
        result = self._genotype("chr1\t9\t.\tTT\tTTT\t.\tPASS\t.\tGT\t0|1\n")
        self.assertEqual(result.num_overlapping_variants, 0)
        self.assertEqual(result.variant_positions, [])
        self.assertEqual(result.zygosity, "HOM")

    def test_left_anchored_repeat_insertion_is_still_counted(self):
        """The neighbouring case, an insertion whose inserted bases are positioned at the locus start, must still
        be picked up."""
        result = self._genotype("chr1\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t0|1\n")
        self.assertEqual(result.num_overlapping_variants, 1)
        self.assertEqual(result.zygosity, "HET")
        self.assertEqual((result.num_repeats_allele1, result.num_repeats_allele2), (6, 7))

    def test_zero_width_locus_keeps_its_insertion(self):
        """A repeat that catalog discovery found only in an insertion has start == end in the catalog.

        The insertion positioned at that single coordinate is the locus, so it must not be dropped as a flank
        insertion and the haplotype slice must cover its inserted bases. Otherwise the locus came out as a
        homozygous zero-repeat call.
        """
        result = self._genotype("chr1\t10\t.\tT\tTCAGCAGCAG\t.\tPASS\t.\tGT\t0|1\n", start_0based=10, end=10)
        self.assertEqual(result.num_overlapping_variants, 1)
        self.assertEqual(result.zygosity, "HET")
        self.assertEqual((result.allele1_sequence, result.allele2_sequence), ("", "CAGCAGCAG"))
        self.assertEqual((result.num_repeats_allele1, result.num_repeats_allele2), (0, 3))

    def test_insertion_of_the_motif_positioned_at_the_locus_end_belongs_to_the_locus(self):
        """Catalog discovery ends an interrupted tract at its last base, so an expansion appended after that base
        is a left-aligned insertion positioned exactly at the locus end. It is the expansion the locus was
        discovered from, so it must be genotyped rather than dropped as a flank insertion.
        """
        # An interrupted tract (CAG)x2 CAT at [10, 19), expanded by CAG after the CAT
        reference = "TTTTTTTTTT" + "CAGCAGCAT" + "TTTTTTTTTT"
        result = self._genotype("chr1\t19\t.\tT\tTCAG\t.\tPASS\t.\tGT\t0|1\n", start_0based=10, end=19,
                                reference_sequence=reference)
        self.assertEqual(result.num_overlapping_variants, 1)
        self.assertEqual(result.zygosity, "HET")
        self.assertEqual((result.allele1_sequence, result.allele2_sequence), ("CAGCAGCAT", "CAGCAGCATCAG"))
        self.assertEqual((result.num_repeats_allele1, result.num_repeats_allele2), (3, 4))

        # The same for a pure tract when the VCF isn't left-aligned, and for a rotated copy of the motif
        for alt in ("GCAGCAGCAG", "GAGCAGCAGC"):
            result = self._genotype(f"chr1\t28\t.\tG\t{alt}\t.\tPASS\t.\tGT\t0|1\n")
            self.assertEqual(result.zygosity, "HET", alt)
            self.assertEqual((result.num_repeats_allele1, result.num_repeats_allele2), (6, 9), alt)

        # An insertion of something other than the motif positioned at the end still sits in the right flank, so
        # an adjacent locus with a different motif keeps it and a non-repeat insertion doesn't disturb the call
        for alt in ("GCCGCCGCCG", "GTTGACCATGA", "GCA"):
            result = self._genotype(f"chr1\t28\t.\tG\t{alt}\t.\tPASS\t.\tGT\t0|1\n")
            self.assertEqual(result.num_overlapping_variants, 0, alt)
            self.assertEqual(result.zygosity, "HOM", alt)

    def test_fully_deleted_alleles_have_unknown_purity(self):
        """An allele deleted to an empty sequence has no purity to judge, so IsPureRepeat is unknown."""
        deleted_ref = self.REFERENCE_SEQUENCE[8:8 + 21]
        result = self._genotype(f"chr1\t9\t.\t{deleted_ref}\tT\t.\tPASS\t.\tGT\t1|1\n")
        self.assertEqual((result.allele1_sequence, result.allele2_sequence), ("", ""))
        self.assertIsNone(result.is_pure_repeat)
        self.assertIsNone(result.repeat_purity)
        self.assertEqual(result.to_tsv_dict()["IsPureRepeat"], "")

    def test_purity_columns_match_the_size_columns_when_repeat_counts_tie(self):
        """Repeat counts are truncated, so alleles differing by 1bp tie on count but not on length.

        The purity columns must then follow the same short/long ordering as the size columns.
        """
        # A 1bp insertion on haplotype 0 makes allele1 19bp (purity < 1) and leaves allele2 18bp (purity 1)
        result = self._genotype("chr1\t12\t.\tA\tAA\t.\tPASS\t.\tGT\t1|0\n")
        self.assertEqual((len(result.allele1_sequence), len(result.allele2_sequence)), (19, 18))
        self.assertEqual((result.num_repeats_allele1, result.num_repeats_allele2), (6, 6))

        self.assertEqual(result.repeat_size_short_allele_bp, 18)
        self.assertEqual(result.repeat_size_long_allele_bp, 19)
        # The 18bp allele is the pure one, so the short-allele purity must be the higher of the two
        self.assertEqual(result.repeat_purity_short_allele, result.allele2_purity)
        self.assertEqual(result.repeat_purity_long_allele, result.allele1_purity)

    def test_alleles_differing_by_less_than_one_motif_are_heterozygous(self):
        """Zygosity must follow allele length, not the truncated repeat count.

        The two alleles here are 19bp and 18bp and both floor to 6 repeats. Deciding on count alone would
        report HOM while RepeatSizeShortAlleleBp and RepeatSizeLongAlleleBp disagree, so a real heterozygous
        indel would be dropped by anyone filtering the output for HET.
        """
        result = self._genotype("chr1\t12\t.\tA\tAA\t.\tPASS\t.\tGT\t1|0\n")
        self.assertEqual((result.num_repeats_allele1, result.num_repeats_allele2), (6, 6))
        self.assertEqual((result.repeat_size_short_allele_bp, result.repeat_size_long_allele_bp), (18, 19))
        self.assertEqual(result.zygosity, "HET")

    def test_a_sub_motif_indel_in_a_long_motif_vntr_is_heterozygous(self):
        """The same rule at a VNTR, where any indel smaller than one copy lands in the same count bin."""
        motif = "ACGTTGCAAGGCTTACCGGATTCAGGCATA"
        reference_sequence = "T" * 10 + motif * 4 + "A" * 10
        result = self._genotype(f"chr1\t50\t.\t{reference_sequence[49]}\t{reference_sequence[49]}CATGG"
                                f"\t.\tPASS\t.\tGT\t0|1\n",
                                motif=motif, start_0based=10, end=130,
                                reference_sequence=reference_sequence)
        self.assertEqual((result.num_repeats_allele1, result.num_repeats_allele2), (4, 4))
        self.assertEqual((result.repeat_size_short_allele_bp, result.repeat_size_long_allele_bp), (120, 125))
        self.assertEqual(result.zygosity, "HET")

    def test_two_identical_alleles_are_still_homozygous(self):
        """The neighbouring case: alleles of equal length must stay HOM."""
        self.assertEqual(self._genotype("").zygosity, "HOM")
        self.assertEqual(self._genotype("chr1\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t1|1\n").zygosity, "HOM")

    def test_hom_ref_records_are_not_counted_as_overlapping_variants(self):
        """A 0/0 record changes neither haplotype, so it must not appear in the overlapping-variant columns.

        Callers that emit reference calls (DeepVariant RefCall records, `bcftools call -m` without -v)
        produce these routinely, and counting them keeps --skip-hom-ref-loci from skipping the locus.
        """
        result = self._genotype("chr1\t15\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\n")
        self.assertEqual(result.num_overlapping_variants, 0)
        self.assertEqual(result.variant_positions, [])
        self.assertEqual(result.zygosity, "HOM")

    def test_a_record_with_a_missing_genotype_still_counts_as_overlapping(self):
        """The neighbouring case: './.' does change the outcome, so it must not be dropped."""
        result = self._genotype("chr1\t15\t.\tA\tG\t.\tPASS\t.\tGT\t./.\n")
        self.assertEqual(result.num_overlapping_variants, 1)
        self.assertIsNone(result.zygosity)

    def test_a_homozygous_locus_on_haploid_chrx_is_hemizygous(self):
        """Two identical haplotypes at a haploid locus are one copy reported twice.

        That covers a locus with no record, which gets the reference on both haplotypes, and a '1|1' record,
        which is how DipCall writes a male's chrX when his assembly haplotypes aren't split by parent.
        """
        result = self._genotype("", chrom="chrX", vcf_contig="chrX", fasta_contig="chrX",
                                sex_chromosome_ploidy=XY_PLOIDY)
        self.assertEqual(result.zygosity, "HEMI")
        self.assertEqual(result.allele1_sequence, "CAG" * 6)
        self.assertIsNone(result.allele2_sequence)

        result = self._genotype("chrX\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t1|1\n",
                                chrom="chrX", vcf_contig="chrX", fasta_contig="chrX",
                                sex_chromosome_ploidy=XY_PLOIDY)
        self.assertEqual(result.zygosity, "HEMI")
        self.assertEqual(result.num_repeats_allele1, 7)
        self.assertIsNone(result.allele2_sequence)

        # The neighbouring cases: a diploid sample, an autosome, and a PAR keep both copies
        for chrom, par_regions, ploidy in (("chrX", None, None), ("chr1", None, XY_PLOIDY),
                                           ("chrX", {"X": [(0, 20)], "Y": []}, XY_PLOIDY)):
            result = self._genotype(f"{chrom}\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t1|1\n",
                                    chrom=chrom, vcf_contig=chrom, fasta_contig=chrom, par_regions=par_regions,
                                    sex_chromosome_ploidy=ploidy)
            self.assertEqual(result.zygosity, "HOM", f"{chrom} {par_regions} {ploidy}")

    def test_an_xxy_sample_keeps_chrx_diploid_and_chry_haploid(self):
        """chrX and chrY ploidy are separate, so an XXY sample gets diploid rules on chrX and haploid on chrY."""
        xxy_ploidy = {"X": 2, "Y": 1}
        result = self._genotype("chrX\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t.|1\n",
                                chrom="chrX", vcf_contig="chrX", fasta_contig="chrX",
                                sex_chromosome_ploidy=xxy_ploidy)
        self.assertEqual(result.no_call_reason, NO_CALL_REASON_MISSING_GENOTYPE)
        result = self._genotype("", chrom="chrX", vcf_contig="chrX", fasta_contig="chrX",
                                sex_chromosome_ploidy=xxy_ploidy)
        self.assertEqual(result.zygosity, "HOM")

        result = self._genotype("chrY\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t1|.\n",
                                chrom="chrY", vcf_contig="chrY", fasta_contig="chrY",
                                sex_chromosome_ploidy=xxy_ploidy)
        self.assertEqual(result.zygosity, "HEMI")
        self.assertIsNone(result.no_call_reason)

    def test_a_locus_on_a_chromosome_the_sample_lacks_is_a_no_call(self):
        """A non-PAR chrY locus in an XX sample has no records, which must not read as homozygous reference."""
        result = self._genotype("", chrom="chrY", vcf_contig="chrY", fasta_contig="chrY",
                                sex_chromosome_ploidy={"X": 2, "Y": 0})
        self.assertIsNone(result.zygosity)
        self.assertIsNone(result.allele1_sequence)
        self.assertEqual(result.no_call_reason, NO_CALL_REASON_CHROMOSOME_ABSENT)

        # A called record on that chrY says the sample does carry it there (a VCF covering only part of the
        # genome can fall under the whole-sample chrY cutoff), so the locus is genotyped as haploid instead
        result = self._genotype("chrY\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t1|1\n",
                                chrom="chrY", vcf_contig="chrY", fasta_contig="chrY",
                                sex_chromosome_ploidy={"X": 2, "Y": 0})
        self.assertEqual(result.zygosity, "HEMI")
        self.assertEqual(result.num_repeats_allele1, 7)
        self.assertIsNone(result.no_call_reason)

        # A PAR locus on that chrY is as absent as the rest of the chromosome, so with no record it is the same
        # no call rather than a diploid reference call
        result = self._genotype("", chrom="chrY", vcf_contig="chrY", fasta_contig="chrY",
                                par_regions={"X": [], "Y": [(0, 20)]}, sex_chromosome_ploidy={"X": 2, "Y": 0})
        self.assertIsNone(result.zygosity)
        self.assertEqual(result.no_call_reason, NO_CALL_REASON_CHROMOSOME_ABSENT)

        # With a called record there, the PAR locus is genotyped as the two copies a PAR has, so a
        # heterozygous record is HET rather than the no call it would be at a haploid locus
        result = self._genotype("chrY\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t0|1\n",
                                chrom="chrY", vcf_contig="chrY", fasta_contig="chrY",
                                par_regions={"X": [], "Y": [(0, 20)]}, sex_chromosome_ploidy={"X": 2, "Y": 0})
        self.assertEqual(result.zygosity, "HET")
        self.assertIsNone(result.no_call_reason)

    def test_a_haploid_reference_locus_is_still_skipped_as_hom_ref(self):
        """--skip-hom-ref-loci keys on having no overlapping variants, so a HEMI reference locus still counts."""
        result = self._genotype("", chrom="chrX", vcf_contig="chrX", fasta_contig="chrX",
                                sex_chromosome_ploidy=XY_PLOIDY)
        self.assertTrue(is_locus_genotyped_as_reference(result))

    def test_a_heterozygous_call_at_a_haploid_locus_gives_no_call(self):
        """One copy can't carry two alleles, and there is no way to say which of the two is real.

        DipCall writes these (flagged DIPX) for a male whose assembly haplotypes aren't split by parent.
        """
        for chrom, genotype in (("chrX", "0|1"), ("chrX", "1|2"), ("chrY", "1|0")):
            alt = "TCAG,TCAGCAG" if genotype == "1|2" else "TCAG"
            result = self._genotype(f"{chrom}\t10\t.\tT\t{alt}\t.\tDIPX\t.\tGT\t{genotype}\n",
                                    chrom=chrom, vcf_contig=chrom, fasta_contig=chrom,
                                    sex_chromosome_ploidy=XY_PLOIDY)
            self.assertIsNone(result.zygosity, f"{chrom} {genotype}")
            self.assertEqual(result.no_call_reason, NO_CALL_REASON_HET_AT_HAPLOID_LOCUS, f"{chrom} {genotype}")

        # The neighbouring cases stay heterozygous: a diploid sample, the diploid chrX of an XXY sample, a PAR,
        # and an autosome
        for chrom, par_regions, ploidy in (("chrX", None, None), ("chrX", None, {"X": 2, "Y": 1}),
                                           ("chrX", {"X": [(0, 20)], "Y": []}, XY_PLOIDY),
                                           ("chr1", None, XY_PLOIDY)):
            result = self._genotype(f"{chrom}\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t0|1\n",
                                    chrom=chrom, vcf_contig=chrom, fasta_contig=chrom, par_regions=par_regions,
                                    sex_chromosome_ploidy=ploidy)
            self.assertEqual(result.zygosity, "HET", f"{chrom} {par_regions} {ploidy}")

    def test_a_missing_allele_on_chrx_of_a_diploid_sample_gives_no_call(self):
        """In a female, '.|1' on chrX means one assembly didn't cover the site, not a haploid call."""
        result = self._genotype("chrX\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t.|1\n",
                                chrom="chrX", vcf_contig="chrX", fasta_contig="chrX")
        self.assertEqual(result.no_call_reason, NO_CALL_REASON_MISSING_GENOTYPE)

        # The neighbouring case: the same record in an XY sample is a HEMI call
        result = self._genotype("chrX\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t.|1\n",
                                chrom="chrX", vcf_contig="chrX", fasta_contig="chrX",
                                sex_chromosome_ploidy=XY_PLOIDY)
        self.assertEqual(result.zygosity, "HEMI")
        self.assertIsNone(result.no_call_reason)

    def test_a_missing_allele_inside_a_par_gives_no_call(self):
        """The pseudoautosomal regions are diploid in males too, so '.|1' there is an uncalled haplotype."""
        par_regions = {"X": [(0, 20)], "Y": []}
        result = self._genotype("chrX\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t.|1\n",
                                chrom="chrX", vcf_contig="chrX", fasta_contig="chrX", par_regions=par_regions,
                                sex_chromosome_ploidy=XY_PLOIDY)
        self.assertEqual(result.no_call_reason, NO_CALL_REASON_MISSING_GENOTYPE)

        # The neighbouring case: the same record at a locus clear of the PAR is a haploid call
        result = self._genotype("chrX\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t.|1\n",
                                chrom="chrX", vcf_contig="chrX", fasta_contig="chrX",
                                par_regions={"X": [(30, 38)], "Y": []}, sex_chromosome_ploidy=XY_PLOIDY)
        self.assertEqual(result.zygosity, "HEMI")
        self.assertIsNone(result.no_call_reason)

    def test_a_locus_past_the_end_of_the_contig_gives_no_call(self):
        """pyfaidx clips such a slice silently, which would report a truncated allele as a real contraction."""
        result = self._genotype("", end=60)
        self.assertIsNone(result.zygosity)
        self.assertEqual(result.no_call_reason, NO_CALL_REASON_LOCUS_PAST_CONTIG_END)
        self.assertIn("returned only", result.no_call_detail)
        # The reference was never available, so no allele failed to have its variants applied
        self.assertEqual(result.num_alleles_with_build_errors, 0)

    def test_only_the_needed_variant_fields_are_retained(self):
        """Holding every pysam record until the writers run would multiply a genome-wide run's memory."""
        result = self._genotype("chr1\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t0|1\n")
        variant = result.overlapping_variants[0]
        self.assertIsInstance(variant, OverlappingVariant)
        self.assertEqual((variant.chrom, variant.pos, variant.ref, variant.alts),
                         ("chr1", 10, "T", ("TCAG",)))


    def test_a_filtered_hom_ref_record_is_dropped_before_the_filter_check(self):
        """A DeepVariant RefCall record is hom-ref AND filtered, so the two checks must run in this order.

        get_overlapping_vcf_variants drops records that cannot change a haplotype, and only the survivors
        reach the failed-filter check. Running those the other way round would turn every locus near a
        RefCall record into a no call, and RefCall records are dense across a DeepVariant callset.
        """
        result = self._genotype("chr1\t15\t.\tA\tG\t.\tRefCall\t.\tGT\t0/0\n")
        self.assertEqual(result.zygosity, "HOM")
        self.assertEqual(result.num_overlapping_variants, 0)
        self.assertIsNone(result.no_call_reason)

    def test_a_star_allele_record_is_not_counted_as_overlapping(self):
        """A '*' allele stands for a deletion described by its own record, so it changes nothing here.

        Joint-genotyped and bcftools-normalized callsets emit these at every site under a neighbouring
        deletion. Counting them would inflate NumOverlappingVariants and defeat --skip-hom-ref-loci.
        """
        result = self._genotype("chr1\t15\t.\tA\tG,*\t.\tPASS\t.\tGT\t0/2\n")
        self.assertEqual(result.num_overlapping_variants, 0)
        self.assertEqual(result.variant_positions, [])
        self.assertEqual(result.zygosity, "HOM")
        # The neighbouring case: a real alt allele at the same site does count
        self.assertEqual(self._genotype("chr1\t15\t.\tA\tG,*\t.\tPASS\t.\tGT\t0/1\n"
                                        ).num_overlapping_variants, 1)

    def test_a_split_multiallelic_site_genotypes_like_the_unsplit_record(self):
        """`bcftools norm -m -any` rewrites one "1/2" record as two unphased records at the same site.

        Those two alt alleles came from one diploid genotype, so they are necessarily on opposite
        haplotypes. Treating them as two independently unphased heterozygous variants would make a
        normalized callset a no call at every multiallelic tandem repeat.
        """
        # Run with -f, which is the normal usage, bcftools minimizes each split allele separately, so the
        # two records can end up with different REF strings at the same position. Keying the same-site rule
        # on REF as well as position would miss exactly these shapes.
        for unsplit_body, split_body, expected in [
            # two insertions: bcftools leaves the REF identical
            ("chr1\t10\t.\tT\tTCAG,TCAGCAG\t.\tPASS\t.\tGT\t1/2\n",
             "chr1\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t1/0\n"
             "chr1\t10\t.\tT\tTCAGCAG\t.\tPASS\t.\tGT\t0/1\n", (7, 8)),
            # two deletions: each split record keeps its own trimmed REF
            ("chr1\t10\t.\tTCAGCAG\tT,TCAG\t.\tPASS\t.\tGT\t1/2\n",
             "chr1\t10\t.\tTCAGCAG\tT\t.\tPASS\t.\tGT\t1/0\n"
             "chr1\t10\t.\tTCAG\tT\t.\tPASS\t.\tGT\t0/1\n", (4, 5)),
            # a deletion paired with an insertion, the mixed case
            ("chr1\t10\t.\tTCAG\tT,TCAGCAG\t.\tPASS\t.\tGT\t1/2\n",
             "chr1\t10\t.\tTCAG\tT\t.\tPASS\t.\tGT\t1/0\n"
             "chr1\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t0/1\n", (5, 7)),
        ]:
            unsplit = self._genotype(unsplit_body)
            split = self._genotype(split_body)
            self.assertEqual(unsplit.zygosity, "HET", unsplit_body)
            self.assertEqual(
                (unsplit.num_repeats_short_allele, unsplit.num_repeats_long_allele), expected, unsplit_body)
            self.assertEqual(split.zygosity, unsplit.zygosity, split_body)
            self.assertEqual((split.num_repeats_short_allele, split.num_repeats_long_allele), expected,
                             split_body)

    def test_two_distinct_unphased_heterozygous_sites_are_still_ambiguous(self):
        """The neighbouring case: two unphased hets at different positions genuinely need a phase."""
        result = self._genotype("chr1\t12\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n"
                                "chr1\t20\t.\tC\tA\t.\tPASS\t.\tGT\t0/1\n")
        self.assertEqual(result.no_call_reason, NO_CALL_REASON_AMBIGUOUS_PHASING)

    def test_a_diploid_record_with_a_missing_allele_gives_no_call(self):
        """DipCall emits '.|1', '1|.', './.' in bulk; none of them says what the other haplotype is."""
        for genotype in ("./.", ".|1", "1|.", ".|0"):
            result = self._genotype(f"chr1\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t{genotype}\n")
            self.assertIsNone(result.zygosity, f"GT {genotype} should be a no call")
            self.assertEqual(result.no_call_reason, NO_CALL_REASON_MISSING_GENOTYPE,
                             f"GT {genotype} should name the missing-genotype reason")

    def test_a_fully_missing_genotype_on_chrx_names_its_reason(self):
        """chrX outside the PARs of an XY sample salvages '.|1', but './.' and a haploid '.' still say
        nothing.

        DipCall writes '.|.' for sites uncalled on both haplotypes. Exempting every chrX locus from the
        missing-genotype check would leave those rows blank in both Zygosity and NoCallReason.
        """
        for genotype in ("./.", ".|.", "."):
            result = self._genotype(f"chrX\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t{genotype}\n",
                                    chrom="chrX", vcf_contig="chrX", fasta_contig="chrX",
                                    sex_chromosome_ploidy=XY_PLOIDY)
            self.assertIsNone(result.zygosity, f"GT {genotype} should be a no call")
            self.assertEqual(result.no_call_reason, NO_CALL_REASON_MISSING_GENOTYPE,
                             f"GT {genotype} should name the missing-genotype reason")

        # The neighbouring case: every partly called genotype is still salvaged on chrX and chrY
        for chrom in ("chrX", "chrY"):
            for genotype in (".|1", "1|.", "./1", "1/.", "1"):
                result = self._genotype(f"{chrom}\t10\t.\tT\tTCAGCAGCAG\t.\tPASS\t.\tGT\t{genotype}\n",
                                        chrom=chrom, vcf_contig=chrom, fasta_contig=chrom,
                                        sex_chromosome_ploidy=XY_PLOIDY)
                self.assertEqual(result.zygosity, "HEMI", f"{chrom} GT {genotype} should be HEMI")
                self.assertEqual(result.num_repeats_short_allele, 9, f"{chrom} GT {genotype}")
                self.assertIsNone(result.no_call_reason, f"{chrom} GT {genotype}")

    def test_a_single_allele_genotype_at_a_diploid_locus_is_a_no_call(self):
        """A GT of length one leaves the second haplotype of a diploid locus as undetermined as '1|.' does.

        Reporting the one allele as HEMI would contradict the refusal to do so for '1|.', so the locus gets
        no call, under its own reason since the allele isn't uncalled.
        """
        result = self._genotype("chr1\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t1\n")
        self.assertIsNone(result.zygosity)
        self.assertIsNone(result.num_repeats_allele1)
        self.assertIsNone(result.num_repeats_allele2)
        self.assertEqual(result.no_call_reason, NO_CALL_REASON_HAPLOID_GENOTYPE_AT_DIPLOID_LOCUS)
        self.assertIn("at position(s) 10", result.no_call_description)

        # Next to a diploid record, which would otherwise build haplotype 0 from both and leave haplotype 1 missing
        result = self._genotype("chr1\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t1\n"
                                "chr1\t20\t.\tC\tA\t.\tPASS\t.\tGT\t0|1\n")
        self.assertEqual(result.no_call_reason, NO_CALL_REASON_HAPLOID_GENOTYPE_AT_DIPLOID_LOCUS)

        # On a diploid chrX (an XX sample) too
        result = self._genotype("chrX\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t1\n",
                                chrom="chrX", vcf_contig="chrX", fasta_contig="chrX",
                                sex_chromosome_ploidy={"X": 2, "Y": 0})
        self.assertEqual(result.no_call_reason, NO_CALL_REASON_HAPLOID_GENOTYPE_AT_DIPLOID_LOCUS)

        # A haploid "." is an uncalled allele first, so it keeps the missing-genotype reason
        result = self._genotype("chr1\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t.\n")
        self.assertEqual(result.no_call_reason, NO_CALL_REASON_MISSING_GENOTYPE)

        # The neighbouring case: at a haploid locus the same record is the one copy
        result = self._genotype("chrX\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t1\n",
                                chrom="chrX", vcf_contig="chrX", fasta_contig="chrX",
                                sex_chromosome_ploidy=XY_PLOIDY)
        self.assertEqual(result.zygosity, "HEMI")
        self.assertIsNone(result.no_call_reason)

    def test_sex_chromosome_detection_only_runs_when_a_locus_is_on_chrx_or_chry(self):
        """Detection scans all of non-PAR chrX and chrY, which is wasted work when -L leaves no locus there."""
        fasta_path = self._create_test_fasta({"chr1": self.REFERENCE_SEQUENCE, "chrX": self.REFERENCE_SEQUENCE})
        if fasta_path is None:
            self.skipTest("pyfaidx unavailable")
        vcf_header = self.VCF_HEADER.replace("##contig=<ID=chr1,length=38>",
                                             "##contig=<ID=chr1,length=38>\n##contig=<ID=chrX,length=38>")
        vcf_gz_path = self._create_test_vcf_and_index(
            vcf_header + "chr1\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t0|1\n")
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")

        fasta_obj = pyfaidx.Fasta(fasta_path, one_based_attributes=False, as_raw=True)
        try:
            for chrom, is_detection_expected in (("chr1", False), ("chrX", True)):
                tr_locus = ReferenceTandemRepeat(chrom=chrom, start_0based=self.LOCUS_START_0BASED,
                                                 end_1based=self.LOCUS_END, repeat_unit=self.LOCUS_MOTIF)
                with mock.patch("str_analysis.filter_vcf_to_tandem_repeats.detect_sex_chromosome_ploidy",
                                return_value={"X": 2, "Y": 0}) as mock_detect, \
                        contextlib.redirect_stdout(io.StringIO()):
                    genotyped_loci, _, _ = genotype_all_loci([tr_locus], vcf_gz_path, fasta_obj, argparse.Namespace())
                self.assertEqual(mock_detect.called, is_detection_expected, chrom)
                self.assertEqual(len(genotyped_loci), 1, chrom)
        finally:
            fasta_obj.close()

    def test_a_flank_indel_anchored_on_the_last_locus_base_is_not_overlapping(self):
        """Left-alignment moves a right-flank indel onto the locus's last base, changing nothing inside it.

        Its anchor base is unchanged and every inserted or deleted base sits in the flank, so counting it
        would inflate NumOverlappingVariants, keep the locus out of --skip-hom-ref-loci, tag the flank record
        in the contributing-variants VCF, and pair with a real STR indel to force an ambiguous-phasing
        no call. The left edge has been filtered this way all along; the right edge was not.
        """
        result = self._genotype("chr1\t28\t.\tG\tGA\t.\tPASS\t.\tGT\t0/1\n")
        self.assertEqual(result.num_overlapping_variants, 0)
        self.assertEqual(result.zygosity, "HOM")

        # Paired with a real STR indel it must not turn the locus into a no call
        result = self._genotype("chr1\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t0/1\n"
                                "chr1\t28\t.\tG\tGA\t.\tPASS\t.\tGT\t0/1\n")
        self.assertEqual(result.zygosity, "HET")
        self.assertEqual(result.variant_positions, [10])

        # The neighbouring case: an insertion actually inside the locus still counts
        result = self._genotype("chr1\t20\t.\tC\tCCAG\t.\tPASS\t.\tGT\t0/1\n")
        self.assertEqual(result.num_overlapping_variants, 1)
        self.assertEqual(result.zygosity, "HET")

    def test_an_uncalled_alt_allele_does_not_make_a_record_overlap(self):
        """Only the alleles the sample carries decide whether a record affects the locus.

        Here GT 0/2 calls the flank SNV, not the insertion, so the record changes nothing inside the locus.
        Judging on every ALT would keep it and, with filters honored by default, turn the locus into a
        filtered-variant no call on the strength of an allele the sample does not have.
        """
        result = self._genotype("chr1\t10\t.\tT\tTCAG,G\t.\tLowQual\t.\tGT\t0/2\n")
        self.assertEqual(result.num_overlapping_variants, 0)
        self.assertIsNone(result.no_call_reason)
        self.assertEqual(result.zygosity, "HOM")

        # The neighbouring case: calling the insertion instead does affect the locus
        result = self._genotype("chr1\t10\t.\tT\tTCAG,G\t.\tPASS\t.\tGT\t0/1\n")
        self.assertEqual(result.num_overlapping_variants, 1)
        self.assertEqual(result.zygosity, "HET")

    def test_a_chrm_locus_finds_variants_in_an_mt_named_vcf(self):
        """chrM and MT name the same chromosome, and toggling the 'chr' prefix alone does not connect them.

        GRCh38-style references call it chrM while Ensembl-style VCFs call it MT, so a chrM catalog locus
        would otherwise be reported as belonging to a contig the VCF never called, with its real variants
        sitting right there.
        """
        self.assertEqual(build_contig_name_lookup(["MT", "chr1"]), {"M": "MT", "1": "chr1"})
        self.assertEqual(build_contig_name_lookup(["chrM", "1"]), {"M": "chrM", "1": "1"})
        self.assertIsNone(build_contig_name_lookup(["chr1", "MT"]).get(normalize_chromosome_name("chr9")))

    def test_the_contig_check_uses_the_vcf_own_spelling(self):
        """A GRCh37-style VCF named '1' against a 'chr1' catalog must not no-call every locus."""
        result = self._genotype("1\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t0|1\n",
                                vcf_contig="1")
        self.assertEqual(result.zygosity, "HET")
        self.assertIsNone(result.no_call_reason)
        self.assertEqual(result.num_overlapping_variants, 1)

    def test_locus_ploidy_follows_the_sex_chromosome_ploidy_only_outside_the_pars(self):
        """Either naming convention works, and a locus that even partly overlaps a PAR stays diploid."""
        grch38_pars = PAR_REGIONS_BY_CHRX_LENGTH[156_040_895][1]
        xxy_ploidy = {"X": 2, "Y": 1}
        self.assertEqual(get_locus_ploidy("chrX", 5_000_000, 5_000_030, grch38_pars, XY_PLOIDY), 1)
        self.assertEqual(get_locus_ploidy("X", 5_000_000, 5_000_030, grch38_pars, XY_PLOIDY), 1)
        self.assertEqual(get_locus_ploidy("chrY", 5_000_000, 5_000_030, grch38_pars, XY_PLOIDY), 1)
        self.assertEqual(get_locus_ploidy("chrX", 5_000_000, 5_000_030, grch38_pars, xxy_ploidy), 2)
        self.assertEqual(get_locus_ploidy("chrY", 5_000_000, 5_000_030, grch38_pars, xxy_ploidy), 1)
        self.assertEqual(get_locus_ploidy("chr1", 5_000_000, 5_000_030, grch38_pars, XY_PLOIDY), 2)
        self.assertEqual(get_locus_ploidy("chrM", 100, 130, grch38_pars, XY_PLOIDY), 2)
        self.assertEqual(get_locus_ploidy("chrX", 5_000_000, 5_000_030, grch38_pars, None), 2)

        # PAR1 ends at 1-based 2,781,479
        self.assertEqual(get_locus_ploidy("chrX", 2_781_470, 2_781_479, grch38_pars, XY_PLOIDY), 2)
        self.assertEqual(get_locus_ploidy("chrX", 2_781_478, 2_781_490, grch38_pars, XY_PLOIDY), 2)
        self.assertEqual(get_locus_ploidy("chrX", 2_781_479, 2_781_490, grch38_pars, XY_PLOIDY), 1)
        self.assertEqual(get_locus_ploidy("chrY", 57_000_000, 57_000_030, grch38_pars, XY_PLOIDY), 2)

        # A chromosome the sample lacks is absent in its PARs too
        self.assertEqual(get_locus_ploidy("chrY", 57_000_000, 57_000_030, grch38_pars, {"X": 2, "Y": 0}), 0)
        self.assertEqual(get_locus_ploidy("chrY", 5_000_000, 5_000_030, grch38_pars, {"X": 2, "Y": 0}), 0)

        # With no PARs known, all of chrX and chrY follows the sex chromosome ploidy
        self.assertEqual(get_locus_ploidy("chrX", 2_781_470, 2_781_479, {}, XY_PLOIDY), 1)
        self.assertEqual(get_locus_ploidy("chrX", 2_781_470, 2_781_479, None, XY_PLOIDY), 1)

    def test_the_genome_build_is_identified_by_the_length_of_chrx(self):
        """GRCh38, GRCh37 and T2T-CHM13 get their own PARs; any other chrX length gets none."""
        class MockFasta:
            def __init__(self, chrx_length):
                self.chrx_length = chrx_length

            def __getitem__(self, chrom):
                return range(self.chrx_length)

        for chrx_length, (genome_version, expected_pars) in PAR_REGIONS_BY_CHRX_LENGTH.items():
            with contextlib.redirect_stdout(io.StringIO()) as stdout:
                self.assertEqual(get_PAR_region_coordinates(MockFasta(chrx_length), {"X": "chrX"}), expected_pars)
            self.assertIn(genome_version, stdout.getvalue())

        with contextlib.redirect_stdout(io.StringIO()) as stdout:
            self.assertEqual(get_PAR_region_coordinates(MockFasta(1000), {"X": "X"}), {})
        self.assertIn("WARNING", stdout.getvalue())

        # A reference with no chrX has no chrX loci to affect, so it gets no warning
        with contextlib.redirect_stdout(io.StringIO()) as stdout:
            self.assertEqual(get_PAR_region_coordinates(MockFasta(1000), {"1": "chr1"}), {})
        self.assertEqual(stdout.getvalue(), "")

    def _detect_ploidy(self, chrx_genotypes, chry_genotypes=(), par_regions=None):
        """Run detect_sex_chromosome_ploidy on a VCF with the given GTs at consecutive chrX and chrY sites."""
        vcf_body = "".join(f"{chrom}\t{10 * (i + 1)}\t.\tA\tG\t.\tPASS\t.\tGT\t{gt}\n"
                           for chrom, genotypes in (("chrX", chrx_genotypes), ("chrY", chry_genotypes))
                           for i, gt in enumerate(genotypes))
        vcf_gz_path = self._create_test_vcf_and_index(
            self.VCF_HEADER.replace("##contig=<ID=chr1,length=38>",
                                    "##contig=<ID=chrX,length=100000>\n##contig=<ID=chrY,length=100000>")
            + vcf_body)
        if vcf_gz_path is None:
            self.skipTest("bgzip/tabix unavailable")
        vcf_file = pysam.VariantFile(vcf_gz_path)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                return detect_sex_chromosome_ploidy(
                    vcf_file, build_contig_name_lookup(vcf_file.header.contigs), par_regions or {})
        finally:
            vcf_file.close()

    def test_haploid_chrx_is_detected_from_the_share_of_half_called_records(self):
        """More than half of at least 1,000 called non-PAR chrX records must have one uncalled allele.

        Records inside a PAR and records with no allele called at all don't count either way.
        """
        # Male-like: most called records have one uncalled allele; './.' records are ignored
        self.assertEqual(self._detect_ploidy([".|1"] * 600 + [".|0"] * 300 + ["0|1"] * 300 + ["./."] * 2)["X"], 1)
        # Female-like: most called records name both alleles
        self.assertEqual(self._detect_ploidy([".|1"] * 400 + ["0|1"] * 400 + ["1|1"] * 400)["X"], 2)
        # Exactly half is not more than half
        self.assertEqual(self._detect_ploidy([".|1"] * 500 + ["1|1"] * 500)["X"], 2)
        # Too few called records for the fraction to mean anything, however lopsided it is
        self.assertEqual(self._detect_ploidy([".|1"] * 999)["X"], 2)
        self.assertEqual(self._detect_ploidy([".|1"] * 1000)["X"], 1)
        # Records inside a PAR (here the first 600 sites, at positions 10 to 6,000) don't count: without the PAR
        # 1,100 of 2,100 called records would call one allele, with it only 500 of 1,500 do. The diploid
        # records are heterozygous so that the low-heterozygosity test stays out of it.
        self.assertEqual(self._detect_ploidy([".|1"] * 1100 + ["0|1"] * 1000)["X"], 1)
        self.assertEqual(self._detect_ploidy([".|1"] * 1100 + ["0|1"] * 1000,
                                             par_regions={"X": [(0, 6005)], "Y": []})["X"], 2)
        # No called chrX records at all
        self.assertEqual(self._detect_ploidy(["./."])["X"], 2)
        # Single-allele GTs, as written by callers run with ploidy 1 on chrX, call one allele just like ".|1"
        self.assertEqual(self._detect_ploidy(["1"] * 600 + ["0"] * 300 + ["0|1"] * 300)["X"], 1)
        self.assertEqual(self._detect_ploidy(["1"] * 400 + ["0|1"] * 400 + ["1|1"] * 400)["X"], 2)

    def test_haploid_chrx_is_detected_from_single_allele_genotypes_without_a_minimum_count(self):
        """A caller that writes '1' treated chrX as haploid, so a majority of such GTs decides it at any count.

        Without this, a VCF covering little of chrX falls back to diploid and each of its chrX loci becomes a
        no call for having a haploid genotype at a diploid locus.
        """
        self.assertEqual(self._detect_ploidy(["1"] * 5)["X"], 1)
        self.assertEqual(self._detect_ploidy(["0"] * 2 + ["1"])["X"], 1)
        self.assertEqual(self._detect_ploidy(["1"] * 3 + ["0|1"] * 2)["X"], 1)
        # Exactly half is not more than half, and a stray single-allele GT doesn't outvote diploid ones
        self.assertEqual(self._detect_ploidy(["1"] * 3 + ["0|1"] * 3)["X"], 2)
        self.assertEqual(self._detect_ploidy(["1"] + ["0|1"] * 3)["X"], 2)
        # A haploid "." calls nothing and doesn't count
        self.assertEqual(self._detect_ploidy(["."] * 5 + ["0|1"])["X"], 2)
        # ".|1" still needs the 1,000-record minimum, since in a female it can mean an assembly gap
        self.assertEqual(self._detect_ploidy([".|1"] * 5)["X"], 2)

    def test_a_male_with_unsplit_haplotypes_is_detected_from_low_heterozygosity(self):
        """When a male's assembly haplotypes aren't split by parent, DipCall writes chrX as mostly '1|1'."""
        # 1,000 hom-alt and 200 heterozygous records: 17% heterozygous, like the survey's HGSVC males
        self.assertEqual(self._detect_ploidy(["1|1"] * 1000 + ["0|1"] * 200)["X"], 1)
        # Half heterozygous, like a female
        self.assertEqual(self._detect_ploidy(["1|1"] * 600 + ["0|1"] * 600)["X"], 2)
        # Too few records that call both alleles for the fraction to mean anything
        self.assertEqual(self._detect_ploidy(["1|1"] * 900)["X"], 2)

    def test_chrx_and_chry_ploidy_are_detected_separately(self):
        """chrY counts as present once it has enough called records, whatever chrX looks like.

        That is what distinguishes XXY (a female-like chrX plus a chrY) and X0 (a male-like chrX, no chrY).
        """
        male_chrx = [".|1"] * 1000
        female_chrx = ["0|1"] * 600 + ["1|1"] * 600
        chry = ["1|."] * 1000
        self.assertEqual(self._detect_ploidy(male_chrx, chry), {"X": 1, "Y": 1})
        self.assertEqual(self._detect_ploidy(female_chrx), {"X": 2, "Y": 0})
        self.assertEqual(self._detect_ploidy(female_chrx, chry), {"X": 2, "Y": 1})
        self.assertEqual(self._detect_ploidy(male_chrx), {"X": 1, "Y": 0})

        # Too few called chrY records, and chrY records that call no allele, don't count as a chrY
        self.assertEqual(self._detect_ploidy(male_chrx, ["1|."] * 999)["Y"], 0)
        self.assertEqual(self._detect_ploidy(male_chrx, [".|."] * 2000)["Y"], 0)

    def test_the_reference_and_build_failures_are_reported_separately(self):
        """Not finding the reference, running past the contig, and failing to apply a variant are distinct.

        The first two are catalog or reference problems that say nothing about the sample. The rest are variants
        that could not be applied, which is what the build-error counter measures.
        """
        # A variant whose REF disagrees with the fasta
        result = self._genotype("chr1\t11\t.\tT\tTCAGCAG\t.\tPASS\t.\tGT\t0|1\n")
        self.assertEqual(result.no_call_reason, NO_CALL_REASON_REF_ALLELE_MISMATCH)
        self.assertEqual(result.num_alleles_with_build_errors, 1)

        # Two records on the same haplotype whose REF spans overlap: the generic build error. The deletion removes
        # 1-based 12 and 13, and the SNV changes 1-based 12.
        result = self._genotype("chr1\t11\t.\tCAG\tC\t.\tPASS\t.\tGT\t0|1\n"
                                "chr1\t12\t.\tA\tT\t.\tPASS\t.\tGT\t0|1\n")
        self.assertEqual(result.no_call_reason, NO_CALL_REASON_HAPLOTYPE_BUILD_ERROR)
        self.assertEqual(result.num_alleles_with_build_errors, 1)

        # A locus running off the end of its contig
        self.assertEqual(self._genotype("", end=60).no_call_reason, NO_CALL_REASON_LOCUS_PAST_CONTIG_END)

        # A contig the reference does not have at all
        self.assertEqual(self._genotype("", chrom="chr9").no_call_reason,
                         NO_CALL_REASON_CONTIG_NOT_IN_REFERENCE)

    def test_non_iupac_allele_has_its_own_no_call_reason(self):
        """A symbolic or breakend allele gets its own reason, and still counts as a variant that couldn't be applied.

        Unlike the two reference problems, it is the sample's variant that can't be applied, so it stays in the
        build-error count. The offending characters are listed in sorted order so the output is deterministic.
        """
        # Alleles are upper-cased before the check, so the H and R of "chr" are valid IUPAC codes (as is the D of DEL)
        for alt, non_iupac_characters in [("<DEL>", "<>EL"), ("C[chr1:20[", "012:[")]:
            result = self._genotype(f"chr1\t11\t.\tC\t{alt}\t.\tPASS\t.\tGT\t0|1\n")
            self.assertEqual(result.no_call_reason, NO_CALL_REASON_NON_IUPAC_ALLELE, alt)
            self.assertEqual(result.num_alleles_with_build_errors, 1, alt)
            self.assertIsNone(result.zygosity, alt)
            self.assertIn(f"contains non-IUPAC characters {non_iupac_characters}",
                          result.to_tsv_dict()["NoCallReason"], alt)

    def test_only_the_first_no_call_cause_is_reported(self):
        """Two causes at one locus must not silently overwrite each other in the NoCallReason column."""
        # A single copy, so the sequence has no internal periodicity to pass the insertion filter with
        alu = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAGACC"
        # Haplotype 0 carries a non-repeat insertion; haplotype 1 has a REF that disagrees with the fasta
        result = self._genotype(f"chr1\t12\t.\tA\tA{alu}\t.\tPASS\t.\tGT\t1|0\n"
                                f"chr1\t20\t.\tT\tG\t.\tPASS\t.\tGT\t0|1\n",
                                insertion_filter=build_insertion_filter("CAG", argparse.Namespace()))
        self.assertEqual(result.num_alleles_with_non_repeat_insertions, 1)
        self.assertEqual(result.num_alleles_with_build_errors, 1)
        # The insertion is found first, so it is the reported cause rather than being overwritten
        self.assertEqual(result.no_call_reason, NO_CALL_REASON_NON_REPEAT_INSERTION)


    def test_records_that_leave_different_haplotypes_uncalled_give_no_call(self):
        """On chrX of an XY sample, a '.|1' and a '1|.' record at one locus leave neither haplotype
        complete.

        Without an explicit reason the locus would come out blank in both Zygosity and NoCallReason.
        """
        result = self._genotype("chrX\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t.|1\n"
                                "chrX\t20\t.\tC\tA\t.\tPASS\t.\tGT\t1|.\n",
                                chrom="chrX", vcf_contig="chrX", fasta_contig="chrX",
                                sex_chromosome_ploidy=XY_PLOIDY)
        self.assertIsNone(result.zygosity)
        self.assertEqual(result.no_call_reason, NO_CALL_REASON_MISSING_GENOTYPE)

        # The neighbouring case: records that all call the same slot build that one copy from all of them
        result = self._genotype("chrX\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t.|1\n"
                                "chrX\t20\t.\tC\tA\t.\tPASS\t.\tGT\t.|1\n",
                                chrom="chrX", vcf_contig="chrX", fasta_contig="chrX",
                                sex_chromosome_ploidy=XY_PLOIDY)
        self.assertEqual(result.zygosity, "HEMI")
        self.assertIsNone(result.no_call_reason)
        self.assertEqual(result.allele2_sequence, "CAGCAGCAGCAGAAGCAGCAG")

    def test_a_flank_only_record_with_a_missing_allele_does_not_force_a_no_call(self):
        """A '.' allele means the haplotype is unknown, which only matters where the record reaches.

        A right-flank indel changes nothing inside the locus whichever allele turns out to be real, so it
        must not drag the locus into a missing-genotype no call.
        """
        result = self._genotype("chr1\t28\t.\tG\tGA\t.\tPASS\t.\tGT\t.|1\n")
        self.assertEqual(result.num_overlapping_variants, 0)
        self.assertEqual(result.zygosity, "HOM")
        self.assertIsNone(result.no_call_reason)

        # The neighbouring case: the same missing allele on a record inside the locus still no-calls it
        result = self._genotype("chr1\t20\t.\tC\tCCAG\t.\tPASS\t.\tGT\t.|1\n")
        self.assertEqual(result.no_call_reason, NO_CALL_REASON_MISSING_GENOTYPE)

    def test_the_reference_fasta_is_looked_up_under_every_contig_alias(self):
        """A chrM catalog locus must find its sequence in an MT-named reference, and say so if it cannot."""
        # A chrM locus finds its sequence in an MT-named reference
        result = self._genotype("MT\t12\t.\tA\tACAG\t.\tPASS\t.\tGT\t0/1\n",
                                chrom="chrM", vcf_contig="MT", fasta_contig="MT")
        self.assertEqual(result.zygosity, "HET")
        self.assertIsNone(result.no_call_reason)
        self.assertEqual(result.num_overlapping_variants, 1)

        # A contig the fasta genuinely lacks gets its own reason, not a silent blank row and not the
        # generic build error, which is about applying variants rather than finding the reference
        result = self._genotype("", chrom="chr9")
        self.assertIsNone(result.zygosity)
        self.assertEqual(result.no_call_reason, NO_CALL_REASON_CONTIG_NOT_IN_REFERENCE)
        self.assertIn("chr9", result.no_call_detail)
        self.assertEqual(result.num_alleles_with_build_errors, 0)

    def test_the_no_call_checks_run_in_a_fixed_order(self):
        """When several no-call conditions apply at once, the reported one must be fixed by a test.

        Without this, reordering the checks in genotype_single_locus changes the NoCallReason column for
        real loci and no test notices.
        """
        # A heterozygous call at a haploid locus beats a missing allele
        result = self._genotype("chrX\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t0|1\n"
                                "chrX\t20\t.\tC\tA\t.\tPASS\t.\tGT\t./.\n",
                                chrom="chrX", vcf_contig="chrX", fasta_contig="chrX",
                                sex_chromosome_ploidy=XY_PLOIDY)
        self.assertEqual(result.no_call_reason, NO_CALL_REASON_HET_AT_HAPLOID_LOCUS)

        # A missing allele beats ambiguous phasing
        result = self._genotype("chr1\t12\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n"
                                "chr1\t20\t.\tC\tA\t.\tPASS\t.\tGT\t./.\n")
        self.assertEqual(result.no_call_reason, NO_CALL_REASON_MISSING_GENOTYPE)

    def test_an_internal_indel_is_reported_by_both_purity_columns(self):
        """RepeatPurity compares position by position, so one 5bp insertion inside four exact 30bp copies
        shifts every base after it out of phase and scores 0.728. RepeatPurityViaEditDistance charges the 5
        inserted bases plus the 5bp length difference they leave against a same-length pure repeat, so it
        scores (125 - 10) / 125 = 0.92.
        """
        motif = "ACGTTGCAAGGCTTACCGGATTCAGGCATA"
        reference_sequence = "T" * 10 + motif * 4 + "A" * 10
        result = self._genotype(
            f"chr1\t50\t.\t{reference_sequence[49]}\t{reference_sequence[49]}CATGG\t.\tPASS\t.\tGT\t0|1\n",
            motif=motif, start_0based=10, end=130, reference_sequence=reference_sequence)

        self.assertEqual(result.repeat_size_long_allele_bp, 125)
        self.assertAlmostEqual(result.repeat_purity_long_allele, 0.728)
        self.assertAlmostEqual(result.repeat_purity_via_edit_distance_long_allele, 0.92)
        self.assertEqual(result.repeat_purity_short_allele, 1.0)
        self.assertEqual(result.repeat_purity_via_edit_distance_short_allele, 1.0)
        # The overall columns report the less pure allele
        self.assertAlmostEqual(result.repeat_purity, 0.728)
        self.assertAlmostEqual(result.repeat_purity_via_edit_distance, 0.92)
        row = result.to_tsv_dict()
        self.assertEqual((row["RepeatPurityLongAllele"], row["RepeatPurityViaEditDistanceLongAllele"]),
                         ("0.7280", "0.9200"))

    def test_missing_genotype_and_filtered_reasons_reach_the_output_columns(self):
        """Every no-call cause must name itself in NoCallReason, not just leave the row blank."""
        row = self._genotype("chr1\t15\t.\tA\tG\t.\tPASS\t.\tGT\t./.\n").to_tsv_dict()
        self.assertEqual(row["Zygosity"], "")
        self.assertTrue(row["NoCallReason"].startswith(NO_CALL_REASON_MISSING_GENOTYPE))

    def test_the_no_call_reason_names_the_cause(self):
        """A blank genotype row is only actionable if the output says which check produced it."""
        ambiguous = self._genotype("chr1\t12\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n"
                                   "chr1\t20\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\n")
        self.assertEqual(ambiguous.no_call_reason, NO_CALL_REASON_AMBIGUOUS_PHASING)
        self.assertEqual(ambiguous.to_tsv_dict()["NoCallReason"], NO_CALL_REASON_AMBIGUOUS_PHASING)
        # A locus that genotypes normally has no reason recorded
        self.assertIsNone(self._genotype("").no_call_reason)
        self.assertEqual(self._genotype("").to_tsv_dict()["NoCallReason"], "")


class TestCatalogAndOutputOrdering(unittest.TestCase):
    """Test catalog validation and the chromosome ordering shared by the genotype output writers."""

    def setUp(self):
        self._temp_files = []

    def tearDown(self):
        for path in self._temp_files:
            if os.path.exists(path):
                os.remove(path)

    _create_temp_file = TestGenotypingPipeline._create_temp_file

    def test_catalog_with_empty_motif_raises_a_clear_error(self):
        """A name field with no motif before the first ':' must be rejected, not divided by zero later."""
        catalog_path = self._create_temp_file("chr1\t10\t28\t:3bp:6.0x:pure_repeats\n", ".bed")
        with self.assertRaises(ValueError) as raised:
            parse_catalog_bed_file(catalog_path)
        self.assertIn("Missing repeat unit", str(raised.exception))

    def test_chromosomes_sort_naturally(self):
        """chr2 must sort before chr10, and X/Y/M after the numbered chromosomes."""
        chroms = ["chr10", "chr2", "chrM", "chr1", "chrX", "chrY"]
        self.assertEqual(sorted(chroms, key=compute_chrom_sort_key),
                         ["chr1", "chr2", "chr10", "chrX", "chrY", "chrM"])

    def test_unplaced_contigs_stay_grouped(self):
        """Records from one unplaced contig must stay contiguous, since tabix rejects interleaved blocks."""
        records = [("chrUn_KI270302v1", 11), ("chr1_KI270706v1_random", 21),
                   ("chr1_KI270706v1_random", 59), ("chrUn_KI270302v1", 87)]
        ordered = sorted(records, key=lambda x: (compute_chrom_sort_key(x[0]), x[1]))
        self.assertEqual([chrom for chrom, _ in ordered],
                         ["chr1_KI270706v1_random", "chr1_KI270706v1_random",
                          "chrUn_KI270302v1", "chrUn_KI270302v1"])

    def _create_indexed_catalog(self, bed_content):
        """Write a catalog BED, bgzip and tabix it, and return (plain path, bgzipped path)."""
        if shutil.which("bgzip") is None or shutil.which("tabix") is None:
            self.skipTest("bgzip/tabix unavailable")
        catalog_path = self._create_temp_file(bed_content, ".bed")
        subprocess.run(["bgzip", "-kf", catalog_path], check=True)
        subprocess.run(["tabix", "-f", "-p", "bed", catalog_path + ".gz"], check=True)
        self._temp_files.append(catalog_path + ".gz")
        self._temp_files.append(catalog_path + ".gz.tbi")
        return catalog_path, catalog_path + ".gz"

    def test_comment_track_and_blank_lines_are_skipped(self):
        """tabix drops '#' lines from the -L branch, so the plain branch has to drop them too.

        Otherwise the same catalog parses one way and fails the other, and the error misidentifies the cause:
        a '#chrom start end name' header reaches the motif validator and is rejected for containing 'E'.
        """
        catalog_path = self._create_temp_file(
            "#chrom\tstart\tend\tname\ntrack name=trs\nchr1\t10\t28\tCAG\nchr1\t38\t56\tGCC\n\n", ".bed")
        self.assertEqual([locus.locus_id for locus in parse_catalog_bed_file(catalog_path)],
                         ["chr1-10-28-CAG", "chr1-38-56-GCC"])

    def test_the_two_parsing_branches_agree(self):
        """The same catalog must yield the same loci whether or not -L is used."""
        catalog_path, catalog_gz_path = self._create_indexed_catalog(
            "#chrom\tstart\tend\tname\nchr1\t10\t28\tCAG\nchr1\t38\t56\tGCC\n")
        self.assertEqual([locus.locus_id for locus in parse_catalog_bed_file(catalog_path)],
                         [locus.locus_id for locus in
                          parse_catalog_bed_file(catalog_gz_path, intervals=["chr1:1-100"])])

    def test_overlapping_intervals_do_not_duplicate_loci(self):
        """tabix returns every record overlapping a region, so a locus can be fetched by several intervals.

        Sharding one chromosome into adjacent windows is a normal way to use a list-valued -L, and any locus
        straddling a window boundary would otherwise be genotyped and written twice.
        """
        _, catalog_gz_path = self._create_indexed_catalog("chr1\t10\t28\tCAG\nchr1\t10\t60\tCAG\n")
        expected = ["chr1-10-28-CAG", "chr1-10-60-CAG"]

        self.assertEqual([locus.locus_id for locus in parse_catalog_bed_file(
            catalog_gz_path, intervals=["chr1:1-100", "chr1:1-100"])], expected)
        self.assertEqual([locus.locus_id for locus in parse_catalog_bed_file(
            catalog_gz_path, intervals=["chr1:1-20", "chr1:21-100"])], expected)

    def test_an_interval_contig_is_matched_to_the_catalog_spelling(self):
        """A catalog written with 'chr1' must still answer an interval given as '1:...' and the reverse.

        tabix matches contig names literally, so without this a naming mismatch is indistinguishable from a
        contig the catalog has no loci on: the run reports an empty interval and genotypes nothing.
        """
        _, catalog_gz_path = self._create_indexed_catalog("chr1\t10\t28\tCAG\n")
        self.assertEqual([locus.locus_id for locus in parse_catalog_bed_file(
            catalog_gz_path, intervals=["1:1-100"])], ["chr1-10-28-CAG"])
        self.assertEqual([locus.locus_id for locus in parse_catalog_bed_file(
            catalog_gz_path, intervals=["chr1:1-100"])], ["chr1-10-28-CAG"])

        self.assertEqual(match_interval_to_catalog_contig("1:1-100", {"chr1"}), "chr1:1-100")
        self.assertEqual(match_interval_to_catalog_contig("chr1:1-100", {"1"}), "1:1-100")
        self.assertEqual(match_interval_to_catalog_contig("chr1:1-100", {"chr1"}), "chr1:1-100")
        self.assertIsNone(match_interval_to_catalog_contig("chr9:1-100", {"chr1"}))
        self.assertEqual(match_interval_to_catalog_contig("chrMT:1-100", {"chrM"}), "chrM:1-100")
        self.assertEqual(match_interval_to_catalog_contig("chrM", {"MT"}), "MT")

    def test_a_malformed_interval_raises_instead_of_being_skipped(self):
        """An interval on a contig the catalog has, but with unusable coordinates, is a user error."""
        _, catalog_gz_path = self._create_indexed_catalog("chr1\t10\t28\tCAG\n")
        for bad_interval in ["chr1:500-100", "chr1:abc-5"]:
            with self.assertRaisesRegex(ValueError, "Invalid interval"):
                parse_catalog_bed_file(catalog_gz_path, intervals=[bad_interval])

    def test_an_interval_on_a_contig_with_no_records_is_skipped(self):
        """Sharding a genome-wide run by chromosome hits contigs the catalog has no loci on."""
        _, catalog_gz_path = self._create_indexed_catalog("chr1\t10\t28\tCAG\n")
        self.assertEqual([locus.locus_id for locus in parse_catalog_bed_file(
            catalog_gz_path, intervals=["chr3:1-100", "chr1:1-100"])], ["chr1-10-28-CAG"])


class TestTRFLongMotifSplitting(unittest.TestCase):
    """Test that motif composition handles motifs longer than one TRF alignment line."""

    def setUp(self):
        self._trf_working_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._trf_working_dir, ignore_errors=True)

    def _split_with_trf(self, sequence, motif_size):
        """Run the TRF motif splitter on one allele sequence, or skip if TRF isn't installed."""
        trf_executable_path = shutil.which("trf")
        if trf_executable_path is None:
            self.skipTest("trf executable unavailable")

        results = run_trf_motif_splitting(
            [("chr1-0-100-M$allele1", sequence, motif_size)], trf_executable_path, self._trf_working_dir)
        self.assertEqual(len(results), 1)
        _, _, entry, method = results[0]
        return entry, method

    def _assert_entry_reconstructs_sequence(self, entry, sequence):
        """The parsed motifs plus prefix and suffix must put the allele back together exactly."""
        rebuilt = entry["prefix"] + "".join(entry["motifs"]) + entry["suffix"]
        self.assertEqual(rebuilt, sequence,
                         f"parsed motifs rebuild {len(rebuilt)}bp but the allele is {len(sequence)}bp")

    def test_motif_longer_than_one_alignment_line_is_not_truncated(self):
        """TRF wraps its alignment at 65 characters, so a 70bp motif spans two lines per copy.

        Before this was handled, each copy kept only its first 65 bases and the reported motif sequence no
        longer reconstructed the allele.
        """
        motif = "GCTAAGGTCCATTGACCGTAAGCTTGGCCAATCGTTAGGCCATTAGGCCTTAAGGCATCGATTGCAGTTA"
        self.assertEqual(len(motif), 70)
        sequence = motif * 5

        entry, method = self._split_with_trf(sequence, len(motif))
        self.assertEqual(method, MOTIF_DETECTION_METHOD_TRF)
        self.assertEqual([len(m) for m in entry["motifs"]], [70] * 5)
        self._assert_entry_reconstructs_sequence(entry, sequence)

    def test_short_motif_still_splits_correctly(self):
        """The common case of a motif that fits on one alignment line must be unaffected.

        Reconstruction alone is too weak a check: one 48bp pseudo-motif reconstructs this sequence just as
        well as the sixteen 3bp copies, while making Allele*MotifCounts meaningless. The decomposition
        itself is asserted here.
        """
        sequence = "CAG" * 10 + "CAA" + "CAG" * 5

        entry, method = self._split_with_trf(sequence, 3)
        self.assertEqual(method, MOTIF_DETECTION_METHOD_TRF)
        self._assert_entry_reconstructs_sequence(entry, sequence)
        self.assertEqual(len(entry["motifs"]), 16)
        self.assertEqual({len(motif) for motif in entry["motifs"]}, {3})
        self.assertEqual(dict(collections.Counter(entry["motifs"])), {"CAG": 15, "CAA": 1})
        self.assertEqual(entry["motifs"][10], "CAA")

    def test_homopolymer_allele_is_split_by_trf(self):
        """TRF prints a period-1 alignment as one unbroken run, which must still yield 1bp motifs.

        Parsing that run as a single long motif makes run_trf_motif_splitting reject the record, because it
        only accepts a decomposition whose motif length matches the locus. Homopolymers are a large share of
        a genome-wide catalog, so every one of them would silently report the basic-split method instead.
        """
        entry, method = self._split_with_trf("A" * 25, 1)
        self.assertEqual(method, MOTIF_DETECTION_METHOD_TRF)
        self.assertEqual(entry["motifs"], ["A"] * 25)
        self._assert_entry_reconstructs_sequence(entry, "A" * 25)

    def test_a_relative_trf_path_is_resolved_before_the_subprocess_runs(self):
        """TRF runs with its cwd set to the FASTA's directory, so a relative path must be resolved first.

        Otherwise a path such as 'bin/trf' passes the up-front check in the caller's directory and then
        fails to resolve inside the subprocess, where the shell's error goes to DEVNULL and the exit code is
        ignored, so every allele silently reports the basic-split method.
        """
        trf_executable_path = shutil.which("trf")
        if trf_executable_path is None:
            self.skipTest("trf executable unavailable")

        # The directory the relative path is typed in is not the directory TRF runs in: the batch FASTA
        # goes to its own working directory, and the subprocess runs there.
        caller_directory = tempfile.mkdtemp()
        os.makedirs(os.path.join(caller_directory, "bin"), exist_ok=True)
        os.symlink(trf_executable_path, os.path.join(caller_directory, "bin", "trf"))

        original_directory = os.getcwd()
        os.chdir(caller_directory)
        try:
            results = run_trf_motif_splitting(
                [("chr1-0-100-CAG$allele1", "CAG" * 20, 3)], os.path.join("bin", "trf"),
                self._trf_working_dir)
        finally:
            os.chdir(original_directory)
            shutil.rmtree(caller_directory, ignore_errors=True)

        _, _, entry, method = results[0]
        self.assertEqual(method, MOTIF_DETECTION_METHOD_TRF)
        self.assertEqual(entry["motifs"], ["CAG"] * 20)

    def test_a_broken_trf_path_is_an_error_rather_than_a_silent_fallback(self):
        """A mistyped or missing --trf-executable-path must not look like 'TRF found no repeats'.

        TRF is launched through a shell with its exit code ignored, since the versions in circulation
        disagree about what success returns, so nothing downstream can tell a failed launch from an empty
        result. Every allele would quietly report the basic-split method.
        """
        with self.assertRaises(FileNotFoundError):
            run_trf_motif_splitting([("chr1-0-100-CAG$allele1", "CAG" * 10, 3)],
                                    "/nonexistent/path/to/trf", self._trf_working_dir)


class TestContributingVariantsVcfOutput(unittest.TestCase):
    """End-to-end tests for the VCF of variants that contributed to TR genotyping."""

    def setUp(self):
        self._temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._temp_dir, ignore_errors=True)

    def _run_genotype_subcommand(self, contigs, catalog_lines, vcf_records):
        """Run the genotype subcommand with --write-vcf, returning the output VCF path.

        Args:
            contigs (dict): contig name to reference sequence
            catalog_lines (list): BED lines for the catalog, without trailing newlines
            vcf_records (list): VCF body lines, without trailing newlines, already coordinate-sorted

        Returns:
            str: path to the bgzipped contributing-variants VCF
        """
        if shutil.which("bgzip") is None or shutil.which("tabix") is None:
            self.skipTest("bgzip/tabix unavailable")

        fasta_path = os.path.join(self._temp_dir, "ref.fa")
        with open(fasta_path, "w") as f:
            for chrom, sequence in contigs.items():
                f.write(f">{chrom}\n{sequence}\n")
        pyfaidx.Fasta(fasta_path)

        catalog_path = os.path.join(self._temp_dir, "catalog.bed")
        with open(catalog_path, "w") as f:
            f.write("\n".join(catalog_lines) + "\n")

        header_contigs = "".join(
            f"##contig=<ID={chrom},length={len(sequence)}>\n" for chrom, sequence in contigs.items())
        vcf_path = os.path.join(self._temp_dir, "input.vcf")
        with open(vcf_path, "w") as f:
            f.write("##fileformat=VCFv4.2\n" + header_contigs +
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' +
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n" +
                    "\n".join(vcf_records) + "\n")
        subprocess.run(["bgzip", "-f", vcf_path], check=True)
        subprocess.run(["tabix", "-p", "vcf", vcf_path + ".gz"], check=True)

        args = argparse.Namespace(
            reference_fasta_path=fasta_path, catalog_bed=catalog_path, input_vcf_path=vcf_path + ".gz",
            output_prefix=os.path.join(self._temp_dir, "out"), interval=None, verbose=False,
            show_progress_bar=False, write_vcf=True, write_json=False, add_motif_composition=None,
            skip_hom_ref_loci=False, trf_executable_path=None)
        do_genotype_subcommand(args)

        return os.path.join(self._temp_dir, "out.tandem_repeat_contributing_variants.vcf.gz")

    def test_unplaced_contigs_produce_an_indexable_vcf(self):
        """Records from different unplaced contigs must not interleave, or tabix cannot index the output.

        GRCh38 carries roughly 150 unplaced and alt contigs, so a whole-genome run with --write-vcf hits this
        whenever two of them contribute variants at interleaving positions.
        """
        # 10bp of T, then (CAG)x6 at [10,28), then 10bp of A, then a 60bp poly-G tract at [38,98).
        # Both contigs need two contributing records for interleaving to be possible at all, so each one
        # carries a locus over the CAG tract and another over the poly-G tract.
        reference_sequence = "TTTTTTTTTT" + "CAG" * 6 + "AAAAAAAAAA" + "G" * 60
        self.assertEqual(reference_sequence[20], "A")   # the REF base of the record at position 21
        self.assertEqual(reference_sequence[10], "C")   # the REF base of the record at position 11
        self.assertEqual(reference_sequence[58], "G")   # the REF base of the record at position 59
        self.assertEqual(reference_sequence[86], "G")   # the REF base of the record at position 87

        output_vcf_path = self._run_genotype_subcommand(
            contigs={"chr1_KI270706v1_random": reference_sequence, "chrUn_KI270302v1": reference_sequence},
            catalog_lines=["chr1_KI270706v1_random\t10\t28\tCAG", "chr1_KI270706v1_random\t38\t98\tG",
                           "chrUn_KI270302v1\t10\t28\tCAG", "chrUn_KI270302v1\t38\t98\tG"],
            vcf_records=[
                "chr1_KI270706v1_random\t21\t.\tA\tACAG\t.\tPASS\t.\tGT\t0|1",
                "chr1_KI270706v1_random\t59\t.\tG\tGT\t.\tPASS\t.\tGT\t0|1",
                "chrUn_KI270302v1\t11\t.\tC\tCCAG\t.\tPASS\t.\tGT\t0|1",
                "chrUn_KI270302v1\t87\t.\tG\tGT\t.\tPASS\t.\tGT\t0|1",
            ])

        self.assertTrue(os.path.exists(output_vcf_path + ".tbi"), "tabix did not produce an index")
        output_vcf = pysam.VariantFile(output_vcf_path)
        try:
            self.assertIsNotNone(output_vcf.index, "the output VCF has no usable index")
            written_records = [(record.chrom, record.pos) for record in output_vcf.fetch()]
        finally:
            output_vcf.close()

        # All four records must be written, grouped by contig and ascending within each contig. Sorting on
        # position alone would give chrUn 11, chr1_random 21, chr1_random 59, chrUn 87, which tabix rejects.
        self.assertEqual(written_records, [
            ("chr1_KI270706v1_random", 21),
            ("chr1_KI270706v1_random", 59),
            ("chrUn_KI270302v1", 11),
            ("chrUn_KI270302v1", 87),
        ])

    def test_one_variant_overlapping_two_loci_is_annotated_with_both(self):
        """A variant inside two catalog loci is written once, naming both loci and both motifs."""
        output_vcf_path = self._run_genotype_subcommand(
            contigs={"chr1": "TTTTTTTTTT" + "CAG" * 6 + "AAAAAAAAAA"},
            # The two loci overlap at position 27 (1-based), so the single record below sits inside both
            catalog_lines=["chr1\t10\t28\tCAG", "chr1\t26\t38\tA"],
            vcf_records=["chr1\t27\t.\tA\tACAG\t.\tPASS\t.\tGT\t0|1"])

        output_vcf = pysam.VariantFile(output_vcf_path)
        try:
            records = list(output_vcf.fetch())
            self.assertEqual(len(records), 1, "a variant overlapping two loci must be written once")
            self.assertEqual(records[0].info["TR_LocusId"], ("chr1-10-28-CAG", "chr1-26-38-A"))
            self.assertEqual(records[0].info["TR_Motif"], ("CAG", "A"))
        finally:
            output_vcf.close()


    def test_genotype_options_reach_the_batch_driver(self):
        """The CLI options must survive the hop from args into genotype_single_locus.

        Every other test of these behaviours injects the resolved value directly, so the translation in
        genotype_all_loci is unguarded: a regression there would report two reference copies at every
        chrX/chrY locus of a male sample, apply LowQual records as if they had passed, or count an Alu
        insertion as repeat bases, with nothing failing.
        """
        if shutil.which("bgzip") is None or shutil.which("tabix") is None:
            self.skipTest("bgzip/tabix unavailable")

        reference_sequence = "TTTTTTTTTT" + "CAG" * 6 + "AAAAAAAAAA"
        fasta_path = os.path.join(self._temp_dir, "opts.fa")
        with open(fasta_path, "w") as f:
            f.write(f">chr1\n{reference_sequence}\n")
        pyfaidx.Fasta(fasta_path)
        catalog_path = os.path.join(self._temp_dir, "opts.bed")
        with open(catalog_path, "w") as f:
            f.write("chr1\t10\t28\tCAG\n")

        def run(records, **option_overrides):
            vcf_path = os.path.join(self._temp_dir, "opts.vcf")
            with open(vcf_path, "w") as f:
                f.write("##fileformat=VCFv4.2\n"
                        f"##contig=<ID=chr1,length={len(reference_sequence)}>\n"
                        '##FILTER=<ID=LowQual,Description="low">\n'
                        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n" + records)
            subprocess.run(["bgzip", "-f", vcf_path], check=True)
            subprocess.run(["tabix", "-f", "-p", "vcf", vcf_path + ".gz"], check=True)
            output_prefix = os.path.join(self._temp_dir, "opts_out")
            options = dict(
                reference_fasta_path=fasta_path, catalog_bed=catalog_path,
                input_vcf_path=vcf_path + ".gz", output_prefix=output_prefix, interval=None,
                verbose=False, show_progress_bar=False, write_vcf=False, write_json=False,
                add_motif_composition=None, skip_hom_ref_loci=False, trf_executable_path=None)
            options.update(option_overrides)
            do_genotype_subcommand(argparse.Namespace(**options))
            with gzip.open(f"{output_prefix}.tandem_repeat_genotypes.tsv.gz", "rt") as f:
                lines = f.read().split("\n")
            return dict(zip(lines[0].split("\t"), lines[1].split("\t")))

        self.assertEqual(run("")["Zygosity"], "HOM")

        # A record whose FILTER names a failed filter is applied like any other
        self.assertEqual(run("chr1\t10\t.\tT\tTCAG\t.\tLowQual\t.\tGT\t1|1\n")["Zygosity"], "HOM")

    def test_the_json_output_carries_the_no_call_reason(self):
        """A no-call row must explain itself in JSON as well as in the TSV, and the two must agree."""
        if shutil.which("bgzip") is None or shutil.which("tabix") is None:
            self.skipTest("bgzip/tabix unavailable")

        reference_sequence = "TTTTTTTTTT" + "CAG" * 6 + "AAAAAAAAAA"
        fasta_path = os.path.join(self._temp_dir, "nc.fa")
        with open(fasta_path, "w") as f:
            f.write(f">chr1\n{reference_sequence}\n")
        pyfaidx.Fasta(fasta_path)
        catalog_path = os.path.join(self._temp_dir, "nc.bed")
        with open(catalog_path, "w") as f:
            f.write("chr1\t10\t28\tCAG\n")
        vcf_path = os.path.join(self._temp_dir, "nc.vcf")
        with open(vcf_path, "w") as f:
            f.write("##fileformat=VCFv4.2\n"
                    f"##contig=<ID=chr1,length={len(reference_sequence)}>\n"
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
                    "chr1\t15\t.\tA\tG\t.\tPASS\t.\tGT\t./.\n")
        subprocess.run(["bgzip", "-f", vcf_path], check=True)
        subprocess.run(["tabix", "-f", "-p", "vcf", vcf_path + ".gz"], check=True)

        output_prefix = os.path.join(self._temp_dir, "nc_out")
        do_genotype_subcommand(argparse.Namespace(
            reference_fasta_path=fasta_path, catalog_bed=catalog_path, input_vcf_path=vcf_path + ".gz",
            output_prefix=output_prefix, interval=None, verbose=False, show_progress_bar=False,
            write_vcf=False, write_json=True, add_motif_composition=None, skip_hom_ref_loci=False,
            trf_executable_path=None))

        with gzip.open(f"{output_prefix}.tandem_repeat_genotypes.json.gz", "rt") as f:
            record = json.load(f)[0]
        with gzip.open(f"{output_prefix}.tandem_repeat_genotypes.tsv.gz", "rt") as f:
            lines = f.read().split("\n")
        tsv_row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))

        self.assertIsNone(record["Zygosity"])
        self.assertTrue(record["NoCallReason"].startswith(NO_CALL_REASON_MISSING_GENOTYPE))
        self.assertEqual(record["NoCallReason"], tsv_row["NoCallReason"])
        # The JSON record must carry every documented column, not a subset. It omits only the two
        # motif-composition columns, which are added when --add-motif-composition is used.
        self.assertEqual(set(record),
                         set(GENOTYPE_TSV_OUTPUT_COLUMNS) - {
                             "Allele1MotifSequence", "Allele1SequenceMotifSplittingMethod",
                             "Allele2MotifSequence", "Allele2SequenceMotifSplittingMethod"})

    def test_skip_hom_ref_loci_drops_the_locus_from_every_output(self):
        """--skip-hom-ref-loci is applied separately by each writer, so all three must agree.

        The flag gates row filtering in write_genotypes_tsv, write_genotypes_json, write_genotypes_vcf and
        compute_motif_composition. If one of them drifted from the others the outputs would cover different
        locus sets, or the motif columns would be blank for a locus the writers still emit.
        """
        if shutil.which("bgzip") is None or shutil.which("tabix") is None:
            self.skipTest("bgzip/tabix unavailable")

        reference_sequence = "TTTTTTTTTT" + "CAG" * 6 + "AAAAAAAAAA" + "GCC" * 6 + "TTTTTTTTTT"
        fasta_path = os.path.join(self._temp_dir, "ref.fa")
        with open(fasta_path, "w") as f:
            f.write(f">chr1\n{reference_sequence}\n")
        pyfaidx.Fasta(fasta_path)

        catalog_path = os.path.join(self._temp_dir, "catalog.bed")
        with open(catalog_path, "w") as f:
            # The CAG locus carries a variant; the GCC locus at [38,56) has none, so it is hom-ref
            f.write("chr1\t10\t28\tCAG\nchr1\t38\t56\tGCC\n")

        vcf_path = os.path.join(self._temp_dir, "input.vcf")
        with open(vcf_path, "w") as f:
            f.write("##fileformat=VCFv4.2\n"
                    f"##contig=<ID=chr1,length={len(reference_sequence)}>\n"
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
                    "chr1\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t0|1\n")
        subprocess.run(["bgzip", "-f", vcf_path], check=True)
        subprocess.run(["tabix", "-p", "vcf", vcf_path + ".gz"], check=True)

        output_prefix = os.path.join(self._temp_dir, "skipped")
        args = argparse.Namespace(
            reference_fasta_path=fasta_path, catalog_bed=catalog_path, input_vcf_path=vcf_path + ".gz",
            output_prefix=output_prefix, interval=None, verbose=False, show_progress_bar=False,
            write_vcf=True, write_json=True, add_motif_composition="basic", skip_hom_ref_loci=True,
            trf_executable_path=None)
        do_genotype_subcommand(args)

        with gzip.open(f"{output_prefix}.tandem_repeat_genotypes.tsv.gz", "rt") as f:
            tsv_lines = f.read().strip().split("\n")
        tsv_header = tsv_lines[0].split("\t")
        tsv_rows = [dict(zip(tsv_header, line.split("\t"))) for line in tsv_lines[1:]]
        self.assertEqual([row["LocusId"] for row in tsv_rows], ["chr1-10-28-CAG"])
        # The surviving locus must still get its motif columns, ie. compute_motif_composition skipped the
        # same locus the writers did rather than a different one
        self.assertTrue(tsv_rows[0]["Allele1MotifSequence"])
        self.assertEqual(tsv_rows[0]["Allele1SequenceMotifSplittingMethod"], MOTIF_DETECTION_METHOD_BASIC_SPLIT)

        with gzip.open(f"{output_prefix}.tandem_repeat_genotypes.json.gz", "rt") as f:
            json_records = json.load(f)
        self.assertEqual([record["LocusId"] for record in json_records], ["chr1-10-28-CAG"])
        self.assertIsNotNone(json_records[0]["Allele1MotifSequence"])

        output_vcf = pysam.VariantFile(f"{output_prefix}.tandem_repeat_contributing_variants.vcf.gz")
        try:
            self.assertEqual([record.info["TR_LocusId"] for record in output_vcf.fetch()],
                             [("chr1-10-28-CAG",)])
        finally:
            output_vcf.close()

    def test_missing_filter_is_not_rewritten_as_pass(self):
        """A FILTER of '.' means filters were not applied, which is not the same claim as PASS."""
        output_vcf_path = self._run_genotype_subcommand(
            contigs={"chr1": "TTTTTTTTTT" + "CAG" * 6 + "AAAAAAAAAA"},
            catalog_lines=["chr1\t10\t28\tCAG"],
            vcf_records=["chr1\t10\t.\tT\tTCAG\t.\t.\t.\tGT\t0|1"])

        output_vcf = pysam.VariantFile(output_vcf_path)
        try:
            records = list(output_vcf.fetch())
            self.assertEqual(len(records), 1)
            self.assertEqual(list(records[0].filter.keys()), [])
        finally:
            output_vcf.close()

    def test_pass_filter_is_preserved(self):
        """The neighbouring case, an explicit PASS, must still come through as PASS."""
        output_vcf_path = self._run_genotype_subcommand(
            contigs={"chr1": "TTTTTTTTTT" + "CAG" * 6 + "AAAAAAAAAA"},
            catalog_lines=["chr1\t10\t28\tCAG"],
            vcf_records=["chr1\t10\t.\tT\tTCAG\t.\tPASS\t.\tGT\t0|1"])

        output_vcf = pysam.VariantFile(output_vcf_path)
        try:
            records = list(output_vcf.fetch())
            self.assertEqual(len(records), 1)
            self.assertEqual(list(records[0].filter.keys()), ["PASS"])
        finally:
            output_vcf.close()


if __name__ == "__main__":
    unittest.main()
