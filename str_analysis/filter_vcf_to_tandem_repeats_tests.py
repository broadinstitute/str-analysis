#!/usr/bin/env python3

"""Test script for the filter_vcf_to_catalog_tandem_repeats.py script."""

import argparse
import collections
import unittest

import pyfaidx

from str_analysis.filter_vcf_to_tandem_repeats import Allele, TandemRepeatAllele, MinimalTandemRepeatAllele, DETECTION_MODE_PURE_REPEATS, \
    merge_overlapping_tandem_repeat_loci, detect_perfect_and_almost_perfect_tandem_repeats


class TestAllele(unittest.TestCase):
    """Test the Allele class."""


    def setUp(self):
        """Set up the test case."""
        self._fasta_obj = pyfaidx.Fasta("str_analysis/data/tests/chr22_11Mb.fa.gz", one_based_attributes=False, as_raw=True)
        self._fasta_obj_with_one_based_attributes = pyfaidx.Fasta("str_analysis/data/tests/chr22_11Mb.fa.gz", one_based_attributes=True, as_raw=True)
        self._poly_a_insertion1 = Allele("chr22", 10513201, "T", "TAAAAA", self._fasta_obj)
        self._poly_a_insertion2 = Allele("chr22", 10513201, "T", "TAAAA", self._fasta_obj)
        self._poly_a_deletion1 = Allele("chr22", 10513201, "TAAAAA", "T", self._fasta_obj)
        self._poly_a_deletion2 = Allele("chr22", 10513201, "TAAAA", "T", self._fasta_obj)

        self._AAGA_insertion = Allele("chr22", 10515040, "T", "TAAGA", self._fasta_obj)
        self._AAGA_deletion = Allele("chr22", 10515040, "TAAGA", "T", self._fasta_obj)

        self._AAGA_tandem_repeat_expansion = TandemRepeatAllele(self._AAGA_insertion, "AAGA", 0, 4, 37, DETECTION_MODE_PURE_REPEATS)
        self._AAGA_tandem_repeat_contraction = TandemRepeatAllele(self._AAGA_deletion, "AAGA", 0, 4, 33, DETECTION_MODE_PURE_REPEATS)


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
                self.assertEqual(len(allele.get_left_flanking_sequence()), expected_length)
                self.assertEqual(len(allele.get_right_flanking_sequence()), expected_length)

                self.assertEqual(allele.get_right_flank_end(), allele.get_right_flank_start_0based() + expected_length)
                self.assertEqual(allele.get_left_flank_end(), allele.get_left_flank_start_0based() + expected_length)
            

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
        

    def test_detect_perfect_and_almost_perfect_tandem_repeats(self):
        """Test the detect_perfect_and_almost_perfect_tandem_repeats function."""
        counters = collections.defaultdict(int)

        args = argparse.Namespace(

            min_repeat_unit_length=1,
            max_repeat_unit_length=1000,
            show_progress_bar=False,
            min_repeats=3,
            min_tandem_repeat_length=9,
            dont_run_trf=False,
            debug=False,

        )
        results, alleles_to_process_next_using_trf = detect_perfect_and_almost_perfect_tandem_repeats( [
            Allele("chr22", 10689286, "A", "ACAGCAGCAGCAGCAG", self._fasta_obj),
        ], counters, args)
        self.assertEqual(alleles_to_process_next_using_trf, [])
        self.assertEqual(len(results), 1, msg=f"Expected 1 result, got {len(results)}: {results}")

        # TODO add more tests here

    def test_merge_overlapping_tandem_repeat_loci(self):
        """Test the merge_overlapping_tandem_repeat_loci function."""

        tandem_repeat_alleles = [
            MinimalTandemRepeatAllele("chr17", 3, 2, "CAG", None),
            MinimalTandemRepeatAllele("chr17", 2, 5, "CAG", None),
            MinimalTandemRepeatAllele("chr17", 5, 14, "CAG", None),
            MinimalTandemRepeatAllele("chr22", 10515040, 10515077, "AAGA", DETECTION_MODE_PURE_REPEATS),
            MinimalTandemRepeatAllele("chr1", 1929384, 1929384, "AGGGTAGGGAGGGAGGGAGAGGAGGGGAGAGGGTAGGGAGGGAGAGGAGGGGGAGGGAGGGAGGGGAGGGAGGGGAG", DETECTION_MODE_PURE_REPEATS),
            MinimalTandemRepeatAllele("chr1", 1929384, 1929490, "AGGGTAGGGAGGGAGGGAGAGGAGGGGAGAGGGTAGGGAGGGAGAGGAGGGGGAGGGAGGGAGGGGAGGGAGGGGAG", DETECTION_MODE_PURE_REPEATS),
        ]
        for n in range(2, 10):
            # [tandem_repeat_alleles[0]] * n, [tandem_repeat_alleles[1]] * n, [tandem_repeat_alleles[2]] * n, 
            for tandem_repeat_alleles_to_merge, exepected_n_results in [
                (tandem_repeat_alleles, 3), 
                (tandem_repeat_alleles * n, 3), 
                ([tandem_repeat_alleles[0]] * n, 1),
                ([tandem_repeat_alleles[1]] * n, 1),
                ([tandem_repeat_alleles[2]] * n, 1),
                ([tandem_repeat_alleles[3]] * n, 1),
                ([tandem_repeat_alleles[4]] * n, 1),
                ([tandem_repeat_alleles[5]] * n, 1),
            ]:
                merged_tandem_repeat_alleles = merge_overlapping_tandem_repeat_loci(tandem_repeat_alleles_to_merge)
                self.assertEqual(len(merged_tandem_repeat_alleles), exepected_n_results, 
                                 msg=f"Expected {exepected_n_results} results, got {len(merged_tandem_repeat_alleles)} when merging "
                                     f"{len(tandem_repeat_alleles_to_merge)} tandem repeat alleles for n={n}: {tandem_repeat_alleles_to_merge}")

    
    def tearDown(self):
        """Tear down the test case."""
        self._fasta_obj.close()


