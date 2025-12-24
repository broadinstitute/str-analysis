#!/usr/bin/env python3

"""Test script for the filter_vcf_to_catalog_tandem_repeats.py script."""

import argparse
import collections
import pkgutil
import tempfile
import unittest

import pyfaidx

from str_analysis.filter_vcf_to_tandem_repeats import Allele, TandemRepeatAllele, MinimalTandemRepeatAllele, \
    DETECTION_MODE_PURE_REPEATS, DETECTION_MODE_ALLOW_INTERRUPTIONS, DETECTION_MODE_TRF, \
    merge_overlapping_tandem_repeat_loci, detect_perfect_and_almost_perfect_tandem_repeats, detect_tandem_repeats_using_trf


class TestAllele(unittest.TestCase):
    """Test the Allele class."""


    def setUp(self):
        """Set up the test case."""
        # create file-like object around a binary string
        fasta_data = pkgutil.get_data("str_analysis", "data/tests/chr22_11Mb.fa.gz")
        # write it to a temporary named file
        with tempfile.NamedTemporaryFile(suffix=".fa.gz") as fasta_file:
            fasta_file.write(fasta_data)
            fasta_file.flush()
            self._fasta_obj = pyfaidx.Fasta(fasta_file.name, one_based_attributes=False, as_raw=True)
            self._fasta_obj_with_one_based_attributes = pyfaidx.Fasta(fasta_file.name, one_based_attributes=True, as_raw=True)
        self._poly_a_insertion1 = Allele("chr22", 10513201, "T", "TAAAAA", self._fasta_obj)
        self._poly_a_insertion2 = Allele("chr22", 10513201, "T", "TAAAA", self._fasta_obj)
        self._poly_a_deletion1 = Allele("chr22", 10513201, "TAAAAA", "T", self._fasta_obj)
        self._poly_a_deletion2 = Allele("chr22", 10513201, "TAAAA", "T", self._fasta_obj)

        self._AAGA_insertion = Allele("chr22", 10515040, "T", "TAAGA", self._fasta_obj)
        self._AAGA_deletion = Allele("chr22", 10515040, "TAAGA", "T", self._fasta_obj)

        self._AAGA_tandem_repeat_expansion = TandemRepeatAllele(self._AAGA_insertion, "AAGA", 0, 4, 37, DETECTION_MODE_PURE_REPEATS)
        self._AAGA_tandem_repeat_contraction = TandemRepeatAllele(self._AAGA_deletion, "AAGA", 0, 4, 33, DETECTION_MODE_PURE_REPEATS)

        self._default_args = argparse.Namespace(

            min_repeat_unit_length=1,
            max_repeat_unit_length=1000,
            show_progress_bar=False,
            min_repeats=3,
            min_tandem_repeat_length=9,
            debug=False,
            trf_working_dir=tempfile.mkdtemp(),
            input_vcf_prefix="test",
            trf_executable_path="trf",
            trf_threads=2,
            verbose=False,
            allow_multiple_trf_results_per_locus=False,
            dont_allow_interruptions=False,
            dont_run_trf=False,
            min_indel_size_to_run_trf=7,
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
        

    def test_detect_tandem_repeats(self):
        """Test the detect_perfect_and_almost_perfect_tandem_repeats function."""
        counters = collections.defaultdict(int)

        for test_interrupted_repeats in [False, True]:
            alleles = [
                Allele("chr22", 10689286, "A", "ACAGCAGCAGCATCAG" if test_interrupted_repeats else "ACAGCAGCAGCAGCAG", self._fasta_obj, allele_order=1),
                Allele("chr22", 10689288, "T", "TATTATTATTATT", self._fasta_obj, allele_order=2),
                Allele("chr22", 10689285, "T", "TATTATTATTATT", self._fasta_obj, allele_order=3),
                Allele("chr22", 10689282, "T", "TATTATTATTATT", self._fasta_obj, allele_order=4),
                Allele("chr22", 10689267, "T", "TATTATTATTATT", self._fasta_obj, allele_order=5),
            ]
            
            for detect_repeats_using_TRF in [False, ]:
                if test_interrupted_repeats and detect_repeats_using_TRF:
                    continue

                if not detect_repeats_using_TRF:
                    results, alleles_to_process_next_using_trf = detect_perfect_and_almost_perfect_tandem_repeats(alleles, counters, self._default_args)
                else:
                    results = detect_tandem_repeats_using_trf(alleles, counters, self._default_args)
                    alleles_to_process_next_using_trf = []
                    
                self.assertEqual(alleles_to_process_next_using_trf, [])
                self.assertEqual(len(results), len(alleles), msg=f"Expected {len(alleles)} results, got {len(results)}: {results}")
        
                if detect_repeats_using_TRF:
                    expected_detection_mode = DETECTION_MODE_TRF
                elif test_interrupted_repeats:
                    expected_detection_mode = DETECTION_MODE_ALLOW_INTERRUPTIONS
                else:
                    expected_detection_mode = DETECTION_MODE_PURE_REPEATS

                msg = ("pure repeats" if not test_interrupted_repeats else "interrupted repeats") + (f" using TRF" if detect_repeats_using_TRF else " not using TRF")

                self.assertEqual(results[0].chrom, "chr22", msg=msg)
                self.assertEqual(results[0].repeat_unit, "CAG", msg=msg)
                self.assertEqual(results[0].start_0based, 10689286, msg=msg)
                self.assertEqual(results[0].end_1based, 10689286, msg=msg)
                self.assertEqual(results[0].num_repeats_ref, 0, msg=msg)
                self.assertEqual(results[0].num_repeats_alt, 5, msg=msg)
                self.assertEqual(results[0].num_repeats_in_left_flank, 0, msg=msg)
                self.assertEqual(results[0].num_repeats_in_right_flank, 0, msg=msg)
                self.assertEqual(results[0].num_repeats_in_variant, 5, msg=msg)
                self.assertEqual(results[0].num_repeats_in_variant_and_flanks, 5, msg=msg)
                self.assertEqual(results[0].detection_mode, expected_detection_mode, msg=msg)

                for results_index, i in enumerate([0, 1, 2, 7], start=1):
                    msg = (f"using TRF" if detect_repeats_using_TRF else f"not using TRF") + f" for i = {i}, allele {results[results_index].allele} result #{results_index + 1}: {results[results_index]}. Left flank: {results[results_index].allele.get_left_flanking_sequence()}. Variant: {results[results_index].allele.variant_bases}. Right flank: {results[results_index].allele.get_right_flanking_sequence()}"
                    self.assertEqual(results[results_index].chrom, "chr22", msg=msg)
                    self.assertTrue(results[results_index].start_0based == 10689265 or results[results_index].start_0based == 10689267, msg=msg)
                    self.assertEqual(results[results_index].end_1based, 10689288, msg=msg)
                    self.assertEqual(results[results_index].repeat_unit, "ATT", msg=msg)
                    self.assertEqual(results[results_index].num_repeats_ref, 7, msg=msg)
                    self.assertEqual(results[results_index].num_repeats_alt, 11, msg=msg)
                    self.assertEqual(results[results_index].num_repeats_in_left_flank, 7 - i, msg=msg)
                    self.assertEqual(results[results_index].num_repeats_in_right_flank, 0 + i, msg=msg)
                    self.assertEqual(results[results_index].num_repeats_in_variant, 4, msg=msg)
                    self.assertEqual(results[results_index].num_repeats_in_variant_and_flanks, 11, msg=msg)


    def test_merge_overlapping_tandem_repeat_loci(self):
        """Test the merge_overlapping_tandem_repeat_loci function."""

        tandem_repeat_alleles = [
            MinimalTandemRepeatAllele("chr17", 3, 3, "CAG", None, allele_order=1),
            MinimalTandemRepeatAllele("chr17", 2, 5, "CAG", None, allele_order=2),
            MinimalTandemRepeatAllele("chr17", 5, 14, "CAG", None, allele_order=3),
            MinimalTandemRepeatAllele("chr22", 10515040, 10515077, "AAGA", DETECTION_MODE_PURE_REPEATS, allele_order=4),
            MinimalTandemRepeatAllele("chr1", 1929384, 1929384, "AGGGTAGGGAGGGAGGGAGAGGAGGGGAGAGGGTAGGGAGGGAGAGGAGGGGGAGGGAGGGAGGGGAGGGAGGGGAG", DETECTION_MODE_PURE_REPEATS, allele_order=5),
            MinimalTandemRepeatAllele("chr1", 1929384, 1929490, "AGGGTAGGGAGGGAGGGAGAGGAGGGGAGAGGGTAGGGAGGGAGAGGAGGGGGAGGGAGGGAGGGGAGGGAGGGGAG", DETECTION_MODE_PURE_REPEATS, allele_order=6),
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


    def test_merge_overlapping_tandem_repeat_loci2(self):
        """Test the merge_overlapping_tandem_repeat_loci function."""

        tandem_repeat_alleles = [
            MinimalTandemRepeatAllele("chr22", 10515040, 10515077, "CAGCCGCAG", DETECTION_MODE_PURE_REPEATS, allele_order=1),
            MinimalTandemRepeatAllele("chr22", 10515040, 10515083, "CAGCCGCAG", DETECTION_MODE_PURE_REPEATS, allele_order=2),
            MinimalTandemRepeatAllele("chr22", 10515040, 10515080, "CAGCCGCAG", DETECTION_MODE_PURE_REPEATS, allele_order=3),
        ]
        for n in range(1, 5):
            # [tandem_repeat_alleles[0]] * n, [tandem_repeat_alleles[1]] * n, [tandem_repeat_alleles[2]] * n, 
            tandem_repeat_alleles_to_merge = tandem_repeat_alleles * n
        
            merged_tandem_repeat_alleles = merge_overlapping_tandem_repeat_loci(tandem_repeat_alleles_to_merge)
            self.assertEqual(len(merged_tandem_repeat_alleles), 1, 
                                msg=f"Expected 1 result, got {len(merged_tandem_repeat_alleles)} when merging "
                                    f"{len(tandem_repeat_alleles_to_merge)} tandem repeat alleles for n={n}: {tandem_repeat_alleles_to_merge}")

            self.assertEqual(merged_tandem_repeat_alleles[0].locus_id, "chr22-10515040-10515083-CAGCCGCAG")

    def test_merge_overlapping_tandem_repeat_loci3(self):
        """Test the merge_overlapping_tandem_repeat_loci function."""

        for n in range(0, 5):
            # [tandem_repeat_alleles[0]] * n, [tandem_repeat_alleles[1]] * n, [tandem_repeat_alleles[2]] * n, 
            tandem_repeat_alleles = [
                MinimalTandemRepeatAllele("chr22", 10515040, 10515077, "CAG"+("CCC"*n), DETECTION_MODE_PURE_REPEATS, allele_order=1),
                MinimalTandemRepeatAllele("chr22", 10515040, 10515083, "CCG"+("CCC"*n), DETECTION_MODE_PURE_REPEATS, allele_order=2),
            ]
        
            merged_tandem_repeat_alleles = merge_overlapping_tandem_repeat_loci(tandem_repeat_alleles)
            expected_results = 2 if n < 2 else 1
            self.assertEqual(len(merged_tandem_repeat_alleles), expected_results, 
                                msg=f"Expected {expected_results} results, got {len(merged_tandem_repeat_alleles)} when merging "
                                    f"{len(tandem_repeat_alleles)} tandem repeat alleles for n={n}: {tandem_repeat_alleles}")

            msg = f"n={n}, input alleles={tandem_repeat_alleles}, merged={merged_tandem_repeat_alleles}"
            if expected_results == 2:
                self.assertEqual(merged_tandem_repeat_alleles[0].locus_id, "chr22-10515040-10515077-CAG" + ("CCC"*n), msg=msg)
                self.assertEqual(merged_tandem_repeat_alleles[1].locus_id, "chr22-10515040-10515083-CCG" + ("CCC"*n), msg=msg)
            else:
                self.assertEqual(merged_tandem_repeat_alleles[0].locus_id, "chr22-10515040-10515083-CCG" + ("CCC"*n), msg=msg)

    def test_merge_overlapping_tandem_repeat_loci4(self):
        """Test the merge_overlapping_tandem_repeat_loci function."""

        # two overlapping
        tandem_repeat_alleles = [
            MinimalTandemRepeatAllele("chr22", 10515040, 10515052, "CAG", DETECTION_MODE_PURE_REPEATS, allele_order=1),
            MinimalTandemRepeatAllele("chr22", 10515046, 10515058, "CAG", DETECTION_MODE_TRF, allele_order=2),
        ]
        merged_tandem_repeat_alleles = merge_overlapping_tandem_repeat_loci(tandem_repeat_alleles)
        self.assertEqual(len(merged_tandem_repeat_alleles), 1)
        self.assertEqual(merged_tandem_repeat_alleles[0].locus_id, "chr22-10515040-10515058-CAG")
        self.assertEqual(merged_tandem_repeat_alleles[0].detection_mode, "merged:pure,trf")

        # four overlapping
        tandem_repeat_alleles = [
            MinimalTandemRepeatAllele("chr22", 10515040, 10515052, "CAG", DETECTION_MODE_PURE_REPEATS, allele_order=1),
            MinimalTandemRepeatAllele("chr22", 10515046, 10515058, "CAG", DETECTION_MODE_TRF, allele_order=2),
            MinimalTandemRepeatAllele("chr22", 10515052, 10515064, "CAG", DETECTION_MODE_ALLOW_INTERRUPTIONS, allele_order=3),
            MinimalTandemRepeatAllele("chr22", 10515064, 10515076, "CAG", DETECTION_MODE_PURE_REPEATS, allele_order=4),
        ]
        merged_tandem_repeat_alleles = merge_overlapping_tandem_repeat_loci(tandem_repeat_alleles)
        self.assertEqual(len(merged_tandem_repeat_alleles), 1)
        self.assertEqual(merged_tandem_repeat_alleles[0].locus_id, "chr22-10515040-10515076-CAG")
        self.assertEqual(merged_tandem_repeat_alleles[0].detection_mode, "merged:interrupted,pure,trf")


        # adjacent
        tandem_repeat_alleles = [
            MinimalTandemRepeatAllele("chr22", 10515040, 10515046, "CAG", DETECTION_MODE_PURE_REPEATS, allele_order=1),
            MinimalTandemRepeatAllele("chr22", 10515047, 10515058, "CAG", DETECTION_MODE_TRF, allele_order=2),
        ]
        merged_tandem_repeat_alleles = merge_overlapping_tandem_repeat_loci(tandem_repeat_alleles)
        self.assertEqual(len(merged_tandem_repeat_alleles), 1)
        self.assertEqual(merged_tandem_repeat_alleles[0].locus_id, "chr22-10515040-10515058-CAG")
        self.assertEqual(merged_tandem_repeat_alleles[0].detection_mode, "merged:pure,trf")

        # A contains B
        tandem_repeat_alleles = [
            MinimalTandemRepeatAllele("chr22", 10515040, 10515058, "CAG", DETECTION_MODE_PURE_REPEATS, allele_order=1),
            MinimalTandemRepeatAllele("chr22", 10515046, 10515052, "CAG", DETECTION_MODE_TRF, allele_order=2),
        ]
        merged_tandem_repeat_alleles = merge_overlapping_tandem_repeat_loci(tandem_repeat_alleles)
        self.assertEqual(len(merged_tandem_repeat_alleles), 1)
        self.assertEqual(merged_tandem_repeat_alleles[0].locus_id, "chr22-10515040-10515058-CAG")
        self.assertEqual(merged_tandem_repeat_alleles[0].detection_mode, "merged:pure,trf")

        # B contains A
        tandem_repeat_alleles = [
            MinimalTandemRepeatAllele("chr22", 10515046, 10515052, "CAG", DETECTION_MODE_PURE_REPEATS, allele_order=1),
            MinimalTandemRepeatAllele("chr22", 10515040, 10515058, "CAG", DETECTION_MODE_TRF, allele_order=2),
        ]
        
        merged_tandem_repeat_alleles = merge_overlapping_tandem_repeat_loci(tandem_repeat_alleles)
        self.assertEqual(len(merged_tandem_repeat_alleles), 1)
        self.assertEqual(merged_tandem_repeat_alleles[0].locus_id, "chr22-10515040-10515058-CAG")
        self.assertEqual(merged_tandem_repeat_alleles[0].detection_mode, "merged:pure,trf")

        # not overlapping
        tandem_repeat_alleles = [
            MinimalTandemRepeatAllele("chr22", 10515040, 10515046, "CAG", DETECTION_MODE_PURE_REPEATS, allele_order=1),
            MinimalTandemRepeatAllele("chr22", 10515048, 10515058, "CAG", DETECTION_MODE_TRF, allele_order=2),
            MinimalTandemRepeatAllele("chr22", 10515048, 10515058, "CAG", DETECTION_MODE_TRF, allele_order=2),
        ]
        merged_tandem_repeat_alleles = merge_overlapping_tandem_repeat_loci(tandem_repeat_alleles)
        self.assertEqual(len(merged_tandem_repeat_alleles), 2)
        self.assertEqual(merged_tandem_repeat_alleles[0].locus_id, "chr22-10515040-10515046-CAG")
        self.assertEqual(merged_tandem_repeat_alleles[0].detection_mode, "pure")
        self.assertEqual(merged_tandem_repeat_alleles[1].locus_id, "chr22-10515048-10515058-CAG")
        self.assertEqual(merged_tandem_repeat_alleles[1].detection_mode, "merged:trf")


    def test_merge_overlapping_tandem_repeat_loci5(self):
        """Test the merge_overlapping_tandem_repeat_loci function."""

        # two overlapping, 4 total
        tandem_repeat_alleles = [
            MinimalTandemRepeatAllele("chr22", 10515040, 10515052, "CAG", DETECTION_MODE_PURE_REPEATS, allele_order=1),
            MinimalTandemRepeatAllele("chr22", 10515042, 10515152, "CAGCGGCAG", DETECTION_MODE_TRF, allele_order=2),
            MinimalTandemRepeatAllele("chr22", 10515043, 10515153, "T", DETECTION_MODE_TRF, allele_order=3),
            MinimalTandemRepeatAllele("chr22", 10515046, 10515058, "AGC", DETECTION_MODE_TRF, allele_order=4),
        ]
        merged_tandem_repeat_alleles = merge_overlapping_tandem_repeat_loci(tandem_repeat_alleles)
        self.assertEqual(len(merged_tandem_repeat_alleles), 3)
        self.assertEqual(merged_tandem_repeat_alleles[0].locus_id, "chr22-10515040-10515058-CAG")
        self.assertEqual(merged_tandem_repeat_alleles[0].detection_mode, "merged:pure,trf")


    def test_merge_overlapping_tandem_repeat_loci6(self):
        # two adjacent, 6 total
        tandem_repeat_alleles = [
            MinimalTandemRepeatAllele("chr22", 10515040, 10515052, "CAG", DETECTION_MODE_PURE_REPEATS, allele_order=1),
            MinimalTandemRepeatAllele("chr22", 10515042, 10515152, "CAGCGGCAG", DETECTION_MODE_TRF, allele_order=2),
            MinimalTandemRepeatAllele("chr22", 10515043, 10515153, "T", DETECTION_MODE_TRF, allele_order=3),
            MinimalTandemRepeatAllele("chr22", 10515053, 10515058, "AGC", DETECTION_MODE_TRF, allele_order=4),
            MinimalTandemRepeatAllele("chr22", 10515055, 10515058, "AGC", DETECTION_MODE_TRF, allele_order=4),
        ]
        merged_tandem_repeat_alleles = merge_overlapping_tandem_repeat_loci(tandem_repeat_alleles)
        self.assertEqual(len(merged_tandem_repeat_alleles), 3)
        self.assertEqual(merged_tandem_repeat_alleles[0].locus_id, "chr22-10515040-10515058-CAG")
        self.assertEqual(merged_tandem_repeat_alleles[0].detection_mode, "merged:pure,trf")


    def test_failed_filtering(self):
        alleles = [
            Allele("chr22", 10689286, "A", "ACAG", self._fasta_obj, allele_order=1),
            Allele("chr22", 10689286, "A", "ACAGCAG", self._fasta_obj, allele_order=2),
            Allele("chr22", 10689286, "A", "ACACACACA", self._fasta_obj, allele_order=3),
            Allele("chr22", 10689286, "A", "ACACACTCAACACACTCA", self._fasta_obj, allele_order=4),
            Allele("chr22", 10689286, "A", "A" + ("CAGTGGCA"*10 + "T")*2, self._fasta_obj, allele_order=5),
        ]
        counters = collections.defaultdict(int)
        results, alleles_to_process_next_using_trf = detect_perfect_and_almost_perfect_tandem_repeats(alleles, counters, self._default_args)
        self.assertEqual(len(results), 0)
        self.assertEqual(alleles_to_process_next_using_trf[0].allele_order, 3)
        self.assertEqual(len(alleles_to_process_next_using_trf), 3)

        results = detect_tandem_repeats_using_trf(alleles, counters, self._default_args)
        self.assertEqual(len(results), 0)
        #self.assertEqual(results[0].allele.allele_order, 5)

    def test_extended_flanks(self):
        pass

    def tearDown(self):
        """Tear down the test case."""
        self._fasta_obj.close()


