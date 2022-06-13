import collections
import os
import pyfaidx
import tempfile
import unittest

from filter_vcf_to_STR_variants import parse_num_alt_from_genotype, compute_short_and_long_allele_size
from str_analysis.filter_vcf_to_STR_variants import get_flanking_reference_sequences, check_if_variant_is_str


class Tests(unittest.TestCase):
    """

        Insertion:

             1234567890--1234567890..
        REF: TGTGTGTGTT--ACACACACTGTGTGTGTGGGGGG

        ALT: TGTGTGTGTTacACACACACTGTGTGTGTGGGGGG

            output-left: TGTGTGTTac
            output-right: acACACACAC

            pos: 10
            ref: T
            alt: TAC


        Deletion:

             12345678901234567890
        REF: TGTGTGTGTTACACACACTGTGTGTGTGGGGGG

        ALT: TGTGTGTGTT--ACACACTGTGTGTGTGGGG

            output-left: TGTGTGTTAC
            output-right: ACACACACTG

            pos: 10
            ref: TAC
            alt: T
    """

    def setUp(self):
        self.temp_fasta_file = tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False)
        self.temp_fasta_file.write(">chrTest1\n")
        self.temp_fasta_file.write("TGTGTGTGTTACACACACTGTGTGTGTGGGGGG\n")

        self.temp_fasta_file.write(">chrTest2\n")
        self.temp_fasta_file.write("TGTGTGTGTTACACACTGTGTGTGTGGGGGG\n")

        self.test3_seq1 = "GTG"*5
        self.test3_seq2 = "ATC"*5
        self.test3_seq3 = "TTC"*5
        self.temp_fasta_file.write(">chrTest3\n")
        self.temp_fasta_file.write(f"{self.test3_seq1}{self.test3_seq2}{self.test3_seq3}\n")

        """
        Test 4:
        
        Insertion

            chrTest4  4  C  CAGCAGCAG
            
            123456789012345678901
            TTTTTTCAGCAGCAGTTTTTT
                    |
                    GCAGCAGCAG
               
            
            left_flank = CAG 
            right_flank = AG
            
            truth: 
                ref:  3 x CAG
                alt:  6 x CAG
            
        ===================
        
        Deletion
            
            chrTest4 9 GCAG G
            
            123456789012345678901
            TTTTTTCAGCAGCAGTTTTTT
                     ---
            
            left_flank = CAG 
            right_flank = AG
            
            ref:  3 x CAG
            alt:  2 x CAG

        """
        self.temp_fasta_file.write(">chrTest4\n")
        self.temp_fasta_file.write("TTTTTTCAGCAGCAGTTTTTT\n")

        """
        ========
        Test 5:  Deletion:
        {
             'ref': 'ACGCCGACGCCGCCGCCGC'
             'alt': 'A',
             'left_flank': '',
             'right_flank': 'CGCCGCCGCCGCCGC'
             'num_repeats_ref': 10,
             'num_repeats_alt': 5,
             'num_repeats_right_flank': 5,
             'num_repeats_within_variant_bases': 5,
             'repeat_unit': 'CGC',
        }
        """

        self.temp_fasta_file.write(">chrTest5\n")
        self.temp_fasta_file.write("TTTTTACGCCGACGCCGCCGCCGCCGCCGCCGCCGCCGCTTTTT\n")

        """
        Test 6:  Insertion:
        {
            'ref': 'A',
            'alt': 'ACAG',
            'left_flank': '',
            'right_flank': 'CAGCAGCAGCAGCAGCAGCAGCAG'            
            'found_by_TRF': False,
            'is_perfect_repeat': True,
            'num_repeats_ref': 8,
            'num_repeats_alt': 9,
            'num_repeats_left_flank': 0,
            'num_repeats_right_flank': 8,
            'num_repeats_within_variant_bases': 1,
            'repeat_unit': 'CAG',
        }
        """
        self.temp_fasta_file.write(">chrTest6\n")
        self.temp_fasta_file.write("TTTTTACAGCAGCAGCAGCAGCAGCAGCAGTTTTT\n")

        self.temp_fasta_file.close()

        self.fasta_obj = pyfaidx.Fasta(self.temp_fasta_file.name, one_based_attributes=False, as_raw=True)

    def test_get_flanking_reference_sequences1(self):
        # insertion 2-STR
        left_sequence, variant_sequence, right_sequence = get_flanking_reference_sequences(self.fasta_obj, "chrTest1", 10, "T", "TAC")
        self.assertEqual(left_sequence, "TGTGTGTGTT")
        self.assertEqual(variant_sequence, "AC")
        self.assertEqual(right_sequence, "ACACACACTG")

        left_sequence, variant_sequence, right_sequence = get_flanking_reference_sequences(self.fasta_obj, "chrTest1", 10, "T", "TACAC")
        self.assertEqual(left_sequence, "TGTGTGTGTT")
        self.assertEqual(variant_sequence, "ACAC")
        self.assertEqual(right_sequence, "ACACACACTGTGTGTGTGGG")

        # deletion
        left_sequence, variant_sequence, right_sequence = get_flanking_reference_sequences(self.fasta_obj, "chrTest2", 10, "TAC", "T")
        self.assertEqual(left_sequence, "TGTGTGTGTT")
        self.assertEqual(variant_sequence, "AC")
        self.assertEqual(right_sequence, "ACACTGTGTG")

        # insertion 3-STR
        left_sequence, variant_sequence, right_sequence = get_flanking_reference_sequences(self.fasta_obj, "chrTest3", 15, self.test3_seq1[-1], self.test3_seq1[-1] + self.test3_seq2)
        self.assertEqual(left_sequence, self.test3_seq1)
        self.assertEqual(variant_sequence, self.test3_seq2)
        self.assertEqual(right_sequence, self.test3_seq2 + self.test3_seq3)

    def test_parse_num_alt_from_genotype(self):
        self.assertEqual(parse_num_alt_from_genotype("1/1"), 2)
        self.assertEqual(parse_num_alt_from_genotype("1|1"), 2)
        self.assertEqual(parse_num_alt_from_genotype("1\\1"), 2)
        self.assertEqual(parse_num_alt_from_genotype(".|1"), 1)
        self.assertEqual(parse_num_alt_from_genotype("1|."), 1)
        self.assertEqual(parse_num_alt_from_genotype("1/."), 1)
        self.assertEqual(parse_num_alt_from_genotype("./1"), 1)
        self.assertEqual(parse_num_alt_from_genotype("0/1"), 1)
        self.assertEqual(parse_num_alt_from_genotype("1/0"), 1)
        self.assertEqual(parse_num_alt_from_genotype("0/0"), 0)
        self.assertEqual(parse_num_alt_from_genotype("0/."), 0)
        self.assertEqual(parse_num_alt_from_genotype("./."), 0)

    def test_check_if_variant_is_str4_insertion(self):
        # insertion 2-STR
        counters = collections.defaultdict(int)

        (
            chrom, start_1based, end_1based, repeat_unit, num_repeats_ref, num_repeats_alt,
            num_repeats_left_flank, num_repeats_right_flank, found_by_TRF, is_perfect_repeat
        ) = check_if_variant_is_str(
            self.fasta_obj,
            "chrTest4", 9, "G", "GCAGCAGCAG",
            min_str_repeats=3, min_str_length=9,
            min_fraction_of_variant_covered_by_repeat=0.9,
            counters=counters,
            use_trf=False)

        self.assertEqual(chrom, "chrTest4")
        self.assertEqual(start_1based, 7)
        self.assertEqual(end_1based, 15)
        self.assertEqual(repeat_unit, "CAG")
        self.assertEqual(num_repeats_ref, 3)
        self.assertEqual(num_repeats_alt, 6)
        self.assertEqual(num_repeats_left_flank, 1)
        self.assertEqual(num_repeats_right_flank, 2)
        self.assertFalse(found_by_TRF)
        self.assertTrue(is_perfect_repeat)

    def test_check_if_variant_is_str4_deletion(self):
        # insertion 2-STR
        counters = collections.defaultdict(int)

        (
            chrom, start_1based, end_1based, repeat_unit, num_repeats_ref, num_repeats_alt,
            num_repeats_left_flank, num_repeats_right_flank, found_by_TRF, is_perfect_repeat
        ) = check_if_variant_is_str(
            self.fasta_obj,
            "chrTest4", 9, "GCAG", "G",
            min_str_repeats=3, min_str_length=9,
            min_fraction_of_variant_covered_by_repeat=0.9,
            counters=counters,
            use_trf=False)

        self.assertEqual(chrom, "chrTest4")
        self.assertEqual(start_1based, 7)
        self.assertEqual(end_1based, 15)
        self.assertEqual(repeat_unit, "CAG")
        self.assertEqual(num_repeats_ref, 3)
        self.assertEqual(num_repeats_alt, 2)
        self.assertEqual(num_repeats_left_flank, 1)
        self.assertEqual(num_repeats_right_flank, 1)
        self.assertFalse(found_by_TRF)
        self.assertTrue(is_perfect_repeat)

    def test_check_if_variant_is_str5_deletion(self):
        # insertion 2-STR
        counters = collections.defaultdict(int)

        (
            chrom, start_1based, end_1based, repeat_unit, num_repeats_ref, num_repeats_alt,
            num_repeats_left_flank, num_repeats_right_flank, found_by_TRF, is_perfect_repeat
        ) = check_if_variant_is_str(
            self.fasta_obj,
            "chrTest5", 6, "ACGCCGACGCCGCCGCCGC", "A",
            min_str_repeats=3, min_str_length=9,
            min_fraction_of_variant_covered_by_repeat=0.8,
            counters=counters,
            use_trf=False)

        self.assertEqual(chrom, "chrTest5")
        self.assertEqual(start_1based, 7)
        self.assertEqual(end_1based, 39)
        self.assertEqual(repeat_unit, "CGC")
        self.assertEqual(num_repeats_ref, 11)
        self.assertEqual(num_repeats_alt, 5)
        self.assertEqual(num_repeats_left_flank, 0)
        self.assertEqual(num_repeats_right_flank, 5)
        self.assertFalse(found_by_TRF)
        self.assertFalse(is_perfect_repeat)

    def test_check_if_variant_is_str6_insertion(self):
        # insertion 2-STR
        counters = collections.defaultdict(int)

        (
            chrom, start_1based, end_1based, repeat_unit, num_repeats_ref, num_repeats_alt,
            num_repeats_left_flank, num_repeats_right_flank, found_by_TRF, is_perfect_repeat
        ) = check_if_variant_is_str(
            self.fasta_obj,
            "chrTest6", 6, "A", "ACAG",
            min_str_repeats=3, min_str_length=9,
            min_fraction_of_variant_covered_by_repeat=0.8,
            counters=counters,
            use_trf=False)

        self.assertEqual(chrom, "chrTest6")
        self.assertEqual(start_1based, 7)
        self.assertEqual(end_1based, 30)
        self.assertEqual(repeat_unit, "CAG")
        self.assertEqual(num_repeats_ref, 8)
        self.assertEqual(num_repeats_alt, 9)
        self.assertEqual(num_repeats_left_flank, 0)
        self.assertEqual(num_repeats_right_flank, 8)
        self.assertFalse(found_by_TRF)
        self.assertTrue(is_perfect_repeat)

    def test_compute_short_and_long_allele_size_str4_insertion(self):
        short_allele_size, long_allele_size = compute_short_and_long_allele_size(7, 15, "G", "GCAGCAGCAG", 3, 1)
        self.assertEqual(short_allele_size, 3)
        self.assertEqual(long_allele_size, 6)

        short_allele_size, long_allele_size = compute_short_and_long_allele_size(7, 15, "G", "GCAGCAGCAG", 3, 2)
        self.assertEqual(short_allele_size, 6)
        self.assertEqual(long_allele_size, 6)

    def test_compute_short_and_long_allele_size_str4_deletion(self):
        short_allele_size, long_allele_size = compute_short_and_long_allele_size(7, 15, "GCAG", "G", 3, 1)
        self.assertEqual(long_allele_size, 3)
        self.assertEqual(short_allele_size, 2)

        short_allele_size, long_allele_size = compute_short_and_long_allele_size(7, 15, "GCAG", "G", 3, 2)
        self.assertEqual(long_allele_size, 2)
        self.assertEqual(short_allele_size, 2)


    def tearDown(self):
        if os.path.isfile(self.temp_fasta_file.name):
            os.remove(self.temp_fasta_file.name)
