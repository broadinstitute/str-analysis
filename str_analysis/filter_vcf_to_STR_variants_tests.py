import collections
import os
from pprint import pprint, pformat
import pyfaidx
import tempfile
import unittest

from str_analysis.filter_vcf_to_STR_variants import get_flanking_reference_sequences, check_if_allele_is_str


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

        # Test7: insert denovo repeat sequence
        self.temp_fasta_file.write(">chrTest7\n")
        self.temp_fasta_file.write("TTTTTCCCCC\n")

        # Test8: delete entire repeat sequence
        self.temp_fasta_file.write(">chrTest8\n")
        self.temp_fasta_file.write("TTTTTCAGCAGCAGCCCCC\n")

        self.temp_fasta_file.close()

        self.fasta_obj = pyfaidx.Fasta(self.temp_fasta_file.name, one_based_attributes=False, as_raw=True)

    def check_results_for_pure_repeats(self, results):
        for key in results:
            key2 = key.replace("Pure", "")
            if key2 in results and "Pure" in key:
                self.assertEqual(results[key], results[key2], f"{key} != {key2} in " + pformat(results))

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

    def test_check_if_allele_is_str4_insertion(self):
        # insertion 2-STR
        counters = collections.defaultdict(int)

        result = check_if_allele_is_str(
            self.fasta_obj,
            "chrTest4", 9, "G", "GCAGCAGCAG",
            min_str_repeats=3, min_str_length=9,
            counters=counters,
            allow_interruptions=True)

        self.assertEqual(result["Chrom"], "chrTest4")
        self.assertEqual(result["Start1Based"], 7)
        self.assertEqual(result["End1Based"], 15)
        self.assertEqual(result["RepeatUnit"], "CAG")
        self.assertEqual(result["NumRepeatsRef"], 3)
        self.assertEqual(result["NumRepeatsAlt"], 6)
        self.assertEqual(result["NumRepeatsLeftFlank"], 1)
        self.assertEqual(result["NumRepeatsRightFlank"], 2)
        self.assertEqual(result["IsPureRepeat"], "Yes")
        self.check_results_for_pure_repeats(result)

    def test_check_if_allele_is_str4_deletion(self):
        # insertion 2-STR
        counters = collections.defaultdict(int)

        result = check_if_allele_is_str(
            self.fasta_obj,
            "chrTest4", 9, "GCAG", "G",
            min_str_repeats=3, min_str_length=9,
            counters=counters,
            allow_interruptions=True)

        self.assertEqual(result["Chrom"], "chrTest4")
        self.assertEqual(result["Start1Based"], 7)
        self.assertEqual(result["End1Based"], 15)
        self.assertEqual(result["RepeatUnit"], "CAG")
        self.assertEqual(result["NumRepeatsRef"], 3)
        self.assertEqual(result["NumRepeatsAlt"], 2)
        self.assertEqual(result["NumRepeatsLeftFlank"], 1)
        self.assertEqual(result["NumRepeatsRightFlank"], 1)
        self.assertEqual(result["IsPureRepeat"], "Yes")
        self.check_results_for_pure_repeats(result)

    def test_check_if_allele_is_str5_deletion(self):
        # insertion 2-STR
        counters = collections.defaultdict(int)

        result = check_if_allele_is_str(
            self.fasta_obj,
            "chrTest5", 6, "ACGCCGACGCCGCCGCCGC", "A",
            min_str_repeats=3, min_str_length=9,
            counters=counters,
            allow_interruptions=True)

        self.assertEqual(result["Chrom"], "chrTest5")
        self.assertEqual(result["Start1Based"], 7)
        self.assertEqual(result["End1Based"], 39)
        self.assertEqual(result["RepeatUnit"], "CGC")
        self.assertEqual(result["NumRepeatsRef"], 11)
        self.assertEqual(result["NumRepeatsAlt"], 5)
        self.assertEqual(result["NumRepeatsLeftFlank"], 0)
        self.assertEqual(result["NumRepeatsRightFlank"], 5)
        self.assertEqual(result["IsPureRepeat"], "No")

        self.assertEqual(result["RepeatUnitInterruptionIndex"], 2)
        self.assertEqual(result["FractionPureRepeats"], 10/11)

        self.assertEqual(result["PureStart1Based"], 7)
        self.assertEqual(result["PureEnd1Based"], 39)
        self.assertEqual(result["NumPureRepeatsInStr"], 10)
        self.assertEqual(result["NumPureRepeatsInVariant"], 5)
        self.assertEqual(result["NumPureRepeatsLeftFlank"], 0)
        self.assertEqual(result["NumPureRepeatsRightFlank"], 5)

    def test_check_if_allele_is_str6_insertion(self):
        # insertion 2-STR
        counters = collections.defaultdict(int)

        result = check_if_allele_is_str(
            self.fasta_obj,
            "chrTest6", 6, "A", "ACAG",
            min_str_repeats=3, min_str_length=9,
            counters=counters,
            allow_interruptions=True)

        self.assertEqual(result["Chrom"], "chrTest6")
        self.assertEqual(result["Start1Based"], 7)
        self.assertEqual(result["End1Based"], 30)
        self.assertEqual(result["RepeatUnit"], "CAG")
        self.assertEqual(result["NumRepeatsRef"], 8)
        self.assertEqual(result["NumRepeatsAlt"], 9)
        self.assertEqual(result["NumRepeatsLeftFlank"], 0)
        self.assertEqual(result["NumRepeatsRightFlank"], 8)
        self.assertEqual(result["IsPureRepeat"], "Yes")
        self.check_results_for_pure_repeats(result)

    def test_check_if_allele_is_str7_insertion(self):
        # insertion 2-STR
        counters = collections.defaultdict(int)

        result = check_if_allele_is_str(
            self.fasta_obj,
            "chrTest6", 6, "A", "ACAG",
            min_str_repeats=3, min_str_length=9,
            counters=counters,
            allow_interruptions=True)

        self.assertEqual(result["Chrom"], "chrTest6")
        self.assertEqual(result["Start1Based"], 7)
        self.assertEqual(result["End1Based"], 30)
        self.assertEqual(result["RepeatUnit"], "CAG")
        self.assertEqual(result["NumRepeatsRef"], 8)
        self.assertEqual(result["NumRepeatsAlt"], 9)
        self.assertEqual(result["NumRepeatsLeftFlank"], 0)
        self.assertEqual(result["NumRepeatsRightFlank"], 8)
        self.assertEqual(result["IsPureRepeat"], "Yes")
        self.check_results_for_pure_repeats(result)

    def test_check_if_allele_is_str8_insertion(self):
        # insertion 2-STR
        counters = collections.defaultdict(int)

        result = check_if_allele_is_str(
            self.fasta_obj,
            "chrTest6", 6, "A", "ACAG",
            min_str_repeats=3, min_str_length=9,
            counters=counters,
            allow_interruptions=True)

        self.assertEqual(result["Chrom"], "chrTest6")
        self.assertEqual(result["Start1Based"], 7)
        self.assertEqual(result["End1Based"], 30)
        self.assertEqual(result["RepeatUnit"], "CAG")
        self.assertEqual(result["NumRepeatsRef"], 8)
        self.assertEqual(result["NumRepeatsAlt"], 9)
        self.assertEqual(result["NumRepeatsLeftFlank"], 0)
        self.assertEqual(result["NumRepeatsRightFlank"], 8)
        self.assertEqual(result["IsPureRepeat"], "Yes")
        self.check_results_for_pure_repeats(result)

    def test_check_if_allele_is_str9_insertion(self):
        # insertion 2-STR
        counters = collections.defaultdict(int)

        # TTTTTCCCCC
        result = check_if_allele_is_str(
            self.fasta_obj,
            "chrTest7", 6, "C", "CCAGCAGCAG",
            min_str_repeats=3, min_str_length=9,
            counters=counters,
            allow_interruptions=True)

        self.assertEqual(result["Chrom"], "chrTest7")
        self.assertEqual(result["Start1Based"], 7)
        self.assertEqual(result["End1Based"], 6)
        self.assertEqual(result["RepeatUnit"], "CAG")
        self.assertEqual(result["NumRepeatsRef"], 0)
        self.assertEqual(result["NumRepeatsAlt"], 3)
        self.assertEqual(result["NumRepeatsLeftFlank"], 0)
        self.assertEqual(result["NumRepeatsRightFlank"], 0)
        self.assertEqual(result["IsPureRepeat"], "Yes")
        self.check_results_for_pure_repeats(result)


    def test_check_if_allele_is_str10_insertion(self):
        # insertion 2-STR
        counters = collections.defaultdict(int)

        # TTTTTCAGCAGCAGCCCCC
        result = check_if_allele_is_str(
            self.fasta_obj,
            "chrTest8", 5, "TCAGCAGCAG", "T",
            min_str_repeats=3, min_str_length=9,
            counters=counters,
            allow_interruptions=True)

        self.assertEqual(result["Chrom"], "chrTest8")
        self.assertEqual(result["Start1Based"], 6)
        self.assertEqual(result["End1Based"], 14)
        self.assertEqual(result["RepeatUnit"], "CAG")
        self.assertEqual(result["NumRepeatsRef"], 3)
        self.assertEqual(result["NumRepeatsAlt"], 0)
        self.assertEqual(result["NumRepeatsLeftFlank"], 0)
        self.assertEqual(result["NumRepeatsRightFlank"], 0)
        self.assertEqual(result["IsPureRepeat"], "Yes")
        self.check_results_for_pure_repeats(result)

    def tearDown(self):
        if os.path.isfile(self.temp_fasta_file.name):
            os.remove(self.temp_fasta_file.name)

