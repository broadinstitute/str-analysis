import collections
import os
from pprint import pformat, pprint
import pyfaidx
import tempfile
import unittest

from str_analysis.filter_vcf_to_STR_variants import get_flanking_reference_sequences, check_if_allele_is_str, \
    get_num_repeats_in_allele, check_if_variant_is_str, compute_indel_variant_bases


class Tests(unittest.TestCase):

    def setUp(self):

        self.maxDiff = None
        self.temp_fasta_file = tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False)

        """
        Test 1 - Insertion:
        
             1234567890--1234567890..
        REF: TGTGTGTGTT--ACACACACTGTGTGTGTGGGGGG
        
        ALT: TGTGTGTGTTacACACACACTGTGTGTGTGGGGGG
        
            output-left: TGTGTGTTac
            output-right: acACACACAC
        
            pos: 10
            ref: T
            alt: TAC
        
        
        Test 1 - Deletion:
        
             12345678901234567890
        REF: TGTGTGTGTTACACACACTGTGTGTGTGGGGGG
        
        ALT: TGTGTGTGTT--ACACACTGTGTGTGTGGGG
        
            output-left: TGTGTGTTAC
            output-right: ACACACACTG
        
            pos: 10
            ref: TAC
            alt: T
        """
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


    def test_get_num_repeats_in_allele(self):
        single_allele_str_spec = [
            {"NumRepeatsRef": 1, "NumRepeatsAlt": 2},
        ]
        self.assertEqual(get_num_repeats_in_allele(single_allele_str_spec, 0), 1)
        self.assertEqual(get_num_repeats_in_allele(single_allele_str_spec, 1), 2)
        self.assertRaises(ValueError, lambda: get_num_repeats_in_allele(single_allele_str_spec, 2))

        multiallelic_str_spec = [
            {"NumRepeatsRef": 1, "NumRepeatsAlt": 2},
            {"NumRepeatsRef": 1, "NumRepeatsAlt": 3},
        ]
        self.assertEqual(get_num_repeats_in_allele(multiallelic_str_spec, 0), 1)
        self.assertEqual(get_num_repeats_in_allele(multiallelic_str_spec, 1), 2)
        self.assertEqual(get_num_repeats_in_allele(multiallelic_str_spec, 2), 3)
        self.assertRaises(ValueError, lambda: get_num_repeats_in_allele(multiallelic_str_spec, 3))

    def test_compute_indel_variant_bases(self):
        # test cases for compute_indel_variant_bases
        self.assertEqual(compute_indel_variant_bases("A", "AC"), "C")
        self.assertEqual(compute_indel_variant_bases("A", "ACG"), "CG")
        self.assertEqual(compute_indel_variant_bases("AC", "A"), "C")
        self.assertEqual(compute_indel_variant_bases("ACAGCAGCAT", "ACAG"), "CAGCAT")
        self.assertEqual(compute_indel_variant_bases("ACAGCAGCAT", "GCAT"), None)
        self.assertEqual(compute_indel_variant_bases("GCAT", "ACAGCAGCAT"), None)
        self.assertEqual(compute_indel_variant_bases("A", "ACAGCAGCAT"), "CAGCAGCAT")
        self.assertEqual(compute_indel_variant_bases("ACAGCAGCAT", "A"), "CAGCAGCAT")


    def test_check_if_variant_is_str(self):

        counters = collections.defaultdict(int)
        str_specs, filter_reason = check_if_variant_is_str(
            self.fasta_obj, "chrTest4", 9, "G", ["GCAGCAGCAG"],
            min_str_repeats=3, min_str_length=9, min_repeat_unit_length=1, max_repeat_unit_length=50,
            counters=counters, interruptions="no")
        self.assertEqual(len(str_specs), 1)
        self.assertDictEqual(str_specs[0], {
            'Chrom': 'chrTest4',
            'Pos': 9,
            'Ref': 'G',
            'Alt': 'GCAGCAGCAG',
            'Start1Based': 7,
            'End1Based': 15,
            'PureStart1Based': 7,
            'PureEnd1Based': 15,
            'FilterReason': None,
            'FractionPureRepeatsRef': 1.0,
            'FractionPureRepeatsAlt': 1.0,
            'IsPureRepeat': True,
            'NumPureRepeatsRef': 3,
            'NumPureRepeatsAlt': 6,
            'NumRepeatsLeftFlank': 1,
            'NumRepeatsRightFlank': 2,

            'NumRepeatsRef': 3,
            'NumRepeatsAlt': 6,
            'NumRepeatsInVariant': 3,

            'NumPureRepeatsInVariant': 3,
            'NumPureRepeatsLeftFlank': 1,
            'NumPureRepeatsRightFlank': 2,
            'RepeatUnit': 'CAG',
            'MotifInterruptionIndex': None,
        })

        # test homopolymer + min_repeat_unit_length
        str_specs, filter_string = check_if_variant_is_str(
            self.fasta_obj, "chrTest4", 9, "G", ["GGGGGGGGGGGGGGG"],
            min_str_repeats=3, min_str_length=9, min_repeat_unit_length=2, max_repeat_unit_length=50,
            counters=counters, interruptions="no")
        self.assertEqual(len(str_specs), 1)
        self.assertDictEqual(str_specs[0], {'FilterReason': 'repeat unit < 2 bp', 'RepeatUnit': None})

        # test dinucleotide + min_repeat_unit_length
        str_specs, filter_string = check_if_variant_is_str(
            self.fasta_obj, "chrTest4", 9, "G", ["GAGAGAGAGAGAGAGAGAGAG"],
            min_str_repeats=3, min_str_length=9, min_repeat_unit_length=3, max_repeat_unit_length=50,
            counters=counters, interruptions="no")
        self.assertEqual(len(str_specs), 1)
        self.assertDictEqual(str_specs[0], {'FilterReason': 'repeat unit < 3 bp', 'RepeatUnit': None})

        # test homopolymer
        str_specs, filter_string = check_if_variant_is_str(
            self.fasta_obj, "chrTest4", 9, "G", ["GAAAAAAAAAA"],
            min_str_repeats=3, min_str_length=9, min_repeat_unit_length=1, max_repeat_unit_length=50,
            counters=counters, interruptions="no")
        self.assertEqual(len(str_specs), 1)
        self.assertEqual(filter_string, None)
        self.assertDictEqual(str_specs[0], {'Alt': 'GAAAAAAAAAA',
                                        'Chrom': 'chrTest4',
                                        'End1Based': 9,
                                        'FilterReason': None,
                                        'FractionPureRepeatsRef': None,
                                        'FractionPureRepeatsAlt': 1.0,
                                        'IsPureRepeat': True,
                                        'NumPureRepeatsAlt': 10,
                                        'NumPureRepeatsInVariant': 10,
                                        'NumPureRepeatsLeftFlank': 0,
                                        'NumPureRepeatsRef': 0,
                                        'NumPureRepeatsRightFlank': 0,
                                        'NumRepeatsAlt': 10,
                                        'NumRepeatsInVariant': 10,
                                        'NumRepeatsLeftFlank': 0,
                                        'NumRepeatsRef': 0,
                                        'NumRepeatsRightFlank': 0,
                                        'Pos': 9,
                                        'PureEnd1Based': 9,
                                        'PureStart1Based': 10,
                                        'Ref': 'G',
                                        'RepeatUnit': 'A',
                                        'MotifInterruptionIndex': None,
                                        'Start1Based': 10})

        # test homopolymer with interruptions, allow_interruptions = False
        str_specs, filter_string = check_if_variant_is_str(
            self.fasta_obj, "chrTest4", 9, "G", ["GAAAAACAAAAAAAC"],
            min_str_repeats=3, min_str_length=9, min_repeat_unit_length=1, max_repeat_unit_length=50,
            counters=counters, interruptions="no")
        self.assertEqual(len(str_specs), 1)
        self.assertEqual(filter_string, str_specs[0]["FilterReason"])
        self.assertEqual(str_specs[0], {'FilterReason': 'INDEL without repeats', 'RepeatUnit': None})

        # test homopolymer with interruptions, allow_interruptions = True
        str_specs, filter_string = check_if_variant_is_str(
            self.fasta_obj, "chrTest4", 9, "G", ["GAAAAACTAAAAAAACAAAAAAAAAAAAAAAAAAAA"],
            min_str_repeats=3, min_str_length=9, min_repeat_unit_length=1, max_repeat_unit_length=50,
            counters=counters, interruptions="only-if-pure-repeats-not-found")
        self.assertEqual(len(str_specs), 1)
        self.assertEqual(filter_string, str_specs[0]["FilterReason"])
        self.assertEqual(str_specs[0], {'FilterReason': 'INDEL without repeats', 'RepeatUnit': None})

        # test homopolymer with interruptions, allow_interruptions = True
        str_specs, filter_string = check_if_variant_is_str(
            self.fasta_obj, "chrTest4", 9, "G", ["GAAAAACAAAAAAACAAAAAAAAAAAAAAAAAAAA"],
            min_str_repeats=3, min_str_length=9, min_repeat_unit_length=1, max_repeat_unit_length=50,
            counters=counters, interruptions="only-if-pure-repeats-not-found")
        self.assertEqual(len(str_specs), 1)
        self.assertEqual(filter_string, None)
        self.assertEqual(str_specs[0], {'Alt': 'GAAAAACAAAAAAACAAAAAAAAAAAAAAAAAAAA',
                                    'Chrom': 'chrTest4',
                                    'End1Based': 9,
                                    'FilterReason': None,
                                    'FractionPureRepeatsRef': None,
                                    'FractionPureRepeatsAlt': 0.9411764705882353,
                                    'IsPureRepeat': False,
                                    'NumPureRepeatsAlt': 32,
                                    'NumPureRepeatsInVariant': 32,
                                    'NumPureRepeatsLeftFlank': 0,
                                    'NumPureRepeatsRef': 0,
                                    'NumPureRepeatsRightFlank': 0,
                                    'NumRepeatsAlt': 34,
                                    'NumRepeatsInVariant': 34,
                                    'NumRepeatsLeftFlank': 0,
                                    'NumRepeatsRef': 0,
                                    'NumRepeatsRightFlank': 0,
                                    'Pos': 9,
                                    'PureEnd1Based': 9,
                                    'PureStart1Based': 10,
                                    'Ref': 'G',
                                    'RepeatUnit': 'A',
                                    'MotifInterruptionIndex': 0,
                                    'Start1Based': 10})

        # test dinucleotides with interruptions, allow_interruptions = True
        str_specs, filter_string = check_if_variant_is_str(
            self.fasta_obj, "chrTest4", 9, "G", ["GAGACACAGAGAG"],
            min_str_repeats=3, min_str_length=9, min_repeat_unit_length=1, max_repeat_unit_length=50,
            counters=counters, interruptions="only-if-pure-repeats-not-found")
        self.assertEqual(len(str_specs), 1)
        self.assertEqual(filter_string, str_specs[0]["FilterReason"])
        self.assertEqual(str_specs[0], {'FilterReason': 'INDEL without repeats', 'RepeatUnit': None})

        str_specs, filter_string = check_if_variant_is_str(
            self.fasta_obj, "chrTest4", 9, "G", ["GAGACACAGAGAGAGAGAGAG"],
            min_str_repeats=3, min_str_length=9, min_repeat_unit_length=1, max_repeat_unit_length=50,
            counters=counters, interruptions="only-if-pure-repeats-not-found")
        self.assertEqual(filter_string, None)
        self.assertEqual(len(str_specs), 1)
        self.assertEqual(str_specs[0], {'Alt': 'GAGACACAGAGAGAGAGAGAG',
                                    'Chrom': 'chrTest4',
                                    'End1Based': 9,
                                    'FilterReason': None,
                                    'FractionPureRepeatsRef': 1,
                                    'FractionPureRepeatsAlt': 0.8181818181818182,
                                    'IsPureRepeat': False,
                                    'NumPureRepeatsAlt': 9,
                                    'NumPureRepeatsInVariant': 8,
                                    'NumPureRepeatsLeftFlank': 1,
                                    'NumPureRepeatsRef': 1,
                                    'NumPureRepeatsRightFlank': 0,
                                    'NumRepeatsAlt': 11,
                                    'NumRepeatsInVariant': 10,
                                    'NumRepeatsLeftFlank': 1,
                                    'NumRepeatsRef': 1,
                                    'NumRepeatsRightFlank': 0,
                                    'Pos': 9,
                                    'PureEnd1Based': 9,
                                    'PureStart1Based': 8,
                                    'Ref': 'G',
                                    'RepeatUnit': 'AG',
                                    'MotifInterruptionIndex': 1,
                                    'Start1Based': 8})

        # trinucleotide repeat
        str_specs, filter_string = check_if_variant_is_str(
            self.fasta_obj, "chrTest4", 9, "G", ["GCAGCAGCATCAACACCAC"],
            min_str_repeats=3, min_str_length=9, min_repeat_unit_length=1, max_repeat_unit_length=50,
            counters=counters, interruptions="only-if-pure-repeats-not-found")
        self.assertEqual(len(str_specs), 1)
        self.assertEqual(filter_string, str_specs[0]["FilterReason"])
        self.assertEqual(str_specs[0], {'Alt': 'GCAGCAGCATCAACACCAC',
                                    'Chrom': 'chrTest4',
                                    'End1Based': 15,
                                    'FilterReason': None,
                                    'FractionPureRepeatsRef': 1,
                                    'FractionPureRepeatsAlt': 0.5555555555555556,
                                    'IsPureRepeat': False,
                                    'NumPureRepeatsAlt': 5,
                                    'NumPureRepeatsInVariant': 2,
                                    'NumPureRepeatsLeftFlank': 1,
                                    'NumPureRepeatsRef': 3,
                                    'NumPureRepeatsRightFlank': 2,
                                    'NumRepeatsAlt': 9,
                                    'NumRepeatsInVariant': 6,
                                    'NumRepeatsLeftFlank': 1,
                                    'NumRepeatsRef': 3,
                                    'NumRepeatsRightFlank': 2,
                                    'Pos': 9,
                                    'PureEnd1Based': 15,
                                    'PureStart1Based': 7,
                                    'Ref': 'G',
                                    'RepeatUnit': 'CAG',
                                    'MotifInterruptionIndex': 2,
                                    'Start1Based': 7})

    def test_check_if_multiallelic_variant_is_str(self):
        counters = collections.defaultdict(int)
        str_specs, filter_string = check_if_variant_is_str(
            self.fasta_obj, "chrTest4", 9, "G", ["GCAGCAGCAG", "GCAGCAGCATCAACACCAC"],
            min_str_repeats=3, min_str_length=9, min_repeat_unit_length=1, max_repeat_unit_length=50,
            counters=counters, interruptions="no")
        self.assertEqual(filter_string, "STR allele;INDEL without repeats")
        self.assertEqual(str_specs[0]["FilterReason"], None)
        self.assertEqual(str_specs[1]["FilterReason"], "INDEL without repeats")
        self.assertDictEqual(counters, {'STR allele counts: INS': 1,
                                        'STR allele counts: TOTAL': 1,
                                        'STR allele delta: 3 repeats': 1,
                                        'STR allele motif size: 3 bp': 1,
                                        'STR allele reference repeats: with both left and right matching ref. repeat': 1,
                                        'STR allele size: 0-25bp': 1,
                                        'allele counts: INS alleles': 2,
                                        'allele filter: INDEL without repeats': 1,
                                        'variant filter: multi-allelic: one STR, one non-STR allele': 1})

        counters = collections.defaultdict(int)
        str_specs, filter_string = check_if_variant_is_str(
            self.fasta_obj, "chrTest4", 9, "G", ["GCAGCAGCAG", "GCAGCAGCATCAACACCAC"],
            min_str_repeats=3, min_str_length=9, min_repeat_unit_length=1, max_repeat_unit_length=50,
            counters=counters, interruptions="only-if-pure-repeats-not-found")

        self.assertDictEqual(str_specs[0], {'Alt': 'GCAGCAGCAG',
                                            'Chrom': 'chrTest4',
                                            'End1Based': 15,
                                            'FilterReason': None,
                                            'FractionPureRepeatsRef': 1.0,
                                            'FractionPureRepeatsAlt': 1.0,
                                            'IsPureRepeat': False,
                                            'NumPureRepeatsAlt': 6,
                                            'NumPureRepeatsInVariant': 3,
                                            'NumPureRepeatsLeftFlank': 1,
                                            'NumPureRepeatsRef': 3,
                                            'NumPureRepeatsRightFlank': 2,
                                            'NumRepeatsAlt': 6,
                                            'NumRepeatsInVariant': 3,
                                            'NumRepeatsLeftFlank': 1,
                                            'NumRepeatsRef': 3,
                                            'NumRepeatsRightFlank': 2,
                                            'Pos': 9,
                                            'PureEnd1Based': 15,
                                            'PureStart1Based': 7,
                                            'Ref': 'G',
                                            'RepeatUnit': 'CAG',
                                            'MotifInterruptionIndex': 2,
                                            'Start1Based': 7})
        self.assertDictEqual(str_specs[1], {'Alt': 'GCAGCAGCATCAACACCAC',
                                            'Chrom': 'chrTest4',
                                            'End1Based': 15,
                                            'FilterReason': None,
                                            'FractionPureRepeatsRef': 1,
                                            'FractionPureRepeatsAlt': 0.5555555555555556,
                                            'IsPureRepeat': False,
                                            'NumPureRepeatsAlt': 5,
                                            'NumPureRepeatsInVariant': 2,
                                            'NumPureRepeatsLeftFlank': 1,
                                            'NumPureRepeatsRef': 3,
                                            'NumPureRepeatsRightFlank': 2,
                                            'NumRepeatsAlt': 9,
                                            'NumRepeatsInVariant': 6,
                                            'NumRepeatsLeftFlank': 1,
                                            'NumRepeatsRef': 3,
                                            'NumRepeatsRightFlank': 2,
                                            'Pos': 9,
                                            'PureEnd1Based': 15,
                                            'PureStart1Based': 7,
                                            'Ref': 'G',
                                            'RepeatUnit': 'CAG',
                                            'MotifInterruptionIndex': 2,
                                            'Start1Based': 7})

        self.assertDictEqual(counters, {'STR allele counts: INS': 2,
                              'STR allele counts: TOTAL': 2,
                              'STR allele delta: 3 repeats': 1,
                              'STR allele delta: 6 repeats': 1,
                              'STR allele motif size: 3 bp': 2,
                              'STR allele reference repeats: with both left and right matching ref. repeat': 2,
                              'STR allele size: 0-25bp': 2,
                              'allele counts: INS alleles': 2})

        counters = collections.defaultdict(int)
        _, filter_string = check_if_variant_is_str(
            self.fasta_obj, "chrTest4", 9, "G", ["GCAGCAGCAG", "GCAGCATCATCATCACCAT"],
            min_str_repeats=3, min_str_length=9, min_repeat_unit_length=1, max_repeat_unit_length=50,
            counters=counters, interruptions="only-if-pure-repeats-not-found")
        self.assertEqual(filter_string, "STR alleles with different motifs;STR alleles with different motifs")

        _, filter_string = check_if_variant_is_str(
            self.fasta_obj, "chrTest4", 9, "G", ["GCAGCAGCAGCAT", "GCTGCAGCAGCAGCAG"],
            min_str_repeats=3, min_str_length=9, min_repeat_unit_length=1, max_repeat_unit_length=50,
            counters=counters, interruptions="only-if-pure-repeats-not-found")
        self.assertEqual(filter_string, "STR alleles with different interruption patterns;STR alleles with different interruption patterns")

    def check_if_multiallelic_variant_is_str_with_allow_interruptions_always(self):
        counters = collections.defaultdict(int)
        str_specs, filter_string = check_if_variant_is_str(
            self.fasta_obj, "chrTest4", 9, "G", ["GCAGCAGCAG", "GCAGCAGCATCAACACCAC"],
            min_str_repeats=3, min_str_length=9, min_repeat_unit_length=1, max_repeat_unit_length=50,
            counters=counters, interruptions="always")

        self.assertDictEqual(str_specs[0], {'Alt': 'GCAGCAGCAG',
                                            'Chrom': 'chrTest4',
                                            'End1Based': 15,
                                            'FilterReason': None,
                                            'FractionPureRepeatsRef': 1.0,
                                            'FractionPureRepeatsAlt': 1.0,
                                            'IsPureRepeat': False,
                                            'NumPureRepeatsAlt': 6,
                                            'NumPureRepeatsInVariant': 3,
                                            'NumPureRepeatsLeftFlank': 1,
                                            'NumPureRepeatsRef': 3,
                                            'NumPureRepeatsRightFlank': 2,
                                            'NumRepeatsAlt': 6,
                                            'NumRepeatsInVariant': 3,
                                            'NumRepeatsLeftFlank': 1,
                                            'NumRepeatsRef': 3,
                                            'NumRepeatsRightFlank': 2,
                                            'Pos': 9,
                                            'PureEnd1Based': 15,
                                            'PureStart1Based': 7,
                                            'Ref': 'G',
                                            'RepeatUnit': 'CAG',
                                            'MotifInterruptionIndex': 2,
                                            'Start1Based': 7})
        self.assertDictEqual(str_specs[1], {'Alt': 'GCAGCAGCATCAACACCAC',
                                            'Chrom': 'chrTest4',
                                            'End1Based': 15,
                                            'FilterReason': None,
                                            'FractionPureRepeatsRef': 1,
                                            'FractionPureRepeatsAlt': 0.5555555555555556,
                                            'IsPureRepeat': False,
                                            'NumPureRepeatsAlt': 5,
                                            'NumPureRepeatsInVariant': 2,
                                            'NumPureRepeatsLeftFlank': 1,
                                            'NumPureRepeatsRef': 3,
                                            'NumPureRepeatsRightFlank': 2,
                                            'NumRepeatsAlt': 9,
                                            'NumRepeatsInVariant': 6,
                                            'NumRepeatsLeftFlank': 1,
                                            'NumRepeatsRef': 3,
                                            'NumRepeatsRightFlank': 2,
                                            'Pos': 9,
                                            'PureEnd1Based': 15,
                                            'PureStart1Based': 7,
                                            'Ref': 'G',
                                            'RepeatUnit': 'CAG',
                                            'MotifInterruptionIndex': 2,
                                            'Start1Based': 7})

        self.assertDictEqual(counters, {'STR allele counts: INS': 2,
                                        'STR allele counts: TOTAL': 2,
                                        'STR allele delta: 3 repeats': 1,
                                        'STR allele delta: 6 repeats': 1,
                                        'STR allele motif size: 3 bp': 2,
                                        'STR allele reference repeats: with both left and right matching ref. repeat': 2,
                                        'STR allele size: 0-25bp': 2,
                                        'allele counts: INS alleles': 2})

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
            min_str_repeats=3, min_str_length=9, min_repeat_unit_length=1, max_repeat_unit_length=50,
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
        self.assertTrue(result["IsPureRepeat"])
        self.check_results_for_pure_repeats(result)

    def test_check_if_allele_is_str4_deletion(self):
        # insertion 2-STR
        counters = collections.defaultdict(int)

        result = check_if_allele_is_str(
            self.fasta_obj,
            "chrTest4", 9, "GCAG", "G",
            min_str_repeats=3, min_str_length=9, min_repeat_unit_length=1, max_repeat_unit_length=50,
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
        self.assertTrue(result["IsPureRepeat"])
        self.check_results_for_pure_repeats(result)

    def test_check_if_allele_is_str5_deletion(self):
        # insertion 2-STR
        counters = collections.defaultdict(int)

        result = check_if_allele_is_str(
            self.fasta_obj,
            "chrTest5", 6, "ACGCCGACGCCGCCGCCGC", "A",
            min_str_repeats=3, min_str_length=9, min_repeat_unit_length=1, max_repeat_unit_length=50,
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
        self.assertFalse(result["IsPureRepeat"])

        self.assertEqual(result["MotifInterruptionIndex"], 2)
        self.assertEqual(result["FractionPureRepeatsRef"], 10/11)
        self.assertEqual(result["FractionPureRepeatsAlt"], 5/5)

        self.assertEqual(result["PureStart1Based"], 7)
        self.assertEqual(result["PureEnd1Based"], 39)
        self.assertEqual(result["NumPureRepeatsInVariant"], 5)
        self.assertEqual(result["NumPureRepeatsLeftFlank"], 0)
        self.assertEqual(result["NumPureRepeatsRightFlank"], 5)

    def test_check_if_allele_is_str6_insertion(self):
        # insertion 2-STR
        counters = collections.defaultdict(int)

        result = check_if_allele_is_str(
            self.fasta_obj,
            "chrTest6", 6, "A", "ACAG",
            min_str_repeats=3, min_str_length=9, min_repeat_unit_length=1, max_repeat_unit_length=50,
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
        self.assertTrue(result["IsPureRepeat"])
        self.check_results_for_pure_repeats(result)

    def test_check_if_allele_is_str7_insertion(self):
        # insertion 2-STR
        counters = collections.defaultdict(int)

        result = check_if_allele_is_str(
            self.fasta_obj,
            "chrTest6", 6, "A", "ACAG",
            min_str_repeats=3, min_str_length=9, min_repeat_unit_length=1, max_repeat_unit_length=50,
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
        self.assertTrue(result["IsPureRepeat"])
        self.check_results_for_pure_repeats(result)

    def test_check_if_allele_is_str8_insertion(self):
        # insertion 2-STR
        counters = collections.defaultdict(int)

        result = check_if_allele_is_str(
            self.fasta_obj,
            "chrTest6", 6, "A", "ACAG",
            min_str_repeats=3, min_str_length=9, min_repeat_unit_length=1, max_repeat_unit_length=50,
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
        self.assertTrue(result["IsPureRepeat"])
        self.check_results_for_pure_repeats(result)

    def test_check_if_allele_is_str9_insertion(self):
        # insertion 2-STR
        counters = collections.defaultdict(int)

        # TTTTTCCCCC
        result = check_if_allele_is_str(
            self.fasta_obj,
            "chrTest7", 6, "C", "CCAGCAGCAG",
            min_str_repeats=3, min_str_length=9, min_repeat_unit_length=1, max_repeat_unit_length=50,
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
        self.assertTrue(result["IsPureRepeat"])
        self.check_results_for_pure_repeats(result)


    def test_check_if_allele_is_str10_insertion(self):
        # insertion 2-STR
        counters = collections.defaultdict(int)

        # TTTTTCAGCAGCAGCCCCC
        result = check_if_allele_is_str(
            self.fasta_obj,
            "chrTest8", 5, "TCAGCAGCAG", "T",
            min_str_repeats=3, min_str_length=9, min_repeat_unit_length=1, max_repeat_unit_length=50,
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
        self.assertTrue(result["IsPureRepeat"])
        self.check_results_for_pure_repeats(result)

    def tearDown(self):
        if os.path.isfile(self.temp_fasta_file.name):
            os.remove(self.temp_fasta_file.name)

