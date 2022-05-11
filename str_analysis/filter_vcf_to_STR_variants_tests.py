import os
import pyfaidx
import tempfile
import unittest

from str_analysis.filter_vcf_to_STR_variants import get_flanking_reference_sequences


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
        self.temp_fasta_file.close()

        self.fasta_obj = pyfaidx.Fasta(self.temp_fasta_file.name, one_based_attributes=False, as_raw=True)

    def test_get_containing_sequences1(self):
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
        self.assertEqual(right_sequence, "ACACACTGTG")

        # insertion 3-STR
        left_sequence, variant_sequence, right_sequence = get_flanking_reference_sequences(self.fasta_obj, "chrTest3", 15, self.test3_seq1[-1], self.test3_seq1[-1] + self.test3_seq2)
        self.assertEqual(left_sequence, self.test3_seq1)
        self.assertEqual(variant_sequence, self.test3_seq2)
        self.assertEqual(right_sequence, self.test3_seq2 + self.test3_seq3)

    def tearDown(self):
        if os.path.isfile(self.temp_fasta_file.name):
            os.remove(self.temp_fasta_file.name)
