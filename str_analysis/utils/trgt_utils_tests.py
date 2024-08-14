import collections
import os
from pprint import pformat, pprint
import pyfaidx
import tempfile
import unittest

from str_analysis.utils.trgt_utils import convert_trgt_locus_to_expansion_hunter_format


class Tests(unittest.TestCase):

    def setUp(self):

        self.maxDiff = None
        self.temp_fasta_file = tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False)

        seq1 = "TG"*5
        spacer1 = "AGACTAGCGAGTAGC"
        seq2 = "CAG"*10
        spacer2 = spacer1[::-1]
        seq3 = "G"*15

        # ---------------------  chr1
        self.temp_fasta_file.write(">chr1\n")
        self.temp_fasta_file.write(f"{spacer1}{seq1}{seq2}{seq3}{spacer2}\n")

        self.chr1_start_0based = [
            len(spacer1),
            len(spacer1) + len(seq1),
            len(spacer1) + len(seq1) + len(seq2),
        ]

        self.chr1_end_1based = [
            len(spacer1) + len(seq1),
            len(spacer1) + len(seq1) + len(seq2),
            len(spacer1) + len(seq1) + len(seq2) + len(seq3),
        ]

        self.chr1_motifs = ["TG", "CAG", "G"]
        self.chr1_locus_ids = [f"chr1-{self.chr1_start_0based[i]}-{self.chr1_end_1based[i]}-{self.chr1_motifs[i]}" for i in range(len(self.chr1_start_0based))]
        self.chr1_reference_regions = [f"chr1:{self.chr1_start_0based[i]}-{self.chr1_end_1based[i]}" for i in range(len(self.chr1_start_0based))]
        self.chr1_trgt_locus_struc = f"{spacer1}({self.chr1_motifs[0]})n({self.chr1_motifs[1]})n({self.chr1_motifs[2]})n{spacer2}"

        # ---------------------  chr2
        self.temp_fasta_file.write(">chr2\n")
        self.temp_fasta_file.write(f"{spacer1}{seq1}{spacer2}{seq2}{seq3}\n")

        self.chr2_start_0based = [
            len(spacer1),
            len(spacer1) + len(seq1) + len(spacer2),
            len(spacer1) + len(seq1) + len(spacer2) + len(seq2),
        ]
        self.chr2_end_1based = [
            len(spacer1) + len(seq1) + len(spacer2),
            len(spacer1) + len(seq1) + len(spacer2) + len(seq2),
            len(spacer1) + len(seq1) + len(spacer2) + len(seq2) + len(seq3),
        ]
        self.chr2_trgt_locus_struc = f"{spacer1}({self.chr1_motifs[0]})n({self.chr1_motifs[1]})n({self.chr1_motifs[2]})n{spacer2}"
        self.chr2_locus_ids = [f"chr2-{self.chr1_start_0based[0]}-{self.chr1_end_1based[0]}-{self.chr1_motifs[0]}"]
        self.chr2_reference_regions = [f"chr2:{self.chr1_start_0based[0]}-{self.chr1_end_1based[0]}"]

        # --------------------- chr3
        seq4 = "GTG"*5
        seq5 = "CAG"*5
        seq6 = "TTC"*5
        self.temp_fasta_file.write(">chr3\n")
        self.temp_fasta_file.write(f"{seq4}{seq5}{seq6}\n")

        self.temp_fasta_file.close()

        self.fasta_obj = pyfaidx.Fasta(self.temp_fasta_file.name, one_based_attributes=False, as_raw=True)


    def test_convert_trgt_locus_to_expansion_hunter_format(self):

        expansion_hunter_specs = list(convert_trgt_locus_to_expansion_hunter_format(
            self.fasta_obj, "chr1", 0, len(self.fasta_obj["chr1"]), self.chr1_trgt_locus_struc, verbose=True))
        self.assertEqual(3, len(expansion_hunter_specs))

        for i in range(len(expansion_hunter_specs)):
            self.assertEqual(expansion_hunter_specs[i], {
                "LocusId": self.chr1_locus_ids[i],
                "ReferenceRegion": self.chr1_reference_regions[i],
                "LocusStructure": f"({self.chr1_motifs[i]})*",
                "VariantType": "Repeat",
            })



    def tearDown(self):
        if os.path.isfile(self.temp_fasta_file.name):
            os.remove(self.temp_fasta_file.name)

