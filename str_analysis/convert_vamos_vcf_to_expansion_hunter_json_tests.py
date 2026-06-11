import tempfile
import unittest

from str_analysis.convert_vamos_vcf_to_expansion_hunter_json import (
    process_vamos_vcf, parse_info, haplotype_num_repeats)

# One header line, then loci covering: homozygous (GT 1/1, only H1), heterozygous (GT 1/2), a hom-ref locus,
# a multi-motif VNTR locus (should be skipped), a het locus whose haplotypes need sorting, and a locus that
# reports ALTANNO but no LEN (exercises the LEN fallback). POS is 1-based, so start_0based = POS - 1.
VAMOS_VCF_CONTENTS1 = """##fileformat=VCFv4.1
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
##INFO=<ID=RU,Number=1,Type=String,Description="Repeat units">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=ALTANNO_H1,Number=1,Type=String,Description="Motif representation for the h1 alternate allele">
##INFO=<ID=ALTANNO_H2,Number=1,Type=String,Description="Motif representation for the h2 alternate allele">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002
chr1\t591734\t.\tN\t<VNTR>\t.\tPASS\tEND=591751;RU=A;SVTYPE=VNTR;ALTANNO_H1=0,0,0;LEN_H1=21;\tGT\t1/1
chr2\t650384\t.\tN\t<VNTR>\t.\tPASS\tEND=650423;RU=AC;SVTYPE=VNTR;ALTANNO_H1=0,0,0;LEN_H1=22;ALTANNO_H2=0,0,0;LEN_H2=20;\tGT\t1/2
chr3\t201\t.\tN\t<VNTR>\t.\tPASS\tEND=218;RU=A;SVTYPE=VNTR;ALTANNO_H1=0,0,0;LEN_H1=18;\tGT\t1/1
chr4\t301\t.\tN\t<VNTR>\t.\tPASS\tEND=340;RU=AC,GT;SVTYPE=VNTR;ALTANNO_H1=0,1;LEN_H1=2;\tGT\t1/1
chr5\t401\t.\tN\t<VNTR>\t.\tPASS\tEND=418;RU=A;SVTYPE=VNTR;ALTANNO_H1=0,0,0;LEN_H1=30;ALTANNO_H2=0,0,0;LEN_H2=10;\tGT\t1/2
chr6\t501\t.\tN\t<VNTR>\t.\tPASS\tEND=520;RU=A;SVTYPE=VNTR;ALTANNO_H1=0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;\tGT\t1/1"""


class Tests(unittest.TestCase):

    def setUp(self):
        self.test_vcf1 = tempfile.NamedTemporaryFile(mode="wt", suffix=".vcf", delete=False)
        self.test_vcf1.write(VAMOS_VCF_CONTENTS1)
        self.test_vcf1.flush()
        self.test_vcf1.seek(0)

    def test_parse_info(self):
        info = parse_info("END=591751;RU=A;SVTYPE=VNTR;ALTANNO_H1=0,0,0;LEN_H1=21;")
        self.assertEqual(info["END"], "591751")
        self.assertEqual(info["RU"], "A")
        self.assertEqual(info["LEN_H1"], "21")
        self.assertNotIn("LEN_H2", info)

    def test_haplotype_num_repeats(self):
        # LEN is used when present
        self.assertEqual(haplotype_num_repeats({"LEN_H1": "21"}, "H1"), 21)
        # an absent haplotype returns None
        self.assertIsNone(haplotype_num_repeats({"LEN_H1": "21"}, "H2"))
        # fall back to counting ALTANNO indices when LEN is missing (',' separator for --read mode)
        self.assertEqual(haplotype_num_repeats({"ALTANNO_H1": "0,0,0,0"}, "H1"), 4)
        # '-' separator (as produced by --contig mode) is also handled
        self.assertEqual(haplotype_num_repeats({"ALTANNO_H1": "0-0-0"}, "H1"), 3)

    def test_discard_hom_ref(self):
        results = process_vamos_vcf(self.test_vcf1.name, discard_hom_ref=True)
        self.assertEqual(results["SampleParameters"]["SampleId"], "HG002")

        # chr3 (hom-ref) and chr4 (multi-motif VNTR) are dropped
        self.assertEqual(set(results["LocusResults"].keys()), {
            "chr1-591733-591751-A",
            "chr2-650383-650423-AC",
            "chr5-400-418-A",
            "chr6-500-520-A",
        })

        v = results["LocusResults"]["chr1-591733-591751-A"]["Variants"]["chr1-591733-591751-A"]
        self.assertEqual(v["Genotype"], "21/21")
        self.assertEqual(v["GenotypeConfidenceInterval"], "21-21/21-21")
        self.assertEqual(v["ReferenceRegion"], "chr1:591733-591751")
        self.assertEqual(v["RepeatUnit"], "A")
        self.assertEqual(results["LocusResults"]["chr1-591733-591751-A"]["AlleleCount"], 2)

        # heterozygous AC locus
        self.assertEqual(
            results["LocusResults"]["chr2-650383-650423-AC"]["Variants"]["chr2-650383-650423-AC"]["Genotype"],
            "20/22")

        # haplotypes are sorted ascending (H2=10 < H1=30)
        self.assertEqual(
            results["LocusResults"]["chr5-400-418-A"]["Variants"]["chr5-400-418-A"]["Genotype"], "10/30")

        # LEN absent -> count is taken from the 23 ALTANNO_H1 indices
        self.assertEqual(
            results["LocusResults"]["chr6-500-520-A"]["Variants"]["chr6-500-520-A"]["Genotype"], "23/23")

    def test_keep_hom_ref(self):
        results = process_vamos_vcf(self.test_vcf1.name, discard_hom_ref=False)
        # the hom-ref locus is kept, the multi-motif VNTR is still skipped
        self.assertIn("chr3-200-218-A", results["LocusResults"])
        self.assertNotIn("chr4-300-340-AC", results["LocusResults"])
        self.assertEqual(
            results["LocusResults"]["chr3-200-218-A"]["Variants"]["chr3-200-218-A"]["Genotype"], "18/18")


if __name__ == "__main__":
    unittest.main()
