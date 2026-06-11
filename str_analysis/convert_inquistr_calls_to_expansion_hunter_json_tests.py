import tempfile
import unittest

from str_analysis.convert_inquistr_calls_to_expansion_hunter_json import (
    process_inquistr_calls, parse_allele_size_in_repeats)

INQUISTR_CONTENTS1 = """# file_type=individual_call
# version=0.27.1
# command=call
# sample=HG002
chromosome\tbegin\tend\tinfo\tHG002_H1\tHG002_H2
chr1\t591733\t591751\tchr1-591733-591751-A\t3\t-2
chr1\t650383\t650423\tchr1-650383-650423-AC\t4\tNaN
chr2\t100\t140\tchr2-100-140-AC\tNaN\tNaN
chr3\t200\t218\tchr3-200-218-A\t0\t0
chr4\t300\t340\tchr4-300-340-AC\t2.5\t2.5
chr5\t400\t418\tchr5-400-418-A\t-50\t-50"""


class Tests(unittest.TestCase):

    def setUp(self):
        self.test_table1 = tempfile.NamedTemporaryFile(mode="wt", suffix=".inq", delete=False)
        self.test_table1.write(INQUISTR_CONTENTS1)
        self.test_table1.flush()
        self.test_table1.seek(0)

    def test_parse_allele_size_in_repeats(self):
        # homopolymer (motif size 1), reference is 18bp
        self.assertEqual(parse_allele_size_in_repeats("3", 18, 1), 21)
        self.assertEqual(parse_allele_size_in_repeats("-2", 18, 1), 16)
        self.assertEqual(parse_allele_size_in_repeats("0", 18, 1), 18)
        # contraction below the reference is clamped at 0 repeats
        self.assertEqual(parse_allele_size_in_repeats("-50", 18, 1), 0)
        # no-call values
        self.assertIsNone(parse_allele_size_in_repeats("NaN", 18, 1))
        self.assertIsNone(parse_allele_size_in_repeats("nan", 18, 1))
        self.assertIsNone(parse_allele_size_in_repeats(".", 18, 1))
        self.assertIsNone(parse_allele_size_in_repeats(None, 18, 1))
        # AC motif (size 2), reference is 40bp; the repeat count is floored to a whole repeat
        self.assertEqual(parse_allele_size_in_repeats("4", 40, 2), 22)     # 44 / 2 = 22
        self.assertEqual(parse_allele_size_in_repeats("2.5", 40, 2), 21)   # 42.5 / 2 = 21.25 -> 21
        self.assertEqual(parse_allele_size_in_repeats("3.5", 40, 2), 21)   # 43.5 / 2 = 21.75 -> 21

    def test_discard_hom_ref(self):
        results = process_inquistr_calls(self.test_table1.name, discard_hom_ref=True)
        self.assertEqual(results["SampleParameters"]["SampleId"], "HG002")

        # chr2 (both NaN, no-call) and chr3 (hom-ref) should be discarded
        self.assertEqual(set(results["LocusResults"].keys()), {
            "chr1-591733-591751-A",
            "chr1-650383-650423-AC",
            "chr4-300-340-AC",
            "chr5-400-418-A",
        })

        v = results["LocusResults"]["chr1-591733-591751-A"]["Variants"]["chr1-591733-591751-A"]
        self.assertEqual(v["Genotype"], "16/21")
        self.assertEqual(v["GenotypeConfidenceInterval"], "16-16/21-21")
        self.assertEqual(v["ReferenceRegion"], "chr1:591733-591751")
        self.assertEqual(v["RepeatUnit"], "A")
        self.assertEqual(results["LocusResults"]["chr1-591733-591751-A"]["AlleleCount"], 2)

        # one haplotype is NaN -> single-allele genotype
        v = results["LocusResults"]["chr1-650383-650423-AC"]["Variants"]["chr1-650383-650423-AC"]
        self.assertEqual(v["Genotype"], "22")
        self.assertEqual(results["LocusResults"]["chr1-650383-650423-AC"]["AlleleCount"], 1)

        # non-integer median bp diff
        self.assertEqual(
            results["LocusResults"]["chr4-300-340-AC"]["Variants"]["chr4-300-340-AC"]["Genotype"], "21/21")

        # contraction clamped at 0; not hom-ref since reference is 18 repeats
        self.assertEqual(
            results["LocusResults"]["chr5-400-418-A"]["Variants"]["chr5-400-418-A"]["Genotype"], "0/0")

    def test_keep_hom_ref(self):
        results = process_inquistr_calls(self.test_table1.name, discard_hom_ref=False)
        # chr3 hom-ref is kept, only the both-NaN no-call (chr2) is dropped
        self.assertIn("chr3-200-218-A", results["LocusResults"])
        self.assertNotIn("chr2-100-140-AC", results["LocusResults"])
        self.assertEqual(
            results["LocusResults"]["chr3-200-218-A"]["Variants"]["chr3-200-218-A"]["Genotype"], "18/18")


if __name__ == "__main__":
    unittest.main()
