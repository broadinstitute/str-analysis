import gzip
import tempfile
import unittest

from str_analysis.convert_atarva_vcf_to_expansion_hunter_json import (
    process_atarva_vcf, parse_info, parse_sample)

# FORMAT is always GT:AL:CN:LPM:AR:SD:DP:SN:SQ:MA:MR:DS:MV; the converter reads only CN from the SAMPLE column and
# MOTIF/START/END from INFO, so the remaining FORMAT fields are filled with "." placeholders below.
_FORMAT = "GT:AL:CN:LPM:AR:SD:DP:SN:SQ:MA:MR:DS:MV"
_PAD = ":.:.:.:.:.:.:.:.:.:."  # the 10 trailing FORMAT fields after GT:AL:CN

ATARVA_VCF1 = "\n".join([
    "##fileformat=VCFv4.2",
    "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "HG002"]),
    # heterozygous, CN given out of order (20,17) -> alleles sorted ascending to "17/20"
    "\t".join(["chr1", "15797", ".", "N", "A", "0", "PASS",
               "AC=1;AN=2;MOTIF=CTT;START=15796;END=15849;ID=.;REFCN=17",
               _FORMAT, "0|1:60,51:20,17" + _PAD]),
    # homozygous expansion (CN 90,90 vs reference 78) -> kept as "90/90"
    "\t".join(["chr1", "11227", ".", "N", "A", "0", "PASS",
               "AC=2;AN=2;MOTIF=AGG;START=11226;END=11462;ID=.;REFCN=78",
               _FORMAT, "1/1:270,270:90,90" + _PAD]),
    # FILTER=LESS_READS no-call (CN ".") -> skipped
    "\t".join(["chr2", "20000", ".", "N", ".", "0", "LESS_READS",
               "AC=0;AN=0;MOTIF=AT;START=19998;END=20040;ID=.;REFCN=21",
               _FORMAT, "." + ":." * 12]),
    # multizygous, three alleles (CN 10,12,14) -> can't be a diploid genotype -> skipped
    "\t".join(["chr3", "500", ".", "N", "A", "0", "PASS",
               "AC=2;AN=2;MOTIF=CAG;START=498;END=540;ID=.;REFCN=14",
               _FORMAT, "1/2:.:10,12,14" + _PAD]),
    # hom-ref (CN 20,20 == reference 20) -> discarded when discard_hom_ref=True, kept otherwise
    "\t".join(["chr4", "101", ".", "N", "A", "0", "PASS",
               "AC=0;AN=2;MOTIF=AT;START=100;END=140;ID=.;REFCN=20",
               _FORMAT, "0/0:40,40:20,20" + _PAD]),
    # hemizygous / haploid call (single CN value, e.g. male chrY) -> single-allele genotype "15"
    "\t".join(["chrY", "301", ".", "N", "A", "0", "PASS",
               "AC=1;AN=1;MOTIF=CAG;START=300;END=330;ID=.;REFCN=10",
               _FORMAT, "1:45:15" + _PAD]),
])


class Tests(unittest.TestCase):

    def setUp(self):
        self.vcf1 = tempfile.NamedTemporaryFile(mode="wt", suffix=".vcf", delete=False)
        self.vcf1.write(ATARVA_VCF1)
        self.vcf1.flush()
        self.vcf1.seek(0)

    def test_parse_info(self):
        info = parse_info("AC=1;AN=2;MOTIF=CTT;START=15796;END=15849;ID=.;REFCN=17")
        self.assertEqual(info["MOTIF"], "CTT")
        self.assertEqual(info["START"], "15796")
        self.assertEqual(info["END"], "15849")
        # a flag field (no '=') maps to True
        self.assertIs(parse_info("SOMEFLAG")["SOMEFLAG"], True)

    def test_parse_sample(self):
        sample = parse_sample("GT:AL:CN", "0|1:60,51:20,17")
        self.assertEqual(sample["GT"], "0|1")
        self.assertEqual(sample["AL"], "60,51")
        self.assertEqual(sample["CN"], "20,17")

    def test_discard_hom_ref(self):
        results = process_atarva_vcf(self.vcf1.name, discard_hom_ref=True)

        # sample id is read from the last column of the #CHROM header
        self.assertEqual(results["SampleParameters"]["SampleId"], "HG002")
        self.assertIsNone(results["SampleParameters"]["Sex"])

        # no-call (chr2), multizygous (chr3), and hom-ref (chr4) are all dropped
        self.assertEqual(set(results["LocusResults"].keys()), {
            "chr1-15796-15849-CTT",
            "chr1-11226-11462-AGG",
            "chrY-300-330-CAG",
        })

        # heterozygous: CN 20,17 -> alleles sorted ascending
        het = results["LocusResults"]["chr1-15796-15849-CTT"]
        v = het["Variants"]["chr1-15796-15849-CTT"]
        self.assertEqual(v["Genotype"], "17/20")
        self.assertEqual(v["GenotypeConfidenceInterval"], "17-17/20-20")
        self.assertEqual(v["ReferenceRegion"], "chr1:15796-15849")
        self.assertEqual(v["RepeatUnit"], "CTT")
        self.assertEqual(v["VariantId"], "chr1-15796-15849-CTT")
        self.assertEqual(v["VariantType"], "Repeat")
        self.assertEqual(het["AlleleCount"], 2)

        # homozygous expansion
        v = results["LocusResults"]["chr1-11226-11462-AGG"]["Variants"]["chr1-11226-11462-AGG"]
        self.assertEqual(v["Genotype"], "90/90")
        self.assertEqual(results["LocusResults"]["chr1-11226-11462-AGG"]["AlleleCount"], 2)

        # hemizygous / haploid: a single CN value yields a single-allele genotype
        hemi = results["LocusResults"]["chrY-300-330-CAG"]
        v = hemi["Variants"]["chrY-300-330-CAG"]
        self.assertEqual(v["Genotype"], "15")
        self.assertEqual(v["GenotypeConfidenceInterval"], "15-15")
        self.assertEqual(hemi["AlleleCount"], 1)

    def test_keep_hom_ref(self):
        results = process_atarva_vcf(self.vcf1.name, discard_hom_ref=False)

        # the hom-ref locus is now kept
        self.assertIn("chr4-100-140-AT", results["LocusResults"])
        self.assertEqual(
            results["LocusResults"]["chr4-100-140-AT"]["Variants"]["chr4-100-140-AT"]["Genotype"], "20/20")

        # no-call and multizygous loci are still dropped regardless of discard_hom_ref
        self.assertNotIn("chr2-19998-20040-AT", results["LocusResults"])
        self.assertNotIn("chr3-498-540-CAG", results["LocusResults"])

    def test_sample_id_argument_overrides_header(self):
        results = process_atarva_vcf(self.vcf1.name, sample_id="OVERRIDE", discard_hom_ref=True)
        self.assertEqual(results["SampleParameters"]["SampleId"], "OVERRIDE")

    def test_reads_gzip_compressed_vcf(self):
        # The pipeline always bgzips the VCF before invoking the converter, so a ".gz" path takes the gzip.open branch.
        # Reading a gzip-compressed fixture must yield the same result as the equivalent plain-text VCF.
        gz = tempfile.NamedTemporaryFile(mode="wb", suffix=".vcf.gz", delete=False)
        gz.write(gzip.compress(ATARVA_VCF1.encode()))
        gz.flush()
        gz.close()

        results = process_atarva_vcf(gz.name, discard_hom_ref=True)
        self.assertEqual(results["SampleParameters"]["SampleId"], "HG002")
        self.assertEqual(set(results["LocusResults"].keys()), {
            "chr1-15796-15849-CTT",
            "chr1-11226-11462-AGG",
            "chrY-300-330-CAG",
        })
        self.assertEqual(
            results["LocusResults"]["chr1-15796-15849-CTT"]["Variants"]["chr1-15796-15849-CTT"]["Genotype"], "17/20")


if __name__ == "__main__":
    unittest.main()
