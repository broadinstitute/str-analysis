import os
import tempfile
import unittest

from str_analysis.convert_hipstr_vcf_to_expansion_hunter_json import process_hipstr_vcf

# The locus in every record below: chr1:101-120 (1-based inclusive), motif ATGGA, so 4 repeats in the reference.
MOTIF = "ATGGA"
LEFT_FLANK = "TTTTT"
RIGHT_FLANK = "CCCCC"
LOCUS_ID = "chr1-100-120-ATGGA"

VCF_HEADER = "\n".join([
    "##fileformat=VCFv4.1",
    "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]),
])


def record(pos_1based, ref, alt, genotype="0/1", start_1based=101, end_1based=120):
    return "\t".join(["chr1", str(pos_1based), LOCUS_ID, ref, alt, ".", ".",
                      f"START={start_1based};END={end_1based};PERIOD={len(MOTIF)};DP=30",
                      "GT:GB:Q:DP", f"{genotype}:0|0:1.00:30"])


class ProcessHipstrVcfTests(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    def variants_of(self, *records):
        """Run the converter over the given records and return {locus_id: variant record} for the loci it kept."""
        vcf_path = os.path.join(self.temp_dir.name, "test.vcf")
        with open(vcf_path, "wt") as f:
            f.write(VCF_HEADER + "\n" + "\n".join(records) + "\n")
        results = process_hipstr_vcf(vcf_path, sample_id="test")["LocusResults"]
        return {locus_id: locus["Variants"][locus_id] for locus_id, locus in results.items()}

    def genotype_of(self, *records):
        """Run the converter over the given records and return {locus_id: Genotype} for the loci it kept."""
        return {locus_id: variant["Genotype"] for locus_id, variant in self.variants_of(*records).items()}

    def test_record_spanning_the_locus_exactly(self):
        # POS == START and REF ends at END, so there is nothing to trim
        self.assertEqual(
            self.genotype_of(record(101, MOTIF * 4, MOTIF * 6)),
            {LOCUS_ID: "4/6"})

    def test_expanded_allele_keeps_its_expansion_when_the_record_is_trimmed(self):
        # HipSTR/LongTR extend the genotyped region 5bp past each end of the locus. Trimming has to remove those
        # 5bp from each end of each allele -- slicing every allele to the locus width instead would report the
        # 6-repeat allele as 4 repeats, i.e. as homozygous reference.
        self.assertEqual(
            self.genotype_of(record(96,
                                    LEFT_FLANK + MOTIF * 4 + RIGHT_FLANK,
                                    LEFT_FLANK + MOTIF * 6 + RIGHT_FLANK)),
            {LOCUS_ID: "4/6"})

    def test_contracted_allele_keeps_its_contraction_when_the_record_is_trimmed(self):
        self.assertEqual(
            self.genotype_of(record(96,
                                    LEFT_FLANK + MOTIF * 4 + RIGHT_FLANK,
                                    LEFT_FLANK + MOTIF * 2 + RIGHT_FLANK)),
            {LOCUS_ID: "2/4"})

    def test_allele_shorter_than_the_flanks_is_skipped(self):
        # a deletion that removes part of the flanks leaves an allele shorter than the 10bp being trimmed away,
        # so it can't be trimmed to the locus interval at all
        self.assertEqual(
            self.genotype_of(record(96, LEFT_FLANK + MOTIF * 4 + RIGHT_FLANK, "TTTTCCCC")),
            {})

    def test_symbolic_allele_is_a_zero_length_allele(self):
        # LongTR writes <DEL> for an allele spanning zero bases of the locus. Measuring the literal text instead
        # would report len("<DEL>") // 5 == 1 repeat.
        self.assertEqual(
            self.genotype_of(record(101, MOTIF * 4, "<DEL>")),
            {LOCUS_ID: "0/4"})

    def test_symbolic_allele_on_both_haplotypes(self):
        # with no allele sequence left, the repeat unit has to come from the reference allele
        variants = self.variants_of(record(101, MOTIF * 4, "<DEL>", genotype="1/1"))
        self.assertEqual({locus_id: v["Genotype"] for locus_id, v in variants.items()}, {LOCUS_ID: "0/0"})
        self.assertEqual(variants[LOCUS_ID]["RepeatUnit"], MOTIF)

    def test_symbolic_allele_is_not_trimmed(self):
        # The symbolic allele already spans zero bases of the locus, so the flanking-base trim doesn't apply to it,
        # while the real allele on the other haplotype keeps its 6 repeats. Pairing the symbolic allele against an
        # EXPANDED allele is what makes this discriminating: slicing every allele to the locus width instead would
        # report 0/4 here, which is also what a symbolic-vs-reference record would produce either way.
        self.assertEqual(
            self.genotype_of(record(96,
                                    LEFT_FLANK + MOTIF * 4 + RIGHT_FLANK,
                                    "<DEL>," + LEFT_FLANK + MOTIF * 6 + RIGHT_FLANK,
                                    genotype="1/2")),
            {LOCUS_ID: "0/6"})

    def test_allele_shorter_than_the_flanks_is_skipped_even_when_paired_with_a_symbolic_allele(self):
        # a symbolic allele is "" and so always sorts first, so checking only the shorter allele would leave the
        # real 8bp allele unchecked and silently trim it away to "" (reporting a spurious 0/0)
        self.assertEqual(
            self.genotype_of(record(96, LEFT_FLANK + MOTIF * 4 + RIGHT_FLANK, "<DEL>,TTTTCCCC", genotype="1/2")),
            {})

    def test_record_that_does_not_span_the_locus_is_skipped(self):
        # POS is inside the locus, so the record is missing the locus's first 5bp
        self.assertEqual(
            self.genotype_of(record(106, MOTIF * 3, MOTIF * 5)),
            {})


if __name__ == "__main__":
    unittest.main()
