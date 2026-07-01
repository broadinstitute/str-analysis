import simplejson as json
import tempfile
import unittest

from str_analysis.convert_ensembletr_vcf_to_expansion_hunter_json import (
    load_catalog_loci, reconcile_locus, process_ensembletr_vcf)

# A catalog with three loci. The second locus' motif (GTGTGGAAACTGCGACACTCACG) is intentionally NOT the alphabetically
# canonical rotation, so the test can confirm the truth-set motif (not EnsembleTR's canonicalized RU) is restored.
CATALOG = [
    {"LocusId": "1-591733-591751-A", "ReferenceRegion": "chr1:591733-591751", "LocusStructure": "(A)*"},
    {"LocusId": "1-597794-598734-GTGTGGAAACTGCGACACTCACG",
     "ReferenceRegion": "chr1:597794-598734", "LocusStructure": "(GTGTGGAAACTGCGACACTCACG)*"},
    {"LocusId": "1-650383-650423-AC", "ReferenceRegion": "chr1:650383-650423", "LocusStructure": "(AC)*"},
]

# EnsembleTR-style records. START/END are 1-based; END = end_1based + 1 (matching observed EnsembleTR output). The first
# record's RU is a canonical rotation that differs from the catalog motif. The 597795 record uses coordinates shifted by
# 1bp to exercise the overlap-based reconciliation fallback. The last record is a no-call (NCOPY ".").
VCF_CONTENTS = """##fileformat=VCFv4.1
##INFO=<ID=START,Number=1,Type=Integer,Description="x">
##INFO=<ID=END,Number=1,Type=Integer,Description="x">
##INFO=<ID=PERIOD,Number=1,Type=Integer,Description="x">
##INFO=<ID=RU,Number=1,Type=String,Description="x">
##INFO=<ID=METHODS,Number=1,Type=String,Description="x">
##FORMAT=<ID=GT,Number=1,Type=String,Description="x">
##FORMAT=<ID=NCOPY,Number=1,Type=String,Description="x">
##FORMAT=<ID=SCORE,Number=1,Type=Float,Description="x">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002
chr1\t591734\t.\tA\t<STR20>\t.\t.\tSTART=591734;END=591752;PERIOD=1;RU=A;METHODS=0|1|0|1\tGT:NCOPY:SCORE\t1/1:97.0,20.0:0.47
chr1\t597796\t.\tG\t<STR8>\t.\t.\tSTART=597796;END=598736;PERIOD=23;RU=AAACTGCGACACTCACGGTGTGG;METHODS=0|1|0|0\tGT:NCOPY:SCORE\t1/1:8.0,8.0:0.95
chr1\t650384\t.\tAC\t<STR18>\t.\t.\tSTART=650384;END=650424;PERIOD=2;RU=AC;METHODS=0|0|1|0\tGT:NCOPY:SCORE\t1/2:20.0,18.0:1.0
chr1\t900000\t.\tA\t<STR5>\t.\t.\tSTART=900000;END=900010;PERIOD=1;RU=A;METHODS=0|1|0|0\tGT:NCOPY:SCORE\t./.:.:."""


class Tests(unittest.TestCase):

    def setUp(self):
        self.catalog_file = tempfile.NamedTemporaryFile(mode="wt", suffix=".json", delete=False)
        json.dump(CATALOG, self.catalog_file)
        self.catalog_file.flush()
        self.vcf_file = tempfile.NamedTemporaryFile(mode="wt", suffix=".vcf", delete=False)
        self.vcf_file.write(VCF_CONTENTS)
        self.vcf_file.flush()
        self.exact_lookup, self.per_chrom_intervals = load_catalog_loci(self.catalog_file.name)

    def test_load_catalog_motif_from_locus_id(self):
        self.assertEqual(self.exact_lookup[("1", 591733, 591751)],
                         ("1-591733-591751-A", "A"))
        self.assertEqual(self.exact_lookup[("1", 597794, 598734)][1], "GTGTGGAAACTGCGACACTCACG")

    def test_reconcile_exact_and_overlap(self):
        # exact match
        self.assertEqual(
            reconcile_locus("1", 591733, 591751, 1, self.exact_lookup, self.per_chrom_intervals)[0],
            "1-591733-591751-A")
        # overlap fallback (record shifted +1bp, same period 23)
        match = reconcile_locus("1", 597795, 598735, 23, self.exact_lookup, self.per_chrom_intervals)
        self.assertIsNotNone(match)
        self.assertEqual(match[0], "1-597794-598734-GTGTGGAAACTGCGACACTCACG")
        # no overlapping locus -> None
        self.assertIsNone(reconcile_locus("1", 900000, 900010, 1, self.exact_lookup, self.per_chrom_intervals))

    def test_process_vcf(self):
        results, counters = process_ensembletr_vcf(
            self.vcf_file.name, self.exact_lookup, self.per_chrom_intervals)
        lr = results["LocusResults"]
        self.assertEqual(results["SampleParameters"]["SampleId"], "HG002")
        # 3 reconciled loci; the 4th record is a no-call and is skipped
        self.assertEqual(len(lr), 3)
        self.assertEqual(counters["no-call records skipped"], 1)

        # canonical-rotation RU is replaced by the truth-set motif from the catalog
        v = lr["1-597794-598734-GTGTGGAAACTGCGACACTCACG"]["Variants"]["1-597794-598734-GTGTGGAAACTGCGACACTCACG"]
        self.assertEqual(v["RepeatUnit"], "GTGTGGAAACTGCGACACTCACG")
        self.assertEqual(v["ReferenceRegion"], "chr1:597794-598734")

        # NCOPY (unsorted "97.0,20.0") becomes sorted short->long repeat counts; SCORE becomes Q; CI is zero-width
        a = lr["1-591733-591751-A"]["Variants"]["1-591733-591751-A"]
        self.assertEqual(a["Genotype"], "20/97")
        self.assertEqual(a["GenotypeConfidenceInterval"], "20-20/97-97")
        self.assertAlmostEqual(a["Q"], 0.47)

        # het call sorted short->long
        self.assertEqual(lr["1-650383-650423-AC"]["Variants"]["1-650383-650423-AC"]["Genotype"], "18/20")


if __name__ == "__main__":
    unittest.main()
