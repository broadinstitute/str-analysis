import gzip
import os
import tempfile
import unittest
from unittest import mock

from str_analysis.convert_multisample_LPS_table_to_allele_frequency_histograms import (
    compute_histograms, compute_row, load_vcf_trid_metadata, main)


def run_converter(lps_table_contents, extra_args=()):
    """Run the converter's main() on an LPS table and return its output rows.

    Args:
        lps_table_contents (str): contents of a headerless LPS table (trid, motif, sample columns).
        extra_args (list): additional command line args to pass to main().

    Returns:
        list: one dict per output row, mapping column name to the unparsed string value.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        input_path = os.path.join(temp_dir, "lps.tsv")
        with open(input_path, "wt") as f:
            f.write(lps_table_contents)

        argv = ["convert_multisample_LPS_table_to_allele_frequency_histograms",
                "--no-header", "--input-table", input_path] + list(extra_args)
        with mock.patch("sys.argv", argv):
            main()

        output_paths = [p for p in os.listdir(temp_dir) if p.endswith(".tsv.gz")]
        if len(output_paths) != 1:
            raise ValueError(f"Expected 1 output file, found {output_paths}")
        with gzip.open(os.path.join(temp_dir, output_paths[0]), "rt") as f:
            header = next(f).rstrip("\n").split("\t")
            return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in f]


def run_converter_with_trid_metadata(lps_table_contents, trid_metadata_contents):
    """Run the converter's main() with --vcf-trid-metadata-tsv and return its output rows.

    Args:
        lps_table_contents (str): contents of a headerless LPS table (trid, motif, sample columns).
        trid_metadata_contents (str): contents of the TRID metadata TSV, including its header line.

    Returns:
        list: one dict per output row, mapping column name to the unparsed string value.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        trid_metadata_path = os.path.join(temp_dir, "trid_metadata.tsv")
        with open(trid_metadata_path, "wt") as f:
            f.write(trid_metadata_contents)

        return run_converter(lps_table_contents, ["--vcf-trid-metadata-tsv", trid_metadata_path])


TRID_METADATA = (
    "trid\tlocus_id\tmotif\tinterval\tvc\n"
    "1-44835-44876-AAAT\t1-44835-44876-AAAT\tAAAT\t1:44835-44876\t\n")

# the chrX:149631602-149631762 cluster as the AoU (Danzi) VCF holds it: two records covering the
# same TCC repeat over the same span, under different TRIDs
TMEM185A_GENE_NAMED_TRID = ("X-149631602-149631617-TCC,X-149631685-149631694-GCT,TMEM185A_CGCCGT,"
                            "X-149631729-149631732-CGC")
TMEM185A_COORDINATE_TRID = ("X-149631602-149631617-TCC,X-149631685-149631694-GCT,"
                            "X-149631723-149631735-CGCCGT")
TMEM185A_TRID_METADATA = (
    "trid\tlocus_id\tmotif\tinterval\tvc\n"
    f"{TMEM185A_GENE_NAMED_TRID}\tX-149631602-149631617-TCC\tTCC\tX:149631602-149631762\tX:149631602-149631762\n"
    f"{TMEM185A_COORDINATE_TRID}\tX-149631602-149631617-TCC\tTCC\tX:149631602-149631762\tX:149631602-149631762\n")


class Tests(unittest.TestCase):

    def test_compute_histograms_excludes_partially_called_samples(self):
        alleles_by_sample_id = {"s1": [8, 8], "s2": [8], "s3": [9, 10]}

        # without the exclusion, s2's single allele is doubled into a 2nd 8/8 genotype
        self.assertEqual(
            compute_histograms([8, 8, 8, 9, 10], alleles_by_sample_id),
            ("8x:3,9x:1,10x:1", "8/8:2,9/10:1"))

        # s2 was "8,." rather than a hemizygous call, so its genotype wasn't fully observed
        self.assertEqual(
            compute_histograms([8, 8, 8, 9, 10], alleles_by_sample_id, {"s2"}),
            ("8x:3,9x:1,10x:1", "8/8:1,9/10:1"))

    def test_compute_row_excludes_partially_called_samples_from_short_alleles(self):
        row = compute_row("1-44835-44876-AAAT", "AAAT", [8, 8, 8, 9, 10],
                          {"s1": [8, 8], "s2": [8], "s3": [9, 10]},
                          partially_called_sample_ids={"s2"})

        self.assertEqual(row["AlleleSizeHistogram"], "8x:3,9x:1,10x:1")
        self.assertEqual(row["BiallelicHistogram"], "8/8:1,9/10:1")
        self.assertEqual(row["NumCalledAlleles"], 5)
        self.assertEqual(row["ShortAlleleMax"], 9)  # s2's 8 is not known to be its shorter allele

    def test_compute_row_with_only_partially_called_samples(self):
        row = compute_row("1-70000-70020-AC", "AC", [8, 9], {"s1": [8], "s2": [9]},
                          partially_called_sample_ids={"s1", "s2"})

        self.assertEqual(row["AlleleSizeHistogram"], "8x:1,9x:1")
        self.assertEqual(row["NumCalledAlleles"], 2)
        self.assertEqual(row["BiallelicHistogram"], "")
        self.assertEqual(row["ShortAlleleMax"], "")
        self.assertEqual(row["ShortAllele99thPercentile"], "")

    def test_hemi_allele_columns_count_only_hemizygous_male_calls(self):
        # s1 is a true hemizygous call, s2 a male partial call, s3 a male diploid call in a
        # pseudoautosomal region, s4 a female diploid call. Only s1 is hemizygous.
        row = compute_row("X-1000-1020-CAG", "CAG", [5, 6, 7, 9, 8, 11],
                          {"s1": [5], "s2": [6], "s3": [7, 9], "s4": [8, 11]},
                          sample_id_to_sex={"s1": "male", "s2": "male", "s3": "male", "s4": "female"},
                          partially_called_sample_ids={"s2"})

        self.assertEqual(row["HemiAlleleMax"], 5)
        self.assertEqual(row["HemiAllele99thPercentile"], 5)

        # every measured allele still counts in the per-allele histogram
        self.assertEqual(row["AlleleSizeHistogram"], "5x:1,6x:1,7x:1,8x:1,9x:1,11x:1")
        self.assertEqual(row["NumCalledAlleles"], 6)

    def test_hemi_allele_columns_are_empty_without_a_hemizygous_call(self):
        row = compute_row("X-1000-1020-CAG", "CAG", [7, 9], {"s1": [7, 9]},
                          sample_id_to_sex={"s1": "male"})

        self.assertEqual(row["HemiAlleleMax"], "")
        self.assertEqual(row["HemiAllele99thPercentile"], "")

    def test_missing_allele_sizes(self):
        rows = run_converter(
            "1-44835-44876-AAAT\tAAAT\t8,8\t.,.\t8,.\t.\t9,10\n"
            "1-50000-50020-AT\tAT\t.,.\t.,.\t.\t.,.\t.,.\n"
            "1-60000-60020-AG\tAG\t5,7\t6,6\t.,6\t7\t.,.\n")

        # the locus where every sample is a no-call has no called alleles, so it produces no row
        self.assertEqual([row["LocusId"] for row in rows], ["1-44835-44876-AAAT", "1-60000-60020-AG"])

        self.assertEqual(rows[0]["AlleleSizeHistogram"], "8x:3,9x:1,10x:1")
        self.assertEqual(rows[0]["BiallelicHistogram"], "8/8:1,9/10:1")
        self.assertEqual(rows[0]["NumCalledAlleles"], "5")
        self.assertEqual(rows[0]["ShortAlleleMax"], "9")

        # the "7" call is hemizygous rather than partial, so it still counts as a 7/7 genotype
        self.assertEqual(rows[1]["AlleleSizeHistogram"], "5x:1,6x:3,7x:2")
        self.assertEqual(rows[1]["BiallelicHistogram"], "5/7:1,6/6:1,7/7:1")
        self.assertEqual(rows[1]["NumCalledAlleles"], "6")

    def test_non_integer_allele_size_still_raises(self):
        with self.assertRaises(ValueError):
            run_converter("1-44835-44876-AAAT\tAAAT\t8,8\tNA\n")

    def test_exact_repeat_of_an_lps_row_is_skipped(self):
        rows = run_converter_with_trid_metadata(
            "1-44835-44876-AAAT\tAAAT\t8,8\t9,10\n"
            "1-44835-44876-AAAT\tAAAT\t8,8\t9,10\n",
            TRID_METADATA)

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["Interval"], "1:44835-44876")
        self.assertEqual(rows[0]["AlleleSizeHistogram"], "8x:2,9x:1,10x:1")

    def test_different_lps_rows_for_the_same_record_raise(self):
        with self.assertRaisesRegex(ValueError, "already used by line #1, and the two rows' values differ"):
            run_converter_with_trid_metadata(
                "1-44835-44876-AAAT\tAAAT\t8,8\t9,10\n"
                "1-44835-44876-AAAT\tAAAT\t8,8\t9,11\n",
                TRID_METADATA)

    def test_trid_metadata_record_without_an_lps_row_raises(self):
        with self.assertRaisesRegex(ValueError, r"1 VCF records in --vcf-trid-metadata-tsv were never consumed"):
            run_converter_with_trid_metadata(
                "1-44835-44876-AAAT\tAAAT\t8,8\t9,10\n",
                TRID_METADATA + "1-50000-50020-AT\t1-50000-50020-AT\tAT\t1:50000-50020\t\n")

    def test_a_trid_on_two_vcf_records_is_skipped(self):
        # records 1 and 3 of the AoU (Danzi) VCF share a TRID, with spans 149631602-149631617 and
        # -149631762, and nothing says which LPS row came from which, so neither is used.
        trid = TMEM185A_COORDINATE_TRID
        rows = run_converter_with_trid_metadata(
            f"{trid}\tTCC\t5,5\t5,6\n"
            f"{trid}\tTCC\t3,3\t3,3\n"
            "1-44835-44876-AAAT\tAAAT\t8,8\t9,10\n",
            "trid\tlocus_id\tmotif\tinterval\tvc\n"
            f"{trid}\tX-149631602-149631617-TCC\tTCC\tX:149631602-149631617\tX:149631602-149631617\n"
            f"{trid}\tX-149631602-149631617-TCC\tTCC\tX:149631602-149631762\tX:149631602-149631762\n"
            "1-44835-44876-AAAT\t1-44835-44876-AAAT\tAAAT\t1:44835-44876\t\n")

        # the unaffected locus is still written
        self.assertEqual([row["LocusId"] for row in rows], ["1-44835-44876-AAAT"])

    def test_two_vcf_records_covering_one_repeat_each_get_a_row(self):
        # the AoU (Danzi) VCF covers this repeat twice, under two TRIDs, so the TRID metadata has a
        # row for each. Only the TRID column tells the two output rows apart.
        rows = run_converter_with_trid_metadata(
            f"{TMEM185A_GENE_NAMED_TRID}\tTCC\t5,5\t5,6\n"
            f"{TMEM185A_COORDINATE_TRID}\tTCC\t5,5\t5,7\n",
            TMEM185A_TRID_METADATA)

        self.assertEqual([row["LocusId"] for row in rows],
                         ["X-149631602-149631617-TCC", "X-149631602-149631617-TCC"])
        self.assertEqual([row["TRID"] for row in rows], [TMEM185A_GENE_NAMED_TRID, TMEM185A_COORDINATE_TRID])
        self.assertEqual([row["AlleleSizeHistogram"] for row in rows], ["5x:3,6x:1", "5x:3,7x:1"])

    def test_the_same_locus_id_twice_in_one_vcf_record_raises(self):
        trid = "1-44835-44876-AAAT,1-44835-44876-AAAT"
        with self.assertRaisesRegex(ValueError, "duplicate output tuple"):
            run_converter_with_trid_metadata(
                f"{trid}\tAAAT\t8,8\t9,10\n",
                "trid\tlocus_id\tmotif\tinterval\tvc\n"
                f"{trid}\t1-44835-44876-AAAT\tAAAT\t1:44835-44876\t\n"
                f"{trid}\t1-44835-44876-AAAT\tAAAT\t1:44835-44876\t\n")

    def test_lps_row_missing_from_the_trid_metadata_is_skipped(self):
        # The AoU record's MOTIFS field listed CGC, which no repeat id in its TRID ends with, so
        # trgt-lps wrote a CGC row that the TRID metadata has no LocusId for.
        trid = TMEM185A_COORDINATE_TRID
        interval = "X:149631602-149631762"
        rows = run_converter_with_trid_metadata(
            f"{trid}\tCGC\t2,2\t2,3\n"
            f"{trid}\tCGCCGT\t2,2\t2,4\n"
            f"{trid}\tTCC\t5,5\t5,6\n"
            f"{trid}\tGCT\t3,3\t3,3\n",
            "trid\tlocus_id\tmotif\tinterval\tvc\n"
            f"{trid}\tX-149631723-149631735-CGCCGT\tCGCCGT\t{interval}\t{interval}\n"
            f"{trid}\tX-149631602-149631617-TCC\tTCC\t{interval}\t{interval}\n"
            f"{trid}\tX-149631685-149631694-GCT\tGCT\t{interval}\t{interval}\n")

        self.assertEqual([row["LocusId"] for row in rows], [
            "X-149631723-149631735-CGCCGT", "X-149631602-149631617-TCC", "X-149631685-149631694-GCT"])
        self.assertEqual(rows[0]["AlleleSizeHistogram"], "2x:3,4x:1")

    def test_load_vcf_trid_metadata_reads_bgz_file(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            tsv_path = os.path.join(temp_dir, "trid_metadata.tsv.bgz")
            with gzip.open(tsv_path, "wt") as f:
                f.write("trid\tlocus_id\tmotif\tinterval\tvc\n"
                        "1-44835-44876-AAAT\t1-44835-44876-AAAT\tAAAT\t1:44835-44876\t\n")

            self.assertEqual(
                load_vcf_trid_metadata(tsv_path),
                ({("1-44835-44876-AAAT", "AAAT"): ("1:44835-44876", "", ["1-44835-44876-AAAT"])}, set()))


if __name__ == "__main__":
    unittest.main()
