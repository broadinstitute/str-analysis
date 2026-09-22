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

    def test_load_vcf_trid_metadata_reads_bgz_file(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            tsv_path = os.path.join(temp_dir, "trid_metadata.tsv.bgz")
            with gzip.open(tsv_path, "wt") as f:
                f.write("trid\tlocus_id\tmotif\tinterval\tvc\n"
                        "1-44835-44876-AAAT\t1-44835-44876-AAAT\tAAAT\t1:44835-44876\t\n")

            self.assertEqual(
                load_vcf_trid_metadata(tsv_path),
                {("1-44835-44876-AAAT", "AAAT"): ("1:44835-44876", "", ["1-44835-44876-AAAT"])})


if __name__ == "__main__":
    unittest.main()
