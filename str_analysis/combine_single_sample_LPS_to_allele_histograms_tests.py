import collections
import gzip
import os
import tempfile
import unittest
from unittest import mock

from str_analysis.combine_single_sample_LPS_to_allele_histograms import (
    get_lps_filename_prefix, main, parse_lps_table)


def parse_lps_table_contents(contents, n_outlier_sample_ids=10):
    """Run parse_lps_table on an LPS table and return the histograms it produced.

    Args:
        contents (str): contents of a single-sample LPS table, including its header line.
        n_outlier_sample_ids (int): passed through to parse_lps_table.

    Returns:
        tuple: (all_allele_histogram, short_allele_histogram, hemizygous_allele_histogram), each a
        dict mapping (trid, motif) to {allele_size: count}, and each holding only the loci that
        got at least one allele.
    """
    histograms = [collections.defaultdict(lambda: collections.defaultdict(int)) for _ in range(3)]
    sample_ids = [collections.defaultdict(lambda: collections.defaultdict(list)) for _ in range(3)]

    with tempfile.TemporaryDirectory() as temp_dir:
        input_path = os.path.join(temp_dir, "sample1.lps.tsv")
        with open(input_path, "wt") as f:
            f.write(contents)

        parse_lps_table(
            input_path,
            None,
            n_outlier_sample_ids,
            histograms[0], sample_ids[0],
            histograms[1], sample_ids[1],
            histograms[2], sample_ids[2],
            use_sample_id_from_header=True,
        )

    return tuple({key: dict(counts) for key, counts in histogram.items()} for histogram in histograms)


def run_main(contents):
    """Run the script end to end on a single LPS table and return its output rows.

    Args:
        contents (str): contents of a single-sample LPS table, including its header line.

    Returns:
        list: one dict per output row, mapping column name to the unparsed string value.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        input_path = os.path.join(temp_dir, "sample1.lps.tsv")
        with open(input_path, "wt") as f:
            f.write(contents)

        output_prefix = os.path.join(temp_dir, "out")
        argv = ["combine_single_sample_LPS_to_allele_histograms",
                "--use-sample-id-from-header", "-o", output_prefix, input_path]
        with mock.patch("sys.argv", argv):
            main()

        with gzip.open(f"{output_prefix}.tsv.gz", "rt") as f:
            header = next(f).rstrip("\n").split("\t")
            return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in f]


class Tests(unittest.TestCase):

    def test_diploid_and_hemizygous_calls(self):
        all_alleles, short_alleles, hemizygous_alleles = parse_lps_table_contents(
            "trid\tmotif\tsample1\n"
            "X-263540-263579-TTTA\tTTTA\t10,11\n"
            "X-264744-264774-AAAG\tAAAG\t5,5\n"
            "X-95117-95209-GTCA\tGTCA\t3\n")

        self.assertEqual(all_alleles[("X-263540-263579-TTTA", "TTTA")], {10: 1, 11: 1})
        self.assertEqual(short_alleles[("X-263540-263579-TTTA", "TTTA")], {10: 1})
        self.assertEqual(all_alleles[("X-264744-264774-AAAG", "AAAG")], {5: 2})
        self.assertEqual(short_alleles[("X-264744-264774-AAAG", "AAAG")], {5: 1})

        # the hemizygous call counts once in every histogram
        self.assertEqual(all_alleles[("X-95117-95209-GTCA", "GTCA")], {3: 1})
        self.assertEqual(short_alleles[("X-95117-95209-GTCA", "GTCA")], {3: 1})
        self.assertEqual(hemizygous_alleles[("X-95117-95209-GTCA", "GTCA")], {3: 1})

        self.assertEqual(len(hemizygous_alleles), 1)

    def test_no_calls_are_skipped(self):
        all_alleles, short_alleles, hemizygous_alleles = parse_lps_table_contents(
            "trid\tmotif\tsample1\n"
            "1-50000-50020-AT\tAT\t.\n"
            "1-60000-60020-AG\tAG\t.,.\n")

        self.assertEqual(all_alleles, {})
        self.assertEqual(short_alleles, {})
        self.assertEqual(hemizygous_alleles, {})

    def test_partial_call_counts_only_in_the_all_allele_histogram(self):
        all_alleles, short_alleles, hemizygous_alleles = parse_lps_table_contents(
            "trid\tmotif\tsample1\n"
            "1-44835-44876-AAAT\tAAAT\t8,.\n"
            "1-70000-70020-AC\tAC\t.,9\n")

        # the measured allele counts once, and isn't doubled the way a homozygous call is
        self.assertEqual(all_alleles[("1-44835-44876-AAAT", "AAAT")], {8: 1})
        self.assertEqual(all_alleles[("1-70000-70020-AC", "AC")], {9: 1})

        # which of the sample's two alleles is shorter is unknown, and it isn't a hemizygous call
        self.assertEqual(short_alleles, {})
        self.assertEqual(hemizygous_alleles, {})

    def test_locus_with_only_partial_calls_still_gets_an_output_row(self):
        rows = run_main(
            "trid\tmotif\tsample1\n"
            "1-44835-44876-AAAT\tAAAT\t8,.\n"
            "1-60000-60020-AG\tAG\t5,7\n")

        self.assertEqual([row["LocusId"] for row in rows], ["1-44835-44876-AAAT", "1-60000-60020-AG"])
        self.assertEqual(rows[0]["AllAlleleHistogram"], "8x:1")
        self.assertEqual(rows[0]["ShortAlleleHistogram"], "")
        self.assertEqual(rows[0]["HemizygousAlleleHistogram"], "")
        self.assertEqual(rows[1]["AllAlleleHistogram"], "5x:1,7x:1")
        self.assertEqual(rows[1]["ShortAlleleHistogram"], "5x:1")

    def test_non_integer_allele_size_raises(self):
        with self.assertRaises(ValueError):
            parse_lps_table_contents("trid\tmotif\tsample1\n1-44835-44876-AAAT\tAAAT\tNA\n")

    def test_get_lps_filename_prefix_strips_gz_and_bgz(self):
        for filename in ["HG00096.lps.tsv", "HG00096.lps.tsv.gz", "HG00096.lps.tsv.bgz",
                         "HG00096.repeat_counts.txt.bgz"]:
            self.assertEqual(get_lps_filename_prefix(f"/data/{filename}"), "HG00096", filename)

    def test_get_lps_filename_prefix_only_strips_suffixes_that_start_with_a_dot(self):
        self.assertEqual(get_lps_filename_prefix("/data/sample_lps.tsv"), "sample_lps")
        self.assertEqual(get_lps_filename_prefix("/data/sample_tsv"), "sample_tsv")


if __name__ == "__main__":
    unittest.main()
