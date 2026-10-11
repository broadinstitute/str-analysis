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

        return combine_tables([input_path], os.path.join(temp_dir, "out"))


def combine_tables(input_paths, output_prefix):
    """Run the script end to end on the given input tables and return its output rows.

    Args:
        input_paths (list): paths of single-sample LPS tables or of combined histogram tables written by this script.
        output_prefix (str): output path prefix, so the output table is written to {output_prefix}.tsv.gz.

    Returns:
        list: one dict per output row, mapping column name to the unparsed string value.
    """
    argv = ["combine_single_sample_LPS_to_allele_histograms",
            "--use-sample-id-from-header", "-o", output_prefix] + input_paths
    with mock.patch("sys.argv", argv):
        main()

    with gzip.open(f"{output_prefix}.tsv.gz", "rt") as f:
        header = next(f).rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in f]


def combine_single_locus_samples_into_shard(temp_dir, shard_name, sample_id_to_lps_value):
    """Write a single-sample LPS table for each sample at one locus, and combine them into a shard table.

    Args:
        temp_dir (str): directory to write the tables to.
        shard_name (str): file name prefix for this shard's tables.
        sample_id_to_lps_value (dict): maps sample id to its LPS value at the locus, e.g. "3,5".

    Returns:
        str: path of the shard's combined histogram table.
    """
    input_paths = []
    for sample_id, lps_value in sample_id_to_lps_value.items():
        input_path = os.path.join(temp_dir, f"{shard_name}.{sample_id}.lps.tsv")
        with open(input_path, "wt") as f:
            f.write(f"trid\tmotif\t{sample_id}\n1-1000-1030-CAG\tCAG\t{lps_value}\n")
        input_paths.append(input_path)

    output_prefix = os.path.join(temp_dir, shard_name)
    combine_tables(input_paths, output_prefix)
    return f"{output_prefix}.tsv.gz"


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

    def test_small_shard_does_not_turn_a_common_allele_size_into_outliers(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            # 3x is common in the large shard, so that shard lists nothing at 3x. The small shard has too few
            # samples for 3x to look common, so it lists all of its samples there.
            large_shard_samples = {f"a{i:02d}": "3,3" for i in range(12)}
            large_shard_samples["a12"] = "3,5"
            small_shard_samples = {f"b{i}": "3,3" for i in range(1, 7)}
            small_shard_samples["b0"] = "3,6"

            large_shard_path = combine_single_locus_samples_into_shard(temp_dir, "large_shard", large_shard_samples)
            small_shard_path = combine_single_locus_samples_into_shard(temp_dir, "small_shard", small_shard_samples)
            rows = combine_tables([large_shard_path, small_shard_path], os.path.join(temp_dir, "combined"))

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["AllAlleleHistogram"], "3x:38,5x:1,6x:1")
        # the two rare alleles, one from each shard, are still recorded, but none of the 3x carriers
        self.assertEqual(rows[0]["OutlierSampleIds_AllAlleles"], "6x:b0,5x:a12")
        self.assertEqual(rows[0]["OutlierSampleIds_ShortAlleles"], "")

    def test_common_allele_size_that_no_shard_lists_still_ends_the_outliers(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            # 10x is common in the large shard, so that shard lists nothing at or below 10x. The small shard has no
            # 10x alleles, so no shard lists 10x, and the small shard's shorter 4x alleles are not outliers.
            large_shard_samples = {f"a{i:02d}": "10,10" for i in range(12)}
            large_shard_samples["a12"] = "4,10"
            small_shard_samples = {f"b{i}": "4,4" for i in range(6)}

            large_shard_path = combine_single_locus_samples_into_shard(temp_dir, "large_shard", large_shard_samples)
            small_shard_path = combine_single_locus_samples_into_shard(temp_dir, "small_shard", small_shard_samples)
            rows = combine_tables([large_shard_path, small_shard_path], os.path.join(temp_dir, "combined"))

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["AllAlleleHistogram"], "4x:13,10x:25")
        self.assertEqual(rows[0]["OutlierSampleIds_AllAlleles"], "")
        self.assertEqual(rows[0]["OutlierSampleIds_ShortAlleles"], "")

    def test_rare_allele_size_carried_by_homozygous_samples_is_recorded(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            # 6 homozygous samples give 12 alleles at 9x. That is more than 10 alleles but fewer than 10 samples,
            # so all 6 are still outliers.
            sample_id_to_lps_value = {f"s{i:02d}": "9,9" for i in range(6)}
            sample_id_to_lps_value.update({f"s{i:02d}": "3,3" for i in range(6, 26)})
            shard_path = combine_single_locus_samples_into_shard(temp_dir, "shard", sample_id_to_lps_value)
            with gzip.open(shard_path, "rt") as f:
                header = next(f).rstrip("\n").split("\t")
                row = dict(zip(header, next(f).rstrip("\n").split("\t")))

        self.assertEqual(row["AllAlleleHistogram"], "3x:40,9x:12")
        self.assertEqual(row["OutlierSampleIds_AllAlleles"], ",".join(f"9x:s{i:02d}" for i in range(6)))
        self.assertEqual(row["OutlierSampleIds_ShortAlleles"], ",".join(f"9x:s{i:02d}" for i in range(6)))

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
