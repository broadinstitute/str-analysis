"""Tests for make_minicram_for_expansion_hunter.py"""

import json
import os
import random
import sys
import tempfile
import unittest
from unittest import mock

import pysam

from str_analysis import make_minicram_for_expansion_hunter
from str_analysis.make_minicram_for_expansion_hunter import merge_nearby_regions

CONTIG_LENGTH = 20_000
READ_LENGTH = 100


def make_read(query_name, flag, reference_start, mate_reference_start, cigar=f"{READ_LENGTH}M"):
    """Returns a read on chr1 whose mate is also on chr1."""
    read = pysam.AlignedSegment()
    read.query_name = query_name
    read.flag = flag
    read.reference_id = 0
    read.reference_start = reference_start
    read.mapping_quality = 0 if flag & 4 else 60
    if not flag & 4:
        read.cigarstring = cigar
    read.next_reference_id = 0
    read.next_reference_start = mate_reference_start
    read.query_sequence = "A" * READ_LENGTH
    read.query_qualities = pysam.qualitystring_to_array("I" * READ_LENGTH)
    return read


def make_read_pair(query_name, read1_start, read2_start):
    """Returns two mapped, paired reads."""
    return [
        make_read(query_name, 1 | 2 | 32 | 64, read1_start, read2_start),
        make_read(query_name, 1 | 2 | 16 | 128, read2_start, read1_start),
    ]


class MergeNearbyRegionsTests(unittest.TestCase):

    def test_merges_regions_at_most_max_distance_apart(self):
        # gaps of 198, 1000 and 1001 bp on chr1
        regions = [("chr1", 2305, 2307), ("chr1", 100, 102), ("chr1", 300, 302), ("chr1", 1302, 1304),
                   ("chr2", 200, 202)]
        self.assertEqual(merge_nearby_regions(regions, max_distance=1000), [
            ("chr1", 100, 1304), ("chr1", 2305, 2307), ("chr2", 200, 202)])

    def test_merges_across_chromosome_naming_conventions_but_not_across_chromosomes(self):
        self.assertEqual(merge_nearby_regions([("chr1", 100, 102), ("1", 150, 152)], max_distance=1000),
                         [("chr1", 100, 152)])
        self.assertEqual(merge_nearby_regions([("chr1", 100, 102), ("chr2", 150, 152)], max_distance=1000),
                         [("chr1", 100, 102), ("chr2", 150, 152)])

    def test_keeps_a_contained_region_inside_the_merged_one(self):
        self.assertEqual(merge_nearby_regions([("chr1", 100, 900), ("chr1", 200, 300)], max_distance=0),
                         [("chr1", 100, 900)])


class MakeMinicramTests(unittest.TestCase):
    """Runs the tool on a small synthetic CRAM and BAM.

    The catalog has one locus at chr1:3000-3030 with an off-target region at chr1:6000-6030, so with the default
    1kb window the reads ExpansionHunter can use are those overlapping chr1:2000-4030 or chr1:5000-7030, plus their
    mates. chr1:15000-15200 is a deep pileup: three of its reads are mates of reads in the windows, the rest aren't.
    """

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.reference_fasta_path = os.path.join(self.temp_dir.name, "reference.fa")
        random_generator = random.Random(1)
        with open(self.reference_fasta_path, "w") as f:
            f.write(">chr1\n" + "".join(random_generator.choice("ACGT") for _ in range(CONTIG_LENGTH)) + "\n")
        pysam.faidx(self.reference_fasta_path)

        self.catalog_path = os.path.join(self.temp_dir.name, "catalog.json")
        with open(self.catalog_path, "w") as f:
            json.dump([{
                "LocusId": "TEST_LOCUS",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr1:3000-3030",
                "VariantType": "RareRepeat",
                "OfftargetRegions": ["chr1:6000-6030"],
            }], f)

        reads = []
        # a pair entirely inside the locus window
        reads += make_read_pair("pair_inside_locus_window", 2500, 2800)
        # a read in the locus window and one in the off-target window whose mates are in the pileup
        reads += make_read_pair("pair_from_locus_window_to_pileup", 2600, 15050)
        reads += make_read_pair("pair_from_offtarget_window_to_pileup", 5500, 15100)
        # a read that overlaps the locus window only with its end (it starts before the window), whose mate is in the
        # pileup
        reads += make_read_pair("pair_from_locus_window_edge_to_pileup", 1950, 15150)
        # a read that ends before the locus window, whose mate is in the pileup: ExpansionHunter needs neither, and the
        # read starts 300bp before the window, more than the filter's widening of 1.5 times the read length
        reads += make_read_pair("pair_from_before_locus_window_to_pileup", 1700, 15120)
        # the rest of the pileup: pairs that have nothing to do with the windows
        for i in range(50):
            reads += make_read_pair(f"unrelated_pair_in_pileup_{i}", 15000 + i, 15100 + i)
        reads.sort(key=lambda read: read.reference_start)

        self.expected_read_names = {
            "pair_inside_locus_window", "pair_from_locus_window_to_pileup", "pair_from_offtarget_window_to_pileup",
            "pair_from_locus_window_edge_to_pileup"}

        header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": CONTIG_LENGTH}]}
        self.input_paths = {}
        for file_type, mode, format_options in ("cram", "wc", [b"no_ref=1"]), ("bam", "wb", None):
            input_path = os.path.join(self.temp_dir.name, f"input.{file_type}")
            with pysam.AlignmentFile(input_path, mode, header=header, reference_filename=self.reference_fasta_path,
                                     format_options=format_options) as f:
                for read in reads:
                    f.write(read)
            pysam.index(input_path)
            self.input_paths[file_type] = input_path

    def tearDown(self):
        self.temp_dir.cleanup()

    def run_make_minicram(self, input_path):
        output_cram_path = os.path.join(self.temp_dir.name, f"{os.path.basename(input_path)}.minicram.cram")
        with mock.patch.object(sys, "argv", [
                "make_minicram_for_expansion_hunter", "-R", self.reference_fasta_path, "-c", self.catalog_path,
                "-o", output_cram_path, input_path]):
            make_minicram_for_expansion_hunter.main()
        with pysam.AlignmentFile(output_cram_path, reference_filename=self.reference_fasta_path) as f:
            return [(read.query_name, read.flag) for read in f]

    def test_keeps_window_reads_and_their_mates_and_drops_unrelated_reads(self):
        for file_type, input_path in self.input_paths.items():
            with self.subTest(file_type=file_type):
                output_reads = self.run_make_minicram(input_path)
                self.assertEqual({query_name for query_name, _ in output_reads}, self.expected_read_names)
                # both reads of every pair, each written once
                self.assertEqual(len(output_reads), 2 * len(self.expected_read_names))
                self.assertEqual(len(set(output_reads)), len(output_reads))


if __name__ == "__main__":
    unittest.main()
