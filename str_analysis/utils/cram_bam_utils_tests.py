import gzip
import hashlib
import os
import pkgutil
import sys
import tempfile
import unittest
import subprocess
import urllib.request
from pathlib import Path

import pysam

from str_analysis.make_bamlet import extract_region
from str_analysis.utils.cram_bam_utils import IntervalReader
from str_analysis.utils.file_utils import set_requester_pays_project


# The test CRAM/BAM files contain reads aligned to chr9 (the FXN locus). Decoding the CRAM requires a
# chr9 reference whose md5 matches the M5 tag in the CRAM header (6c198acf68b5af7b9d676dfdd531b5de), which
# is the standard GRCh38 chr9 distributed in the Broad public reference bucket. Rather than committing a
# ~36Mb FASTA or relying on htslib's md5-registry fallback (which fails on GitHub Actions), download just
# chr9 from the Broad reference via an HTTP range request on first use and cache it locally.
_BROAD_HG38_FASTA_URL = (
    "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta")
_CHR9_EXPECTED_MD5 = "6c198acf68b5af7b9d676dfdd531b5de"
# chr9 entry from Homo_sapiens_assembly38.fasta.fai: length, byte offset, bases per line, bytes per line.
_CHR9_LENGTH, _CHR9_OFFSET, _CHR9_LINEBASES, _CHR9_LINEWIDTH = 138394717, 1551854835, 100, 101
_TEST_REFERENCE_CACHE_DIR = Path(
    os.environ.get("STR_ANALYSIS_TEST_REF_DIR", "~/.cache/str_analysis_test_refs")).expanduser()


def get_chr9_reference_fasta():
    """Returns the path to a local bgzipped chr9 FASTA for decoding the test CRAM files.

    On first use, downloads the chr9 sequence (only) from the Broad public hg38 reference via an HTTP range
    request, bgzips and indexes it, then caches it under _TEST_REFERENCE_CACHE_DIR. Subsequent calls reuse
    the cached file, so the reference is downloaded at most once per machine (or once per CI cache).

    Returns:
        str: Path to the bgzipped, faidx-indexed chr9 reference FASTA.
    """
    chr9_fasta_gz = _TEST_REFERENCE_CACHE_DIR / "chr9.fa.gz"
    if chr9_fasta_gz.is_file() and (_TEST_REFERENCE_CACHE_DIR / "chr9.fa.gz.fai").is_file():
        return str(chr9_fasta_gz)

    _TEST_REFERENCE_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    sequence_byte_length = ((_CHR9_LENGTH // _CHR9_LINEBASES) * _CHR9_LINEWIDTH
                            + (_CHR9_LENGTH % _CHR9_LINEBASES)
                            + (1 if _CHR9_LENGTH % _CHR9_LINEBASES else 0))
    request = urllib.request.Request(_BROAD_HG38_FASTA_URL, headers={
        "Range": f"bytes={_CHR9_OFFSET}-{_CHR9_OFFSET + sequence_byte_length - 1}"})
    sequence_bytes = urllib.request.urlopen(request, timeout=300).read()
    if hashlib.md5(sequence_bytes.replace(b"\n", b"").upper()).hexdigest() != _CHR9_EXPECTED_MD5:
        raise ValueError(f"Downloaded chr9 sequence md5 doesn't match expected {_CHR9_EXPECTED_MD5}")

    chr9_fasta = _TEST_REFERENCE_CACHE_DIR / "chr9.fa"
    with open(chr9_fasta, "wb") as f:
        f.write(b">chr9\n")
        f.write(sequence_bytes if sequence_bytes.endswith(b"\n") else sequence_bytes + b"\n")
    pysam.tabix_compress(str(chr9_fasta), str(chr9_fasta_gz), force=True)
    pysam.faidx(str(chr9_fasta_gz))
    chr9_fasta.unlink()
    return str(chr9_fasta_gz)


class TestCramBamUtils(unittest.TestCase):

	def setUp(self):
		self._FXN_intervals = [
			("chr9", 69037287, 69037304),
		]
		set_requester_pays_project("cmg-analysis")

		# Write a local copy of a small test CRAM so the interval-merging tests don't require GCS access.
		self._temp_dir = tempfile.TemporaryDirectory()
		self._local_cram_path = os.path.join(self._temp_dir.name, "FXN.wgsim_HET_250xGAA.cram")
		with open(self._local_cram_path, "wb") as f:
			f.write(pkgutil.get_data("str_analysis", "data/tests/FXN.wgsim_HET_250xGAA.cram"))
		with open(self._local_cram_path + ".crai", "wb") as f:
			f.write(pkgutil.get_data("str_analysis", "data/tests/FXN.wgsim_HET_250xGAA.cram.crai"))

	def tearDown(self):
		self._temp_dir.cleanup()

	def test_interval_reader_non_overlapping_intervals(self):
		reader = IntervalReader(self._local_cram_path)

		reader.add_interval("chr1", 1, 5)
		reader.add_interval("chr1", 6, 15)
		reader.add_interval("chr2", 16, 30)
		reader.add_interval("chr2", 1, 15)

		self.assertEqual(reader._get_merged_intervals(), [
			("chr1", 1, 5),
			("chr1", 6, 15),
			("chr2", 1, 15),
			("chr2", 16, 30),
		])

	def test_interval_reader_overlapping_intervals(self):
		reader = IntervalReader(self._local_cram_path)

		reader.add_interval("chr1", 1, 10)
		reader.add_interval("chr1", 5, 15)
		reader.add_interval("chr1", 20, 30)
		reader.add_interval("chr2", 1, 10)
		reader.add_interval("chr2", 5, 15)
		reader.add_interval("chr2", 20, 30)

		self.assertEqual(reader._get_merged_intervals(), [
			("chr1", 1, 15),
			("chr1", 20, 30),
			("chr2", 1, 15),
			("chr2", 20, 30),
		])

	def test_interval_reader_adjacent_intervals(self):
		reader = IntervalReader(self._local_cram_path)

		reader.add_interval("chr1", 1, 5)
		reader.add_interval("chr1", 5, 15)
		reader.add_interval("chr2", 15, 30)
		reader.add_interval("chr2", 1, 15)

		self.assertEqual(reader._get_merged_intervals(), [
			("chr1", 1, 15),
			("chr2", 1, 30),
		])


	def test_cram_reader_on_local_files(self):
		chr9_reference_fasta = get_chr9_reference_fasta()
		with tempfile.NamedTemporaryFile(suffix=".cram") as input_cram_file, \
			  tempfile.NamedTemporaryFile(suffix=".crai") as input_crai_file, \
			  tempfile.NamedTemporaryFile(suffix=".bam") as input_bam_file, \
			  tempfile.NamedTemporaryFile(suffix=".bai") as input_bai_file, \
			  tempfile.NamedTemporaryFile(suffix=".cram") as output_cram_file, \
			  tempfile.NamedTemporaryFile(suffix=".bam") as output_bam_file:

			input_cram_file.write(pkgutil.get_data("str_analysis", "data/tests/FXN.wgsim_HET_250xGAA.cram"))
			input_cram_file.flush()
			input_crai_file.write(pkgutil.get_data("str_analysis", "data/tests/FXN.wgsim_HET_250xGAA.cram.crai"))
			input_crai_file.flush()

			input_bam_file.write(pkgutil.get_data("str_analysis", "data/tests/FXN.wgsim_HET_250xGAA.bam"))
			input_bam_file.flush()
			input_bai_file.write(pkgutil.get_data("str_analysis", "data/tests/FXN.wgsim_HET_250xGAA.bam.bai"))
			input_bai_file.flush()

			cram_reader = IntervalReader(input_cram_file.name, input_crai_file.name, reference_fasta_path=chr9_reference_fasta)
			for interval in self._FXN_intervals:
				cram_reader.add_interval(*interval)
			cram_reads_counter = cram_reader.save_to_file(output_cram_file.name)

			bam_reader = IntervalReader(input_bam_file.name, input_bai_file.name)
			for interval in self._FXN_intervals:
				bam_reader.add_interval(*interval)
			bam_reads_counter = bam_reader.save_to_file(output_bam_file.name)

			#print(f"Retrieved {cram_reads_counter} reads from CRAM and {bam_reads_counter} reads from BAM")
			self.assertEqual(cram_reads_counter, bam_reads_counter)

	def test_cram_reader_skips_interval_with_no_crai_entry(self):
		# Regression test: an interval on a contig that is present in the CRAM header but has no CRAI entries
		# (e.g. a mate landing on a decoy contig) must be skipped, not treated as fatal. Previously
		# _load_cram_containers did `return 0` on the first such interval, so save_to_file wrote no output at
		# all and every other valid interval's reads were silently dropped.
		chr9_reference_fasta = get_chr9_reference_fasta()
		with tempfile.NamedTemporaryFile(suffix=".cram") as baseline_output, \
			  tempfile.NamedTemporaryFile(suffix=".cram") as combined_output:

			baseline_reader = IntervalReader(
				self._local_cram_path, self._local_cram_path + ".crai", reference_fasta_path=chr9_reference_fasta)
			for interval in self._FXN_intervals:
				baseline_reader.add_interval(*interval)
			baseline_read_count = baseline_reader.save_to_file(baseline_output.name)

			combined_reader = IntervalReader(
				self._local_cram_path, self._local_cram_path + ".crai", reference_fasta_path=chr9_reference_fasta)
			for interval in self._FXN_intervals:
				combined_reader.add_interval(*interval)
			# a decoy contig present in the CRAM header but absent from the CRAI index
			combined_reader.add_interval("chrUn_JTFH01000963v1_decoy", 1, 2)
			combined_read_count = combined_reader.save_to_file(combined_output.name)

			# the CRAI-less decoy interval is skipped, not fatal: the output is still written with the same reads
			self.assertTrue(os.path.isfile(combined_output.name))
			self.assertGreater(combined_read_count, 0)
			self.assertEqual(combined_read_count, baseline_read_count)

	def test_cram_reader_does_not_double_count_cached_bytes(self):
		# Regression test: _get_byte_range must count each downloaded byte range exactly once. Previously it
		# incremented the byte counter before the cache check AND again after each read, so cache-miss ranges
		# were counted twice and cache hits (which perform no I/O) were counted too. Re-saving the same
		# intervals from the in-memory cache must therefore not change the reported total_bytes.
		chr9_reference_fasta = get_chr9_reference_fasta()
		with tempfile.NamedTemporaryFile(suffix=".cram") as first_output, \
			  tempfile.NamedTemporaryFile(suffix=".cram") as second_output:

			reader = IntervalReader(
				self._local_cram_path, self._local_cram_path + ".crai",
				reference_fasta_path=chr9_reference_fasta, cache_byte_ranges=True)
			for interval in self._FXN_intervals:
				reader.add_interval(*interval)

			reader.save_to_file(first_output.name)
			bytes_after_first = reader.get_total_bytes_loaded_from_cram()
			ranges_after_first = reader.get_total_byte_ranges_loaded_from_cram()
			self.assertGreater(bytes_after_first, 0)

			# the second save over identical intervals is served entirely from the byte-range cache, so it
			# must not increase either the byte total or the byte-range (container) count
			reader.save_to_file(second_output.name)
			self.assertEqual(reader.get_total_bytes_loaded_from_cram(), bytes_after_first)
			self.assertEqual(reader.get_total_byte_ranges_loaded_from_cram(), ranges_after_first)

	def test_crai_index_intervals_are_0based(self):
		# Regression test: CRAI alignment_start is 1-based, so parse_crai_index must store it as 0-based
		# (alignment_start - 1) to match the 0-based half-open overlap queries. Otherwise a request whose
		# 0-based half-open end equals a container's 1-based start would miss that container by one base.
		crai_bytes = pkgutil.get_data("str_analysis", "data/tests/FXN.wgsim_HET_250xGAA.cram.crai")
		try:
			crai_text = gzip.decompress(crai_bytes).decode()
		except (OSError, gzip.BadGzipFile):
			crai_text = crai_bytes.decode()
		# first record with a non-negative alignment_span: ref_id, alignment_start(1-based), alignment_span, ...
		ref_id, alignment_start_1based, alignment_span = next(
			(r[0], r[1], r[2]) for r in (list(map(int, line.split("\t"))) for line in crai_text.strip().splitlines())
			if r[2] >= 0)

		reader = IntervalReader(self._local_cram_path, self._local_cram_path + ".crai")
		crai_interval_tree = reader._crai_interval_trees[ref_id]

		# a 0-based half-open request ending exactly at the container's 1-based alignment_start must overlap it
		self.assertTrue(crai_interval_tree.overlap(alignment_start_1based - 2, alignment_start_1based))
		# and the stored 0-based start is present in the tree
		self.assertIn(alignment_start_1based - 1, {interval.begin for interval in crai_interval_tree})

	def test_save_to_file_deduplicates_reads_across_intervals(self):
		# Regression test: a read overlapping two non-overlapping requested intervals is returned once per
		# interval by fetch(), and must be written to the output only once. Duplicate copies would otherwise
		# be misread downstream as a complete read pair and suppress real mate retrieval.
		chr9_reference_fasta = get_chr9_reference_fasta()
		with tempfile.NamedTemporaryFile(suffix=".cram") as output_cram_file:
			reader = IntervalReader(
				self._local_cram_path, self._local_cram_path + ".crai", reference_fasta_path=chr9_reference_fasta)
			# two non-overlapping 1bp intervals ~33bp apart; a single ~150bp read spans both (26 such reads exist)
			reader.add_interval("chr9", 69037287, 69037288)
			reader.add_interval("chr9", 69037320, 69037321)
			written_read_count = reader.save_to_file(output_cram_file.name)

			with pysam.AlignmentFile(output_cram_file.name, reference_filename=chr9_reference_fasta) as f:
				read_keys = [(read.query_name, read.flag, read.reference_start) for read in f]

			self.assertGreater(written_read_count, 0)
			self.assertEqual(len(read_keys), len(set(read_keys)))  # no duplicate alignments in the output
			self.assertEqual(written_read_count, len(read_keys))   # returned count matches unique reads written

	def test_save_to_file_emits_coordinate_order_without_sorting(self):
		# Guard the assumption that lets callers index save_to_file's stream directly without sorting: source
		# CRAM is coordinate-sorted, intervals are traversed in coordinate order, and a read spanning two
		# intervals is emitted by the earlier interval before de-duplication suppresses its later copy.
		with tempfile.TemporaryDirectory() as temp_dir:
			input_cram_path = os.path.join(temp_dir, "coordinate_sorted.cram")
			with pysam.AlignmentFile(
					input_cram_path, "wc",
					header={
						"HD": {"VN": "1.6", "SO": "coordinate"},
						"SQ": [
							{"SN": "chr1", "LN": 2000},
							{"SN": "chr2", "LN": 2000},
							{"SN": "chr10", "LN": 2000},
						],
					},
					format_options=[b"no_ref=1"]) as input_cram:
				for query_name, reference_id, reference_start, flag, mate_reference_id, mate_start, read_length in (
						# read2 precedes read1, so read1's mate lies at an earlier coordinate
						("split_pair", 0, 100, 1 | 2 | 128, 0, 900, 50),
						# spans two non-adjacent requested intervals and must be written only by the earlier one
						("boundary", 0, 195, 0, -1, -1, 120),
						("split_pair", 0, 900, 1 | 2 | 64, 0, 100, 50),
						("chr2_read", 1, 400, 0, -1, -1, 50),
						# chr10 follows chr2 by reference id despite preceding it lexicographically
						("chr10_read", 2, 300, 0, -1, -1, 50),
						("unmapped", -1, -1, 4, -1, -1, 50),
				):
					read = pysam.AlignedSegment()
					read.query_name = query_name
					read.flag = flag
					read.reference_id = reference_id
					read.reference_start = reference_start
					read.mapping_quality = 60 if reference_id >= 0 else 0
					if reference_id >= 0:
						read.cigarstring = f"{read_length}M"
					read.next_reference_id = mate_reference_id
					read.next_reference_start = mate_start
					read.query_sequence = "A" * read_length
					read.query_qualities = pysam.qualitystring_to_array("I" * read_length)
					input_cram.write(read)
			pysam.index(input_cram_path)

			reader = IntervalReader(
				input_cram_path, input_cram_path + ".crai", include_unmapped_read_pairs=True)
			# Add at least three non-adjacent intervals in deliberately non-coordinate order. The 200 and 300
			# intervals both return "boundary"; the 900 interval contains a read whose mate is at position 100.
			for interval in (
					("chr10", 300, 310),
					("chr1", 900, 910),
					("chr2", 400, 410),
					("chr1", 300, 301),
					("chr1", 200, 201),
					("chr1", 100, 110),
			):
				reader.add_interval(*interval)

			output_cram_path = os.path.join(temp_dir, "output.cram")
			try:
				# create_index=True exercises the changed branch that indexes the output directly with
				# pysam.index instead of re-sorting it first. pysam.index accepts an out-of-order CRAM without
				# error, so the coordinate order itself is verified by the assertions below, not by this call.
				self.assertEqual(reader.save_to_file(output_cram_path, create_index=True), 6)
			finally:
				reader.close()
			with pysam.AlignmentFile(output_cram_path) as output_cram:
				output_reads = list(output_cram)
				# SAM coordinate order places unmapped reads after all reference ids.
				coordinate_keys = [
					(read.reference_id if read.reference_id >= 0 else len(output_cram.references),
					 read.reference_start if read.reference_start >= 0 else sys.maxsize)
					for read in output_reads
				]

			self.assertEqual(coordinate_keys, sorted(coordinate_keys))
			self.assertEqual(sum(read.query_name == "boundary" for read in output_reads), 1)
			self.assertEqual(
				[(read.reference_start, read.next_reference_start)
				 for read in output_reads if read.query_name == "split_pair"],
				[(100, 900), (900, 100)])
			self.assertTrue(output_reads[-1].is_unmapped)

	def test_cram_reader_handles_cram_v21_input(self):
		# Regression test: the CRAM EOF marker length is version-dependent (38 bytes for CRAM 3, 30 for CRAM 2.x).
		# Sizing the last container with a hard-coded 38-byte EOF truncates a CRAM 2.1 file's last container by
		# 8 bytes, corrupting the reconstructed subset. The FXN interval lands in the last container here.
		chr9_reference_fasta = get_chr9_reference_fasta()
		with tempfile.TemporaryDirectory() as temp_dir:
			input_bam_path = os.path.join(temp_dir, "FXN.bam")
			with open(input_bam_path, "wb") as f:
				f.write(pkgutil.get_data("str_analysis", "data/tests/FXN.wgsim_HET_250xGAA.bam"))
			cram_v21_path = os.path.join(temp_dir, "FXN.v2.1.cram")
			pysam.view("-O", "cram,version=2.1,no_ref=1", "-o", cram_v21_path, input_bam_path, catch_stdout=False)
			pysam.index(cram_v21_path)

			output_cram_path = os.path.join(temp_dir, "out.cram")
			reader = IntervalReader(cram_v21_path, cram_v21_path + ".crai", reference_fasta_path=chr9_reference_fasta)
			for interval in self._FXN_intervals:
				reader.add_interval(*interval)
			written_read_count = reader.save_to_file(output_cram_path)

			self.assertGreater(written_read_count, 0)
			self.assertTrue(os.path.isfile(output_cram_path))

	def test_cli_exits_cleanly_when_no_reads_in_region(self):
		# Regression test (covers the negative-window clamp and the empty-first-pass guard): a region with no
		# overlapping reads must produce a clear error + exit code 1, not a pysam ValueError. chr9:1-2 both
		# exercises window_start = max(0, 1 - 1000) (no negative-coordinate fetch) and the no-reads exit path.
		chr9_reference_fasta = get_chr9_reference_fasta()
		with tempfile.TemporaryDirectory() as temp_dir:
			output_cram_path = os.path.join(temp_dir, "out.cram")
			result = subprocess.run(
				[sys.executable, "-m", "str_analysis.make_minicram_for_expansion_hunter",
				 "-R", chr9_reference_fasta, "-L", "chr9:1-2", "-o", output_cram_path,
				 "-i", self._local_cram_path + ".crai", self._local_cram_path],
				capture_output=True, text=True)

			self.assertEqual(result.returncode, 1)
			self.assertFalse(os.path.isfile(output_cram_path))
			self.assertIn("No reads were found", result.stdout + result.stderr)

	def test_save_to_file_resolves_chromosome_naming_convention(self):
		# Regression test: a requested region whose chromosome naming convention differs from the CRAM header's
		# (e.g. "9" vs a header listing "chr9") must resolve to the header's actual reference name for fetch(),
		# rather than being passed through raw and crashing with "invalid contig".
		chr9_reference_fasta = get_chr9_reference_fasta()
		_, fxn_start, fxn_end = self._FXN_intervals[0]
		with tempfile.NamedTemporaryFile(suffix=".cram") as chr_prefixed_output, \
			  tempfile.NamedTemporaryFile(suffix=".cram") as unprefixed_output:
			chr_prefixed_reader = IntervalReader(
				self._local_cram_path, self._local_cram_path + ".crai", reference_fasta_path=chr9_reference_fasta)
			chr_prefixed_reader.add_interval("chr9", fxn_start, fxn_end)
			chr_prefixed_read_count = chr_prefixed_reader.save_to_file(chr_prefixed_output.name)

			unprefixed_reader = IntervalReader(
				self._local_cram_path, self._local_cram_path + ".crai", reference_fasta_path=chr9_reference_fasta)
			unprefixed_reader.add_interval("9", fxn_start, fxn_end)  # header uses "chr9"
			unprefixed_read_count = unprefixed_reader.save_to_file(unprefixed_output.name)

			self.assertGreater(unprefixed_read_count, 0)
			self.assertEqual(unprefixed_read_count, chr_prefixed_read_count)

	def test_load_cram_containers_skips_contig_absent_from_header(self):
		# Regression test: an interval on a contig that is not present in the CRAM header must be skipped with a
		# warning (like the CRAI-less case), not raise an uncaught KeyError that aborts the whole export.
		chr9_reference_fasta = get_chr9_reference_fasta()
		with tempfile.NamedTemporaryFile(suffix=".cram") as output_cram_file:
			reader = IntervalReader(
				self._local_cram_path, self._local_cram_path + ".crai", reference_fasta_path=chr9_reference_fasta)
			reader.add_interval(*self._FXN_intervals[0])          # valid chr9 interval
			reader.add_interval("chrNonexistentContig", 1, 100)   # contig absent from the CRAM header
			written_read_count = reader.save_to_file(output_cram_file.name)  # must not raise KeyError

			self.assertGreater(written_read_count, 0)
			self.assertTrue(os.path.isfile(output_cram_file.name))

	def test_cram_reader_handles_empty_cram(self):
		# Regression test: a CRAM with zero mapped reads has an empty CRAI (no data-container records), which
		# leaves the end-of-header offset unset; the reader must still construct and return a clean zero-read
		# result instead of crashing with a TypeError in the constructor.
		chr9_reference_fasta = get_chr9_reference_fasta()
		with tempfile.TemporaryDirectory() as temp_dir:
			with pysam.AlignmentFile(self._local_cram_path, check_sq=False) as source:
				header = source.header.to_dict()
			empty_cram_path = os.path.join(temp_dir, "empty.cram")
			with pysam.AlignmentFile(empty_cram_path, "wc", header=header, format_options=[b"no_ref=1"]):
				pass  # write a header-only CRAM with no reads
			pysam.index(empty_cram_path)

			reader = IntervalReader(  # must not raise despite the empty CRAI
				empty_cram_path, empty_cram_path + ".crai", reference_fasta_path=chr9_reference_fasta)
			for interval in self._FXN_intervals:
				reader.add_interval(*interval)
			output_cram_path = os.path.join(temp_dir, "out.cram")
			self.assertEqual(reader.save_to_file(output_cram_path), 0)

	def test_extract_region_skips_unpaired_reads(self):
		# Regression test: an unpaired (single-end / long-read) alignment has no mate and its next_reference_name
		# is None; extract_region must not emit a (None, ...) mate region, which would crash add_interval.
		with tempfile.TemporaryDirectory() as temp_dir:
			unpaired_bam_path = os.path.join(temp_dir, "unpaired.bam")
			header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
			with pysam.AlignmentFile(unpaired_bam_path, "wb", header=header) as out:
				read = pysam.AlignedSegment()
				read.query_name = "unpaired_read"
				read.flag = 0  # unpaired
				read.reference_id = 0
				read.reference_start = 100
				read.mapping_quality = 60
				read.cigartuples = [(0, 50)]
				read.query_sequence = "A" * 50
				read.query_qualities = pysam.qualitystring_to_array("I" * 50)
				out.write(read)
			pysam.index(unpaired_bam_path)

			with pysam.AlignmentFile(unpaired_bam_path) as input_bam:
				genomic_regions = extract_region("chr1", 90, 160, input_bam=input_bam, bamlet=None)

			self.assertFalse(any(region[0] is None for region in genomic_regions))

	def test_cram_reader_on_google_storage_files(self):
		try:
			os.environ["GCS_OAUTH_TOKEN"] = subprocess.check_output("gcloud auth application-default print-access-token", shell=True).decode("utf-8").strip()
		except Exception as e:
			print(f"WARNING: Unable to set GCS_OAUTH_TOKEN: {e}. Skipping gs:// tests..")
			return

		#with tempfile.NamedTemporaryFile(suffix=".cram") as output_cram_file, \
		#	tempfile.NamedTemporaryFile(suffix=".bam") as output_bam_file:
		#
		#	# check cloud paths
		#	cram_reader = IntervalReader(
		#		"gs://str-analysis/tests/FXN.wgsim_HET_250xGAA.cram",
		#		"gs://str-analysis/tests/FXN.wgsim_HET_250xGAA.cram.crai",
		#		verbose=True)
		#	for interval in self._FXN_intervals:
		#		cram_reader.add_interval(*interval)
#
		#	cram_reads_counter = cram_reader.save_to_file(output_cram_file.name)
		#	#print(f"Retrieved {cram_reads_counter} reads from CRAM using default mode")
#
		#	cram_reader2 = IntervalReader(
		#		"gs://str-analysis/tests/FXN.wgsim_HET_250xGAA.cram",
		#		"gs://str-analysis/tests/FXN.wgsim_HET_250xGAA.cram.crai")
		#	for interval in self._FXN_intervals:
		#		cram_reader2.add_interval(*interval)
		#	cram_reads_counter2 = cram_reader2.save_to_file(output_bam_file.name)
		#	#print(f"Retrieved {cram_reads_counter2} reads from CRAM")
		#	self.assertEqual(cram_reads_counter, cram_reads_counter2)
#
		#	bam_reader = IntervalReader(
		#		"gs://str-analysis/tests/FXN.wgsim_HET_250xGAA.bam",
		#		"gs://str-analysis/tests/FXN.wgsim_HET_250xGAA.bam.bai")
		#	for interval in self._FXN_intervals:
		#		bam_reader.add_interval(*interval)
		#	bam_reads_counter = bam_reader.save_to_file(output_bam_file.name)
#
		#	#print(f"Retrieved {bam_reads_counter} reads from BAM")
		#	self.assertEqual(cram_reads_counter, bam_reads_counter)
