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
from unittest import mock

import pysam

from str_analysis.make_bamlet import extract_region
from str_analysis.utils import cram_bam_utils
from str_analysis.utils.cram_bam_utils import IntervalReader, merge_adjacent_byte_ranges, get_total_mapped_reads, parse_crai_index
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

	def test_interval_reader_normalizes_and_validates_intervals(self):
		reader = IntervalReader(self._local_cram_path)
		try:
			reader.add_interval("chr9", -1, 1)
			self.assertEqual(reader._get_merged_intervals(), [("chr9", 0, 1)])
			with self.assertRaisesRegex(ValueError, "chromosome must be a non-empty string"):
				reader.add_interval(None, 0, 1)
			with self.assertRaisesRegex(ValueError, "end must be greater than start"):
				reader.add_interval("chr9", -1, 0)
		finally:
			reader.close()


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

	def test_save_to_file_reference_compression_option(self):
		chr9_reference_fasta = get_chr9_reference_fasta()
		with tempfile.NamedTemporaryFile(suffix=".cram") as reference_compressed_output, \
			  tempfile.NamedTemporaryFile(suffix=".cram") as self_contained_output:
			reader = IntervalReader(
				self._local_cram_path, self._local_cram_path + ".crai", reference_fasta_path=chr9_reference_fasta)
			for interval in self._FXN_intervals:
				reader.add_interval(*interval)

			observed_format_options = {}
			real_alignment_file = pysam.AlignmentFile

			def record_output_format_options(*args, **kwargs):
				"""Records output format options while delegating to pysam.AlignmentFile.

				save_to_file writes to a sibling "<output>.partial.cram" and only renames it onto the requested
				path once indexing has succeeded, so match on that prefix and record under the requested name.
				"""
				# args[0] is not always a path -- the CRAM branch opens the temp container file object itself
				for output_name in (reference_compressed_output.name, self_contained_output.name):
					if isinstance(args[0], str) and (
							args[0] == output_name or args[0].startswith(f"{output_name}.partial.")):
						observed_format_options[output_name] = kwargs.get("format_options")
				return real_alignment_file(*args, **kwargs)

			try:
				with mock.patch(
						"str_analysis.utils.cram_bam_utils.pysam.AlignmentFile",
						side_effect=record_output_format_options):
					reader.save_to_file(
						reference_compressed_output.name,
						create_index=False,
						disable_reference_compression=False)
					reader.save_to_file(
						self_contained_output.name,
						create_index=False,
						disable_reference_compression=True)
			finally:
				reader.close()

			self.assertIsNone(observed_format_options[reference_compressed_output.name])
			self.assertEqual(observed_format_options[self_contained_output.name], [b"no_ref=1"])

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
		# intervals from the disk-backed cache must therefore not change the reported total_bytes.
		chr9_reference_fasta = get_chr9_reference_fasta()
		with tempfile.NamedTemporaryFile(suffix=".cram") as first_output, \
			  tempfile.NamedTemporaryFile(suffix=".cram") as second_output:

			reader = IntervalReader(
				self._local_cram_path, self._local_cram_path + ".crai",
				reference_fasta_path=chr9_reference_fasta, cache_byte_ranges=True)
			cache_file = reader._byte_ranges_cache._file
			try:
				for interval in self._FXN_intervals:
					reader.add_interval(*interval)

				reader.save_to_file(first_output.name)
				bytes_after_first = reader.get_total_bytes_loaded_from_cram()
				ranges_after_first = reader.get_total_byte_ranges_loaded_from_cram()
				self.assertGreater(bytes_after_first, 0)
				self.assertGreater(os.fstat(cache_file.fileno()).st_size, 0)
				self.assertTrue(all(
					isinstance(location, tuple) for location in reader._byte_ranges_cache._index.values()))

				# the second save over identical intervals is served entirely from the byte-range cache, so it
				# must not increase either the byte total or the byte-range (container) count
				reader.save_to_file(second_output.name)
				self.assertEqual(reader.get_total_bytes_loaded_from_cram(), bytes_after_first)
				self.assertEqual(reader.get_total_byte_ranges_loaded_from_cram(), ranges_after_first)
			finally:
				reader.close()
			self.assertTrue(cache_file.closed)

	def test_cram_reader_does_not_refetch_containers_when_intervals_grow(self):
		# Regression test: the byte-range cache must be keyed on canonical per-container ranges, not on the merged
		# runs those containers get coalesced into. make_minicram calls save_to_file twice -- once for the primary
		# intervals, then again after mate intervals have been added -- and the second call's larger interval set
		# produces runs with different (start, end) boundaries. When the cache was keyed on those run boundaries,
		# every run missed on the second call and re-downloaded bytes already fetched: on the 50000-locus benchmark
		# total network reads reached 133% of the full CRAM size. Growing the interval set must only ever fetch
		# containers that are genuinely new.
		reader = IntervalReader(
			self._local_cram_path, self._local_cram_path + ".crai", cache_byte_ranges=True)
		# Record every range that actually comes off the wire/disk, so overlap between reads is detectable.
		# Asserting on the byte COUNTER alone is not enough: it stays self-consistent under either cache keying,
		# because whatever gets fetched is also what gets stored.
		physical_reads = []
		original_read_bytes = reader._read_bytes

		def recording_read_bytes(start, end):
			physical_reads.append((start, end))
			return original_read_bytes(start, end)

		reader._read_bytes = recording_read_bytes
		try:
			# _load_cram_containers is exercised directly rather than through save_to_file: the container fetching
			# is where the fix lives, and going through save_to_file would additionally require pysam to DECODE the
			# exported containers against a matching reference (the intervals below deliberately span two different
			# contigs, so no single-contig test reference can decode them).
			#
			# The two intervals sit in ADJACENT containers of the test CRAM -- chr5 in the container at byte offset
			# 99943 and chr9 in the one at 101430 -- so adding the second makes those containers merge into a single
			# larger run. That shifting run boundary is precisely what used to invalidate the cache key and
			# re-download the first container. Two intervals inside one container would not exercise it.
			with tempfile.NamedTemporaryFile(suffix=".cram") as first_output, \
				  tempfile.NamedTemporaryFile(suffix=".cram") as second_output:
				reader.add_interval("chr5", 127247318, 127247487)
				reader._load_cram_containers(first_output)
				self.assertGreater(reader.get_total_bytes_loaded_from_cram(), 0)
				self.assertGreater(len(physical_reads), 0)

				reader.add_interval(*self._FXN_intervals[0])
				reader._load_cram_containers(second_output)

			# no byte may be read twice: the total read equals the size of the union of the ranges read
			bytes_read = sum(end - start for start, end in physical_reads)
			bytes_covered = sum(end - start for start, end in merge_adjacent_byte_ranges(physical_reads))
			self.assertEqual(bytes_read, bytes_covered)
		finally:
			reader.close()

	def test_save_to_file_dedup_window_does_not_grow_with_output_size(self):
		# Regression test: save_to_file used to keep one (query_name, flag, reference_start) key for EVERY read it
		# wrote, for the whole call. At ~219 bytes per entry that made it the dominant memory consumer -- an
		# estimated ~8GB of the 13.7GB peak on the 50000-locus benchmark, and the reason the 6.5GB and 13GB VM tiers
		# were OOM-killed. Reads are traversed in coordinate order, so only reads overlapping the current position
		# need to be remembered; the window must therefore track local depth, NOT the number of reads written.
		chr9_reference_fasta = get_chr9_reference_fasta()

		def peak_window_and_reads_written(region_end):
			with tempfile.NamedTemporaryFile(suffix=".cram") as output_cram_file:
				reader = IntervalReader(
					self._local_cram_path, self._local_cram_path + ".crai",
					reference_fasta_path=chr9_reference_fasta)
				try:
					for interval_start in range(69035900, region_end, 25):
						reader.add_interval("chr9", interval_start, interval_start + 10)
					reads_written = reader.save_to_file(output_cram_file.name)
					return reader.get_peak_dedup_window_size(), reads_written
				finally:
					reader.close()

		narrow_window, narrow_reads = peak_window_and_reads_written(69037300)
		wide_window, wide_reads = peak_window_and_reads_written(69038700)

		# the wider region really does write substantially more reads ...
		self.assertGreater(wide_reads, narrow_reads * 1.5)
		# ... but the dedup window tracks read depth, so it must not grow in proportion
		self.assertLess(wide_window, narrow_window * 1.5)
		# and it stays far below the total number of reads written, which is what used to be retained
		self.assertLess(wide_window, wide_reads / 2)

	def test_interval_reader_debug_mode_runs(self):
		# Regression test: a debug print referenced a variable that had been renamed, so constructing an
		# IntervalReader with debug=True raised NameError before doing any work. Nothing else in the suite
		# exercises debug=True, so the whole --debug code path was unprotected.
		reader = IntervalReader(
			self._local_cram_path, self._local_cram_path + ".crai", debug=True, cache_byte_ranges=True)
		try:
			reader.add_interval(*self._FXN_intervals[0])
			with tempfile.NamedTemporaryFile(suffix=".cram") as output_cram_file:
				reader._load_cram_containers(output_cram_file)
		finally:
			reader.close()

	def test_container_and_request_counters_are_distinct(self):
		# Regression test: one counter was used for both "containers loaded" and "read requests issued", but
		# adjacent containers are coalesced into a single request, so the value reported as a container count
		# (including in the --output-data-transfer-stats TSV) undercounted containers.
		reader = IntervalReader(
			self._local_cram_path, self._local_cram_path + ".crai", cache_byte_ranges=True)
		try:
			# these two intervals live in ADJACENT containers, so they coalesce into one request
			reader.add_interval("chr5", 127247318, 127247487)
			reader.add_interval(*self._FXN_intervals[0])
			with tempfile.NamedTemporaryFile(suffix=".cram") as output_cram_file:
				reader._load_cram_containers(output_cram_file)

			self.assertEqual(reader.get_total_containers_loaded_from_cram(), 2)
			# both containers came from a single coalesced request, so requests < containers for this fetch;
			# the reader also issues the version probe and header reads, which are requests but not containers
			self.assertLess(
				reader.get_total_byte_ranges_loaded_from_cram() - 2, reader.get_total_containers_loaded_from_cram())
		finally:
			reader.close()

	def test_cram_version_probe_is_counted_in_transfer_totals(self):
		# Regression test: _detect_cram_eof_container read the first 6 bytes directly instead of through
		# _read_bytes, so those bytes and that request were missing from the reported transfer totals even
		# though _read_bytes is documented as the only place bytes come off the wire.
		reader = IntervalReader(self._local_cram_path, self._local_cram_path + ".crai")
		try:
			# the probe and the CRAM header read both happen during construction, before any interval is added
			self.assertGreaterEqual(reader.get_total_byte_ranges_loaded_from_cram(), 2)
			self.assertGreaterEqual(reader.get_total_bytes_loaded_from_cram(), 6)
		finally:
			reader.close()

	def _run_make_bamlet(self, output_name, regions):
		"""Runs the make_bamlet CLI on the local BAM fixture and returns (output_path, CompletedProcess)."""
		local_bam_path = os.path.join(self._temp_dir.name, "FXN.wgsim_HET_250xGAA.bam")
		if not os.path.isfile(local_bam_path):
			with open(local_bam_path, "wb") as f:
				f.write(pkgutil.get_data("str_analysis", "data/tests/FXN.wgsim_HET_250xGAA.bam"))
			with open(local_bam_path + ".bai", "wb") as f:
				f.write(pkgutil.get_data("str_analysis", "data/tests/FXN.wgsim_HET_250xGAA.bam.bai"))

		output_path = os.path.join(self._temp_dir.name, output_name)
		return output_path, subprocess.run(
			[sys.executable, "-m", "str_analysis.make_bamlet", "-R", get_chr9_reference_fasta(),
			 "-o", output_path, local_bam_path] + regions,
			capture_output=True, text=True)

	def test_make_bamlet_handles_locus_near_contig_start(self):
		# Regression test: main() padded every region by 2000bp without clamping, so any locus within 2000bp of
		# a contig start passed a negative coordinate to pysam.fetch and raised "start out of range".
		output_path, result = self._run_make_bamlet("near_origin.bam", ["chr9:1-2"])
		self.assertNotIn("start out of range", result.stdout + result.stderr)
		self.assertEqual(result.returncode, 0, msg=result.stderr[-500:])

	def test_make_bamlet_does_not_duplicate_reads_across_regions(self):
		# Regression test: extract_region was called once per requested region with its own local read_pairs and
		# no cross-region deduplication, so a read overlapping two regions was written once per region.
		# Duplicates are misread downstream as a complete read pair and suppress real mate retrieval.
		output_path, result = self._run_make_bamlet(
			"two_regions.bam", ["chr9:69037287-69037288", "chr9:69037320-69037321"])
		self.assertEqual(result.returncode, 0, msg=result.stderr[-500:])
		with pysam.AlignmentFile(output_path) as f:
			read_keys = [(read.query_name, read.flag, read.reference_start) for read in f]
		self.assertGreater(len(read_keys), 0)
		self.assertEqual(len(read_keys), len(set(read_keys)))

	def test_make_bamlet_writes_cram_when_cram_output_requested(self):
		# Regression test: a .cram output path was opened in CRAM mode but then unconditionally sorted into a
		# .sorted.bam and renamed over the .cram path, so the tool exited 0 having produced BAM bytes (and a
		# .bai) behind a .cram filename.
		output_path, result = self._run_make_bamlet("out.cram", ["chr9:69037287-69037304"])
		self.assertEqual(result.returncode, 0, msg=result.stderr[-500:])
		with open(output_path, "rb") as f:
			self.assertEqual(f.read(4), b"CRAM")
		self.assertTrue(os.path.isfile(output_path + ".crai"))
		self.assertFalse(os.path.isfile(output_path + ".bai"))

	def test_make_minicram_accepts_a_different_chromosome_naming_convention(self):
		# Regression test: the mate-discovery phase passed the user's raw -L chromosome name straight to
		# pysam.fetch against the extracted temp CRAM, whose header carries the INPUT CRAM's naming convention.
		# Requesting "9:..." against a CRAM whose header says "chr9" therefore raised "invalid contig", even
		# though IntervalReader already normalizes naming for its own fetch calls.
		output_path = os.path.join(self._temp_dir.name, "naming.cram")
		result = subprocess.run(
			[sys.executable, "-m", "str_analysis.make_minicram_for_expansion_hunter",
			 "-R", get_chr9_reference_fasta(), "-o", output_path,
			 "-i", self._local_cram_path + ".crai", "-L", "9:69037287-69037304", self._local_cram_path],
			capture_output=True, text=True)
		# Only the naming failure is asserted, not overall success: once the contig resolves, mate discovery finds
		# mates on OTHER contigs, and the final export then decodes them against this test's chr9-ONLY reference,
		# which some pysam/htslib builds reject. That is a limitation of the single-contig test reference, not of
		# the naming fix, so asserting returncode == 0 here would make this test environment-dependent.
		self.assertNotIn("invalid contig", result.stdout + result.stderr)

	def test_get_total_mapped_reads_counts_cram_records(self):
		# Regression test: the CRAM and BAM paths shared get_index_statistics(), but a CRAI carries no per-contig
		# alignment counts, so every nonempty CRAM was reported as having 0 mapped reads.
		local_bam_path = os.path.join(self._temp_dir.name, "FXN.wgsim_HET_250xGAA.bam")
		with open(local_bam_path, "wb") as f:
			f.write(pkgutil.get_data("str_analysis", "data/tests/FXN.wgsim_HET_250xGAA.bam"))
		with open(local_bam_path + ".bai", "wb") as f:
			f.write(pkgutil.get_data("str_analysis", "data/tests/FXN.wgsim_HET_250xGAA.bam.bai"))

		mapped_in_bam = get_total_mapped_reads(local_bam_path)
		self.assertGreater(mapped_in_bam, 0)
		self.assertEqual(get_total_mapped_reads(self._local_cram_path, scan_cram=True), mapped_in_bam)
		# the scan is opt-in because it decodes the whole file, so the default reports "unavailable" rather than
		# the 0 that index statistics would wrongly yield for a CRAM
		self.assertIsNone(get_total_mapped_reads(self._local_cram_path))

	def test_transfer_counters_available_for_bam_backed_reader(self):
		# Regression test: the transfer counters were initialized only inside __init__'s CRAM branch, so calling
		# the (public, unguarded) getters or save_data_transfer_stats on a BAM-backed reader raised AttributeError.
		local_bam_path = os.path.join(self._temp_dir.name, "FXN.wgsim_HET_250xGAA.bam")
		with open(local_bam_path, "wb") as f:
			f.write(pkgutil.get_data("str_analysis", "data/tests/FXN.wgsim_HET_250xGAA.bam"))
		with open(local_bam_path + ".bai", "wb") as f:
			f.write(pkgutil.get_data("str_analysis", "data/tests/FXN.wgsim_HET_250xGAA.bam.bai"))

		reader = IntervalReader(local_bam_path, local_bam_path + ".bai")
		try:
			self.assertEqual(reader.get_total_bytes_loaded_from_cram(), 0)
			self.assertEqual(reader.get_total_containers_loaded_from_cram(), 0)
			self.assertEqual(reader.get_total_byte_ranges_loaded_from_cram(), 0)
		finally:
			reader.close()

	def test_constructing_reader_does_not_load_containers(self):
		# Regression test: the constructor built its chromosome lookup by calling _load_cram_containers, which
		# appends every CRAI -1 (unmapped) container when include_unmapped_read_pairs is set -- so merely creating
		# a reader downloaded all the unmapped-read containers before a single interval had been added. Only the
		# reference names are needed there, and those come from the already-fetched CRAM header.
		#
		# Asserted by patching _load_cram_containers rather than by counting loaded containers: this test CRAM's
		# CRAI has no -1 (unmapped) entries at all, so a container-count assertion would pass either way.
		with mock.patch.object(IntervalReader, "_load_cram_containers") as mock_load_containers:
			reader = IntervalReader(
				self._local_cram_path, self._local_cram_path + ".crai",
				include_unmapped_read_pairs=True, cache_byte_ranges=True)
			try:
				mock_load_containers.assert_not_called()
				# the header still has to yield a usable chromosome lookup
				self.assertIn("9", reader._chrom_index_lookup)
			finally:
				reader.close()

	def test_get_total_mapped_reads_accepts_an_explicit_index_path(self):
		# Regression test: get_total_mapped_reads reopened the input with no way to supply a nonstandard index,
		# so callers that track their index path separately (parse_motif_composition passes one everywhere else)
		# could not use it here.
		self.assertGreater(
			get_total_mapped_reads(self._local_cram_path, index_filename=self._local_cram_path + ".crai",
								   scan_cram=True), 0)

	def test_save_to_file_does_not_leak_temp_index_when_it_fails(self):
		# Regression test: save_to_file indexes its temp container CRAM early, but the os.remove of that sidecar
		# .crai sat ~100 lines later. Anything raising in between (opening the output, the fetch/write loop) left
		# the sidecar behind, since NamedTemporaryFile only tracks the .cram itself.
		# tempfile.tempdir is set directly rather than via the TMPDIR env var: tempfile caches the resolved temp
		# directory on first use, so setting TMPDIR here would be ignored and the temp files would land elsewhere,
		# making this assertion pass no matter what.
		temp_root = tempfile.mkdtemp()
		before = set(os.listdir(temp_root))
		original_tempdir = tempfile.tempdir
		tempfile.tempdir = temp_root
		try:
			reader = IntervalReader(self._local_cram_path, self._local_cram_path + ".crai")
			try:
				reader.add_interval(*self._FXN_intervals[0])
				# an unwritable output path makes save_to_file raise after the temp index has been created
				with self.assertRaises(Exception):
					reader.save_to_file(os.path.join(temp_root, "no_such_subdir", "out.cram"))
			finally:
				reader.close()
		finally:
			tempfile.tempdir = original_tempdir

		leaked = [name for name in set(os.listdir(temp_root)) - before if name.endswith(".crai")]
		self.assertEqual(leaked, [])

	def test_save_to_file_closes_alignment_handles_when_it_fails(self):
		# Regression test: save_to_file closed its pysam input/output handles only on the success path, so an
		# exception while opening or writing the output (the finally added for the temp file covered only that
		# temp file) left the input AlignmentFile open.
		opened_files = []
		original_alignment_file = pysam.AlignmentFile

		def recording_alignment_file(*args, **kwargs):
			alignment_file = original_alignment_file(*args, **kwargs)
			opened_files.append(alignment_file)
			return alignment_file

		reader = IntervalReader(self._local_cram_path, self._local_cram_path + ".crai")
		with mock.patch.object(pysam, "AlignmentFile", side_effect=recording_alignment_file):
			try:
				reader.add_interval(*self._FXN_intervals[0])
				with self.assertRaises(Exception):
					reader.save_to_file(os.path.join(self._temp_dir.name, "no_such_subdir", "out.cram"))
			finally:
				reader.close()

		self.assertGreater(len(opened_files), 0)
		self.assertEqual([f for f in opened_files if not f.closed], [])

	def test_fetch_uncached_containers_caps_run_length(self):
		# Regression test: CRAI-derived containers are always exactly adjacent (a container's size is
		# next-offset minus this-offset), so every consecutive stretch of requested containers used to coalesce
		# into ONE _read_bytes call with no upper bound -- a genome-wide catalog requests essentially all of
		# them, which meant holding the whole input CRAM in memory as a single bytes object. MAX_RUN_BYTES is
		# patched down here so the cap can be exercised with a small fixture; the assertion is that the run is
		# split into several reads AND that the per-container bytes cached are byte-for-byte what one
		# unrestricted read produced.
		containers = [(0, 400), (400, 800), (800, 1200), (1200, 1600), (1600, 2000)]

		def cached_bytes_with_cap(cap):
			reader = IntervalReader(self._local_cram_path, self._local_cram_path + ".crai",
									cache_byte_ranges=True)
			try:
				read_calls = []
				real_read_bytes = reader._read_bytes

				def counting_read_bytes(start, end):
					read_calls.append((start, end))
					return real_read_bytes(start, end)

				with mock.patch.object(cram_bam_utils, "MAX_RUN_BYTES", cap), \
						mock.patch.object(reader, "_read_bytes", counting_read_bytes):
					reader._fetch_uncached_containers(containers)
				return read_calls, {key: reader._byte_ranges_cache[key] for key in containers}
			finally:
				reader.close()

		uncapped_calls, uncapped_bytes = cached_bytes_with_cap(10**9)
		capped_calls, capped_bytes = cached_bytes_with_cap(800)

		# the whole adjacent stretch is one read when the cap is out of reach, several once it bites
		self.assertEqual(len(uncapped_calls), 1)
		self.assertEqual(uncapped_calls[0], (0, 2000))
		self.assertGreater(len(capped_calls), 1)
		self.assertTrue(all(end - start <= 800 for start, end in capped_calls))

		# capping changes only how the bytes are requested, never which bytes each container ends up with
		self.assertEqual(capped_bytes, uncapped_bytes)
		# and no byte of the range is skipped or fetched twice
		self.assertEqual(sum(end - start for start, end in capped_calls), 2000)

	def test_parse_crai_index_skips_zero_span_records(self):
		# Regression test: the record filter was 'alignment_span < 0', so a zero-span slice got through and
		# produced start == end, which intervaltree rejects as a null interval -- parse_crai_index runs from
		# IntervalReader.__init__, so such a CRAI made construction raise.
		# A zero-span slice must not become a null interval (intervaltree rejects end <= start), but it must
		# ALSO NOT be dropped outright: its container offset participates in container sizing (a container's
		# size is derived as next-offset minus this-offset, so a missing offset inflates its predecessor), and
		# record 0 establishes end_of_cram_header_byte_offset. Unmapped slices (ref id -1) conventionally carry
		# alignment_span 0, so the zero-span-first layout below is the common real-world case.
		cram_path = os.path.join(self._temp_dir.name, "zero_span.cram")
		with open(cram_path, "wb") as f:
			f.write(b"\0" * 100000)
		crai_path = cram_path + ".crai"
		with gzip.open(crai_path, "wt") as f:
			f.write("-1\t1\t0\t1000\t100\t200\n")       # zero-span unmapped slice, first physical record
			f.write("0\t2000\t500\t5000\t100\t200\n")
			f.write("0\t3000\t500\t7000\t100\t200\n")

		end_of_header_offset, interval_trees = parse_crai_index(crai_path, cram_path)

		# the header ends at the FIRST container, even though that record is zero-span
		self.assertEqual(end_of_header_offset, 1000)
		# every container is represented, including the unmapped one
		self.assertEqual(sum(len(tree) for tree in interval_trees.values()), 3)
		self.assertIn(-1, interval_trees)
		# the zero-span record's offset still bounds the next container, so sizes stay correct
		byte_ranges = sorted(
			(interval.data.start, interval.data.end)
			for tree in interval_trees.values() for interval in tree)
		self.assertEqual(byte_ranges[0], (1000, 5000))
		self.assertEqual(byte_ranges[1], (5000, 7000))

	def test_make_bamlet_accepts_a_different_chromosome_naming_convention(self):
		# Regression test: make_bamlet passed the region's chromosome name straight to pysam.fetch without
		# resolving it against the input file's header, so "9:..." against a file whose header says "chr9"
		# raised "invalid contig". The sibling make_minicram tool already normalized names this way.
		output_path, result = self._run_make_bamlet("naming.bam", ["9:69037287-69037304"])
		self.assertNotIn("invalid contig", result.stdout + result.stderr)
		self.assertEqual(result.returncode, 0, msg=result.stderr[-500:])

		# a genuinely absent contig should still be rejected, with a clear message rather than a pysam traceback
		_, missing_contig_result = self._run_make_bamlet("missing.bam", ["chrZZ:1-2"])
		self.assertNotEqual(missing_contig_result.returncode, 0)
		self.assertIn("is not present in", missing_contig_result.stdout + missing_contig_result.stderr)

	def test_merge_adjacent_byte_ranges(self):
		self.assertEqual(merge_adjacent_byte_ranges([]), [])
		# overlapping and exactly-adjacent ranges collapse; disjoint ones stay separate
		self.assertEqual(merge_adjacent_byte_ranges([(0, 10), (10, 20)]), [(0, 20)])
		self.assertEqual(merge_adjacent_byte_ranges([(0, 15), (10, 20)]), [(0, 20)])
		self.assertEqual(merge_adjacent_byte_ranges([(30, 40), (0, 10)]), [(0, 10), (30, 40)])
		# a range fully contained in an earlier one must not shrink the merged run
		self.assertEqual(merge_adjacent_byte_ranges([(0, 100), (10, 20)]), [(0, 100)])

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

	def test_extract_region_clamps_mate_interval_at_contig_start(self):
		# A mate aligned at 0 previously produced (-1, 1), which survived interval-tree lookup and crashed only
		# after all final-pass containers had been downloaded and indexed.
		with tempfile.TemporaryDirectory() as temp_dir:
			input_cram_path = os.path.join(temp_dir, "mate_at_contig_start.cram")
			with pysam.AlignmentFile(
					input_cram_path, "wc",
					header={"HD": {"VN": "1.6", "SO": "coordinate"},
							"SQ": [{"SN": "chr1", "LN": 5000}]},
					format_options=[b"no_ref=1"]) as input_cram:
				for reference_start, flag, mate_start in (
						(0, 1 | 2 | 128, 2000),
						(2000, 1 | 2 | 64, 0),
				):
					read = pysam.AlignedSegment()
					read.query_name = "mate_at_contig_start"
					read.flag = flag
					read.reference_id = 0
					read.reference_start = reference_start
					read.mapping_quality = 60
					read.cigarstring = "50M"
					read.next_reference_id = 0
					read.next_reference_start = mate_start
					read.query_sequence = "A" * 50
					read.query_qualities = pysam.qualitystring_to_array("I" * 50)
					input_cram.write(read)
			pysam.index(input_cram_path)

			with pysam.AlignmentFile(input_cram_path) as input_cram:
				genomic_regions = extract_region("chr1", 1000, 3000, input_bam=input_cram, bamlet=None)
			self.assertEqual(genomic_regions, [("chr1", 1000, 3000), ("chr1", 0, 1)])

			reader = IntervalReader(input_cram_path, input_cram_path + ".crai")
			try:
				for genomic_region in genomic_regions:
					reader.add_interval(*genomic_region)
				output_cram_path = os.path.join(temp_dir, "output.cram")
				self.assertEqual(reader.save_to_file(output_cram_path), 2)
			finally:
				reader.close()

			with pysam.AlignmentFile(output_cram_path) as output_cram:
				self.assertEqual([read.reference_start for read in output_cram], [0, 2000])

	def test_extract_region_skips_unmapped_mates(self):
		with tempfile.TemporaryDirectory() as temp_dir:
			input_bam_path = os.path.join(temp_dir, "unmapped_mate.bam")
			with pysam.AlignmentFile(
					input_bam_path, "wb",
					header={"HD": {"VN": "1.6", "SO": "coordinate"},
							"SQ": [{"SN": "chr1", "LN": 1000}]}) as input_bam:
				read = pysam.AlignedSegment()
				read.query_name = "mapped_read_with_unmapped_mate"
				read.flag = 1 | 8 | 64
				read.reference_id = 0
				read.reference_start = 100
				read.mapping_quality = 60
				read.cigarstring = "50M"
				read.next_reference_id = -1
				read.next_reference_start = -1
				read.query_sequence = "A" * 50
				read.query_qualities = pysam.qualitystring_to_array("I" * 50)
				input_bam.write(read)
			pysam.index(input_bam_path)

			with pysam.AlignmentFile(input_bam_path) as input_bam:
				self.assertEqual(
					extract_region("chr1", 90, 160, input_bam=input_bam, bamlet=None),
					[("chr1", 90, 160)])

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
