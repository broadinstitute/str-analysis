import binascii
import collections
import hashlib
import heapq
import intervaltree
import os
import pysam
import re
import tempfile

from google.cloud import storage

from str_analysis.utils import file_utils
from str_analysis.utils.file_utils import open_file, get_file_size, get_byte_range_from_google_storage


CRAI_FILE_HEADER = [
	"reference_sequence_id", # reference sequence identifier, or -1 for unmapped reads, -2 for multiple reference sequences. All slices in this container must have a reference sequence id matching this value.
	"alignment_start",  # ignored for unmapped slices
	"alignment_span",   # ignored for unmapped slices
	"absolute_container_header_byte_offset",
	"relative_slice_header_byte_offset",
	"slice_size_in_bytes",
]

CRAM_EOF_CONTAINER = binascii.unhexlify("0f000000ffffffff0fe0454f4600000000010005bdd94f0001000606010001000100ee63014b")
# CRAM 2.x uses a shorter (30-byte) EOF marker than the 38-byte CRAM 3 one above
CRAM_EOF_CONTAINER_V2 = binascii.unhexlify("0b000000ffffffff0fe0454f460000000001000001000606010001000100")

# Upper bound on how many bytes _fetch_uncached_containers will pull in a single read. CRAI-derived containers are
# always exactly adjacent (each container's size is next_offset - this_offset), so without a cap every consecutive
# stretch of requested containers coalesces into one read -- and a genome-wide catalog requests essentially all of
# them, which would have meant holding the entire input CRAM in memory as one bytes object. Capping the run bounds
# peak memory without changing the total bytes read; only the number of requests goes up, and only for runs this
# long. Typical jobs never reach it (a 5-locus run's longest run was ~12Mb).
MAX_RUN_BYTES = 500 * 10**6


class ByteRange:
	def __init__(self, start, end):
		self.start = start
		self.end = end
	def __repr__(self):
		return f"ByteRange[{self.start}, {self.end}]"
	def __str__(self):
		return self.__repr__()


class _DiskBackedByteRangesCache:
	"""Stores CRAM byte-range contents in a temporary file instead of memory."""

	def __init__(self):
		self._file = tempfile.TemporaryFile(buffering=0)
		self._index = {}

	def __contains__(self, byte_range):
		return byte_range in self._index

	def __getitem__(self, byte_range):
		cache_location = self._index[byte_range]
		self._file.seek(cache_location[0])
		return self._file.read(cache_location[1])

	def __setitem__(self, byte_range, value):
		self._file.seek(0, os.SEEK_END)
		self._index[byte_range] = (self._file.tell(), len(value))
		self._file.write(value)

	def close(self):
		self._file.close()


def merge_adjacent_byte_ranges(byte_ranges):
	"""Merges overlapping or exactly-adjacent (start, end) byte ranges into the smallest set of covering runs.

	Args:
		byte_ranges: iterable of (start, end) tuples.

	Returns:
		list: sorted, non-overlapping (start, end) tuples covering the same bytes.
	"""
	merged = []
	for start, end in sorted(byte_ranges):
		if merged and merged[-1][1] >= start:
			merged[-1] = (min(merged[-1][0], start), max(merged[-1][1], end))
		else:
			merged.append((start, end))
	return merged


def parse_crai_index(crai_path, cram_path, eof_container_length=len(CRAM_EOF_CONTAINER)):
	"""Takes a .crai and .cram file path (either local or on gs://) and returns a 2-tuple. The first
	value in the tuple is the byte offset where the CRAM header ends in the input CRAM file.
	The second value is a dictionary of interval trees that, for each reference sequence id in the CRAM file,
	provides a way to quickly look up the CRAM file byte range for a given genomic interval.

	Args:
		eof_container_length: byte length of the input CRAM's EOF marker (38 for CRAM 3, 30 for CRAM 2.x),
			subtracted when computing the size of the last container.
	"""
	container_byte_offsets = set()  # collect unique offsets
	end_of_cram_header_byte_offset = None
	interval_trees = {}

	with open_file(crai_path, download_local_copy_before_opening=True, gunzip=True, is_text_file=True) as crai_file:
		for i, line in enumerate(crai_file):
			fields = line.strip().split("\t")
			if len(fields) != len(CRAI_FILE_HEADER):
				raise ValueError(f"Expected {len(CRAI_FILE_HEADER)} columns but found {len(fields)} in line #{i} of {crai_path}: {line}")

			crai_record = dict(zip(CRAI_FILE_HEADER, map(int, fields)))
			reference_sequence_id = crai_record["reference_sequence_id"]
			# CRAI alignment_start is 1-based; convert to 0-based so it matches the 0-based half-open
			# genomic intervals used in the overlap queries (add_interval / _load_cram_containers)
			start = crai_record["alignment_start"] - 1
			end = start + crai_record["alignment_span"]
			absolute_container_header_byte_offset = crai_record["absolute_container_header_byte_offset"]

			# The header boundary and the container offset are recorded for EVERY record, before any filtering.
			# Skipping a record entirely would drop its offset from container_byte_offsets, and since a
			# container's size is derived as (next offset - this offset), a missing offset silently inflates the
			# preceding container's size; skipping record 0 would also push end_of_cram_header_byte_offset out to
			# a later container, making the "header" span most of the file. Unmapped slices
			# (reference_sequence_id -1) conventionally carry alignment_span 0, so this is the common case.
			if i == 0:
				end_of_cram_header_byte_offset = absolute_container_header_byte_offset

			container_byte_offsets.add(absolute_container_header_byte_offset)

			if reference_sequence_id not in interval_trees:
				interval_trees[reference_sequence_id] = intervaltree.IntervalTree()

			if crai_record["alignment_span"] > 0:
				interval_trees[reference_sequence_id].add(
					intervaltree.Interval(start, end, data=ByteRange(absolute_container_header_byte_offset, None)))
			else:
				# A zero/negative-span slice covers no reference bases, so it can never satisfy an overlap query,
				# and intervaltree rejects a null interval (end <= start). It still needs an entry for the -1
				# (unmapped) tree, which _load_cram_containers iterates wholesale rather than querying by overlap,
				# so store it as a minimal 1-base placeholder.
				interval_trees[reference_sequence_id].add(
					intervaltree.Interval(start, start + 1,
										  data=ByteRange(absolute_container_header_byte_offset, None)))

	cram_file_size = get_file_size(cram_path)

	if end_of_cram_header_byte_offset is None:
		# the CRAI has no data-container records (e.g. a CRAM with zero mapped reads); the header container
		# then spans everything up to the EOF marker, so use that as the end-of-header offset
		end_of_cram_header_byte_offset = cram_file_size - eof_container_length

	# compute the container_sizes map which maps container byte offset to container size in bytes
	container_sizes = {}
	previous_offset = None
	for container_byte_offset in list(sorted(container_byte_offsets)):
		if previous_offset is not None:
			container_sizes[previous_offset] = container_byte_offset - previous_offset
		previous_offset = container_byte_offset

	if previous_offset is not None:
		# add last container
		container_sizes[previous_offset] = cram_file_size - previous_offset - eof_container_length

	for interval_tree in interval_trees.values():
		for interval in interval_tree:
			interval.data.end = interval.data.start + container_sizes[interval.data.start]

	return end_of_cram_header_byte_offset, interval_trees


def normalize_chromosome_name(chrom):
	chrom = re.sub("^chr", "", chrom).upper()
	if chrom == "MT":
		chrom = "M"
	return chrom


class IntervalReader:
	"""This class retrieves genomic regions from a CRAM or BAM file and saves them to a local mini-cram or mini-bam.
	The input file can be local or on gs://.

	CRAM inputs ALWAYS use the i/o-efficient path: the CRAI is parsed up front and only the containers overlapping
	the requested interval(s) are read, rather than the larger spans htslib/htsjdk fetch. This matters when reading
	from Google Cloud Nearline storage, where cost is proportional to bytes read. gs:// CRAM byte ranges are fetched
	with the google-cloud-storage client (see get_byte_range_from_google_storage), so no GCS_OAUTH_TOKEN is involved.

	BAM inputs are read through pysam/htslib directly, since the byte-range algorithm is CRAM-specific. A gs:// BAM
	therefore goes through htslib's own Google Storage support, which requires the GCS_OAUTH_TOKEN environment
	variable to be set first, for example by running:

	export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token);
	"""

	def __init__(self,
				 cram_or_bam_path,
				 crai_or_bai_path=None,
				 reference_fasta_path=None,
				 include_unmapped_read_pairs=False,
				 cache_byte_ranges=False,
				 verbose=False,
				 debug=False):
		"""Constructor.

		Args:
			cram_or_bam_path: CRAM or BAM file path. This can be a local or a gs:// path.
			crai_or_bai_path: Optional CRAI  or BAI index file path. This can be a local or a gs:// path.
				If not specified, it will be inferred based on the cram_or_bam_path.
			reference_fasta_path: Optional reference genome FASTA path to use when reading CRAM files using pysam.
			include_unmapped_read_pairs: If True, also output any umapped read pairs stored at the end of the cram file.
			cache_byte_ranges: This option is only relevant for CRAM files, and is ignored for BAMs.
				Cache byte ranges read from the CRAM file in a temporary disk file. This can be useful when loading a
				relatively small number of intervals and then reusing the same IntervalReader instance across multiple
				rounds of CRAM access.
			verbose: If True, will print more detailed information.
			debug: If True, will print even more detailed information.
		"""
		self._cram_or_bam_path = cram_or_bam_path
		self._crai_or_bai_path = crai_or_bai_path
		self._is_cram_file = self._cram_or_bam_path.endswith(".cram")
		self._is_bam_file = self._cram_or_bam_path.endswith(".bam")
		self._is_file_in_google_storage = self._cram_or_bam_path.startswith("gs://")

		if self._is_cram_file:
			if self._crai_or_bai_path is None:
				self._crai_or_bai_path = f"{self._cram_or_bam_path}.crai"
		elif self._is_bam_file:
			if self._crai_or_bai_path is None:
				self._crai_or_bai_path = f"{self._cram_or_bam_path}.bai"
		else:
			raise ValueError(f"Input file {cram_or_bam_path} must end with .cram or .bam")

		self._reference_fasta_path = reference_fasta_path
		self._include_unmapped_read_pairs = include_unmapped_read_pairs
		self._verbose = verbose or debug
		self._debug = debug

		self._byte_ranges_cache = _DiskBackedByteRangesCache() if cache_byte_ranges and self._is_cram_file else None
		self._genomic_intervals = collections.defaultdict(intervaltree.IntervalTree)
		self._peak_dedup_window_size = 0


		# The raw handle below is only ever read by _read_bytes, which serves the CRAM byte-range path; the BAM
		# path opens its own pysam.AlignmentFile. Opening it for BAM too left an unused fd open for the reader's
		# whole lifetime, so both attributes default to None and only the one actually needed is acquired.
		self._storage_client = None
		self._cram_or_bam_file = None
		if self._is_file_in_google_storage:
			# guarded on _is_cram_file for the same reason as the local handle below: only the CRAM byte-range
			# path reads through this client, so a gs:// BAM-backed reader was building one it never used
			if self._is_cram_file:
				# project is passed EXPLICITLY, matching every storage.Client in file_utils. storage.Client
				# distinguishes an omitted project from an explicit None via a private _marker default: an
				# explicit None means "no project" and never raises, while omitting it falls through to
				# _ClientProjectMixin, which raises EnvironmentError when no default project can be resolved
				# from the environment. Omitting it therefore crashed reader construction on any machine
				# without a configured default project, even though the file_exists() checks that run first
				# succeed there. Read off the module (not imported by value) because
				# set_requester_pays_project() assigns it after this module is imported.
				self._storage_client = storage.Client(project=file_utils.gcloud_requester_pays_project)
		else:
			if not os.path.isfile(self._cram_or_bam_path):
				raise ValueError(f"{self._cram_or_bam_path} not found")

			if self._is_cram_file:
				self._cram_or_bam_file = open(self._cram_or_bam_path, "rb")

		# Two distinct quantities, previously conflated: _total_byte_ranges_loaded_from_cram counts PHYSICAL read
		# requests, while _total_containers_loaded_from_cram counts CRAM containers. They differ because adjacent
		# containers are coalesced into one request, so a single request can load many containers.
		# Initialized for BAM inputs too (they simply stay at 0): the getters and save_data_transfer_stats are
		# public and unguarded, so leaving them undefined made a BAM-backed reader raise AttributeError.
		self._total_byte_ranges_loaded_from_cram = 0
		self._total_containers_loaded_from_cram = 0
		self._total_bytes_loaded_from_cram = 0

		if self._is_cram_file:
			# initialize objects used by the self._load_cram_containers(..) method

			# detect the input CRAM's EOF marker (version-dependent length) so the last container is sized
			# correctly and the reconstructed file gets a compatible footer
			self._cram_eof_container = self._detect_cram_eof_container()
			self._end_of_cram_header_byte_offset, self._crai_interval_trees = parse_crai_index(
				self._crai_or_bai_path, self._cram_or_bam_path, eof_container_length=len(self._cram_eof_container))

			# load the CRAM header
			self._cram_header_bytes = self._get_byte_range(0, self._end_of_cram_header_byte_offset)
			self._init_chrom_index_lookup()

	def add_interval(self, chrom, start_0based, end):
		if not isinstance(chrom, str) or not chrom:
			raise ValueError(f"Interval chromosome must be a non-empty string: {chrom}")
		start_0based = max(0, start_0based)
		if end <= start_0based:
			raise ValueError(f"Interval end must be greater than start: {chrom}:{start_0based}-{end}")

		interval_tree = self._genomic_intervals[normalize_chromosome_name(chrom)]
		new_interval = intervaltree.Interval(start_0based, end + 0.1, data=(chrom, start_0based, end))

		overlapping_intervals = interval_tree.overlap(new_interval)
		if not overlapping_intervals:
			interval_tree.add(new_interval)
			return

		merged_start_0based = min([i.data[1] for i in overlapping_intervals] + [start_0based])
		merged_end = max([i.data[2] for i in overlapping_intervals] + [end])
		merged_interval = intervaltree.Interval(
			merged_start_0based, merged_end + 0.1, data=(chrom, merged_start_0based, merged_end))
		for i in overlapping_intervals:
			interval_tree.remove(i)

		interval_tree.add(merged_interval)

	def get_total_byte_ranges_loaded_from_cram(self):
		"""Returns the total number of physical read requests issued against the input CRAM so far.

		This is NOT the number of containers loaded: adjacent containers are coalesced into a single request,
		so one request can cover many containers. Use get_total_containers_loaded_from_cram() for that.
		This only works for CRAM and not for BAM files.
		"""
		return self._total_byte_ranges_loaded_from_cram

	def get_total_containers_loaded_from_cram(self):
		"""Returns the total number of CRAM containers that were loaded from the input CRAM so far.

		Counts each container once, regardless of how many were coalesced into the same read request, and
		excludes containers served from the byte-range cache (those cost no I/O).
		This only works for CRAM and not for BAM files.
		"""
		return self._total_containers_loaded_from_cram

	def get_total_bytes_loaded_from_cram(self):
		"""Returns the total number of bytes that were loaded from the input CRAM so far.
		This is only works for CRAM and not for BAM files.
		"""
		return self._total_bytes_loaded_from_cram

	def get_peak_dedup_window_size(self):
		"""Returns the largest number of read keys save_to_file held at once while deduplicating its output.

		This is bounded by the local read depth rather than by the total number of reads written, and is the
		high-water mark of what was previously this class's dominant memory consumer.
		"""
		return self._peak_dedup_window_size

	def save_to_file(self, local_path, create_index=True, disable_reference_compression=True):
		"""Load and save the added genomic intervals to a local file at the given local_path.

		The output format follows the input file type (a CRAM-backed reader writes CRAM, a BAM-backed reader
		writes BAM), regardless of local_path's extension, so local_path should be given a matching extension.

		Args:
			local_path: Output file path (written as CRAM when the input is a CRAM, as BAM when the input is a BAM).
			create_index: Whether to create an index for the output.
			disable_reference_compression: For CRAM output, use htslib no_ref=1 to avoid external-reference-based
				compression. This increases output size but avoids reference MD5 validation. Ignored when the input
				is a BAM.
		"""
		# temp_cram_container_file auto-deletes on close, but the .crai pysam.index writes beside it is a separate
		# path NamedTemporaryFile does not track. The removal below is in a finally: everything between indexing and
		# cleanup (opening the output, the fetch/write loop, closing files) can raise, and on that path the sidecar
		# used to be left behind in the temp dir.
		# the pysam handles are closed in the same finally for the same reason: they were closed only on the
		# success path, so an exception while opening or writing the output leaked them
		temp_cram_container_file = None
		pysam_input_file = None
		pysam_output_file = None
		partial_path = None
		try:
			if self._is_cram_file:
				temp_cram_container_file = tempfile.NamedTemporaryFile(suffix=".cram")
				if self._debug:
					print(f"DEBUG: Writing containers to {temp_cram_container_file.name}")
				bytes_written = self._load_cram_containers(temp_cram_container_file)
				if bytes_written == 0:
					print(f"WARNING: No CRAM containers were loaded from {self._cram_or_bam_path}")
					return 0

				if self._debug:
					print(f"DEBUG: {temp_cram_container_file.name} has file size "
						  f"{os.path.getsize(temp_cram_container_file.name):,d} bytes")
					print(f"DEBUG: Using pysam to generate a CRAM index for {temp_cram_container_file.name}")
				temp_cram_container_file.seek(0)
				pysam.index(temp_cram_container_file.name)
				if self._debug:
					print(f"DEBUG: Generated CRAM index {temp_cram_container_file.name}.crai (size: {os.path.getsize(temp_cram_container_file.name + '.crai'):,d} bytes)")
				pysam_input_file = pysam.AlignmentFile(
					temp_cram_container_file, require_index=True, reference_filename=self._reference_fasta_path)
				pysam_input_filename = temp_cram_container_file.name
			elif self._is_bam_file:
				pysam_input_file = pysam.AlignmentFile(
					self._cram_or_bam_path, index_filename=self._crai_or_bai_path,
					reference_filename=self._reference_fasta_path, require_index=True,
				)
				# the file actually being fetched from is the BAM, not its index (the CRAM branch above names the
				# container file it opened); this feeds the "DEBUG: Fetching ... from {pysam_input_filename}" print
				pysam_input_filename = self._cram_or_bam_path
			# no else: __init__ already rejects any input that is neither .cram nor .bam, so one of the two
			# branches above always applies

			chom_order = [normalize_chromosome_name(c) for c in pysam_input_file.references]
			# map normalized chrom name -> the exact reference name in this file's header, so fetch() always uses a
			# name valid for the header even when the requested region used a different naming convention (e.g. "9" vs "chr9")
			normalized_to_reference_name = {
				normalize_chromosome_name(name): name for name in pysam_input_file.references
			}
			# no_ref=1 stores read sequences verbatim (like BAM) instead of reference-compressing them. This avoids htslib
			# validating each contig's reference md5 against the header @SQ M5 tag, which fails when the reference build
			# differs from the input alignment or a read's contig cannot be populated from the supplied reference.
			# Everything is written to a sibling ".partial" path and moved into place only once the reads are
			# written AND the index has been built. Writing local_path directly truncated it up front, so an
			# exception anywhere downstream destroyed a previously valid output and could leave its stale index
			# beside a half-written file. The partial lives in the same directory so the final os.replace is an
			# atomic same-filesystem rename.
			partial_path = f"{local_path}.partial.cram" if self._is_cram_file else f"{local_path}.partial.bam"
			pysam_output_file = pysam.AlignmentFile(partial_path, mode="wc" if self._is_cram_file else "wb",
												template=pysam_input_file, reference_filename=self._reference_fasta_path,
												format_options=[b"no_ref=1"]
												if self._is_cram_file and disable_reference_compression else None)
			read_counter = 0
			# A read overlapping two non-overlapping requested intervals is returned by fetch() once per interval;
			# track written alignments so each is emitted at most once (duplicates would otherwise be interpreted
			# as a complete read pair downstream and suppress real mate retrieval).
			#
			# Only a SLIDING WINDOW of keys is retained, not every read written. Intervals are traversed in header
			# reference order then ascending position, and fetch() returns reads in coordinate order, so once the scan
			# passes a read's reference_end that read cannot be returned again -- not by the rest of this interval, nor
			# by any later interval, which starts at an even higher position. Keeping the full history instead made this
			# set the single largest memory consumer: at ~219 bytes per entry it reached an estimated ~8GB of the 13.7GB
			# peak on the 50000-locus benchmark, and OOM-killed the 6.5GB and 13GB VM tiers outright. The window holds
			# only reads overlapping the current position, i.e. on the order of the local coverage depth.
			written_read_keys = set()
			# min-heap of (reference_end, read_key), used to drop keys once the scan has moved past them
			eviction_heap = []
			current_fetch_chrom = None
			if self._verbose:
				print("Writing reads to", local_path)
			# _get_merged_intervals already sorts chromosomes by the supplied header order and intervals by position.
			# Wrapping its result in sorted() would incorrectly restore lexicographic chromosome order (chr10 before chr2).
			for chrom, start, end in self._get_merged_intervals(
					chrom_sort_order=lambda ch: chom_order.index(normalize_chromosome_name(ch))
						if normalize_chromosome_name(ch) in chom_order else len(chom_order)):
				normalized_chrom = normalize_chromosome_name(chrom)
				if normalized_chrom not in normalized_to_reference_name:
					# contig absent from this file's header (e.g. an off-header interval that _load_cram_containers
					# already skipped) — there is nothing to fetch for it
					continue
				fetch_chrom = normalized_to_reference_name[normalized_chrom]
				if self._debug:
					print(f"DEBUG: Fetching {fetch_chrom}:{start}-{end} from {pysam_input_filename}")

				if fetch_chrom != current_fetch_chrom:
					# positions restart on a new contig, so nothing already seen can recur -- reset the whole window
					current_fetch_chrom = fetch_chrom
					written_read_keys.clear()
					eviction_heap.clear()

				for read in pysam_input_file.fetch(fetch_chrom, start, end):
					# drop keys for reads that end before this one starts; coordinate-ordered traversal guarantees
					# they cannot be returned again
					while eviction_heap and eviction_heap[0][0] < read.reference_start:
						written_read_keys.discard(heapq.heappop(eviction_heap)[1])

					read_key = (read.query_name, read.flag, read.reference_start)
					if read_key in written_read_keys:
						continue
					written_read_keys.add(read_key)
					# reference_end is None for an unmapped read placed at its mate's position; such a read spans no
					# reference bases, so key it off reference_start to keep the heap ordering consistent
					heapq.heappush(
						eviction_heap,
						(read.reference_end if read.reference_end is not None else read.reference_start, read_key))
					# high-water mark of the dedup window, exposed via get_peak_dedup_window_size() so the bound can be
					# asserted in tests and eyeballed in production -- this set is what used to grow to multiple GB
					self._peak_dedup_window_size = max(self._peak_dedup_window_size, len(written_read_keys))
					read_counter += 1
					pysam_output_file.write(read)

			if self._include_unmapped_read_pairs:
				# fetch("*") is a single linear pass over the unplaced reads, so it cannot return the same read twice
				# and no deduplication is needed here. Keeping a key per unmapped read would reintroduce exactly the
				# unbounded-set memory growth the sliding window above exists to avoid.
				for read in pysam_input_file.fetch("*"):
					read_counter += 1
					pysam_output_file.write(read)

			# closed here (not only in the finally) because create_index below reads the finished output file
			pysam_input_file.close()
			pysam_output_file.close()
			pysam_input_file = pysam_output_file = None

			index_suffix = ".crai" if self._is_cram_file else ".bai"
			if create_index:
				if self._debug:
					print(f"DEBUG: Using pysam to generate a CRAM index for {local_path}")

				try:
					# The output is already coordinate-sorted: this code assumes the input is coordinate-sorted (the
					# byte-range CRAM algorithm and require_index=True both rely on that), intervals are traversed in
					# header reference order and then coordinate order, and written_read_keys keeps a read spanning two
					# intervals at its first (earlier) occurrence. A regression test exercises this invariant with
					# non-adjacent intervals, an earlier mate, a duplicated boundary read, and unmapped reads.
					if self._is_cram_file:
						# A CRAI is a flat per-container record list, so pysam.index does not require a sort here;
						# re-sorting would only write another full-size CRAM.
						pysam.index(partial_path)
					else:
						# BAI construction (unlike CRAI) hard-fails on out-of-order positions, so keep an explicit
						# sort here as insurance even though the stream above is already coordinate-sorted.
						pysam.sort("-o", f"{partial_path}.sorted.bam", partial_path)
						os.rename(f"{partial_path}.sorted.bam", partial_path)
						pysam.index(partial_path)
					if self._debug:
						# the BAM branch above writes a .bai, not a .crai -- stat'ing the wrong suffix raised
						# FileNotFoundError, which the except below swallowed into a false "Failed to prepare and
						# index" warning even though indexing had succeeded
						print(f"DEBUG: Generated index {partial_path}{index_suffix} "
							  f"(size: {os.path.getsize(partial_path + index_suffix):,d} bytes)")
				except Exception as e:
					# raised, not warned: the index is not optional. Downstream tools (ExpansionHunter, pysam
					# fetch) cannot use an unindexed output, and a stale .crai left beside a freshly written
					# CRAM is worse than none. Swallowing this returned a nonzero read count, so callers -- and
					# make_minicram's "did the file get written" check -- reported success for an unusable file.
					raise RuntimeError(f"Failed to prepare and index {local_path}: {e}") from e

			# the output and its index are only now moved into place, together. os.replace is atomic within a
			# filesystem, so a reader either sees the previous output or the new one, never a half-written file.
			os.replace(partial_path, local_path)
			if create_index:
				os.replace(f"{partial_path}{index_suffix}", f"{local_path}{index_suffix}")
			elif os.path.isfile(f"{local_path}{index_suffix}"):
				# no index was built, so any index already sitting here belongs to the file just replaced and
				# would silently mis-describe the new one
				os.remove(f"{local_path}{index_suffix}")

			if self._verbose:
				print(f"Wrote {read_counter:,d} reads to {local_path}")

			return read_counter

		finally:
			for pysam_file in (pysam_input_file, pysam_output_file):
				if pysam_file is not None:
					pysam_file.close()
			# on the success path these were already renamed away, so this only fires when something failed
			if partial_path is not None:
				for leftover in (partial_path, f"{partial_path}.crai", f"{partial_path}.bai",
								 f"{partial_path}.sorted.bam"):
					if os.path.isfile(leftover):
						os.remove(leftover)
			if temp_cram_container_file is not None:
				temp_crai_path = f"{temp_cram_container_file.name}.crai"
				if os.path.isfile(temp_crai_path):
					os.remove(temp_crai_path)
				temp_cram_container_file.close()

	def save_data_transfer_stats(self, stats_tsv_path=None):
		# the counters below are only incremented by the CRAM byte-range path, so for a BAM input they are all
		# still 0 and a report would read as "nothing was transferred" even though the BAM was read in full
		# through pysam/htslib. Say so instead of writing a file of zeros.
		if not self._is_cram_file:
			print(f"Data transfer stats are only tracked for CRAM inputs; skipping them for "
				  f"{self._cram_or_bam_path}")
			return

		total_bytes = self.get_total_bytes_loaded_from_cram()
		total_containers = self.get_total_containers_loaded_from_cram()
		total_requests = self.get_total_byte_ranges_loaded_from_cram()
		print(f"Downloaded {total_containers:,d} CRAM containers ({total_requests:,d} read requests) "
			  f"and {total_bytes/10**6:0,.1f}Mb total")
		if not stats_tsv_path:
			stats_tsv_path = re.sub(".cram$", "", os.path.basename(self._cram_or_bam_path))
			stats_tsv_path += ".data_transfer_stats.tsv"

		# total_read_requests is appended last so the first two numeric columns keep their existing positions,
		# for readers that parse this file positionally
		with open(stats_tsv_path, "wt") as stats_file:
			stats_file.write("\t".join([
				"file_path", "total_bytes_loaded", "total_containers_loaded_from_cram", "total_read_requests"]) + "\n")
			stats_file.write("\t".join(map(str, [
				self._cram_or_bam_path, total_bytes, total_containers, total_requests])) + "\n")

		print(f"Wrote stats to {stats_tsv_path}")


	def close(self):
		"""Release any system resources opened by IntervalReader"""
		if self._byte_ranges_cache is not None:
			self._byte_ranges_cache.close()
			self._byte_ranges_cache = None
		if self._storage_client is not None:
			self._storage_client.close()
			self._storage_client = None
		if self._cram_or_bam_file is not None:
			self._cram_or_bam_file.close()
			self._cram_or_bam_file = None


	def _load_cram_containers(self, cram_container_file):
		"""Load all CRAM containers that overlap the genomic intervals added so far to the given file handle.
		Args:
			cram_container_file (file): A file handle which is already open for writing in binary mode

		Returns:
			int: Number of bytes written to the given file handle (not counting the CRAM header and footer)
		"""
		if self._verbose:
			print(f"Copying requested CRAM containers to {cram_container_file.name}")

		# use the CRAI index to compute which byte ranges to load from the CRAM file
		byte_ranges_to_load = []
		for chrom, start_0based, end in self._get_merged_intervals():
			reference_sequence_id = self._chrom_index_lookup.get(normalize_chromosome_name(chrom))
			if reference_sequence_id is None or reference_sequence_id not in self._crai_interval_trees:
				# skip this interval (a contig absent from the CRAM header, or a mate on a decoy contig with
				# no CRAI entries) rather than aborting the whole export, which would discard all other
				# valid intervals' containers
				print(f"WARNING: No CRAI entries found for {chrom} (reference_sequence_id={reference_sequence_id}); skipping {chrom}:{start_0based}-{end}")
				continue

			crai_interval_tree = self._crai_interval_trees[reference_sequence_id]
			overlapping_crai_intervals = list(crai_interval_tree.overlap(start_0based, end))
			if not overlapping_crai_intervals:
				print(f"WARNING: None of the {len(crai_interval_tree)} CRAI entries on {chrom} overlap {chrom}:{start_0based}-{end}; skipping this interval")
				continue

			if self._verbose:
				print(f"Found {len(overlapping_crai_intervals):4,d} CRAI entries that overlapped {chrom}:{start_0based}-{end}")

			for crai_interval in overlapping_crai_intervals:
				byte_ranges_to_load.append((crai_interval.data.start, crai_interval.data.end))

		if self._include_unmapped_read_pairs and self._crai_interval_trees.get(-1):
			crai_interval_tree = self._crai_interval_trees[-1]
			crai_intervals = list(crai_interval_tree)
			if self._verbose:
				print(f"Found {len(crai_intervals):4,d} CRAI entries that contained only unmapped read pairs")

			for crai_interval in crai_intervals:
				byte_ranges_to_load.append((crai_interval.data.start, crai_interval.data.end))

		# Deduplicate to canonical per-container byte ranges. These come straight from the CRAI, so a given container
		# always has the SAME (start, end) key no matter which intervals selected it -- which is what makes the cache
		# reusable across save_to_file calls. Merging containers into larger runs BEFORE caching (as this used to do)
		# keyed the cache on run boundaries instead, and those boundaries shift when the interval set changes between
		# calls, so the second call missed on every run and re-downloaded bytes it already had. On the 50000-locus
		# benchmark that pushed total network reads to 133% of the full CRAM size -- worse than downloading it once,
		# which defeats the entire point of building a minicram.
		containers_to_load = sorted(set(byte_ranges_to_load))
		# snapshot taken before any bytes are read, so the debug print below can report the number of PHYSICAL reads
		# this call issued (every read from here on goes through _read_bytes, which increments this counter)
		byte_ranges_loaded_before = self._total_byte_ranges_loaded_from_cram
		if self._byte_ranges_cache is not None:
			# pre-fetch the containers this call needs but does not already have, coalescing contiguous misses so the
			# request count stays close to what merging used to give, then serve every container from the cache below
			self._fetch_uncached_containers(containers_to_load)
			byte_ranges_to_write = containers_to_load
		else:
			# no cache to reuse across calls, so there is nothing to keep canonical -- merge into runs to keep the
			# number of network requests down, exactly as before. Every container is (re-)loaded on each call here,
			# since without a cache nothing carries over.
			self._total_containers_loaded_from_cram += len(containers_to_load)
			byte_ranges_to_write = merge_adjacent_byte_ranges(containers_to_load)

		# write all CRAM containers that overlap the requested intervals to the given CRAM file
		if self._debug:
			print(f"DEBUG: Wrote CRAM header [0-{len(self._cram_header_bytes):,d}] (hash: {hashlib.md5(self._cram_header_bytes).hexdigest()[:20]})")
		cram_container_file.write(self._cram_header_bytes)
		bytes_counter = 0
		for bytes_start, bytes_end in byte_ranges_to_write:
			cram_container_bytes = self._get_byte_range(bytes_start, bytes_end)
			cram_container_file.write(cram_container_bytes)
			bytes_counter += len(cram_container_bytes)
			if self._debug:
				print(f"DEBUG: Wrote byte range [{bytes_start:,d}-{bytes_end:,d}] (hash: {hashlib.md5(cram_container_bytes).hexdigest()[:20]})")

		cram_container_file.write(self._cram_eof_container)
		if self._debug:
			print(f"DEBUG: Wrote EOF container [0-{len(self._cram_eof_container):,d}] (hash: {hashlib.md5(self._cram_eof_container).hexdigest()[:20]})")

		cram_container_file.flush()
		if self._debug:
			# the request count is the measured number of _read_bytes calls, not len(byte_ranges_to_write): on the
			# cache branch byte_ranges_to_write IS containers_to_load, so its length is the container count over
			# again, and containers served from the cache cost no read at all
			print(f"DEBUG: Wrote {len(containers_to_load):,d} CRAM containers "
				  f"({self._total_byte_ranges_loaded_from_cram - byte_ranges_loaded_before:,d} read requests) "
				  f"to {cram_container_file.name}")

		if self._verbose and self._total_byte_ranges_loaded_from_cram > 0:
			print(f"Copied {self._total_bytes_loaded_from_cram/10**6:5.1f}Mb total from {self._cram_or_bam_path}  to  {cram_container_file.name}")

		return bytes_counter

	def _init_chrom_index_lookup(self):
		"""Initialize a chromosome name to chromosome index lookup dictionary by reading the CRAM or BAM header and
		parsing the list of reference sequence names and in order.
		"""

		# Write just the CRAM header + EOF marker to a temp file. This used to call _load_cram_containers, which
		# does far more than needed here: with include_unmapped_read_pairs set it appends EVERY CRAI -1 (unmapped)
		# container, so merely constructing a reader downloaded all the unmapped-read containers before a single
		# interval had been added. Only the reference names are needed, and those live in the header.
		with tempfile.NamedTemporaryFile(suffix=".cram") as cram_header_file:
			cram_header_file.write(self._cram_header_bytes)
			cram_header_file.write(self._cram_eof_container)
			cram_header_file.flush()

			# By default, if a file is opened in mode 'r', it is checked for a valid header (`check_header` = True)
			# and a definition of chromosome names (`check_sq` = True).
			with pysam.AlignmentFile(cram_header_file.name) as file:
				references = file.references

		# create a lookup table for looking up the chromosome index of a given chromosome name
		self._chrom_index_lookup = {
			normalize_chromosome_name(name): idx for idx, name in enumerate(references)
		}

		if len(self._chrom_index_lookup) == 0:
			raise ValueError(f"Unable to load chromosome list from the header of {self._cram_or_bam_path}")

	def _get_merged_intervals(self, chrom_sort_order=lambda chrom: chrom):
		"""Return a list of genomic intervals that reperesent the union of the intervals provided by the user so far
		via the add_interval(..) method. The intervals are represented as 3-tuples (chrom, start_0based, end),
		and the sorted order is by genomic position, with the chromosome sort order specified by the optional
		chrom_sort_order argument.

		Args:
			chrom_sort_order: A function that takes a chromosome name and returns a value that can be used for sorting
				the chromosome names. The default is to sort the chromosomes by their lexicographic order.

		Return:
			list: A list of 3-tuples (chrom, start_0based, end) representing the merged genomic intervals.
		"""

		return [
			i.data
			for chrom, genomic_intervals in sorted(self._genomic_intervals.items(), key=lambda t: chrom_sort_order(t[0]))
			for i in sorted(genomic_intervals)
		]

	def _detect_cram_eof_container(self):
		"""Reads the input CRAM's magic bytes and returns the EOF marker matching its major version
		(CRAM 2.x uses the 30-byte marker, CRAM 3 the 38-byte one). Defaults to the CRAM 3 marker for
		anything that doesn't look like a CRAM 2 file."""
		# routed through _read_bytes so this probe is included in the transfer counters -- it is a real (if tiny)
		# network read for gs:// inputs, and _read_bytes is documented as the single place bytes come off the wire
		magic_and_version = self._read_bytes(0, 6)

		if magic_and_version[:4] == b"CRAM" and len(magic_and_version) >= 5 and magic_and_version[4] == 2:
			return CRAM_EOF_CONTAINER_V2
		return CRAM_EOF_CONTAINER

	def _fetch_uncached_containers(self, containers):
		"""Fetches every container not already in the byte-range cache, and stores each one under its own key.

		Contiguous runs of missing containers are fetched in a single request, so this keeps the request count close
		to what fetching merged runs gave, while still caching at the per-container granularity that lets a later
		save_to_file call reuse them regardless of how its intervals merge.

		Args:
			containers: sorted list of canonical per-container (start, end) byte ranges needed by this call.
		"""
		missing = [byte_range for byte_range in containers if byte_range not in self._byte_ranges_cache]
		if not missing:
			return
		# containers served from the cache cost no I/O, so only the missing ones count as loaded
		self._total_containers_loaded_from_cram += len(missing)

		# group the missing containers into runs of exactly-adjacent containers (end of one == start of the next),
		# keeping each run's member containers so the fetched bytes can be split back up and cached individually.
		# A run is also cut off once it would exceed MAX_RUN_BYTES, since the whole run is held in memory below.
		# A single container larger than the cap still gets its own run and is read whole -- a container cannot be
		# split without breaking the per-container cache keys -- but CRAM containers are far smaller than the cap.
		runs = []
		for start, end in missing:
			if runs and runs[-1][-1][1] == start and end - runs[-1][0][0] <= MAX_RUN_BYTES:
				runs[-1].append((start, end))
			else:
				runs.append([(start, end)])

		for run in runs:
			run_start, run_end = run[0][0], run[-1][1]
			# read the run directly rather than via _get_byte_range, so the run itself is never stored under its own
			# (non-canonical) key -- that would both defeat the point of this fix and store the bytes twice on disk
			run_bytes = self._read_bytes(run_start, run_end)
			for start, end in run:
				self._byte_ranges_cache[(start, end)] = run_bytes[start - run_start:end - run_start]

	def _read_bytes(self, start, end):
		"""Reads [start, end) from the input file (over the network for gs:// paths) and updates the transfer counters.

		This is the only place bytes actually come off the wire, so the counters it maintains are what
		--output-data-transfer-stats reports; cache hits deliberately never reach here.
		"""
		if self._is_file_in_google_storage:
			byte_range = get_byte_range_from_google_storage(
				self._cram_or_bam_path, start, end, client=self._storage_client)
		else:
			self._cram_or_bam_file.seek(start)
			byte_range = self._cram_or_bam_file.read(end - start)

		self._total_byte_ranges_loaded_from_cram += 1
		self._total_bytes_loaded_from_cram += len(byte_range)

		if len(byte_range) != end - start:
			raise ValueError(f"Expected to read {end - start} bytes (from {start} to {end}) but read {len(byte_range)}")

		return byte_range

	def _get_byte_range(self, start, end):
		if self._byte_ranges_cache is not None and (start, end) in self._byte_ranges_cache:
			return self._byte_ranges_cache[(start, end)]

		byte_range = self._read_bytes(start, end)

		if self._byte_ranges_cache is not None:
			self._byte_ranges_cache[(start, end)] = byte_range

		return byte_range



def get_total_mapped_reads(input_bam_or_cram_path, reference_fasta_path=None, index_filename=None,
                           scan_cram=False):
    """Returns the total number of mapped reads in a BAM/CRAM file, or None when a CRAM count is unavailable.

    For BAM this is the fast index-statistics path. A CRAI carries no per-contig alignment counts, so
    get_index_statistics() reports mapped=0 for every contig and a shared index-statistics path would silently
    return 0 for any nonempty CRAM. The only way to get a real CRAM count is to decode every record, which is a
    whole-file scan -- far more expensive than the locus-scoped work its callers otherwise do. So CRAM returns
    None (unavailable) unless the caller explicitly asks for the scan.

    Args:
        input_bam_or_cram_path: BAM or CRAM path.
        reference_fasta_path: Reference FASTA. Required for a reference-compressed CRAM, which cannot be decoded
            without it; ignored for BAM.
        index_filename: Optional index path, for inputs whose index does not sit at the default location.
        scan_cram: If True, decode the whole CRAM to count mapped records. Off by default because of the cost.

    Returns:
        int: number of mapped reads, or None for a CRAM when scan_cram is False.
    """
    if input_bam_or_cram_path.endswith(".cram"):
        if not scan_cram:
            return None
        with pysam.AlignmentFile(input_bam_or_cram_path, "rc", index_filename=index_filename,
                                 reference_filename=reference_fasta_path) as cram:
            return sum(1 for read in cram.fetch(until_eof=True) if not read.is_unmapped)

    mapped = 0
    with pysam.AlignmentFile(input_bam_or_cram_path, "rb", index_filename=index_filename) as bam:
        for contig in bam.get_index_statistics():
            mapped += contig.mapped

    return mapped
