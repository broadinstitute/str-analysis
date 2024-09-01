import binascii
import collections
import hashlib
import intervaltree
import os
import pysam
import re
import tempfile

from google.cloud import storage

from str_analysis.utils.file_utils import file_exists, open_file, get_file_size, get_byte_range_from_google_storage


CRAI_FILE_HEADER = [
	"reference_sequence_id", # reference sequence identifier, or -1 for unmapped reads, -2 for multiple reference sequences. All slices in this container must have a reference sequence id matching this value.
	"alignment_start",  # ignored for unmapped slices
	"alignment_span",   # ignored for unmapped slices
	"absolute_container_header_byte_offset",
	"relative_slice_header_byte_offset",
	"slice_size_in_bytes",
]

CRAM_EOF_CONTAINER = binascii.unhexlify("0f000000ffffffff0fe0454f4600000000010005bdd94f0001000606010001000100ee63014b")


class ByteRange:
	def __init__(self, start, end):
		self.start = start
		self.end = end
	def __repr__(self):
		return f"ByteRange[{self.start}, {self.end}]"
	def __str__(self):
		return self.__repr__()


def parse_crai_index(crai_path, cram_path):
	"""Takes a .crai and .cram file path (either local or on gs://) and returns a 2-tuple. The first
	value in the tuple is the byte offset where the CRAM header ends in the input CRAM file.
	The second value is a dictionary of interval trees that, for each reference sequence id in the CRAM file,
	provides a way to quickly look up the CRAM file byte range for a given genomic interval.
	"""
	container_byte_offsets = set()  # collect unique offsets
	end_of_cram_header_byte_offset = None
	interval_trees = {}

	with open_file(crai_path, download_local_copy_before_opening=True, gunzip=True, is_text_file=True) as crai_file:
		for i, line in enumerate(crai_file):
			fields = line.strip().split("\t")
			if len(fields) != len(CRAI_FILE_HEADER):
				p.error(f"Expected {len(CRAI_FILE_HEADER)} columns but found {len(fields)} in line #{i} of {crai_path}: {line}")

			crai_record = dict(zip(CRAI_FILE_HEADER, map(int, fields)))
			if crai_record["alignment_span"] < 0:
				continue

			reference_sequence_id = crai_record["reference_sequence_id"]
			start = crai_record["alignment_start"]
			end = start + crai_record["alignment_span"]
			absolute_container_header_byte_offset = crai_record["absolute_container_header_byte_offset"]

			if i == 0:
				end_of_cram_header_byte_offset = absolute_container_header_byte_offset

			interval = intervaltree.Interval(start, end, data=ByteRange(absolute_container_header_byte_offset, None))

			if reference_sequence_id not in interval_trees:
				interval_trees[reference_sequence_id] = intervaltree.IntervalTree()
			interval_trees[reference_sequence_id].add(interval)

			container_byte_offsets.add(absolute_container_header_byte_offset)

	cram_file_size = get_file_size(cram_path)

	# compute the container_sizes map which maps container byte offset to container size in bytes
	container_sizes = {}
	previous_offset = None
	for container_byte_offset in list(sorted(container_byte_offsets)):
		if previous_offset is not None:
			container_sizes[previous_offset] = container_byte_offset - previous_offset
		previous_offset = container_byte_offset

	if previous_offset is not None:
		# add last container
		container_sizes[previous_offset] = cram_file_size - previous_offset - len(CRAM_EOF_CONTAINER)

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
	The original file can be local or remove. For CRAM files, this class optionally uses a more i/o-efficient algorithm
	than htslib/htsjdk/samtools, reading only containers that overlap the requesteed interval(s), while htslib/htsjdk
	read more than that for some reason. This becomes important when reading from Google Cloud Nearline storage which
	currently costs $0.01 per gigabyte to read.

	In default mode, when not using the i/o-efficient algorithm, this class can still read CRAMs or BAMs directly from
	Google Storage by using the pysam/htslib library. In that case, the GCS_OAUTH_TOKEN environment variable must be
	set first, for example by running:

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
				Cache the byte ranges that are read from the CRAM file in memory. This can be useful when loading a
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

		self._byte_ranges_cache = {} if cache_byte_ranges else None
		self._genomic_intervals = collections.defaultdict(intervaltree.IntervalTree)


		if self._is_file_in_google_storage:
			self._storage_client = storage.Client()
		else:
			if not os.path.isfile(self._cram_or_bam_path):
				raise ValueError(f"{self._cram_or_bam_path} not found")

			self._cram_or_bam_file = open(self._cram_or_bam_path, "rb")

		if self._is_cram_file:
			# initialize objects used by the self._load_cram_containers(..) method
			self._total_byte_ranges_loaded_from_cram = 0
			self._total_bytes_loaded_from_cram = 0

			self._end_of_cram_header_byte_offset, self._crai_interval_trees = parse_crai_index(
				self._crai_or_bai_path, self._cram_or_bam_path)

			# load the CRAM header
			self._cram_header_bytes = self._get_byte_range(0, self._end_of_cram_header_byte_offset)
			self._init_chrom_index_lookup()

	def set_verbose(self, verbose):
		self._verbose = verbose

	def add_interval(self, chrom, start_0based, end):
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
		"""Returns the total number of CRAM containers that were loaded from the input CRAM so far.
		This only works for CRAM and not for BAM files.
		"""
		return self._total_byte_ranges_loaded_from_cram

	def get_total_bytes_loaded_from_cram(self):
		"""Returns the total number of bytes that were loaded from the input CRAM so far.
		This is only works for CRAM and not for BAM files.
		"""
		return self._total_bytes_loaded_from_cram

	def clear_intervals(self):
		"""Clears all the genomic intervals that were added so far."""
		self._genomic_intervals = collections.defaultdict(intervaltree.IntervalTree)

	def reset_total_byte_ranges_loaded_from_cram_counter(self):
		"""Resets the counter of the total number of byte ranges that were retrieved from the input CRAM so far."""
		self._total_byte_ranges_loaded_from_cram = 0

	def reset_total_bytes_loaded_from_cram(self):
		"""Resets the counter of the total number of bytes loaded from the input CRAM so far.
		This is only works for CRAM and not for BAM files
		"""
		self._total_bytes_loaded_from_cram = 0

	def save_to_file(self, local_path, create_index=True):
		"""Load and save the added genomic intervals to a local CRAM file at the given local_path, or alternatively,
		write them into the given fileobj, which must be a file object already open for writing.
		"""
		if self._is_cram_file:
			local_path_suffix = ".cram"
			temp_cram_container_file = tempfile.NamedTemporaryFile(suffix=local_path_suffix)
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
			local_path_suffix = ".bam"
			pysam_input_file = pysam.AlignmentFile(
				self._cram_or_bam_path, index_filename=self._crai_or_bai_path,
				reference_filename=self._reference_fasta_path, require_index=True,
			)
			pysam_input_filename = os.path.basename(self._crai_or_bai_path)
		else:
			raise ValueError(f"Output path {local_path} must end with .cram or .bam")

		chom_order = [normalize_chromosome_name(c) for c in pysam_input_file.references]
		pysam_output_file = pysam.AlignmentFile(local_path, mode="wc" if self._is_cram_file else "wb",
												template=pysam_input_file, reference_filename=self._reference_fasta_path)
		read_counter = 0
		if self._verbose:
			print("Writing reads to", local_path)
		for chrom, start, end in sorted(self._get_merged_intervals(chrom_sort_order=lambda chrom: chom_order.index(normalize_chromosome_name(chrom)))):
			if self._debug:
				print(f"DEBUG: Fetching {chrom}:{start}-{end} from {pysam_input_filename}")

			for read in pysam_input_file.fetch(chrom, start, end):
				read_counter += 1
				pysam_output_file.write(read)

		if self._include_unmapped_read_pairs:
			for read in pysam_input_file.fetch("*"):
				read_counter += 1
				pysam_output_file.write(read)

		pysam_input_file.close()
		pysam_output_file.close()

		if create_index:
			if self._debug:
				print(f"DEBUG: Using pysam to generate a CRAM index for {local_path}")

			try:
				pysam.sort("-o", f"{local_path}.sorted.{local_path_suffix}", local_path)
				os.rename(f"{local_path}.sorted.{local_path_suffix}", local_path)
				pysam.index(local_path)
				if self._debug:
					print(f"DEBUG: Generated CRAM index {local_path}.crai (size: {os.path.getsize(local_path + '.crai'):,d} bytes)")
			except Exception as e:
				print(f"WARNING: Failed to sort and index {local_path}: {e}")


		if self._verbose:
			print(f"Wrote {read_counter:,d} reads to {local_path}")

		return read_counter

	def count_reads(self):
		"""Returns the number of reads that overlap the added genomic intervals"""
		with tempfile.NamedTemporaryFile(suffix=".cram") as temp_cram_file:
			read_count = self.save_to_file(temp_cram_file.name, create_index=False)
		return read_count

	def get_header(self):
		"""Returns the CRAM or BAM file header as a dictionary"""
		saved_genomic_intervals = self._genomic_intervals  # save the genomic intervals
		self.clear_intervals()
		with tempfile.NamedTemporaryFile(suffix=".cram") as temp_cram_file:
			self.save_to_file(temp_cram_file.name, create_index=True)
			self._genomic_intervals = saved_genomic_intervals # restore the genomic intervals
			with pysam.AlignmentFile(temp_cram_file.name) as input_file:
				return input_file.header.to_dict()

	def compute_read_stats(self):
		"""Returns a dictionary of read stats for the added genomic intervals. The dictionary contains the following
		keys:
			"coverage": The average coverage across all the added intervals.
			"read_length": The read length.
			"read_count": The total number of reads across all the added intervals.
		"""

		stats = {}
		read_length = None
		read_count = 0
		coverage = []

		with tempfile.NamedTemporaryFile(suffix=".cram") as temp_cram_file:
			self.save_to_file(temp_cram_file.name, create_index=True)

			with pysam.AlignmentFile(temp_cram_file.name) as input_file:
				for chrom, start_0based, end in self._get_merged_intervals():
					bases_in_interval = 0
					for read in input_file.fetch(chrom, start_0based, end):
						if read.is_unmapped or read.is_secondary:
							continue

						current_read_length = read.infer_query_length()
						if current_read_length is None:
							continue

						if read_length is None or read_length < current_read_length:
							read_length = current_read_length

						read_start_1based = read.reference_start + 1
						read_end_1based = read.reference_start + current_read_length

						if read_end_1based < start_0based + 1 or read_start_1based > end:
							continue

						read_count += 1
						bases_in_interval += min(read_end_1based, end) - max(read_start_1based, start_0based + 1) + 1

					if end - start_0based > 0:
						coverage.append(bases_in_interval/float(end - start_0based))

		stats["coverage"] = sum(coverage)/len(coverage) if coverage else 0
		stats["read_length"] = read_length
		stats["read_count"] = read_count

		return stats

	def close(self):
		"""Release any system resources opened by IntervalReader"""
		if self._is_file_in_google_storage:
			self._storage_client.close()
		else:
			self._cram_or_bam_file.close()


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
			reference_sequence_id = self._chrom_index_lookup[normalize_chromosome_name(chrom)]
			if reference_sequence_id not in self._crai_interval_trees:
				print(f"WARNING: No CRAI entries found for {chrom} (reference_sequence_id={reference_sequence_id})")
				return 0

			crai_interval_tree = self._crai_interval_trees[reference_sequence_id]
			overlapping_crai_intervals = list(crai_interval_tree.overlap(start_0based, end))
			if not overlapping_crai_intervals:
				print(f"WARNING: None of the {len(crai_interval_tree)} CRAI entries on {chrom} overlap {chrom}:{start_0based}-{end}")
				return 0

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

		# compute merged byte-range intervals
		merged_byte_ranges = []
		for bytes_start, bytes_end in sorted(byte_ranges_to_load):
			if not merged_byte_ranges:
				merged_byte_ranges.append((bytes_start, bytes_end))
			else:
				previous_interval = merged_byte_ranges[-1]
				if previous_interval[1] >= bytes_start:
					merged_byte_ranges[-1] = (min(previous_interval[0], bytes_start), max(previous_interval[1], bytes_end))
				else:
					merged_byte_ranges.append((bytes_start, bytes_end))

		# write all CRAM containers that overlap the requested intervals to the given CRAM file
		if self._debug:
			print(f"DEBUG: Wrote CRAM header [0-{len(self._cram_header_bytes):,d}] (hash: {hashlib.md5(self._cram_header_bytes).hexdigest()[:20]})")
		cram_container_file.write(self._cram_header_bytes)
		bytes_counter = 0
		for bytes_start, bytes_end in merged_byte_ranges:
			cram_container_bytes = self._get_byte_range(bytes_start, bytes_end)
			cram_container_file.write(cram_container_bytes)
			bytes_counter += len(cram_container_bytes)
			if self._debug:
				print(f"DEBUG: Wrote byte range [{bytes_start:,d}-{bytes_end:,d}] (hash: {hashlib.md5(cram_container_bytes).hexdigest()[:20]})")

		cram_container_file.write(CRAM_EOF_CONTAINER)
		if self._debug:
			print(f"DEBUG: Wrote EOF container [0-{len(CRAM_EOF_CONTAINER):,d}] (hash: {hashlib.md5(CRAM_EOF_CONTAINER).hexdigest()[:20]})")

		cram_container_file.flush()
		if self._debug:
			print(f"DEBUG: Wrote {len(merged_byte_ranges):,d} CRAM containers to {cram_container_file.name}")

		if self._verbose and self._total_byte_ranges_loaded_from_cram > 0:
			print(f"Copied {self._total_bytes_loaded_from_cram/10**6:5.1f}Mb total from {self._cram_or_bam_path}  to  {cram_container_file.name}")

		return bytes_counter

	def _init_chrom_index_lookup(self):
		"""Initialize a chromosome name to chromosome index lookup dictionary by reading the CRAM or BAM header and
		parsing the list of reference sequence names and in order.
		"""

		# load the CRAM header of the input CRAM file into a local temp file
		with tempfile.NamedTemporaryFile(suffix=".cram") as cram_header_file:
			# load the CRAM header from the temp file using pysam to get the list of reference sequence names
			was_verbose = self._verbose
			self._verbose = False
			bytes_written = self._load_cram_containers(cram_header_file)
			self._verbose = was_verbose

			if bytes_written > 0:
				pysam.index(cram_header_file.name)

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

	def _get_byte_range(self, start, end):
		self._total_bytes_loaded_from_cram += end - start

		if self._byte_ranges_cache is not None and (start, end) in self._byte_ranges_cache:
			return self._byte_ranges_cache[(start, end)]

		if self._is_file_in_google_storage:
			byte_range = get_byte_range_from_google_storage(self._cram_or_bam_path, start, end)
		else:
			self._cram_or_bam_file.seek(start)
			byte_range = self._cram_or_bam_file.read(end - start)

		self._total_byte_ranges_loaded_from_cram += 1
		self._total_bytes_loaded_from_cram += len(byte_range)

		if len(byte_range) != end - start:
			raise ValueError(f"Expected to read {end - start} bytes (from {start} to {end}) but read {len(byte_range)}")

		if self._byte_ranges_cache is not None:
			self._byte_ranges_cache[(start, end)] = byte_range

		return byte_range


