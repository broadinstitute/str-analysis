import binascii
import intervaltree
import pysam
import re
import tempfile
from google.cloud import storage
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
	The second value as a dictionary of interval trees that, for each reference sequence id in the CRAM file,
	provides a way to quickly look up the CRAM file byte range for a given genomic interval.
	"""
	container_byte_offsets = set()  # collect unique offsets
	end_of_cram_header_byte_offset = None
	interval_trees = {}

	with open_file(crai_path, gunzip=True) as crai_file:
		for i, line in enumerate(crai_file):
			fields = line.strip().split("\t")
			if len(fields) != len(CRAI_FILE_HEADER):
				p.error(f"Expected {len(CRAI_FILE_HEADER)} columns but found {len(fields)} in line #{i} of {crai_path}: {line}")

			crai_record = dict(zip(CRAI_FILE_HEADER, map(int, fields)))
			if crai_record["reference_sequence_id"] < 0:
				continue
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


class CramIntervalReader:
	"""This class implements a way to retrieve genomic intervals from CRAM files in a way that's more i/o efficient than
	htslib, in the sense that it reads only the containers that overlap the given interval(s), while htslib and
	htsjdk appear to read adjacent containers. This is important for minimizing costs when reading from Google Cloud
	Nearline storage (which costs $0.01 per gigabyte to read).
	"""

	def __init__(self, cram_path, crai_path, verbose=False, cache_byte_ranges=False):
		"""Construct a new CramIntervalReader instance.

		Args:
			cram_path: The path to the CRAM file. This can be a local or a gs:// path.
			crai_path: The path to the CRAI index file. This can be a local or a gs:// path.
			verbose: If True, will print more detailed information about the byte ranges being loaded.
			cache_byte_ranges: If True, will cache the byte ranges that are loaded from the CRAM file in memory. This
				can be useful when loading a relatively small number of intervals and then reusing the same
				CramIntervalReader instance across multiple rounds of CRAM access.
		"""
		self.total_containers_loaded_from_cram = 0
		self.total_bytes_loaded_from_cram = 0

		self._cram_path = cram_path
		self._crai_path = crai_path
		self._client = storage.Client()
		self._verbose = verbose

		self._byte_ranges_to_load = set()
		self._byte_ranges_cache = None if not cache_byte_ranges else {}

		self._cram = open(self._cram_path, "rb") if not self._cram_path.startswith("gs://") else None

		# parse the CRAM index
		self._end_of_cram_header_byte_offset, self._crai_interval_trees = parse_crai_index(
			self._crai_path, self._cram_path)

		# load the CRAM header
		self._cram_header_bytes = self._get_byte_range(0, self._end_of_cram_header_byte_offset)

		self._init_chrom_index_lookup()

	def add_interval(self, chrom, start, end):
		"""Add a genomic interval to the list of intervals to be loaded from the CRAM file."""

		# use the CRAI index to convert this genomic interval to a byte range
		reference_sequence_id = self._chrom_index_lookup[normalize_chromosome_name(chrom)]
		if reference_sequence_id not in self._crai_interval_trees:
			print(f"WARNING: No CRAM containers found for {chrom} (reference_sequence_id={reference_sequence_id})")
			return

		interval_tree = self._crai_interval_trees[reference_sequence_id]
		overlapping_intervals = list(interval_tree.overlap(start, end))
		if not overlapping_intervals:
			print(f"WARNING: None of the {len(interval_tree)} CRAM containers on {chrom} overlap {chrom}:{start}-{end}")
			return

		if self._verbose:
			print(f"Found {len(overlapping_intervals):4,d} CRAM container(s) that overlapped {chrom}:{start}-{end}")

		for interval in overlapping_intervals:
			self._byte_ranges_to_load.add((interval.data.start, interval.data.end))

	def save_to_file(self, local_path=None, fileobj=None):
		"""Load and save the added genomic intervals to a local CRAM file at the given local_path, or alternatively,
		write them into the given fileobj, which must be a file object already open for writing.
		"""

		if (local_path is None and fileobj is None) or (local_path is not None and fileobj is not None):
			raise ValueError("Must specify either local_path or fileobj")

		if local_path:
			fileobj = open(local_path, "wb")

		merged_byte_ranges = []
		for start, end in sorted(self._byte_ranges_to_load):
			if not merged_byte_ranges:
				merged_byte_ranges.append((start, end))
			else:
				previous_interval = merged_byte_ranges[-1]
				if previous_interval[1] >= start:
					merged_byte_ranges[-1] = (min(previous_interval[0], start), max(previous_interval[1], end))
				else:
					merged_byte_ranges.append((start, end))

		fileobj.write(self._cram_header_bytes)
		for bytes_start, bytes_end in merged_byte_ranges:
			if self._verbose:
				print(f"Loading CRAM container with {bytes_end - bytes_start:,d} bytes")
			cram_container_bytes = self._get_byte_range(bytes_start, bytes_end)
			self.total_containers_loaded_from_cram += 1
			self.total_bytes_loaded_from_cram += len(cram_container_bytes)

			fileobj.write(cram_container_bytes)

		fileobj.write(CRAM_EOF_CONTAINER)
		if local_path:
			fileobj.close()

		pysam.index(fileobj.name)

		if self.total_containers_loaded_from_cram > 0:
			total_megabytes_loaded = self.total_bytes_loaded_from_cram/10**6
			print(f"Retrieved {self.total_containers_loaded_from_cram:3,d} CRAM containers .. "
				  f"({total_megabytes_loaded:5.1f}Mb total @ "
				  f"{total_megabytes_loaded/self.total_containers_loaded_from_cram:5.2f}Mb/container) "
				  f"from {self._cram_path}  to  {fileobj.name}")


	def _init_chrom_index_lookup(self):

		# write the CRAM header of the input CRAM file into a local temp file
		with tempfile.NamedTemporaryFile(suffix=".cram") as temp_file:
			self.save_to_file(fileobj=temp_file)
			temp_file.flush()
			temp_file.seek(0)

			# load the CRAM header from the temp file using pysam to get the list of reference sequence names
			with pysam.AlignmentFile(temp_file.name, mode="rb", check_sq=False, require_index=False) as file:
				references = file.references

		# create a lookup table for looking up the chromosome index of a given chromosome name
		self._chrom_index_lookup = {
			normalize_chromosome_name(name): idx for idx, name in enumerate(references)
		}

		if len(self._chrom_index_lookup) == 0:
			raise ValueError(f"Unable to load chromosome list from the header of {self._cram_path}")

	def _get_byte_range(self, start, end):
		self.total_bytes_loaded_from_cram += end - start

		if self._byte_ranges_cache is not None and (start, end) in self._byte_ranges_cache:
			return self._byte_ranges_cache[(start, end)]

		if self._cram_path.startswith("gs://"):
			byte_range = get_byte_range_from_google_storage(self._cram_path, start, end)
		else:
			self._cram.seek(start)
			byte_range = self._cram.read(end - start)

		if len(byte_range) != end - start:
			raise ValueError(f"Expected to read {end - start} bytes (from {start} to {end}) but read {len(byte_range)}")

		if self._byte_ranges_cache is not None:
			self._byte_ranges_cache[(start, end)] = byte_range

		return byte_range


class BamIntervalReader:
	def __init__(self, bam_path, bai_path, verbose=False):
		self._bam_path = bam_path
		self._bai_path = bai_path
		self._client = storage.Client()
		self._verbose = verbose

		self._bam = pysam.AlignmentFile(bam_path, mode="rb", require_index=True)

		self.total_reads_loaded_from_bam = 0
		self._intervals = set()

	def add_interval(self, chrom, start, end):
		self._intervals.add((chrom, start, end))

	def save_to_file(self, local_path=None, fileobj=None):
		if (local_path is None and fileobj is None) or (local_path is not None and fileobj is not None):
			raise ValueError("Must specify either local_path or fileobj")

		if local_path:
			fileobj = pysam.AlignmentFile(local_path, mode="wb", template=self._bam)

		for chrom, start, end in sorted(self._intervals):
			for read in self._bam.fetch(chrom, start, end):
				self.total_reads_loaded_from_bam += 1
				fileobj.write(read)

		if local_path:
			fileobj.close()

		pysam.index(local_path)
		print(f"Copied {self.total_reads_loaded_from_bam:,d} reads from {self._bam_path}  to  {local_path or fileobj.name}")

