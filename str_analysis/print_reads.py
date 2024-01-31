"""This script is a lighter-weight alternative to GATK PrintReads
(https://gatk.broadinstitute.org/hc/en-us/articles/360036883571-PrintReads).
It exracts data for genomic regions of a CRAM or BAM file.
"""

import argparse
import binascii
import collections
import gzip
import intervaltree
import os
import re
import sys
import pysam
import tempfile

from google.cloud import storage
import hailtop.fs as hfs

from str_analysis.utils.file_utils import set_requester_pays_project, file_exists, open_file, get_file_size
from str_analysis.utils.misc_utils import parse_interval

pysam.set_verbosity(0)

CRAM_EOF_CONTAINER = binascii.unhexlify("0f000000ffffffff0fe0454f4600000000010005bdd94f0001000606010001000100ee63014b")

CRAI_FILE_HEADER = [
	"reference_sequence_id", # reference sequence identifier, or -1 for unmapped reads, -2 for multiple reference sequences. All slices in this container must have a reference sequence id matching this value.
	"alignment_start",  # ignored for unmapped slices
	"alignment_span",   # ignored for unmapped slices
	"absolute_container_header_byte_offset",
	"relative_slice_header_byte_offset",
	"slice_size_in_bytes",
]

class ByteRange:
	def __init__(self, start, end):
		self.start = start
		self.end = end
	def __repr__(self):
		return f"ByteRange[{self.start}, {self.end}]"
	def __str__(self):
		return self.__repr__()


def normalize_chromosome_name(chrom):
	chrom = re.sub("^chr", "", chrom).upper()
	if chrom == "MT":
		chrom = "M"
	return chrom


def get_byte_range_from_google_storage(client, google_storage_path, start_bytes, end_bytes, gcloud_project):
	"""Downloads a byte range from the bucket."""
	if not google_storage_path.startswith("gs://"):
		raise ValueError(f"Path {google_storage_path} must start with gs://")

	google_storage_path_match = re.match("^gs://([^/]+)/(.+)", google_storage_path)
	if not google_storage_path_match:
		raise ValueError(f"Path {google_storage_path} must be of the form gs://bucket/path/to/file")

	bucket_name, object_name = google_storage_path_match.groups()
	bucket = client.bucket(bucket_name, user_project=gcloud_project)
	blob = storage.Blob(object_name, bucket)
	if not blob.exists():
		raise ValueError(f"{google_storage_path} not found")

	#print(f"Downloading {google_storage_path} [{start_bytes}-{end_bytes-1}]")
	return blob.download_as_bytes(start=start_bytes, end=end_bytes-1, raw_download=True)

def parse_crai_index(crai_path, cram_path):
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

	if len(interval_trees) == 0:
		raise ValueError(f"No CRAI records loaded from {crai_path}")

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


class BamReader:
	def __init__(self, bam_path, bai_path, gcloud_project=None, verbose=False):
		self._bam_path = bam_path
		self._bai_path = bai_path
		self._gcloud_project = gcloud_project
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


class CramReader:

	def __init__(self, cram_path, crai_path, gcloud_project=None, verbose=False, cache_byte_ranges=False):
		self.total_containers_loaded_from_cram = 0
		self.total_bytes_loaded_from_cram = 0

		self._cram_path = cram_path
		self._crai_path = crai_path
		self._gcloud_project = gcloud_project
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
		if (local_path is None and fileobj is None) or (local_path is not None and fileobj is not None):
			raise ValueError("Must specify either local_path or fileobj")

		if local_path:
			fileobj = open(local_path, "wb")

		merged_byte_ranges = []
		for start, end in sorted(self._byte_ranges_to_load):
			if merged_byte_ranges and merged_byte_ranges[-1][1] >= start:
				merged_byte_ranges[-1][1] = max(merged_byte_ranges[-1][1], end)
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
			print(f"Copied {self.total_containers_loaded_from_cram:,d} containers @ {total_megabytes_loaded:0.1f}Mb total "
				  f"({total_megabytes_loaded/self.total_containers_loaded_from_cram:0.2f}Mb/container) "
				  f"from {self._cram_path}  to  {fileobj.name}")


	def _init_chrom_index_lookup(self):

		# write a temp file with just the CRAM header
		with tempfile.NamedTemporaryFile(suffix=".cram") as temp_file:
			self.save_to_file(fileobj=temp_file)
			temp_file.flush()
			temp_file.seek(0)

			# load the CRAM header using pysam to get the list of reference sequence names
			with pysam.AlignmentFile(temp_file.name, mode="rb", check_sq=False, require_index=False) as file:
				references = file.references

		# create a lookup table for looking up the index of a given chromosome name
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
			byte_range = get_byte_range_from_google_storage(self._client, self._cram_path, start, end, self._gcloud_project)
		else:
			self._cram.seek(start)
			byte_range = self._cram.read(end - start)

		if len(byte_range) != end - start:
			raise ValueError(f"Expected to read {end - start} bytes (from {start} to {end}) but read {len(byte_range)}")

		if self._byte_ranges_cache is not None:
			self._byte_ranges_cache[(start, end)] = byte_range

		return byte_range


def main():
	parser = argparse.ArgumentParser(description="A script to generate BAMlets")
	parser.add_argument("-u", "--gcloud-project", help="Google Cloud project to use for GCS requester pays buckets.")
	parser.add_argument("-o", "--output", help="Output file path. If not specified, it will be based on the input filename.")
	parser.add_argument("--read-index", help="Optional path of the input BAM or CRAM index file. This can be a local "
											 "or a gs:// path")
	parser.add_argument("-L", "--interval", action="append", help="Region(s) to extract reads. This can be a .bed file "
																  "or an interval specified as \"chr:start-end\"")
	parser.add_argument("--verbose", action="store_true")
	parser.add_argument("input_bam_or_cram", help="Input BAM or CRAM file. This can a local or a gs:// path")
	args = parser.parse_args()

	# validate args
	set_requester_pays_project(args.gcloud_project)
	for path in args.input_bam_or_cram, args.read_index:
		if path and not file_exists(path):
			parser.error(f"{path} not found")

	output_path = args.output
	if args.input_bam_or_cram.endswith(".bam"):
		if args.read_index is None:
			args.read_index = f"{args.input_bam_or_cram}.bai"
		if output_path is None:
			output_path = re.sub(".bam$", "", os.path.basename(args.input_bam_or_cram))
			output_path += ".print_reads.bam"
		elif not output_path.endswith(".bam"):
			parser.error(f"Output path {output_path} must end with .bam")

		reader = BamReader(args.input_bam_or_cram, args.read_index, args.gcloud_project, verbose=args.verbose)

	elif args.input_bam_or_cram.endswith(".cram"):
		if args.read_index is None:
			args.read_index = f"{args.input_bam_or_cram}.crai"
		if output_path is None:
			output_path = re.sub(".cram$", "", os.path.basename(args.input_bam_or_cram))
			output_path += ".print_reads.cram"
		elif not output_path.endswith(".cram"):
			parser.error(f"Output path {output_path} must end with .cram")

		reader = CramReader(args.input_bam_or_cram, args.read_index, args.gcloud_project, verbose=args.verbose)
	else:
		parser.error(f"Input {args.input_bam_or_cram} must be a .bam or a .cram file")

	# write out a temp CRAM file with just the cram header and region intervals
	for interval in args.interval:
		chrom, start, end = parse_interval(interval)
		reader.add_interval(chrom, start, end)

	reader.save_to_file(output_path)

if __name__ == "__main__":
	main()