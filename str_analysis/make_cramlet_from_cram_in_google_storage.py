"""
This is a modified version of the https://github.com/Illumina/ExpansionHunterDenovo/blob/master/scripts/make-bamlet.py
script from ExpansionHunterDenovo. It takes a CRAM file as input, and extracts all relevant reads needed to genotype
the given locus using ExpansionHunter. It writes these reads to a much smaller CRAMlet file which can then be given to
ExpansionHunter and should yield the same genotype as the full BAM or CRAM.

It minimizes the number of bytes read from
the input file. This is useful for CRAM files stored in Nearline or other cloud storage types where the cost is
proportional to the number of bytes read (as is currently the case for AllOfUs CRAMS).

For a given STR region (for example the HTT repeat @ chr4:3074877-3074933), this script will
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

from str_analysis.make_bamlet import extract_region
from str_analysis.utils.file_utils import open_file, get_file_size
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


def get_byte_range(client, google_storage_path, start_bytes, end_bytes, gcloud_project):
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

def parse_crai_index(crai_index_path, cram_path):
	container_byte_offsets = set()  # collect unique offsets
	end_of_cram_header_byte_offset = None
	interval_trees = {}

	with open_file(crai_index_path, gunzip=True) as crai_file:
		for i, line in enumerate(crai_file):
			fields = line.strip().split("\t")
			if len(fields) != len(CRAI_FILE_HEADER):
				p.error(f"Expected {len(CRAI_FILE_HEADER)} columns but found {len(fields)} in line #{i} of {crai_index_path}: {line}")

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



class GoogleStorageCramReader:

	def __init__(self, local_cram_header_path, cram_path, crai_path, gcloud_project):
		self.total_cram_containers_loaded_from_google_storage = 0
		self.total_bytes_read_from_google_storage = 0

		self._local_cram_header_path = local_cram_header_path
		self._cram_path = cram_path
		self._crai_path = crai_path
		self._gcloud_project = gcloud_project
		self._client = storage.Client()

		self._cram_data_containers_cache = {}

		# parse the CRAM index
		self._end_of_cram_header_byte_offset, self._crai_interval_trees = parse_crai_index(
			self._crai_path, self._cram_path)

		if len(self._crai_interval_trees) == 0:
			raise ValueError(f"No CRAI records loaded from {self._crai_path}")

		if local_cram_header_path:
			# load the cram header from the local file
			with open(local_cram_header_path, "rb") as f:
				self._cram_header_bytes = f.read()
		else:
			# get the header from the input cram
			self._cram_header_bytes = self._get_byte_range(0, self._end_of_cram_header_byte_offset)

		self._init_chrom_index_lookup()

	def load_cram_containers_for_interval(self, chrom, start, end, max_containers_to_load=None, dry_run=False):
		reference_sequence_id = self._chrom_index_lookup[normalize_chromosome_name(chrom)]
		if reference_sequence_id not in self._crai_interval_trees:
			print(f"WARNING: No CRAM containers found for {chrom} (reference_sequence_id={reference_sequence_id})")
			return

		interval_tree = self._crai_interval_trees[reference_sequence_id]
		overlapping_intervals = list(interval_tree.overlap(start, end))
		if not overlapping_intervals:
			print(f"WARNING: None of the {len(interval_tree)} CRAM containers on {chrom} overlap {chrom}:{start}-{end}")
			return

		print(f"Found {len(overlapping_intervals):4,d} CRAM containers that overlap {chrom}:{start}-{end}")
		containers_loaded_counter = 0
		for interval in overlapping_intervals:
			bytes_start, bytes_end = interval.data.start, interval.data.end
			if bytes_start in self._cram_data_containers_cache:
				cram_container_bytes = self._cram_data_containers_cache[bytes_start]
			else:
				if not dry_run:
					print(f"Loading CRAM container with {bytes_end - bytes_start:,d} bytes for {chrom}:{start}-{end}")
				cram_container_bytes = self._get_byte_range(bytes_start, bytes_end, dry_run=dry_run)
				self._cram_data_containers_cache[bytes_start] = cram_container_bytes
				self.total_cram_containers_loaded_from_google_storage += 1
				
				containers_loaded_counter += 1
				if max_containers_to_load is not None and containers_loaded_counter >= max_containers_to_load:
					print(f"WARNING: Only loaded {max_containers_to_load} out of {len(overlapping_intervals)} overlapping CRAM containers on {chrom} for {chrom}:{start}-{end}")
					break
			#yield cram_container_bytes  # WARNING: this causes lazy loading

	def save_cram_to_file(self, local_path=None, fileobj=None):
		if (local_path is None and fileobj is None) or (local_path is not None and fileobj is not None):
			raise ValueError("Must specify either local_path or fileobj")

		if local_path:
			fileobj = open(local_path, "wb")

		fileobj.write(self._cram_header_bytes)
		for _, cram_container_bytes in sorted(self._cram_data_containers_cache.items()):
			fileobj.write(cram_container_bytes)

		fileobj.write(CRAM_EOF_CONTAINER)
		if local_path:
			fileobj.close()

	def _init_chrom_index_lookup(self):

		# write a temp file with just the cram header
		with tempfile.NamedTemporaryFile(suffix=".cram") as temp_file:
			self.save_cram_to_file(fileobj=temp_file)
			temp_file.flush()
			temp_file.seek(0)
			with pysam.AlignmentFile(temp_file.name, mode="rb", check_sq=False, require_index=False) as file:
				references = file.references

		self._chrom_index_lookup = {
			normalize_chromosome_name(name): idx for idx, name in enumerate(references)
		}

		if len(self._chrom_index_lookup) == 0:
			raise ValueError(f"Unable to load chromosome list from the header of {self._cram_path}")

	def _get_byte_range(self, start, end, dry_run=False):
		if dry_run:
			byte_range = b""
		else:
			byte_range = get_byte_range(self._client, self._cram_path, start, end, self._gcloud_project)
			if len(byte_range) != end - start:
				raise ValueError(f"Expected to read {end - start} bytes (from {start} to {end}) but read {len(byte_range)}")

		self.total_bytes_read_from_google_storage += end - start

		return byte_range


def main():
	parser = argparse.ArgumentParser(description="A script to generate BAMlets")
	parser.add_argument("-d", "--merge-regions-distance", type=int, default=1000, help="Region merge distance. "
		"When retrieving mates, regions that are within this distance of each other will be merged "
		"and retrieved using a single disk read operation. To reduce number of the disk reads, increase "
		"this parameter, or decrease it to reduce the total number of bytes read.")
	parser.add_argument("-u", "--gcloud-project", required=True, help="Google Cloud project name to use when reading the input cram.")
	parser.add_argument("--cram-header", help="Optional path of a local CRAM file that only contains the CRAM header to "
						"use for the output cramlet. Providing this further reduces the number of bytes that need to be "
						"transferred from Google Storage by making it unnecessary to read the header from the input CRAM")
	parser.add_argument("-R", "--reference-fasta", required=True, help="Reference genome FASTA file used for reading the CRAM file")
	parser.add_argument("-o", "--cramlet", help="Output file path prefix")
	parser.add_argument("-i", "--crai-index-path", help="Optional path of the input CRAM index file. This can be a local or a gs:// path")
	parser.add_argument("-e", "--save-cram-header-and-exit", action="store_true", help="Save the CRAM header from the "
					    "given input cram and exit. This output file can then be used as the argument to --cram-header.")
	parser.add_argument("--dry-run", action="store_true", help="Only compute stats for mate regions without actually loading them")
	parser.add_argument("--max-containers-to-load", type=int, help="Maximum number of CRAM containers to load from Google Storage")
	parser.add_argument("input_cram", help="Input CRAM file gs:// path")
	parser.add_argument("region", nargs="+", help="Region(s) for which to extract reads (chr:start-end). For example, "
												  "for the HTT repeat locus on hg38, specify chr4:3074877-3074933")

	args = parser.parse_args()

	# validate args
	if not args.input_cram.endswith(".cram"):
		parser.error(f"Input {args.input_cram} must be a .cram file")
	if not args.input_cram.startswith("gs://"):
		parser.error(f"CRAM path {args.input_cram} must be on Google Storage (ie. start with 'gs://')")

	crai_index_path = args.crai_index_path if args.crai_index_path else f"{args.input_cram}.crai"
	fexists = hfs.exists if crai_index_path.startswith("gs://") else os.path.isfile
	if not fexists(crai_index_path):
		parser.error(f"CRAM index path {crai_index_path} not found")

	for path in args.reference_fasta, args.cram_header:
		if path and not os.path.isfile(path):
			parser.error(f"File not found {path}")

	if not args.cramlet:
		args.cramlet = re.sub(".cram$", "", os.path.basename(args.input_cram))
		args.cramlet += ".cramlet.cram" if not args.save_cram_header_and_exit else ".header.cram"

	# initialize the CRAM reader
	cram_reader = GoogleStorageCramReader(args.cram_header, args.input_cram, crai_index_path, args.gcloud_project)

	# write a temp file with just the cram header
	if args.save_cram_header_and_exit:
		cram_reader.save_cram_to_file(args.cramlet)
		return

	# write out a temp CRAM file with just the cram header and region intervals
	#temporary_cram_file = tempfile.NamedTemporaryFile(suffix=".cram")
	temporary_cram_file = open("temp_file.cram", "wb")
	for region in args.region:
		chrom, start, end = parse_interval(region)
		window_start = start - 2000
		window_end = end + 2000

		cram_reader.load_cram_containers_for_interval(chrom, window_start, window_end)

	cram_reader.save_cram_to_file(fileobj=temporary_cram_file)
	temporary_cram_file.flush()
	pysam.index(temporary_cram_file.name)

	# parse the temp CRAM file and get byte ranges for mates
	temporary_cram_file.seek(0)
	input_bam_file = pysam.AlignmentFile(temporary_cram_file.name, "rc", reference_filename=args.reference_fasta)

	for region in args.region:
		print("-"*100)
		chrom, start, end = parse_interval(region)
		window_start = start - 2000
		window_end = end + 2000

		genomic_regions = extract_region(
			chrom, window_start, window_end,
			input_bam=input_bam_file,
			bamlet=None,
			merge_regions_distance=args.merge_regions_distance)

		for genomic_region in genomic_regions:
			cram_reader.load_cram_containers_for_interval(
				*genomic_region, max_containers_to_load=args.max_containers_to_load, dry_run=args.dry_run)

	temporary_cram_file.close()

	cram_reader.save_cram_to_file(args.cramlet)

	total_megabytes_read = cram_reader.total_bytes_read_from_google_storage/10**6
	print(f"Read {cram_reader.total_cram_containers_loaded_from_google_storage:,d} containers ({total_megabytes_read:0.1f}Mb) from Google Storage, with {total_megabytes_read/cram_reader.total_cram_containers_loaded_from_google_storage:0.2f}Mb/container")


if __name__ == "__main__":
	main()