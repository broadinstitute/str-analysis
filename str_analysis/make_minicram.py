"""Takes a CRAM file and one or more genomic intervals, and outputs overlapping reads in SAM or CRAM format."""

import argparse
import binascii
import collections
import intervaltree
import gzip
import hail as hl
import os
import re
import pysam
import tempfile

from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.file_utils import file_exists, open_file

CRAM_EOF_CONTAINER = binascii.unhexlify("0f000000ffffffff0fe0454f4600000000010005bdd94f0001000606010001000100ee63014b")

CRAI_FILE_HEADER = [
	"reference_sequence_id", # reference sequence identifier, or -1 for unmapped reads, -2 for multiple reference sequences. All slices in this container must have a reference sequence id matching this value.
	"alignment_start",  # ignored for unmapped slices
	"alignment_span",   # ignored for unmapped slices
	"absolute_container_header_byte_offset",
	"relative_slice_header_byte_offset",
	"slice_size_in_bytes",
]


def main():
	p = argparse.ArgumentParser(description="Convert Tandem Repeat Finder output .dat format to .bed")
	#p.add_argument("-R", "--reference-path", required=True, help="Reference FASTA path")
	p.add_argument("-L", "--interval", action="append", help="Interval filter (eg. 'chr1:12345-54321'). "
															 "Can be specified multiple times.")
	p.add_argument("-o", "--output-path", required=True, help="Output CRAM file path.")
	p.add_argument("cram_path", help="Input CRAM file path. This can be a local path or a gs:// path.")
	args = p.parse_args()

	if not args.interval:
		p.error("At least one interval must be specified via -L")

	# parse intervals
	user_intervals = collections.defaultdict(intervaltree.IntervalTree)
	for i in args.interval:
		chrom, start, end = parse_interval(i)
		chrom = chrom.replace("chr", "").upper()
		user_intervals[chrom].add(intervaltree.Interval(start, end))

	if not file_exists(args.cram_path):
		p.error(f"{args.cram_path} not found")

	crai_path = f"{args.cram_path}.crai"
	if not file_exists(crai_path):
		crai_path2 = re.sub(".cram$", ".crai", args.cram_path)
		if not file_exists(crai_path2):
			p.error(f"Index file not found at {crai_path}.crai or {crai_path2}")
		crai_path = crai_path2

	container_byte_offsets = set()  # collect unique offsets
	with open_file(crai_path, gunzip=True) as crai_file:
		for i, line in enumerate(crai_file):
			fields = line.strip().split("\t")
			if len(fields) != len(CRAI_FILE_HEADER):
				p.error(f"Expected {len(CRAI_FILE_HEADER)} columns but found {len(fields)} in line #{i} of {crai_index_path}: {line}")

			crai_record = dict(zip(CRAI_FILE_HEADER, fields))
			if int(crai_record["reference_sequence_id"]) < 0:
				continue
			if float(crai_record['alignment_span']) < 0:
				continue

			container_byte_offsets.add(int(crai_record["absolute_container_header_byte_offset"]))

	container_sizes = {}  # maps container byte offset to container length
	previous_offset = None
	for container_byte_offset in list(sorted(container_byte_offsets)):
		if previous_offset is not None:
			container_sizes[previous_offset] = container_byte_offset - previous_offset
		previous_offset = container_byte_offset
		# TODO add last container

	containers_written_counter = 0
	bytes_written_counter = 0
	with (
		open_file(args.cram_path) as cram_file,
		open_file(crai_path, gunzip=True) as crai_file,
		open(args.output_path, "wb") as output_file,
	):
		for i, line in enumerate(crai_file):
			fields = line.strip().split("\t")
			if len(fields) != len(CRAI_FILE_HEADER):
				p.error(f"Expected {len(CRAI_FILE_HEADER)} columns but found {len(fields)} in line #{i} of {crai_index_path}: {line}")

			crai_record = dict(zip(CRAI_FILE_HEADER, fields))
			if int(crai_record["reference_sequence_id"]) < 0:
				continue

			if float(crai_record['alignment_span']) < 0:
				continue

			if i == 0:
				end_of_cram_header_byte_offset = int(crai_record["absolute_container_header_byte_offset"])
				bytes_before_first_data_container = cram_file.read(end_of_cram_header_byte_offset)

				# write a temp file with just the cram header
				with tempfile.NamedTemporaryFile() as temp_file:
					temp_file.write(bytes_before_first_data_container)
					temp_file.write(CRAM_EOF_CONTAINER)
					temp_file.flush()
					temp_file.seek(0)

					verbosity = pysam.set_verbosity(0)
					cram_references = pysam.AlignmentFile(temp_file.name, mode="rb", check_sq=False, require_index=False).references
					verbosity = pysam.set_verbosity(verbosity)

				output_file.write(bytes_before_first_data_container)

			crai_record_chrom = cram_references[int(crai_record["reference_sequence_id"])].replace("chr", "").upper()
			#if crai_record_chrom not in intervals:
			#	continue

			crai_record_start = int(crai_record["alignment_start"])
			crai_record_end = crai_record_start + int(crai_record["alignment_span"])
			if not user_intervals[crai_record_chrom].overlaps(crai_record_start, crai_record_end):
				continue

			print(f"Found overlapping interval for {crai_record_chrom}:{crai_record_start:,d}-{crai_record_end:,d}")
			absolute_container_header_byte_offset = int(crai_record["absolute_container_header_byte_offset"])
			cram_file.seek(absolute_container_header_byte_offset)

			container_size = container_sizes[absolute_container_header_byte_offset]
			print("Writing container of size", container_size)
			output_file.write(cram_file.read(container_size))

			containers_written_counter += 1
			bytes_written_counter += container_size

		output_file.write(CRAM_EOF_CONTAINER)

	print(f"Wrote {containers_written_counter:,d} containers to {args.output_path}")
	print(f"Wrote {bytes_written_counter:,d} bytes to {args.output_path}")
	print("Done")


if __name__ == "__main__":
	main()
