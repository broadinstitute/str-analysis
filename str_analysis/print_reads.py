"""This script is a lighter-weight alternative to GATK PrintReads
(https://gatk.broadinstitute.org/hc/en-us/articles/360036883571-PrintReads).
It exracts data for genomic regions of a CRAM or BAM file.
"""

import argparse
import binascii
import intervaltree
import os
import re
import pysam

from google.cloud import storage

from str_analysis.utils.file_utils import set_requester_pays_project, file_exists
from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.cram_bam_utils import IntervalReader
pysam.set_verbosity(0)


def main():
	parser = argparse.ArgumentParser(description="A script to generate BAMlets")
	parser.add_argument("-u", "--gcloud-project", help="Google Cloud project to use for GCS requester pays buckets.")
	parser.add_argument("-o", "--output", help="Output file path. If not specified, it will be based on the input filename.")
	parser.add_argument("--read-index", help="Optional path of the input BAM or CRAM index file. This can be a local "
											 "or a gs:// path")
	parser.add_argument("-L", "--interval", action="append", required=True, help="Region(s) to extract reads. This can "
						"be a .bed file or an interval specified as \"chr:start-end\"")
	parser.add_argument("--verbose", action="store_true")
	parser.add_argument("input_bam_or_cram", help="Input BAM or CRAM file. This can a local or a gs:// path")
	args = parser.parse_args()

	# validate args
	set_requester_pays_project(args.gcloud_project)
	for path in args.input_bam_or_cram, args.read_index:
		if path and not file_exists(path):
			parser.error(f"{path} not found")

	if args.output is None:
		if args.input_bam_or_cram.endswith(".cram"):
			args.output = re.sub(".cram$", "", os.path.basename(args.input_bam_or_cram))
			args.output += ".print_reads.cram"
		else:
			args.output = re.sub(".bam$", "", os.path.basename(args.input_bam_or_cram))
			args.output += ".print_reads.bam"

	reader = IntervalReader(
		args.input_bam_or_cram,
		crai_or_bai_path=args.read_index,
		retrieve_cram_containers_from_google_storage=True,
		verbose=args.verbose)

	# write out a temp CRAM file with just the cram header and region intervals
	for interval in args.interval:
		try:
			chrom, start, end = parse_interval(interval)
			reader.add_interval(chrom, start, end)
		except ValueError:
			parser.error(f"Invalid interval {interval}")

	read_counts = reader.save_to_file(args.output)

if __name__ == "__main__":
	main()