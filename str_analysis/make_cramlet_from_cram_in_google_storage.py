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
from str_analysis.utils.cram_bam_utils import CramIntervalReader
from str_analysis.utils.file_utils import set_requester_pays_project, file_exists, open_file, get_file_size
from str_analysis.utils.misc_utils import parse_interval

pysam.set_verbosity(0)


def main():
    parser = argparse.ArgumentParser(description="A script to generate BAMlets")
    parser.add_argument("-d", "--merge-regions-distance", type=int, default=1000, help="Region merge distance. "
        "When retrieving mates, regions that are within this distance of each other will be merged "
        "and retrieved using a single disk read operation. To reduce number of the disk reads, increase "
        "this parameter, or decrease it to reduce the total number of bytes read.")
    parser.add_argument("-u", "--gcloud-project", required=True, help="Google Cloud project name to use when reading the input cram.")
    parser.add_argument("-R", "--reference-fasta", required=True, help="Reference genome FASTA file used for reading the CRAM file")
    parser.add_argument("-o", "--cramlet", help="Output file path prefix")
    parser.add_argument("-i", "--crai-index-path", help="Optional path of the input CRAM index file. This can be a local or a gs:// path")
    parser.add_argument("--dry-run", action="store_true", help="Only compute stats for mate regions without actually loading them")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("input_cram", help="Input CRAM file gs:// path")
    parser.add_argument("region", nargs="+", help="Region(s) for which to extract reads (chr:start-end). For example, "
                                                  "for the HTT repeat locus on hg38, specify chr4:3074877-3074933")

    args = parser.parse_args()

    # validate args
    if not args.input_cram.endswith(".cram"):
        parser.error(f"Input {args.input_cram} must be a .cram file")
    if not args.input_cram.startswith("gs://"):
        parser.error(f"CRAM path {args.input_cram} must be on Google Storage (ie. start with 'gs://'). For local BAM " 
                     f"or CRAM files, use the make_bamlet.py script.")

    set_requester_pays_project(args.gcloud_project)

    input_crai_path = args.crai_index_path if args.crai_index_path else f"{args.input_cram}.crai"
    if not file_exists(input_crai_path):
        parser.error(f"CRAM index path {input_crai_path} not found")

    if not file_exists(args.reference_fasta):
        parser.error(f"Reference file not found {path}")

    if not args.cramlet:
        args.cramlet = re.sub(".cram$", "", os.path.basename(args.input_cram))
        args.cramlet += ".cramlet.cram"

    # create a CramIntervalRreader and use it to generate a temp CRAM file containing the CRAM header and any reads
    # overlapping the user-specified region interval(s)
    window_size = 2000
    print(f"Retrieving reads that are within {window_size:,d}bp of", ", ".join(args.region))
    cram_reader = CramIntervalReader(args.input_cram, input_crai_path, verbose=args.verbose, cache_byte_ranges=True)

    temporary_cram_file = tempfile.NamedTemporaryFile(suffix=".cram")
    #temporary_cram_file = open("temp_file.cram", "wb")
    for region in args.region:
        chrom, start, end = parse_interval(region)
        window_start = start - window_size
        window_end = end + window_size

        cram_reader.add_interval(chrom, window_start, window_end)

    cram_reader.save_to_file(fileobj=temporary_cram_file)
    temporary_cram_file.flush()

    # parse the temp CRAM file and get byte ranges for mates
    print("Retrieving regions that contain far-away mates.")
    temporary_cram_file.seek(0)
    input_bam_file = pysam.AlignmentFile(temporary_cram_file.name, "rc", reference_filename=args.reference_fasta)

    for region in args.region:
        if args.verbose:
            print("-"*100)
        chrom, start, end = parse_interval(region)
        window_start = start - window_size
        window_end = end + window_size

        genomic_regions = extract_region(
            chrom, window_start, window_end,
            input_bam=input_bam_file,
            bamlet=None,
            merge_regions_distance=args.merge_regions_distance,
            verbose=args.verbose)

        for genomic_region in genomic_regions:
            cram_reader.add_interval(*genomic_region)

    temporary_cram_file.close()

    cram_reader.save_to_file(args.cramlet)

if __name__ == "__main__":
    main()