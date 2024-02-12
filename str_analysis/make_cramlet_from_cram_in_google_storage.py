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
import time

from google.cloud import storage
import hailtop.fs as hfs

from str_analysis.make_bamlet import extract_region
from str_analysis.utils.cram_bam_utils import IntervalReader
from str_analysis.utils.file_utils import set_requester_pays_project, file_exists, open_file, get_file_size
from str_analysis.utils.misc_utils import parse_interval

pysam.set_verbosity(0)


def main():
    parser = argparse.ArgumentParser(description="A script to generate BAMlets")
    parser.add_argument("-d", "--merge-regions-distance", type=int, default=1000, help="Region merge distance. "
        "When retrieving mates, regions that are within this distance of each other will be merged "
        "and retrieved using a single disk read operation. To reduce number of the disk reads, increase "
        "this parameter, or decrease it to reduce the total number of bytes read.")
    parser.add_argument("-u", "--gcloud-project", help="Google Cloud project name to use when reading the input cram.")
    parser.add_argument("-R", "--reference-fasta", help="Reference genome FASTA file used for reading the CRAM file")
    parser.add_argument("-o", "--cramlet", help="Output file path prefix")
    parser.add_argument("-i", "--crai-index-path", help="Optional path of the input CRAM index file. This can be a "
                        "local or a gs:// path")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("-t", "--output-download-stats", action="store_true", help="Write out a TSV file with stats "
                        "about the total number of bytes and containers downloaded from the CRAM")
    parser.add_argument("input_cram", help="Input CRAM file path. This can a local or a gs:// path")
    parser.add_argument("region", nargs="+", help="Region(s) for which to extract reads (chr:start-end). For example, "
                                                  "for the HTT repeat locus on hg38, specify chr4:3074877-3074933")

    args = parser.parse_args()

    window_size = 2000

    start_time = time.time()

    # validate args
    intervals = []
    for region in args.region:
        try:
            chrom, start, end = parse_interval(region)
        except ValueError as e:
            parser.error(f"Unable to parse region {region}: {e}")

        window_start = start - window_size
        window_end = end + window_size
        intervals.append((chrom, window_start, window_end))

    if not args.input_cram.endswith(".cram"):
        parser.error(f"Input {args.input_cram} must be a .cram file")

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
    print(f"Retrieving reads within {window_size:,d}bp of", ", ".join(args.region))
    cram_reader = IntervalReader(args.input_cram, input_crai_path, verbose=args.verbose,
                                 reference_fasta_path=args.reference_fasta,
                                 retrieve_cram_containers_from_google_storage=True,
                                 cache_byte_ranges=True)

    for chrom, start, end in intervals:
        cram_reader.add_interval(chrom, start, end)

    temporary_cram_file = tempfile.NamedTemporaryFile(suffix=".cram", delete=False)
    cram_reader.save_to_file(temporary_cram_file.name)
    temporary_cram_file.seek(0)

    # parse the temp CRAM file and get byte ranges for mates
    input_bam_file = pysam.AlignmentFile(
        temporary_cram_file.name, reference_filename=args.reference_fasta)

    for chrom, start, end in intervals:
        if args.verbose and len(intervals) > 1:
            print("-"*100)

        genomic_regions = extract_region(
            chrom, window_start, window_end,
            input_bam=input_bam_file,
            bamlet=None,
            merge_regions_distance=args.merge_regions_distance,
            verbose=args.verbose)

        for genomic_region in genomic_regions:
            cram_reader.add_interval(*genomic_region)

    print(f"Exporting data for {len(intervals)} intervals to {args.cramlet}")
    cram_reader.save_to_file(args.cramlet)

    total_duration_seconds = time.time() - start_time
    if args.output_download_stats:
        total_containers = cram_reader.get_total_containers_downloaded_from_cram()
        total_bytes = cram_reader.get_total_bytes_downloaded_from_cram()
        stats_tsv_path = re.sub("(.cramlet)?.cram$", "", args.cramlet) + ".stats.tsv"
        add_header = not os.path.isfile(stats_tsv_path) or os.path.getsize(stats_tsv_path) == 0

        with open(stats_tsv_path, "a") as stats_file:
            if add_header:
                stats_file.write("\t".join([
                    "input_cram",
                    "total_containers_downloaded",
                    "total_bytes_downloaded",
                    "total_duration_seconds",
                ]) + "\n")
            stats_file.write("\t".join(map(str, [
                args.input_cram,
                total_containers,
                total_bytes,
                round(total_duration_seconds, 2),
            ])) + "\n")

        if args.verbose:
            print(f"Wrote stats to {stats_tsv_path}: Downloaded {total_containers:,d} containers, "
                  f"{total_bytes:,d} bytes in {round(total_duration_seconds, 2)} seconds")
            
    input_bam_file.close()

if __name__ == "__main__":
    main()