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
import json
import os
import re
import sys
import pysam
import tempfile
import time

from google.cloud import storage

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
    parser.add_argument("-w", "--window-size", type=int, default=1000, help="Window size in bp to include around the "
                        "user-specified region(s). This is useful for including read pairs that may overlap the region(s)")
    parser.add_argument("-u", "--gcloud-project", help="Google Cloud project name to use when reading the input cram.")
    parser.add_argument("-R", "--reference-fasta", required=True, help="Reference genome FASTA file used for reading the CRAM file")
    parser.add_argument("-o", "--output-cram", help="Output file path prefix")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-L", "--region", action="append", help="Region(s) for which to extract reads (chr:start-end). "
                       "For example, for the HTT repeat locus on hg38, specify chr4:3074877-3074933")
    group.add_argument("-c", "--variant-catalog", help="Variant catalog JSON path")

    parser.add_argument("-i", "--crai-index-path", help="Optional path of the input CRAM index file. This can be a "
                        "local or a gs:// path")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--output-data-transfer-stats", action="store_true", help="Write out a TSV file with stats "
                        "about the total number of bytes and containers downloaded from the CRAM")
    parser.add_argument("input_cram", help="Input CRAM file path. This can a local or a gs:// path")

    args = parser.parse_args()

    if not args.input_cram.endswith(".cram"):
        parser.error(f"Input CRAM file must have a .cram extension: {args.input_cram}")

    if not args.output_cram:
        args.output_cram = re.sub(".cram$", "", os.path.basename(args.input_cram))
        args.output_cram += ".subset.cram"
    elif not args.output_cram.endswith(".cram"):
        parser.error(f"Output CRAM file must have a .cram extension: {args.output_cram}")

    if args.variant_catalog:
        args.region = []
        with open_file(args.variant_catalog, download_local_copy_before_opening=True) as f:
            variant_catalog_records = json.load(f)
            for i, record in enumerate(variant_catalog_records):
                if "ReferenceRegion" not in record:
                    parser.error(f"Record #{i+1} does not have a ReferenceRegion field: {record}")

                reference_regions = record["ReferenceRegion"]
                if not isinstance(reference_regions, list):
                    reference_regions = [reference_regions]
                for region in reference_regions:
                    if args.debug:
                        print(f"DEBUG: Adding", record["LocusId"], "region", region)
                    args.region.append(region)

    if args.debug:
        args.verbose=True

    start_time = time.time()

    # validate args
    intervals = []
    for region in args.region:
        try:
            chrom, start, end = parse_interval(region)
        except ValueError as e:
            parser.error(f"Unable to parse region {region}: {e}")

        window_start = start - args.window_size
        window_end = end + args.window_size
        intervals.append((chrom, window_start, window_end))

    if not args.input_cram.endswith(".cram"):
        parser.error(f"Input {args.input_cram} must be a .cram file")

    set_requester_pays_project(args.gcloud_project)

    input_crai_path = args.crai_index_path if args.crai_index_path else f"{args.input_cram}.crai"
    if not file_exists(input_crai_path):
        parser.error(f"CRAM index path {input_crai_path} not found")

    if not file_exists(args.reference_fasta):
        parser.error(f"Reference file not found {args.reference_fasta}")

    # create a CramIntervalRreader and use it to generate a temp CRAM file containing the CRAM header and any reads
    # overlapping the user-specified region interval(s)
    print(f"Retrieving reads within {args.window_size:,d}bp of", ", ".join(args.region))
    cram_reader = IntervalReader(args.input_cram, input_crai_path, verbose=args.verbose, debug=args.debug,
                                 reference_fasta_path=args.reference_fasta, cache_byte_ranges=True)

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
            chrom, start, end,
            input_bam=input_bam_file,
            bamlet=None,
            merge_regions_distance=args.merge_regions_distance,
            verbose=args.verbose)

        for genomic_region in genomic_regions:
            cram_reader.add_interval(*genomic_region)

    print(f"Exporting data for {len(intervals)} intervals to {args.output_cram}")
    cram_reader.save_to_file(args.output_cram)

    total_bytes = cram_reader.get_total_bytes_loaded_from_cram()
    total_containers = cram_reader.get_total_byte_ranges_loaded_from_cram()
    total_duration_seconds = time.time() - start_time
    print(f"Downloaded {total_containers:,d} containers, {total_bytes/10**6:0,.1f}Mb in {round(total_duration_seconds, 2)} seconds")
    if args.output_data_transfer_stats:
        cram_reader.save_data_transfer_stats()

    input_bam_file.close()

if __name__ == "__main__":
    main()
