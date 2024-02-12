"""
This is an modified version of the https://github.com/Illumina/ExpansionHunterDenovo/blob/master/scripts/make-bamlet.py
script from ExpansionHunterDenovo. It has been optimized to reduce the total number of i/o read operations.

For a given STR region (for example the HTT repeat @ chr4:3074877-3074933), this script will extract all
relevant reads from the input BAM or CRAM file into a much smaller BAMlet which can then be used as the input
to ExpansionHunter instead of the full BAM or CRAM but yield the same genotype.
"""

import argparse
import binascii
import collections
import gzip
import sys
import os
import pandas as pd
import pysam
import re
 
from str_analysis.utils.misc_utils import parse_interval

def is_close(chrom, pos, region, max_dist=1000):
    reg_chrom, reg_start, reg_end = region

    if chrom != reg_chrom:
        return False

    if reg_start < pos < reg_end:
        return True

    dist = min(abs(pos - reg_start), abs(pos - reg_end))
    if dist > max_dist:
        return False

    return True


def jump_for_mates(bam, chrom, start, end, read_names_set):
    read_pairs = collections.defaultdict(list)
    sys.stdout.write(f"Jumping to {chrom}:{start}-{end} ({end-start}bp window) to retrieve {len(read_names_set)} mates.. ")
    alignment_counter = 0
    for alignment in bam.fetch(chrom, start, end + 1):
        if alignment.is_secondary:
            continue
        alignment_counter += 1
        read_name = alignment.query_name
        if read_name in read_names_set and alignment.reference_start >= start and alignment.reference_start <= end + 1:
            if len(read_pairs[read_name]) > 0:
                print(f"[WARNING: Multiple reads found for {read_name}]")
                continue
            read_pairs[read_name].append(alignment)

        if len(read_pairs) == len(read_names_set):
            print(f"found all {len(read_pairs)} read names after processing {alignment_counter} reads")
            break
    else:
        print(f"[WARNING: Could not locate {len(read_names_set) - len(read_pairs)} read names: {set(read_names_set) - set(read_pairs)}]")

    return [alignment for read_pair_list in read_pairs.values() for alignment in read_pair_list]


def extract_region(chrom, start, end, input_bam, bamlet, merge_regions_distance=1000, verbose=False):
    genomic_regions_to_fetch = [
        (chrom, start, end)
    ]

    # cache all read pairs that overlap the region
    read_pairs = collections.defaultdict(list)
    for alignment in input_bam.fetch(chrom, start, end):
        if alignment.is_secondary or alignment.is_supplementary:
            # skip secondary and supplementary alignments as in
            # https://github.com/bw2/ExpansionHunter/blob/master/ehunter/sample/MateExtractor.cpp#L143
            continue

        read_pairs[alignment.query_name].append(alignment)

    if verbose and len(read_pairs) > 0:
        print(f"Extracted {len(read_pairs)} read pairs cotaining {sum(len(read_pair) for read_pair in read_pairs.values())} reads from region {chrom}:{start}-{end}")

    # compute a dictionary that maps (chrom, start, end) to a set of read names that need to be fetched from that region
    mate_regions = collections.defaultdict(set)
    for read_name, read_pair in read_pairs.items():
        if len(read_pair) >= 2:
            continue

        # see if mate is close to other mates that need to be fetched
        mate_chrom = read_pair[0].next_reference_name
        mate_pos = int(read_pair[0].next_reference_start)

        #if mate_chrom == chrom and min(abs(mate_pos - start), abs(mate_pos - end)) <= 1000:
        #    # skip mates that are close to the ends of the region
        #    continue

        for mate_region in mate_regions:
            if is_close(mate_chrom, mate_pos, mate_region, max_dist=merge_regions_distance):
                previous_read_names_set = mate_regions.pop(mate_region)
                previous_read_names_set.add(read_name)
                key = (mate_chrom, min(mate_pos - 1, mate_region[1]), max(mate_pos + 1, mate_region[2]))
                mate_regions[key] = previous_read_names_set
                break
        else:
            key = mate_chrom, mate_pos - 1, mate_pos + 1
            mate_regions[key].add(read_name)

    genomic_regions_to_fetch.extend(mate_regions.keys())

    print(f"{chrom}:{start}-{end}: Need to fetch {len(read_pairs)} mates from {len(mate_regions)} regions")

    if bamlet is not None:
        for (mate_chrom, mate_start, mate_end), read_names_set in mate_regions.items():
            for alignment in jump_for_mates(input_bam, mate_chrom, mate_start, mate_end, read_names_set):
                read_pairs[alignment.query_name].append(alignment)

        read_counter = 0
        for read_name, read_pair in read_pairs.items():
            for read in read_pair:
                read_counter += 1
                bamlet.write(read)
        print(f"Wrote {read_counter:,d} reads to {getattr(bamlet, 'filename', b'bamlet').decode()}")

    return genomic_regions_to_fetch



def main():
    parser = argparse.ArgumentParser(description="A script to generate BAMlets")
    parser.add_argument("-d", "--merge-regions-distance", type=int, default=1000, help="Region merge distance. "
                        "When retrieving mates, regions that are within this distance of each other will be merged "
                        "and retrieved using a single disk read operation. To reduce number of the disk reads, increase "
                        "this parameter, or decrease it to reduce the total number of bytes read.")
    parser.add_argument("-R", "--reference-fasta", required=True, help="Reference genome FASTA file to use when reading from CRAM")
    parser.add_argument("-o", "--bamlet", help="Output file path prefix")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("--read-index", help="Optional path of the input BAM or CRAM index file. This can be a local "
                                             "or a gs:// path")
    parser.add_argument("input_bam_or_cram", help="Input BAM or CRAM file")
    parser.add_argument("region", nargs="+", help="Region(s) for which to extract reads (chr:start-end). For example, "
                                                  "for the HTT repeat locus on hg38, specify chr4:3074877-3074933")

    args = parser.parse_args()

    if args.bamlet is None:
        args.bamlet = re.sub("(.bam|.cram)$", "", os.path.basename(args.input_bam_or_cram)) + ".bamlet.bam"


    input_bam_file = pysam.AlignmentFile(args.input_bam_or_cram, "r", index_filename=args.read_index, reference_filename=args.reference_fasta)
    bamlet_file = pysam.AlignmentFile(args.bamlet, "wc" if args.bamlet.endswith(".cram") else "wb", template=input_bam_file)

    for region in args.region:
        chrom, start, end = parse_interval(region)

        # get the genomic regions first
        extract_region(
            chrom, start - 2000, end + 2000,
            input_bam=input_bam_file,
            bamlet=bamlet_file,
            merge_regions_distance=args.merge_regions_distance,
            verbose=args.verbose)

    bamlet_file.close()
    input_bam_file.close()

    try:
        pysam.sort("-o", f"{args.bamlet}.sorted.bam", args.bamlet)
        os.rename(f"{args.bamlet}.sorted.bam", args.bamlet)
        pysam.index(args.bamlet)
    except Exception as e:
        print(f"WARNING: Failed to sort and index {args.bamlet}: {e}")

if __name__ == "__main__":
    main()