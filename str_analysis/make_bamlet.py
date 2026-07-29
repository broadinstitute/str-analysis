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
from str_analysis.utils.cram_bam_utils import normalize_chromosome_name

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
        # same filter extract_region applies when collecting primary reads, matching
        # https://github.com/bw2/ExpansionHunter/blob/master/ehunter/sample/MateExtractor.cpp#L143 -- the two
        # mate-selection paths had drifted apart, letting supplementary alignments through only on this one
        if alignment.is_secondary or alignment.is_supplementary:
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


def extract_region(chrom, start, end, input_bam, bamlet, merge_regions_distance=1000, verbose=False,
                   written_read_keys=None):
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
        alignment = read_pair[0]
        # skip unpaired reads (single-end / long-read alignments): they have no mate to fetch, and their
        # next_reference_name is None, which would crash add_interval downstream
        if not alignment.is_paired:
            continue

        # a pair is complete only when both mates (read1 and read2) are present; testing len(read_pair) >= 2
        # would wrongly treat two copies of the same end as complete and skip fetching the real mate
        if any(a.is_read1 for a in read_pair) and any(a.is_read2 for a in read_pair):
            continue

        # An unmapped mate has no genomic interval to retrieve. htslib represents its reference id and start as -1;
        # letting either value reach add_interval would produce an invalid chromosome or negative fetch coordinate.
        if alignment.mate_is_unmapped or alignment.next_reference_id < 0 or alignment.next_reference_start < 0:
            continue

        # see if mate is close to other mates that need to be fetched
        mate_chrom = alignment.next_reference_name
        mate_pos = int(alignment.next_reference_start)
        mate_start = max(0, mate_pos - 1)
        mate_end = mate_pos + 1

        #if mate_chrom == chrom and min(abs(mate_pos - start), abs(mate_pos - end)) <= 1000:
        #    # skip mates that are close to the ends of the region
        #    continue

        for mate_region in mate_regions:
            if is_close(mate_chrom, mate_pos, mate_region, max_dist=merge_regions_distance):
                previous_read_names_set = mate_regions.pop(mate_region)
                previous_read_names_set.add(read_name)
                key = (mate_chrom, min(mate_start, mate_region[1]), max(mate_end, mate_region[2]))
                mate_regions[key] = previous_read_names_set
                break
        else:
            key = mate_chrom, mate_start, mate_end
            mate_regions[key].add(read_name)

    genomic_regions_to_fetch.extend(mate_regions.keys())

    if verbose:
        # read_pairs holds every pair collected from the region; the mates still to fetch are the read names
        # accumulated in mate_regions, which is a strictly smaller set
        mate_count = sum(len(read_names) for read_names in mate_regions.values())
        print(f"{chrom}:{start}-{end}: Need to fetch {mate_count} mates from {len(mate_regions)} regions")

    if bamlet is not None:
        for (mate_chrom, mate_start, mate_end), read_names_set in mate_regions.items():
            for alignment in jump_for_mates(input_bam, mate_chrom, mate_start, mate_end, read_names_set):
                read_pairs[alignment.query_name].append(alignment)

        read_counter = 0
        for read_name, read_pair in read_pairs.items():
            for read in read_pair:
                if written_read_keys is not None:
                    # read_pairs is local to this call, so a read overlapping two requested regions would otherwise
                    # be written once per region; duplicates are misread downstream as a complete read pair and
                    # suppress real mate retrieval
                    read_key = (read.query_name, read.flag, read.reference_start)
                    if read_key in written_read_keys:
                        continue
                    written_read_keys.add(read_key)
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
    parser.add_argument("-o", "--bamlet", help="Output BAMlet path; a .cram extension writes CRAM, otherwise BAM")
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

    # map normalized chrom name -> the exact reference name in this file's header, so fetch() always uses a name
    # valid for the header even when the region was given with a different naming convention ("9" vs "chr9").
    # Without this, pysam raises "invalid contig"; make_minicram_for_expansion_hunter resolves names the same way.
    normalized_to_reference_name = {
        normalize_chromosome_name(name): name for name in input_bam_file.references
    }

    # every region is resolved BEFORE the output file is opened. parser.error raises SystemExit, and running this
    # check inside the write loop below meant an unknown contig in the 2nd of several regions aborted after the
    # 1st region's reads had already been written, leaving an unsorted, unindexed partial BAMlet at the -o path.
    regions = []
    for region in args.region:
        chrom, start, end = parse_interval(region)
        fetch_chrom = normalized_to_reference_name.get(normalize_chromosome_name(chrom))
        if fetch_chrom is None:
            parser.error(f"Chromosome {chrom} from region {region} is not present in "
                         f"{args.input_bam_or_cram}")
        regions.append((fetch_chrom, start, end))

    # CRAM output needs the reference to encode records; without it htslib falls back to resolving the header's
    # M5/UR tags. no_ref=1 stores sequences verbatim instead, avoiding reference-md5 validation against a
    # reference build that may not match the input alignment (same choice save_to_file makes for CRAM output).
    write_cram = args.bamlet.endswith(".cram")
    bamlet_file = pysam.AlignmentFile(
        args.bamlet, "wc" if write_cram else "wb", template=input_bam_file,
        reference_filename=args.reference_fasta if write_cram else None,
        format_options=[b"no_ref=1"] if write_cram else None)

    # shared across regions so a read overlapping more than one requested region is written only once
    written_read_keys = set()
    for fetch_chrom, start, end in regions:
        # get the genomic regions first. max(0, ...) keeps the padded start on the contig -- pysam's fetch()
        # raises "start out of range" for a negative coordinate, which any locus within 2000bp of a contig start
        # would otherwise produce.
        extract_region(
            fetch_chrom, max(0, start - 2000), end + 2000,
            input_bam=input_bam_file,
            bamlet=bamlet_file,
            merge_regions_distance=args.merge_regions_distance,
            verbose=args.verbose,
            written_read_keys=written_read_keys)

    bamlet_file.close()
    input_bam_file.close()

    try:
        # sort into a file of the SAME format as the requested output: samtools infers the format from the
        # extension, so sorting a .cram request into a .sorted.bam and renaming it over the .cram path would
        # leave BAM bytes (and a .bai) behind a .cram name
        is_cram_output = args.bamlet.endswith(".cram")
        sorted_path = f"{args.bamlet}.sorted.cram" if is_cram_output else f"{args.bamlet}.sorted.bam"
        sort_args = ["-o", sorted_path]
        if is_cram_output:
            # no_ref=1 must be repeated here: the file was written self-contained to avoid reference-MD5
            # validation, and re-encoding it during the sort WITHOUT that option would reference-compress it
            # and reintroduce exactly the M5 mismatch the initial write avoided
            sort_args += ["--output-fmt", "CRAM", "--output-fmt-option", "no_ref=1",
                          "--reference", args.reference_fasta]
        pysam.sort(*sort_args, args.bamlet)
        os.rename(sorted_path, args.bamlet)
        pysam.index(args.bamlet)
    except Exception as e:
        print(f"WARNING: Failed to sort and index {args.bamlet}: {e}")

if __name__ == "__main__":
    main()
