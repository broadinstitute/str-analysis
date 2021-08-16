#!/usr/bin/env python3

"""
This script takes a bam/cram file and outputs a .json file with informatino
"""

import argparse
import collections
import json
import os
import pysam
from pprint import pprint
import re

from str_analysis.utils.canonical_repeat_unit import compute_canonical_repeat_unit
from str_analysis.utils.ehdn_info_for_locus import parse_ehdn_info_for_locus
from str_analysis.utils.most_frequent_repeat_unit import compute_most_frequent_repeat_unit

CANVAS_REPEAT_UNIT_SIZE = 5
MARGIN = 7   # base pairs
FLANK_SIZE = 2000   # base pairs
MIN_MAPQ = 3
MIN_FRACTION_OF_BASES_COVERED = 0.7


RFC1_LOCUS_COORDS_0BASED = {
    "37": ("4", 39350044, 39350099),
    "38": ("chr4", 39348424, 39348479),   # chr4:39348424-39348479
}

RFC1_LOCUS_KNOWN_ALLELES_BY_CATEGORY = {
    "BENIGN": {"AAAAG", "AAAGG"},
    "PATHOGENIC": {"AAGGG", "ACAGG"}, # ACAGG is from "A Māori specific RFC1 pathogenic repeat..." [Beecroft 2021]
    #"UNCERTAIN": {"AAGAG", "AGAGG",}, # from [Akcimen 2019)]
}

ALL_RFC1_LOCUS_KNOWN_ALLELES = {a for allele_set in RFC1_LOCUS_KNOWN_ALLELES_BY_CATEGORY.values() for a in allele_set}

GENOME_VERSION_ALIASES = {
    "GRCh37": "37", "hg19": "37", "hg37": "37", "37": "37",
    "GRCh38": "38", "hg38": "38", "38": "38",
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-g", "--genome-version", choices=GENOME_VERSION_ALIASES.keys(), required=True)
    p.add_argument("-R", "--reference", help="Reference fasta path. The reference fasta is sometimes necessary for decoding cram files.")
    p.add_argument("-e", "--ehdn-profile", help="If specified, information relevant to the RFC1 locus will be "
                                                "transferred from this ExpansionHunterDenovo profile to the output")
    p.add_argument("-s", "--sample-id", help="The sample id to put in the output json file. If not specified, it "
                                             "will be retrieved from the bam/cram header or filename prefix.")
    p.add_argument("-o", "--output-prefix", help="Output filename prefix")
    p.add_argument("-v", "--verbose", action="store_true", help="Print detailed log messages")
    p.add_argument("bam_or_cram_path", help="bam or cram path")

    args = p.parse_args()
    args.genome_version = GENOME_VERSION_ALIASES[args.genome_version]

    if not os.path.isfile(args.bam_or_cram_path):
        p.error(f"{args.bam_or_cram_path} not found")

    if args.reference and not os.path.isfile(args.reference):
        p.error(f"{args.reference} not found")

    return args


def main():
    args = parse_args()

    bam_cram_prefix = re.sub(".bam$|.cram$", "", os.path.basename(args.bam_or_cram_path))
    args.output_prefix = args.output_prefix or bam_cram_prefix

    chrom, locus_start_0based, locus_end = RFC1_LOCUS_COORDS_0BASED[args.genome_version]

    result = {}

    # process bam/cram
    print(f"Processing {args.bam_or_cram_path}")
    with pysam.Samfile(args.bam_or_cram_path, reference_filename=args.reference) as f:

        # try to get sample id from bam/cram header
        if not args.sample_id:
            read_groups = f.header.as_dict().get("RG", [])
            if read_groups:
                args.sample_id = read_groups[0].get("SM")
                if args.sample_id:
                    print(f"Using sample id '{args.sample_id}' from the bam/cram header")

            if not args.sample_id:
                args.sample_id = bam_cram_prefix
                print(f"Using sample id '{args.sample_id}' based on the input filename prefix")

        read_length = None
        while read_length is None:
            try:
                r = next(f)
                read_length = r.infer_read_length()
            except StopIteration:
                raise ValueError(f"{args.bam_or_cram_path} contains too few reads")

        # count reads in the left & right flanks to estimate read depth
        # NOTE: f.fetch retrieves all reads that *overlap* the given interval
        left_flank_n_well_aligned_bases = sum((r.query_alignment_length for r in f.fetch(
                chrom,
                locus_start_0based - MARGIN - FLANK_SIZE,
                locus_start_0based - MARGIN,
            ) if r.mapq >= MIN_MAPQ))
        right_flank_n_well_aligned_bases = sum((r.query_alignment_length for r in f.fetch(
                chrom,
                locus_end + MARGIN,
                locus_end + MARGIN + FLANK_SIZE,
            ) if r.mapq >= MIN_MAPQ))

        # get all sequences that overlap the RFC1 locus (regardless of whether they're soft-clipped)
        overlapping_sequences = []
        for r in f.fetch(chrom, locus_start_0based - MARGIN, locus_end + MARGIN):
            # see https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.query_alignment_sequence
            if r.mapq < MIN_MAPQ:
                continue

            read_sequence = r.seq
            #has_soft_clipped_bases_on_left = r.query_alignment_start > 0
            #has_soft_clipped_bases_on_right = r.query_alignment_end < len(read_sequence)
            read_start_pos_including_soft_clips = r.reference_start - r.query_alignment_start
            read_end_pos_including_soft_clips = read_start_pos_including_soft_clips + len(read_sequence)
            start_offset = 0
            if read_start_pos_including_soft_clips < locus_start_0based:
                start_offset = locus_start_0based - read_start_pos_including_soft_clips

            end_offset = len(read_sequence)
            if read_end_pos_including_soft_clips > locus_end:
                end_offset = locus_end - read_end_pos_including_soft_clips

            relevant_bases = read_sequence[start_offset:end_offset]
            if len(relevant_bases) >= CANVAS_REPEAT_UNIT_SIZE:
                overlapping_sequences.append(relevant_bases)

    left_flank_coverage = left_flank_n_well_aligned_bases / FLANK_SIZE
    right_flank_coverage = right_flank_n_well_aligned_bases / FLANK_SIZE
    result.update({
        "sample_id": args.sample_id,
        "left_flank_coverage": left_flank_coverage,
        "right_flank_coverage": right_flank_coverage,
    })

    # compute the repeat unit(s) found in the soft-clipped reads, and how many times each one occurs
    repeat_unit_to_n_occurrences = collections.defaultdict(int)
    repeat_unit_to_read_count = collections.defaultdict(int)
    for overlapping_sequence in overlapping_sequences:
        # in gnomAD, EHdn sometimes finds 6bp repeat units (eg. AAAGGG), so check for those as well
        for repeat_unit_size in (CANVAS_REPEAT_UNIT_SIZE, CANVAS_REPEAT_UNIT_SIZE + 1):
            if len(overlapping_sequence) < repeat_unit_size:
                continue

            repeat_unit, count = compute_most_frequent_repeat_unit(
                overlapping_sequence,
                repeat_unit_size=repeat_unit_size,
                min_occurrences=3,
                min_fraction_bases_covered=MIN_FRACTION_OF_BASES_COVERED)

            if args.verbose:
                if repeat_unit:
                    print(f"Found {compute_canonical_repeat_unit(repeat_unit)} occurs {count}x in read bases that "
                          f"overlap the RFC1 locus: {overlapping_sequence}")
                else:
                    print(f"Didn't find a consistent {repeat_unit_size}bp repeat unit in read bases "
                          f"that overlap the RFC1 locus: {overlapping_sequence}")

            if repeat_unit is not None:
                canonical_repeat_unit = compute_canonical_repeat_unit(repeat_unit)
                repeat_unit_to_read_count[canonical_repeat_unit] += 1
                repeat_unit_to_n_occurrences[canonical_repeat_unit] += count

    result.update({
        "found_n_reads_overlap_rfc1_locus": len(overlapping_sequences),
        "found_repeats_in_n_reads": sum(repeat_unit_to_read_count.values()),
        "found_repeats_in_fraction_of_reads": sum(repeat_unit_to_read_count.values())/len(overlapping_sequences),
    })

    # evaluate the repeat units
    well_supported_repeat_units = []
    for repeat_unit, read_count in repeat_unit_to_read_count.items():
        if "N" in repeat_unit:
            continue
        if read_count < 3:
            continue
        well_supported_repeat_units.append(repeat_unit)

    well_supported_repeat_units.sort(key=lambda repeat_unit: repeat_unit_to_n_occurrences[repeat_unit], reverse=True)

    n_pathogenic_alleles = 0
    n_benign_alleles = 0
    n_total_well_supported_alleles = 0
    for i in 0, 1:
        allele_number = i + 1
        result.update({
            f"allele{allele_number}_repeat_unit": None,
            f"allele{allele_number}_read_count": None,
            f"allele{allele_number}_n_occurrences": None,
        })

        if len(well_supported_repeat_units) <= i:
            continue
        n_total_well_supported_alleles += 1
        repeat_unit = well_supported_repeat_units[i]
        read_count = repeat_unit_to_read_count.get(repeat_unit)
        n_occurrences = repeat_unit_to_n_occurrences.get(repeat_unit)
        if repeat_unit in RFC1_LOCUS_KNOWN_ALLELES_BY_CATEGORY["PATHOGENIC"]:
            n_pathogenic_alleles += 1
        elif repeat_unit in RFC1_LOCUS_KNOWN_ALLELES_BY_CATEGORY["BENIGN"]:
            n_benign_alleles += 1

        result.update({
            f"allele{allele_number}_repeat_unit": repeat_unit,
            f"allele{allele_number}_read_count": read_count,
            f"allele{allele_number}_n_occurrences": n_occurrences,
        })

    # decide which combination of alleles is supported by the data
    # NOTE: there's no attempt to determine the size of the expansion and whether it's in the pathogenic range
    final_call = None
    if n_total_well_supported_alleles > 0 and left_flank_coverage >= 5 and right_flank_coverage >= 5:
        if n_pathogenic_alleles == n_total_well_supported_alleles:
            # reads support only known pathogenic allele(s)
            final_call = "PATHOGENIC MOTIF / PATHOGENIC MOTIF"
        elif n_benign_alleles == n_total_well_supported_alleles:
            # reads support only known benign allele(s)
            final_call = "BENIGN MOTIF / BENIGN MOTIF"
        elif n_benign_alleles == 0 and n_pathogenic_alleles == 0:
            # reads support one or more non-reference alleles of unknown significance
            final_call = "MOTIF OF UNCERTAIN SIGNIFICANCE / MOTIF OF UNCERTAIN SIGNIFICANCE"
        elif n_benign_alleles > 0 and n_pathogenic_alleles > 0:
            # reads support one known benign allele and one pathogenic allele
            final_call = "BENIGN MOTIF / PATHOGENIC MOTIF"
        elif n_pathogenic_alleles > 0:
            # reads support one pathogenic allele and at least one other allele of unknown significance
            final_call = "PATHOGENIC MOTIF / MOTIF OF UNCERTAIN SIGNIFICANCE"
        elif n_benign_alleles > 0:
            # reads support one known benign allele and at least one other allele of unknown significance
            final_call = "BENIGN MOTIF / MOTIF OF UNCERTAIN SIGNIFICANCE"

    result.update({
        "n_total_well_supported_alleles": n_total_well_supported_alleles,
        "n_benign_alleles": n_benign_alleles,
        "n_pathogenic_alleles": n_pathogenic_alleles,
    })

    result["call"] = final_call

    print(f"Final call: {final_call}")

    # process EHdn profile if one was provided
    if args.ehdn_profile:
        print(f"Parsing {args.ehdn_profile}")
        with open(args.ehdn_profile, "rt") as f:
            data = json.load(f)

        records, sample_read_depth, _ = parse_ehdn_info_for_locus(data, chrom, locus_start_0based, locus_end)
        result["ehdn_sample_read_depth"] = sample_read_depth

        # get the 2 alleles with the most read support
        records.sort(key=lambda r: (
            -r["anchored_irr_count_for_this_repeat_unit_and_region"],
            -r["total_irr_count_for_this_repeat_unit_and_region"],
            r["repeat_unit"],
        ))

        for i in 0, 1:
            record = records[i] if len(records) > i else {}
            allele_number = i + 1
            result.update({
                f"ehdn_allele{allele_number}_repeat_unit": record.get("repeat_unit"),
                f"ehdn_allele{allele_number}_anchored_irr_count": record.get("anchored_irr_count_for_this_repeat_unit_and_region"),
                f"ehdn_allele{allele_number}_n_anchored_regions": record.get("n_anchored_regions_for_this_repeat_unit"),
                f"ehdn_allele{allele_number}_paired_irr_count": record.get("paired_irr_count_for_this_repeat_unit"),
                f"ehdn_allele{allele_number}_total_irr_count": record.get("total_irr_count_for_this_repeat_unit_and_region"),
            })

    # generate output
    output_filename = f"{args.output_prefix}.rfc1_canvas_alleles.json"
    with open(output_filename, "wt") as f:
        json.dump(result, f, indent=2)
    print(f"Wrote results to {output_filename}")
    if args.verbose:
        pprint(result)


if __name__ == "__main__":
    main()