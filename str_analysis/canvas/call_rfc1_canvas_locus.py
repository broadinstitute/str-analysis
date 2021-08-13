"""
This script takes a WGS bam or cram file and outputs a .json file with info related to the RFC1/CANVAS STR locus,
including the following fields:

"rfc1_locus_n_well_aligned_reads": the number of reads with well-aligned (non-soft-clipped) bases within the RFC1 locus,
    implying read support for the presence of the benign reference allele
"call": this has a format analogous to a VCF genotype, and can be:
    "PATHOGENIC MOTIF / PATHOGENIC MOTIF"
    "BENIGN MOTIF / BENIGN MOTIF"
    "MOTIF OF UNCERTAIN SIGNIFICANCE / MOTIF OF UNCERTAIN SIGNIFICANCE"
    "BENIGN MOTIF / PATHOGENIC MOTIF"
    "PATHOGENIC MOTIF / MOTIF OF UNCERTAIN SIGNIFICANCE"
    "BENIGN MOTIF / MOTIF OF UNCERTAIN SIGNIFICANCE"
    or it will be None if there's not enough evidence in the read data to support any of the above options.

Also, this script optionally takes the ExpansionHunterDenovo profile for this sample and adds fields to the
output json based on values in the ExpansionHunterDenovo profile. The ExpansionHunterDenovo profile is optional -
the values are simply transfered to the output file, and aren't used for any calculations.
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

MARGIN = 7   # base pairs
FLANK_SIZE = 2000   # base pairs
MIN_MAPQ = 3
CANVAS_REPEAT_UNIT_SIZE = 5

RFC1_LOCUS_COORDS_0BASED = {
    "37": ("4", 39350044, 39350099),
    "38": ("chr4", 39348427, 39348475),   # chr4:39,348,427-39,348,475
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
    p.add_argument("-e", "--ehdn-profile", help="If specified, the relevant information about repeat motifs and read counts at the RFC1 locus "
                                                "will be retrieved from this ExpansionHunterDenovo profile and added to the output")
    p.add_argument("-s", "--sample-id", help="The sample id to put in the output json file. If not specified, it will be retrieved from the bam/cram file.")
    p.add_argument("-o", "--output-prefix", help="Output filename prefix")
    p.add_argument("-v", "--verbose", action="store_true", help="Print detailed log messages")
    p.add_argument("reads", help="bam or cram path")

    args = p.parse_args()
    args.genome_version = GENOME_VERSION_ALIASES[args.genome_version]

    if not os.path.isfile(args.reads):
        p.error(f"{args.reads} not found")

    if args.reference and not os.path.isfile(args.reference):
        p.error(f"{args.reference} not found")

    return args


def main():
    args = parse_args()

    bam_cram_prefix = re.sub(".bam$|.cram$", "", os.path.basename(args.reads))
    args.output_prefix = args.output_prefix or bam_cram_prefix

    chrom, start, end = RFC1_LOCUS_COORDS_0BASED[args.genome_version]

    result = {}

    # process bam/cram
    print(f"Processing {args.reads}")
    with pysam.Samfile(args.reads, reference_filename=args.reference) as f:

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
                raise ValueError(f"{args.reads} contains too few reads")

        # count reads in the left & right flanks to estimate read depth
        # NOTE: f.fetch retrieves all reads that *overlap* the given interval
        left_flank_n_well_aligned_bases = sum(
            (r.query_alignment_end - r.query_alignment_start
             for r in f.fetch(chrom, start - MARGIN - FLANK_SIZE, start - MARGIN) if r.mapq >= MIN_MAPQ))
        right_flank_n_well_aligned_bases = sum(
            (r.query_alignment_end - r.query_alignment_start
             for r in f.fetch(chrom, end + MARGIN, end + MARGIN + FLANK_SIZE) if r.mapq >= MIN_MAPQ))

        # count the number of reads that have MAPQ above the threshold, and have non-soft-clipped bases aligning
        # within the RFC1 AAAAG repeat in the reference genome. Add MARGIN to avoid counting reads where a few
        # bases match AAAAG by chance.
        rfc1_locus_n_well_aligned_reads = sum(
            (1 for r in f.fetch(chrom, start + MARGIN, end - MARGIN) if r.mapq >= MIN_MAPQ
             and (r.query_alignment_start < end - MARGIN or r.query_alignment_end > start + MARGIN))
        )

        # get soft-clipped sequences
        soft_clipped_sequences = []
        for r in f.fetch(chrom, start - MARGIN, end + MARGIN):
            # see https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.query_alignment_sequence
            if r.mapq < MIN_MAPQ:
                continue
            read_sequence = r.seq
            left_softclipped_bases = read_sequence[:r.query_alignment_start]
            if len(left_softclipped_bases) >= CANVAS_REPEAT_UNIT_SIZE:
                soft_clipped_sequences.append(left_softclipped_bases)
            right_softclipped_bases = read_sequence[r.query_alignment_end:]
            if len(right_softclipped_bases) >= CANVAS_REPEAT_UNIT_SIZE:
                soft_clipped_sequences.append(right_softclipped_bases)

    left_flank_coverage = left_flank_n_well_aligned_bases / FLANK_SIZE
    right_flank_coverage = right_flank_n_well_aligned_bases / FLANK_SIZE
    result.update({
        "sample_id": args.sample_id,
        "rfc1_locus_n_well_aligned_reads": rfc1_locus_n_well_aligned_reads,
        "left_flank_coverage": left_flank_coverage,
        "right_flank_coverage": right_flank_coverage,
    })

    # compute the repeat unit(s) found in the soft-clipped reads, and how many times each one occurs
    repeat_unit_to_n_occurrences = collections.defaultdict(int)
    repeat_unit_to_read_count = collections.defaultdict(int)
    for soft_clipped_sequence in soft_clipped_sequences:
        # in gnomAD, EHdn sometimes finds 6bp repeat units (eg. AAAGGG), so check for those as well
        for repeat_unit_size in (CANVAS_REPEAT_UNIT_SIZE, CANVAS_REPEAT_UNIT_SIZE + 1):
            if len(soft_clipped_sequence) < repeat_unit_size:
                continue
            repeat_unit, count = compute_most_frequent_repeat_unit(
                soft_clipped_sequence, repeat_unit_size=repeat_unit_size)
            if args.verbose:
                if repeat_unit:
                    print(f"Found {compute_canonical_repeat_unit(repeat_unit)} occurs {count}x in these "
                          f"soft-clipped bases within the RFC1 locus: {soft_clipped_sequence}")
                else:
                    print(f"Didn't find a consistent {repeat_unit_size}bp repeat unit in these softed-clipped bases "
                          f"within the RFC1 locus: {soft_clipped_sequence}")

            if repeat_unit is not None:
                canonical_repeat_unit = compute_canonical_repeat_unit(repeat_unit)
                repeat_unit_to_read_count[canonical_repeat_unit] += 1
                repeat_unit_to_n_occurrences[canonical_repeat_unit] += count

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
    # absence of reads without soft-clips at the RFC1 locus suggests absence of the benign reference allele (AAAAG)
    if rfc1_locus_n_well_aligned_reads >= 5:
        # well-aligned reads within the RFC1 locus imply presence of the benign reference allele (AAAAG)
        n_benign_alleles += 1
        n_total_well_supported_alleles += 1

    # NOTE: we don't try to determine the size of the expansion and whether it's in the pathogenic range
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

        records, sample_read_depth, _ = parse_ehdn_info_for_locus(data, chrom, start, end, motifs_of_interest=None)
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
    output_filename = f"{args.output_prefix}.rfc1_canvas_locus.json"
    with open(output_filename, "wt") as f:
        json.dump(result, f)
    print(f"Wrote results to {output_filename}")
    if args.verbose:
        pprint(result)


if __name__ == "__main__":
    main()