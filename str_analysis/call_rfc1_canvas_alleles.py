#!/usr/bin/env python3

"""
This script takes a bam/cram file and outputs a .json file with information about repeat motifs it detects at the
RFC1/CANVAS locus.
"""

import argparse
import ast
import collections
import gzip
import json
import os
import pkgutil
import pysam
from pprint import pprint, pformat
import re
import subprocess

from str_analysis.utils.canonical_repeat_unit import compute_canonical_repeat_unit
from str_analysis.utils.ehdn_info_for_locus import parse_ehdn_info_for_locus
from str_analysis.utils.most_frequent_repeat_unit import compute_most_frequent_repeat_unit
from str_analysis.utils.parse_interval import parse_interval


CANVAS_MOTIF_SIZE = 5
MARGIN = 7   # base pairs
FLANK_SIZE = 2000   # base pairs
MIN_MAPQ = 3
MIN_FRACTION_OF_BASES_COVERED = 0.7

NORMALIZE_TO_COVERAGE = 40

RFC1_LOCUS_COORDS_0BASED = {
    "37": ("4", 39350044, 39350099),      # chr4:39350044-39350099
    "38": ("chr4", 39348424, 39348479),   # chr4:39348425-39348479
}

RFC1_LOCUS_KNOWN_MOTIFS_BY_CATEGORY = {
        "BENIGN": {"AAAAG", "AAAGG"},
    "PATHOGENIC": {"AAGGG", "ACAGG"}, # ACAGG is from "A Māori specific RFC1 pathogenic repeat..." [Beecroft 2021]
    #"UNCERTAIN": {"AAGAG", "AGAGG",}, # from [Akcimen 2019]
}

ALL_RFC1_LOCUS_KNOWN_MOTIFS = {a for repeat_unit_set in RFC1_LOCUS_KNOWN_MOTIFS_BY_CATEGORY.values() for a in repeat_unit_set}

GENOME_VERSION_ALIASES = {
    "GRCh37": "37", "hg19": "37", "hg37": "37", "37": "37",
    "GRCh38": "38", "hg38": "38", "38": "38",
}


OFFTARGET_REGIONS = json.loads(gzip.decompress(pkgutil.get_data(__name__, "data/offtarget_regions.json.gz")))


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-g", "--genome-version", choices=GENOME_VERSION_ALIASES.keys(), required=True)
    p.add_argument("-R", "--reference-fasta", help="Reference fasta path. The reference fasta is sometimes necessary "
                                                   "for decoding cram files.")
    p.add_argument("-o", "--output-prefix", help="Output filename prefix")
    p.add_argument("-s", "--sample-id", help="The sample id to put in the output json file. If not specified, it "
        "will be retrieved from the bam/cram header or filename prefix.")

    p.add_argument("-e", "--ehdn-profile", help="If specified, information relevant to the RFC1 locus will be "
        "transferred from this ExpansionHunterDenovo profile to the output json")

    p.add_argument("-r", "--run-expansion-hunter", action="store_true", help="If this option is specified, this "
         "script will run ExpansionHunter once for each of the motif(s) it detects at the RFC1 locus. "
         "ExpansionHunter doesn't currently support genotyping multiallelic repeats such as RFC1 where "
         "an individual may have 2 alleles with motifs that differ from eachother (and from the reference motif). "
         "Running ExpansionHunter separately for each motif provides a work-around.")
    p.add_argument("--expansion-hunter-path", help="The path of the ExpansionHunter executable to use if -r is "
        "specified. This must be ExpansionHunter v3 or v4.", default="ExpansionHunter")

    p.add_argument("--ignore-offtarget-regions", action="store_true", help="Don't compute read support in off-target "
        "regions. The output will be the same, except that it won't contain *_read_count_with_offtargets fields. "
        "Also, if --run-expansion-hunter is used, ExpansionHunter will be run without off-target regions.")

    p.add_argument("-v", "--verbose", action="store_true", help="Print detailed log messages")
    p.add_argument("bam_or_cram_path", help="bam or cram path")

    args = p.parse_args()
    args.genome_version = GENOME_VERSION_ALIASES[args.genome_version]

    if not os.path.isfile(args.bam_or_cram_path):
        p.error(f"{args.bam_or_cram_path} not found")

    if args.reference_fasta and not os.path.isfile(args.reference_fasta):
        p.error(f"{args.reference_fasta} not found")

    if args.run_expansion_hunter and not args.reference_fasta:
        p.error("--reference-fasta is required when --run-expansion-hunter is used")

    return args


def generate_variant_catalog(locus_id, repeat_unit, chrom, start_1based, end_1based, offtarget_regions=None):
    return {
        "LocusId": locus_id,
        "LocusStructure": f"({repeat_unit})*",
        "ReferenceRegion": f"{chrom}:{start_1based - 1}-{end_1based}",
        "VariantType": "RareRepeat",
        "OfftargetRegions": [] if not offtarget_regions else offtarget_regions,
    }


def run_expansion_hunter(
        sample_id,
        expansion_hunter_path,
        genome_version,
        reference_fasta_path,
        bam_or_cram_path,
        repeat_units,
        result,
        use_offtarget_regions=True,
        output_dir=".",
        verbose=False,
):
    """
    :param well_supported_repeat_units:
    :param expansion_hunter_path:
    :param result:
    :param verbose:
    :return:
    """

    if genome_version not in ("37", "38"):
        raise ValueError(f"Unexpected genome version: {genome_version}. Must be '37' or '38'")

    chrom, start_1based, end_1based = RFC1_LOCUS_COORDS_0BASED[genome_version]

    for repeat_unit_number, repeat_unit in enumerate(repeat_units):
        # generate variant catalog
        variant_catalog_locus_label = f"RFC1_{repeat_unit}"
        variant_catalog = generate_variant_catalog(
            variant_catalog_locus_label,
            repeat_unit, chrom, start_1based, end_1based,
            offtarget_regions=OFFTARGET_REGIONS[genome_version][repeat_unit] if use_offtarget_regions else [])

        variant_catalog_path = os.path.join(output_dir, f"{repeat_unit}.variant_catalog.json")
        with open(variant_catalog_path, "wt") as f:
            json.dump([variant_catalog], f)

        # run expansion hunter
        print("--"*10)
        print(f"Running ExpansionHunter on {sample_id} for repeat unit {repeat_unit}")
        if verbose:
            print("Using variant catalog: ")
            pprint(variant_catalog)

        filename_prefix = f"{sample_id}.{repeat_unit}"
        output_prefix = f"{os.path.join(output_dir, filename_prefix)}.expansion_hunter4"
        expansion_hunter_command = f"""{expansion_hunter_path} \
--sex male \
--reference {reference_fasta_path} \
--reads {bam_or_cram_path} \
--variant-catalog {variant_catalog_path} \
--output-prefix {output_prefix} \
--log-level debug
"""

        if verbose:
            print(expansion_hunter_command)

        subprocess.run(expansion_hunter_command, shell=True, stderr=subprocess.STDOUT, check=False)

        # parse result
        if os.path.isfile(f"{output_prefix}.json"):
            with open(f"{output_prefix}.json", "rt") as f:
                expansion_hunter_output_json = json.load(f)
        else:
            print(f"ERROR: ExpansionHunter didn't produce a {output_prefix}.json file. Skipping...")
            continue

        if verbose:
            print(f"ExpansionHunter output: {pformat(expansion_hunter_output_json)}")


        eh_result = expansion_hunter_output_json.get("LocusResults", {}).get(variant_catalog_locus_label, {}).get(
            "Variants", {}).get(variant_catalog_locus_label, {})

        if not eh_result:
            (result[f"motif{repeat_unit_number}_expansion_hunter_short_allele_genotype"],
             result[f"motif{repeat_unit_number}_expansion_hunter_long_allele_genotype"]) = None, None
            continue


        if eh_result.get("Genotype"):
            (
                result[f"motif{repeat_unit_number}_expansion_hunter_short_allele_genotype"],
                result[f"motif{repeat_unit_number}_expansion_hunter_long_allele_genotype"]
            ) = [
                int(g) for g in eh_result["Genotype"].split("/")]

        if eh_result.get("GenotypeConfidenceInterval"):
            (
                result[f"motif{repeat_unit_number}_expansion_hunter_short_allele_CI_start"],
                result[f"motif{repeat_unit_number}_expansion_hunter_short_allele_CI_end"],
                result[f"motif{repeat_unit_number}_expansion_hunter_long_allele_CI_start"],
                result[f"motif{repeat_unit_number}_expansion_hunter_long_allele_CI_end"],
            ) = [
                int(b) for ci in eh_result["GenotypeConfidenceInterval"].split("/") for b in ci.split("-")
            ]

            result[f"motif{repeat_unit_number}_expansion_hunter_short_allele_CI_size"] = (
                    result[f"motif{repeat_unit_number}_expansion_hunter_short_allele_CI_end"] -
                    result[f"motif{repeat_unit_number}_expansion_hunter_short_allele_CI_start"]
            )

            result[f"motif{repeat_unit_number}_expansion_hunter_long_allele_CI_size"] = (
                    result[f"motif{repeat_unit_number}_expansion_hunter_long_allele_CI_end"] -
                    result[f"motif{repeat_unit_number}_expansion_hunter_long_allele_CI_start"]
            )

        for output_label in "spanning_reads", "flanking_reads", "inrepeat_reads":
            read_count_label = "CountsOf" + "".join(word.title() for word in output_label.split('_'))

            if not eh_result.get(read_count_label) or eh_result[read_count_label] == "()":
                # eg. 'CountsOfSpanningReads': '()'
                total = 0
            else:
                # eg .'CountsOfInrepeatReads': '(30, 4), (31, 6)',
                try:
                    read_count_tuples = ast.literal_eval(f"[{eh_result[read_count_label]}]")
                    total = sum(t[1] for t in read_count_tuples)
                except Exception as e:
                    print(f"ERROR: unable to parse {read_count_label}: {read_count_tuples}. {e}")
                    continue

            result[f"motif{repeat_unit_number}_expansion_hunter_total_{output_label}"] = total


def main():
    args = parse_args()

    bam_cram_prefix = re.sub(".bam$|.cram$", "", os.path.basename(args.bam_or_cram_path))
    args.output_prefix = args.output_prefix or bam_cram_prefix

    locus_chrom, locus_start_0based, locus_end = RFC1_LOCUS_COORDS_0BASED[args.genome_version]

    result = {}

    # process bam/cram
    print(f"Processing {args.bam_or_cram_path}")
    with pysam.Samfile(args.bam_or_cram_path, reference_filename=args.reference_fasta) as f:

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
                locus_chrom,
                locus_start_0based - MARGIN - FLANK_SIZE,
                locus_start_0based - MARGIN,
            ) if r.mapq >= MIN_MAPQ))
        right_flank_n_well_aligned_bases = sum((r.query_alignment_length for r in f.fetch(
                locus_chrom,
                locus_end + MARGIN,
                locus_end + MARGIN + FLANK_SIZE,
            ) if r.mapq >= MIN_MAPQ))

        # get all sequences that overlap the RFC1 locus (regardless of whether they're soft-clipped)
        overlapping_sequences = []
        for r in f.fetch(locus_chrom, locus_start_0based - MARGIN, locus_end + MARGIN):
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
            if len(relevant_bases) >= CANVAS_MOTIF_SIZE:
                overlapping_sequences.append(relevant_bases)

    left_flank_coverage = left_flank_n_well_aligned_bases / FLANK_SIZE
    right_flank_coverage = right_flank_n_well_aligned_bases / FLANK_SIZE
    result.update({
        "sample_id": args.sample_id,
        "locus": f"{locus_chrom}:{locus_start_0based}-{locus_end}",
        "genome_version": args.genome_version,
        "left_flank_coverage": left_flank_coverage,
        "right_flank_coverage": right_flank_coverage,
    })

    # compute the motif(s) found in the soft-clipped reads, and how many times each one occurs
    motif_to_n_occurrences = collections.defaultdict(int)
    motif_to_read_count = collections.defaultdict(int)
    for overlapping_sequence in overlapping_sequences:
        # in gnomAD, EHdn sometimes finds 6bp repeat units (eg. AAAGGG), so check for those as well
        for motif_size in (CANVAS_MOTIF_SIZE, CANVAS_MOTIF_SIZE + 1):
            if len(overlapping_sequence) < motif_size:
                continue

            motif, count = compute_most_frequent_repeat_unit(
                overlapping_sequence,
                repeat_unit_size=motif_size,
                min_occurrences=3,
                min_fraction_bases_covered=MIN_FRACTION_OF_BASES_COVERED)

            if args.verbose:
                if motif:
                    print(f"Found {compute_canonical_repeat_unit(motif)} occurs {count}x in read bases that "
                          f"overlap the RFC1 locus: {overlapping_sequence}")
                else:
                    if motif_size == CANVAS_MOTIF_SIZE:
                        print(f"Didn't find a consistent {motif_size}bp repeat unit in read bases "
                            f"that overlap the RFC1 locus: {overlapping_sequence}")

            if motif is not None:
                canonical_motif = compute_canonical_repeat_unit(motif)
                motif_to_read_count[canonical_motif] += 1
                motif_to_n_occurrences[canonical_motif] += count

    result.update({
        "found_n_reads_overlap_rfc1_locus": len(overlapping_sequences),
        "found_repeats_in_n_reads": sum(motif_to_read_count.values()),
        "found_repeats_in_fraction_of_reads": sum(motif_to_read_count.values())/len(overlapping_sequences) if overlapping_sequences else 0,
    })

    # evaluate the repeat units
    well_supported_motifs = []
    for motif, read_count in motif_to_read_count.items():
        if "N" in motif:
            continue
        if read_count < 3:
            continue
        well_supported_motifs.append(motif)

    # select the repeat unit(s) with the most read support
    well_supported_motifs.sort(key=lambda motif: motif_to_n_occurrences[motif], reverse=True)
    selected_motifs = well_supported_motifs[:2]

    # sort then into BENIGN .. PATHOGENIC .. UNCERTAIN SIGNIFICANCE to match the order in the "call" output field
    selected_motifs = sorted(selected_motifs, key=lambda motif:
        1 if motif in RFC1_LOCUS_KNOWN_MOTIFS_BY_CATEGORY["BENIGN"] else
        2 if motif in RFC1_LOCUS_KNOWN_MOTIFS_BY_CATEGORY["PATHOGENIC"] else
        3)

    if args.run_expansion_hunter:
        run_expansion_hunter(
            args.sample_id,
            args.expansion_hunter_path,
            args.genome_version,
            args.reference_fasta,
            args.bam_or_cram_path,
            selected_motifs,
            result,
            use_offtarget_regions=not args.ignore_offtarget_regions,
            verbose=args.verbose)

    flank_coverage_mean = (left_flank_coverage + right_flank_coverage) / 2.0
    n_pathogenic_motifs = 0
    n_benign_motifs = 0
    n_total_well_supported_motifs = 0
    for i in 0, 1:
        motif_number = i + 1
        if len(selected_motifs) <= i:
            continue

        n_total_well_supported_motifs += 1
        motif = selected_motifs[i]
        read_count = motif_to_read_count.get(motif)
        n_occurrences = motif_to_n_occurrences.get(motif)
        if motif in RFC1_LOCUS_KNOWN_MOTIFS_BY_CATEGORY["PATHOGENIC"]:
            n_pathogenic_motifs += 1
        elif motif in RFC1_LOCUS_KNOWN_MOTIFS_BY_CATEGORY["BENIGN"]:
            n_benign_motifs += 1

        result.update({
            f"motif{motif_number}_repeat_unit": motif,
            f"motif{motif_number}_read_count": read_count,
            f"motif{motif_number}_normalized_read_count":
                read_count * NORMALIZE_TO_COVERAGE / flank_coverage_mean if flank_coverage_mean > 0 else 0,
            f"motif{motif_number}_n_occurrences": n_occurrences,
        })

        if not args.ignore_offtarget_regions:
            read_count_with_offtargets = read_count
            with pysam.Samfile(args.bam_or_cram_path, reference_filename=args.reference_fasta) as f:
                for offtarget_region in OFFTARGET_REGIONS[args.genome_version][motif]:
                    offtarget_chrom, offtarget_start, offtarget_end = parse_interval(offtarget_region)
                    sequences = (r.seq for r in f.fetch(offtarget_chrom, offtarget_start, offtarget_end) if r.mapq >= MIN_MAPQ)
                    c, t = count_repeat_in_sequences(
                        sequences,
                        motif,
                        min_occurrences=3,
                        min_fraction_bases_covered=MIN_FRACTION_OF_BASES_COVERED)

                    read_count_with_offtargets += c
                    if args.verbose:
                        print(f"{c} out of {t} reads contained {motif} in off-target region {offtarget_region}")

            result.update({
                f"motif{motif_number}_read_count_with_offtargets": read_count_with_offtargets,
                f"motif{motif_number}_normalized_read_count_with_offtargets":
                    read_count_with_offtargets * NORMALIZE_TO_COVERAGE / flank_coverage_mean if flank_coverage_mean > 0 else 0,
            })

    # decide which combination of motifs is supported by the data
    # NOTE: there's no attempt to determine the size of the expansion and whether it's in the pathogenic range
    if n_total_well_supported_motifs == 0:
        final_call = f"NO CALL (no motif has sufficient read support)"
    elif n_pathogenic_motifs == n_total_well_supported_motifs:
        # reads support only known pathogenic motif(s)
        final_call = "PATHOGENIC MOTIF / PATHOGENIC MOTIF"
    elif n_benign_motifs == n_total_well_supported_motifs:
        # reads support only known benign motif(s)
        final_call = "BENIGN MOTIF / BENIGN MOTIF"
    elif n_benign_motifs == 0 and n_pathogenic_motifs == 0:
        # reads support one or more non-reference motifs of unknown significance
        final_call = "MOTIF OF UNCERTAIN SIGNIFICANCE / MOTIF OF UNCERTAIN SIGNIFICANCE"
    elif n_benign_motifs > 0 and n_pathogenic_motifs > 0:
        # reads support one known benign motif and one pathogenic motif
        final_call = "BENIGN MOTIF / PATHOGENIC MOTIF"
    elif n_pathogenic_motifs > 0:
        # reads support one pathogenic motif and at least one other motif of unknown significance
        final_call = "PATHOGENIC MOTIF / MOTIF OF UNCERTAIN SIGNIFICANCE"
    elif n_benign_motifs > 0:
        # reads support one known benign motif and at least one other motif of unknown significance
        final_call = "BENIGN MOTIF / MOTIF OF UNCERTAIN SIGNIFICANCE"

    result.update({
        "n_total_well_supported_motifs": n_total_well_supported_motifs,
        "n_benign_motifs": n_benign_motifs,
        "n_pathogenic_motifs": n_pathogenic_motifs,
    })

    result["call"] = final_call

    print(f"Final call: {final_call}")

    # process EHdn profile if one was provided
    if args.ehdn_profile:
        print(f"Parsing {args.ehdn_profile}")
        with open(args.ehdn_profile, "rt") as f:
            data = json.load(f)

        records, sample_read_depth, _ = parse_ehdn_info_for_locus(data, locus_chrom, locus_start_0based, locus_end)
        result["ehdn_sample_read_depth"] = sample_read_depth

        # get the 2 motifs with the most read support
        records.sort(key=lambda r: (
            -r["anchored_irr_count_for_this_repeat_unit_and_region"],
            -r["total_irr_count_for_this_repeat_unit_and_region"],
            r["repeat_unit"],
        ))

        for i in 0, 1:
            record = records[i] if len(records) > i else {}
            motif_number = i + 1
            result.update({
                f"ehdn_motif{motif_number}_repeat_unit": record.get("repeat_unit"),
                f"ehdn_motif{motif_number}_anchored_irr_count": record.get("anchored_irr_count_for_this_repeat_unit_and_region"),
                f"ehdn_motif{motif_number}_n_anchored_regions": record.get("n_anchored_regions_for_this_repeat_unit"),
                f"ehdn_motif{motif_number}_paired_irr_count": record.get("paired_irr_count_for_this_repeat_unit"),
                f"ehdn_motif{motif_number}_total_irr_count": record.get("total_irr_count_for_this_repeat_unit_and_region"),
            })

    # generate output
    output_filename = f"{args.output_prefix}.rfc1_canvas_alleles.json"
    with open(output_filename, "wt") as f:
        json.dump(result, f, indent=2)
    print(f"Wrote results to {output_filename}")
    pprint(result)


def count_repeat_in_sequences(sequences, repeat_unit, min_occurrences=3, min_fraction_bases_covered=0.8):
    """Count how many of the given sequences support a specific repeat unit.

    Args:
        sequences (str): iterator over read sequences
        repeat_unit (str): the repeat unit to search for within the sequences
        min_occurrences (int): the repeat unit must occur in the sequence at least this many times for the sequence to be counted
        min_fraction_bases_covered (float): the repeat unit must cover this fraction of the sequence for the sequence to be counted

    Returns:
        2-tuple: (the number of sequences that contain the repeat unit, the total number of sequences in the input iterator)
    """
    read_count = 0
    total = 0
    for sequence in sequences:
        total += 1
        count = sequence.count(repeat_unit)
        if count >= min_occurrences and count * len(repeat_unit)/len(sequence) >= min_fraction_bases_covered:
            read_count += 1

    return read_count, total


if __name__ == "__main__":
    main()