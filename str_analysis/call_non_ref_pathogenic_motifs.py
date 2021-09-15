#!/usr/bin/env python3

DESCRIPTION = """This script takes a bam or cram file and determines which motifs are present at known pathogenic loci 
(such as RFC1, BEAN1, DAB1, etc.) where several motifs are known to segregate in the population. It then optionally runs
 ExpansionHunterDenovo and/or ExpansionHunter on the detected motifs and gathers relevant fields from their outputs. 
Finally it outputs a json file per locus with collected information as well as a "call" field indicating whether 
pathogenic motifs were detected.
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
from str_analysis.utils.misc_utils import parse_interval

LOCUS_INFO = json.loads(pkgutil.get_data(__name__, "data/locus_info.json"))
OFFTARGET_REGIONS = json.loads(gzip.decompress(pkgutil.get_data(__name__, "data/offtarget_regions.json.gz")))

MARGIN = 7   # base pairs
FLANK_SIZE = 2000   # base pairs
MIN_MAPQ = 3
MIN_FRACTION_OF_BASES_COVERED = 0.7

NORMALIZE_TO_COVERAGE = 40

GENOME_VERSION_ALIASES = {
    "GRCh37": "37", "hg19": "37", "hg37": "37", "37": "37",
    "GRCh38": "38", "hg38": "38", "38": "38",
}


PATHOGENIC_MOTIF_CALL = "PATHOGENIC MOTIF"
BENIGN_MOTIF_CALL = "BENIGN MOTIF"
UNCERTAIN_SIG_CALL = "MOTIF OF UNCERTAIN SIGNIFICANCE"

NO_CALL1 = f"NO CALL (no motif has sufficient read support)"

PATHOGENIC_PATHOGENIC_CALL = f"{PATHOGENIC_MOTIF_CALL} / {PATHOGENIC_MOTIF_CALL}"
BENIGN_BENIGN_CALL = f"{BENIGN_MOTIF_CALL} / {BENIGN_MOTIF_CALL}"
UNCERTAIN_SIG_UNCERTAIN_SIG_CALL = f"{UNCERTAIN_SIG_CALL} / {UNCERTAIN_SIG_CALL}"

BENIGN_PATHOGENIC_CALL = f"{BENIGN_MOTIF_CALL} / {PATHOGENIC_MOTIF_CALL}"
PATHOGENIC_UNCERTAIN_SIG_CALL = f"{PATHOGENIC_MOTIF_CALL} / {UNCERTAIN_SIG_CALL}"
BENIGN_UNCERTAIN_SIG_CALL = f"{BENIGN_MOTIF_CALL} / {UNCERTAIN_SIG_CALL}"


def parse_args():
    p = argparse.ArgumentParser(description=DESCRIPTION)
    p.add_argument("-g", "--genome-version", choices=GENOME_VERSION_ALIASES.keys(), required=True)
    p.add_argument("-R", "--reference-fasta", help="Reference fasta path. The reference fasta is sometimes necessary "
        "for decoding cram files.", required=True)
    p.add_argument("-o", "--output-prefix", help="Output filename prefix")
    p.add_argument("-s", "--sample-id", help="The sample id to put in the output json file. If not specified, it "
        "will be retrieved from the bam/cram header or filename prefix.")

    p.add_argument("--run-expansion-hunter-denovo", action="store_true", help="Optionally run ExpansionHunterDenovo "
        "and copy information relevant to the locus from ExpansionHunterDenovo results to the output json.")
    p.add_argument("--expansion-hunter-denovo-path", help="The path of the ExpansionHunterDenovo executable to use "
        "when --expansion-hunter-denovo-path is specified.", default="ExpansionHunterDenovo")
    p.add_argument("--expansion-hunter-denovo-profile", help="Optionally copy information relevant to the locus "
        "from this ExpansionHunterDenovo profile to the output json. This is instead of --run-expansion-hunter-denovo")

    p.add_argument("-r", "--run-expansion-hunter", action="store_true", help="If this option is specified, this "
         "script will run ExpansionHunter once for each of the motif(s) it detects at the locus. "
         "ExpansionHunter doesn't currently support genotyping multiallelic repeats such as RFC1 where "
         "an individual may have 2 alleles with motifs that differ from eachother (and from the reference motif). "
         "Running ExpansionHunter separately for each motif provides a work-around.")
    p.add_argument("--expansion-hunter-path", help="The path of the ExpansionHunter executable to use if -r is "
        "specified. This must be ExpansionHunter v3 or v4.", default="ExpansionHunter")

    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument("--all-loci", action="store_true", help="Generate calls for all these loci: " + ", ".join(LOCUS_INFO.keys()))
    grp.add_argument("-l", "--locus", action="append", help="Call a subset of the known pathogenic loci with alternate motifs.",
        choices=LOCUS_INFO.keys())

    group = p.add_mutually_exclusive_group()
    group.add_argument("--run-reviewer", action="store_true", help="Run the REViewer tool to visualize "
        "ExpansionHunter output. --run-expansion-hunter must also be specified.")
    group.add_argument("--run-reviewer-for-pathogenic-calls", action="store_true", help="Run the REViewer tool "
        f"to visualize ExpansionHunter output when this script calls a sample as having {PATHOGENIC_PATHOGENIC_CALL}. "
        f"--run-expansion-hunter must also be specified.")

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

    if args.run_reviewer and not args.run_expansion_hunter:
        p.error("--run-expansion-hunter is required when --run-reviewer is used since REViewer depends on the output of ExpansionHunter")

    return args


def generate_variant_catalog(locus_id, repeat_unit, chrom, start_1based, end_1based, offtarget_regions=None):
    """Generate the ExpansionHunter variant catalog contents for a particular locus"""

    return {
        "LocusId": locus_id,
        "LocusStructure": f"({repeat_unit})*",
        "ReferenceRegion": f"{chrom}:{start_1based - 1}-{end_1based}",
        "VariantType": "RareRepeat" if offtarget_regions else "Repeat",
        "OfftargetRegions": offtarget_regions or [],
    }


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


def run_expansion_hunter(
    locus_id,
    locus_coords_1based,
    repeat_units,
    args,
    locus_results_json,
    run_reviewer=False,
    use_offtarget_regions=False,
):
    """Run ExpansionHunter and parse relevant output fields + add them to results.
    Then optionally run REViewer to generate read visualizations.

    Args:
        locus_id (str): locus id (eg. "RFC1")
        locus_coords_1based (str): locus coordinates (eg. "chr1:12345-54321")
        repeat_units (list): Repeat units to process at this locus.
        args (object): command-line arguments from argparse
        locus_results_json (dict): Results will be added to this dictionary
        run_reviewer (bool): Whether to run REViewer to generate read visualizations
    """

    if args.genome_version not in ("37", "38"):
        raise ValueError(f"Unexpected genome version: {args.genome_version}. Must be '37' or '38'")

    chrom, start_1based, end_1based = parse_interval(locus_coords_1based)

    for repeat_unit_number, repeat_unit in enumerate(repeat_units):
        motif_number = repeat_unit_number + 1 
        # generate variant catalog
        variant_catalog_locus_label = f"{locus_id}_{repeat_unit}"
        offtarget_regions = []
        if use_offtarget_regions:
            offtarget_regions = OFFTARGET_REGIONS[args.genome_version][repeat_unit]

        variant_catalog = generate_variant_catalog(
            variant_catalog_locus_label,
            repeat_unit, chrom, start_1based, end_1based,
            offtarget_regions=offtarget_regions)

        variant_catalog_path = f"{repeat_unit}.variant_catalog.json"
        with open(variant_catalog_path, "wt") as f:
            json.dump([variant_catalog], f)

        # run expansion hunter
        print("--"*10)
        print(f"Running ExpansionHunter on {args.sample_id} for repeat unit {repeat_unit}")
        if args.verbose:
            print("Using variant catalog: ")
            pprint(variant_catalog)

        output_prefix = f"{args.sample_id}.{locus_id}_{repeat_unit}.expansion_hunter"
        expansion_hunter_command = f"""{args.expansion_hunter_path} \
--sex male \
--reference {args.reference_fasta} \
--reads {args.bam_or_cram_path} \
--variant-catalog {variant_catalog_path} \
--output-prefix {output_prefix} \
--log-level debug
"""

        if args.verbose:
            print(f"Running: {expansion_hunter_command}")

        subprocess.run(expansion_hunter_command, shell=True, stderr=subprocess.STDOUT, check=False)

        # parse result
        if os.path.isfile(f"{output_prefix}.json"):
            with open(f"{output_prefix}.json", "rt") as f:
                expansion_hunter_output_json = json.load(f)
        else:
            print(f"ERROR: ExpansionHunter didn't produce a {output_prefix}.json file. Skipping...")
            continue

        if args.verbose:
            print(f"ExpansionHunter output: {pformat(expansion_hunter_output_json)}")

        eh_result = expansion_hunter_output_json.get("LocusResults", {}).get(variant_catalog_locus_label, {}).get(
            "Variants", {}).get(variant_catalog_locus_label, {})

        if not eh_result:
            continue

        locus_results_json[f"expansion_hunter_motif{motif_number}_json_output_file"] = f"{output_prefix}.json"
        locus_results_json[f"expansion_hunter_motif{motif_number}_repeat_unit"] = repeat_unit

        if eh_result.get("Genotype"):
            (
                locus_results_json[f"expansion_hunter_motif{motif_number}_short_allele_genotype"],
                locus_results_json[f"expansion_hunter_motif{motif_number}_long_allele_genotype"]
            ) = [
                int(g) for g in eh_result["Genotype"].split("/")]

        if eh_result.get("GenotypeConfidenceInterval"):
            (
                locus_results_json[f"expansion_hunter_motif{motif_number}_short_allele_CI_start"],
                locus_results_json[f"expansion_hunter_motif{motif_number}_short_allele_CI_end"],
                locus_results_json[f"expansion_hunter_motif{motif_number}_long_allele_CI_start"],
                locus_results_json[f"expansion_hunter_motif{motif_number}_long_allele_CI_end"],
            ) = [
                int(b) for ci in eh_result["GenotypeConfidenceInterval"].split("/") for b in ci.split("-")
            ]

            locus_results_json[f"expansion_hunter_motif{motif_number}_short_allele_CI_size"] = (
                    locus_results_json[f"expansion_hunter_motif{motif_number}_short_allele_CI_end"] -
                    locus_results_json[f"expansion_hunter_motif{motif_number}_short_allele_CI_start"]
            )

            locus_results_json[f"expansion_hunter_motif{motif_number}_long_allele_CI_size"] = (
                    locus_results_json[f"expansion_hunter_motif{motif_number}_long_allele_CI_end"] -
                    locus_results_json[f"expansion_hunter_motif{motif_number}_long_allele_CI_start"]
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

            locus_results_json[f"expansion_hunter_motif{motif_number}_total_{output_label}"] = total

        if run_reviewer:
            reviewer_command = f"""samtools sort {output_prefix}_realigned.bam -o {output_prefix}.sorted.bam \
&& samtools index {output_prefix}.sorted.bam \
&& REViewer --reads {output_prefix}.sorted.bam  \
    --vcf {output_prefix}.vcf \
    --reference {args.reference_fasta} \
    --catalog {variant_catalog_path} \
    --locus {variant_catalog_locus_label} \
    --output-prefix {output_prefix}_reviewer
"""
            if args.verbose:
                print(f"Running: {reviewer_command}")

            subprocess.run(reviewer_command, shell=True, stderr=subprocess.STDOUT, check=False)

            reviewer_output_filename = f"{output_prefix}_reviewer.{variant_catalog_locus_label}.svg"
            if os.path.isfile(reviewer_output_filename):
                locus_results_json[f"expansion_hunter_motif{motif_number}_reviewer_svg"] = reviewer_output_filename


def run_expansion_hunter_denovo(args):
    """Run ExpansionHunterDenovo.

    Arguments:
        args (object): command-line arguments from argparse

    Return:
        Path of ExpansionHunterDenovo *.str_profile.json output file.
    """

    print("--"*10)
    print(f"Running ExpansionHunterDenovo on {args.sample_id}")
    output_prefix = f"{args.sample_id}.expansion_hunter_denovo"
    expansion_hunter_denovo_command = f"""{args.expansion_hunter_denovo_path} profile \
--reference {args.reference_fasta} \
--reads {args.bam_or_cram_path} \
--output-prefix {output_prefix} \
"""

    if args.verbose:
        print(f"Running: {expansion_hunter_denovo_command}")

    subprocess.run(expansion_hunter_denovo_command, shell=True, stderr=subprocess.STDOUT, check=False)

    # parse result
    output_path = f"{output_prefix}.str_profile.json"
    if not os.path.isfile(output_path):
        print(f"ERROR: ExpansionHunterDenovo didn't produce a {output_path} file.")
        return None

    return output_path


def process_reads_in_locus(
        bam_or_cram_path,
        reference_fasta,
        locus_chrom,
        locus_start_0based,
        locus_end,
        motif_size):
    """
    Parses reads from the given bam/cram file that are within a window around the given STR locus. Then it computes
    a list that, for each well-aligned read that overlaps the locus, contains the subsequence of bases from that read
    that fall between the locus start and end coordinates. Also, it returns the total number of well-aligned read bases
    to the left and to the right of the STR locus.

    Args:
        bam_or_cram_path (str): bam or cram path
        reference_fasta (str): reference fasta path
        locus_chrom (str): chromosome name of STR locus
        locus_start_0based (int):  0-based start coordinate of STR locus
        locus_end (int): end coordinate of STR locus
        motif_size (int): motif size

    Return:
          overlapping_sequences, left_flank_n_well_aligned_bases, right_flank_n_well_aligned_bases
    """
    with pysam.Samfile(bam_or_cram_path, reference_filename=reference_fasta) as f:

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

        # get all sequences that overlap the locus (regardless of whether they're soft-clipped)
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
            if len(relevant_bases) >= motif_size:
                overlapping_sequences.append(relevant_bases)

    return overlapping_sequences, left_flank_n_well_aligned_bases, right_flank_n_well_aligned_bases


def compute_final_call(n_total_well_supported_motifs, n_pathogenic_motifs, n_benign_motifs):
    """Decide which combination of motifs is supported by the data.

    NOTE: there's no attempt to determine the size of the expansion and whether it's in the pathogenic range.

    Return:
        str: a string such as "BENIGN MOTIF / PATHOGENIC MOTIF" indicating which motifs were detected
    """

    if n_total_well_supported_motifs == 0:
        final_call = NO_CALL1
    elif n_pathogenic_motifs == n_total_well_supported_motifs:
        # reads support only known pathogenic motif(s)
        final_call = PATHOGENIC_PATHOGENIC_CALL
    elif n_benign_motifs == n_total_well_supported_motifs:
        # reads support only known benign motif(s)
        final_call = BENIGN_BENIGN_CALL
    elif n_benign_motifs == 0 and n_pathogenic_motifs == 0:
        # reads support one or more non-reference motifs of unknown significance
        final_call = UNCERTAIN_SIG_UNCERTAIN_SIG_CALL
    elif n_benign_motifs > 0 and n_pathogenic_motifs > 0:
        # reads support one known benign motif and one pathogenic motif
        final_call = BENIGN_PATHOGENIC_CALL
    elif n_pathogenic_motifs > 0:
        # reads support one pathogenic motif and at least one other motif of unknown significance
        final_call = PATHOGENIC_UNCERTAIN_SIG_CALL
    elif n_benign_motifs > 0:
        # reads support one known benign motif and at least one other motif of unknown significance
        final_call = BENIGN_UNCERTAIN_SIG_CALL
    else:
        raise Exception(f"Unexpected state when n_total_well_supported_motifs={n_total_well_supported_motifs}, "
                        f"n_pathogenic_motifs={n_pathogenic_motifs}, n_benign_motifs={n_benign_motifs}")

    return final_call


def process_offtarget_regions(motif, motif_number, read_count, flank_coverage_mean, args, locus_results_json):
    """Look for off-target reads for a specific motif.

    Args:
        motif (str): Repeat unit
        motif_number (int): Motif number relative to other selected motifs.
        read_count (int): Not-offtarget read count.
        flank_coverage_mean (float): mean coverage in flanking regions
        args (object): command-line arguments from argparse.
        locus_results_json (dict): Results will be added to this dictionary.
    """
    offtarget_regions = OFFTARGET_REGIONS[args.genome_version].get(motif)
    if offtarget_regions is None:
        print(f"WARNING: off-target regions not available for {motif}")
        return

    read_count_with_offtargets = read_count
    with pysam.Samfile(args.bam_or_cram_path, reference_filename=args.reference_fasta) as f:
        for offtarget_region in offtarget_regions:
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

    locus_results_json.update({
        f"motif{motif_number}_read_count_with_offtargets": read_count_with_offtargets,
        f"motif{motif_number}_normalized_read_count_with_offtargets":
            read_count_with_offtargets * NORMALIZE_TO_COVERAGE / flank_coverage_mean if flank_coverage_mean > 0 else 0,
    })


def process_locus(locus_id, args):
    """Compute results for a single locus and write them to a json file.

    Args:
        locus_id (str): locus id
        args (object): command-line args from argparse
    """
    locus_coords_1based = LOCUS_INFO[locus_id]["LocusCoords_1based"][args.genome_version]
    locus_chrom, start_1based, end_1based = parse_interval(locus_coords_1based)
    locus_start_0based = start_1based - 1
    locus_end = end_1based
    use_offtarget_regions = LOCUS_INFO[locus_id]["UseOfftargetRegions"]

    known_pathogenic_motifs = list(map(compute_canonical_repeat_unit, LOCUS_INFO[locus_id]["Motifs"]["PATHOGENIC"]))
    known_benign_motifs = list(map(compute_canonical_repeat_unit, LOCUS_INFO[locus_id]["Motifs"]["BENIGN"]))

    pathogenic_motif_size = len(known_pathogenic_motifs[0])

    # process bam/cram
    overlapping_sequences, left_flank_n_well_aligned_bases, right_flank_n_well_aligned_bases = process_reads_in_locus(
        args.bam_or_cram_path, args.reference_fasta, locus_chrom, locus_start_0based, locus_end, pathogenic_motif_size)

    locus_results_json = {}
    left_flank_coverage = left_flank_n_well_aligned_bases / FLANK_SIZE
    right_flank_coverage = right_flank_n_well_aligned_bases / FLANK_SIZE
    locus_results_json.update({
        "sample_id": args.sample_id,
        "locus_id": locus_id,
        "locus_coords": f"{locus_chrom}:{start_1based}-{end_1based}",
        "genome_version": args.genome_version,
        "left_flank_coverage": left_flank_coverage,
        "right_flank_coverage": right_flank_coverage,
    })

    # compute the motif(s) found in the soft-clipped reads, and how many times each one occurs
    motif_to_n_occurrences = collections.defaultdict(int)
    motif_to_read_count = collections.defaultdict(int)
    if locus_id == "RFC1":
        # in gnomAD, EHdn sometimes finds 6bp repeat units (eg. AAAGGG), so check for those as well
        motif_sizes_to_check = [pathogenic_motif_size, pathogenic_motif_size + 1]
    else:
        motif_sizes_to_check = [pathogenic_motif_size]

    for overlapping_sequence in overlapping_sequences:
        for motif_size in motif_sizes_to_check:
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
                          f"overlap the {locus_id} locus: {overlapping_sequence}")
                else:
                    if motif_size == pathogenic_motif_size:
                        print(f"Didn't find a consistent {motif_size}bp repeat unit in read bases "
                              f"that overlap the {locus_id} locus: {overlapping_sequence}")

            if motif is not None:
                canonical_motif = compute_canonical_repeat_unit(motif)
                motif_to_read_count[canonical_motif] += 1
                motif_to_n_occurrences[canonical_motif] += count

    locus_results_json.update({
        "found_n_reads_overlap_the_locus": len(overlapping_sequences),
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
        1 if motif in known_benign_motifs else
        2 if motif in known_pathogenic_motifs else
        3)

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
        if motif in known_pathogenic_motifs:
            n_pathogenic_motifs += 1
        elif motif in known_benign_motifs:
            n_benign_motifs += 1

        read_count = motif_to_read_count.get(motif)
        n_occurrences = motif_to_n_occurrences.get(motif)

        locus_results_json.update({
            f"motif{motif_number}_repeat_unit": motif,
            f"motif{motif_number}_read_count": read_count,
            f"motif{motif_number}_normalized_read_count":
                read_count * NORMALIZE_TO_COVERAGE / flank_coverage_mean if flank_coverage_mean > 0 else 0,
            f"motif{motif_number}_n_occurrences": n_occurrences,
        })

        if use_offtarget_regions:
            process_offtarget_regions(motif, motif_number, read_count, flank_coverage_mean, args, locus_results_json)

    final_call = compute_final_call(n_total_well_supported_motifs, n_pathogenic_motifs, n_benign_motifs)

    locus_results_json.update({
        "n_total_well_supported_motifs": n_total_well_supported_motifs,
        "n_benign_motifs": n_benign_motifs,
        "n_pathogenic_motifs": n_pathogenic_motifs,
        "call": final_call,
    })

    print(f"Final call: {final_call}")

    # run ExpansionHunter if requested
    if args.run_expansion_hunter:
        run_expansion_hunter(
            locus_id,
            locus_coords_1based,
            selected_motifs,
            args,
            locus_results_json,
            run_reviewer=args.run_reviewer or (
                args.run_reviewer_for_pathogenic_calls and final_call == PATHOGENIC_PATHOGENIC_CALL
            ),
            use_offtarget_regions=use_offtarget_regions,
        )

    # process EHdn profile if one was provided
    if args.expansion_hunter_denovo_profile:
        print(f"Parsing {args.expansion_hunter_denovo_profile}")
        open_func = gzip.open if args.expansion_hunter_denovo_profile.endswith("gz") else open
        with open_func(args.expansion_hunter_denovo_profile, "rt") as f:
            expansion_hunter_denovo_json = json.load(f)

        records, sample_read_depth, _ = parse_ehdn_info_for_locus(
            expansion_hunter_denovo_json, locus_chrom, locus_start_0based, locus_end)

        #locus_results_json[f"expansion_hunter_denovo_profile"] = args.expansion_hunter_denovo_profile
        locus_results_json["ehdn_sample_read_depth"] = sample_read_depth

        # get the 2 motifs with the most read support
        records.sort(key=lambda r: (
            -r["anchored_irr_count_for_this_repeat_unit_and_region"],
            -r["total_irr_count_for_this_repeat_unit_and_region"],
            r["repeat_unit"],
        ))

        for i in 0, 1:
            if i >= len(records):
                continue

            record = records[i]
            motif_number = i + 1
            locus_results_json.update({
                f"ehdn_motif{motif_number}_repeat_unit": record.get("repeat_unit"),
                f"ehdn_motif{motif_number}_anchored_irr_count": record.get("anchored_irr_count_for_this_repeat_unit_and_region"),
                f"ehdn_motif{motif_number}_n_anchored_regions": record.get("n_anchored_regions_for_this_repeat_unit"),
                f"ehdn_motif{motif_number}_paired_irr_count": record.get("paired_irr_count_for_this_repeat_unit"),
                f"ehdn_motif{motif_number}_total_irr_count": record.get("total_irr_count_for_this_repeat_unit_and_region"),
            })

    # generate output json
    output_filename = f"{args.output_prefix}.{locus_id}_motifs.json"
    with open(output_filename, "wt") as f:
        json.dump(locus_results_json, f, indent=2)
    print(f"Wrote results to {output_filename}")
    pprint(locus_results_json)


def main():
    args = parse_args()

    bam_cram_prefix = re.sub(".bam$|.cram$", "", os.path.basename(args.bam_or_cram_path))
    args.output_prefix = args.output_prefix or bam_cram_prefix

    # process bam/cram
    if not args.sample_id:
        # try to get sample id from bam/cram header
        with pysam.Samfile(args.bam_or_cram_path, reference_filename=args.reference_fasta) as f:

            read_groups = f.header.as_dict().get("RG", [])
            if read_groups:
                args.sample_id = read_groups[0].get("SM")
                if args.sample_id:
                    print(f"Using sample id '{args.sample_id}' from the bam/cram header")

            if not args.sample_id:
                args.sample_id = bam_cram_prefix
                print(f"Using sample id '{args.sample_id}' based on the input filename prefix")

    if args.run_expansion_hunter_denovo:
        args.expansion_hunter_denovo_profile = run_expansion_hunter_denovo(args)

    if args.locus:
        loci = args.locus
    elif args.all_loci:
        loci = LOCUS_INFO.keys()
    else:
        raise ValueError("Must specify --locus or --all-loci")

    for locus_id in loci:
        print(f"Processing {locus_id} in {args.bam_or_cram_path}")
        process_locus(locus_id, args)


if __name__ == "__main__":
    main()