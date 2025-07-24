#!/usr/bin/env python3

"""
This script takes a VCF (either single-sample or multi-sample) and filters it to the subset of insertions and deletions
that represent tandem repeat (TR) expansions or contractions. It then extends these repeats into their immediate
left and right flanking regions in the reference. This yields a set of tandem repeat loci that can be used for
further analysis, such as genotyping.

This is the next iteration of the filter_vcf_to_STR_variants.py script. It adds the ability to run TandemRepeatFinder to
discover more imperfect repeats (particularly VNTRs) and separates tandem repeat locus discovery from genotyping
(which is now done by the filter_vcf_to_genotype_tandem_repeats.py script).
"""

import argparse
import collections
import datetime
import gzip
import intervaltree
import math
import os
from pprint import pformat
import pyfaidx
import pysam
import re
import shutil
import threading
import tqdm

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif, reverse_complement
from str_analysis.utils.find_repeat_unit import find_repeat_unit_allowing_interruptions
from str_analysis.utils.find_repeat_unit import find_repeat_unit_without_allowing_interruptions
from str_analysis.utils.find_repeat_unit import extend_repeat_into_sequence_allowing_interruptions
from str_analysis.utils.find_repeat_unit import extend_repeat_into_sequence_without_allowing_interruptions
from str_analysis.utils.file_utils import open_file
from str_analysis.utils.trf_runner import TRFRunner


DETECTION_MODE_PURE_REPEATS = "pure_repeats"
DETECTION_MODE_ALLOW_INTERRUPTIONS = "allow_interruptions"
DETECTION_MODE_TRF = "trf"

DETECTION_MODE_ORDER = [
    DETECTION_MODE_PURE_REPEATS,
    DETECTION_MODE_ALLOW_INTERRUPTIONS,
    DETECTION_MODE_TRF,
]

CURRENT_TIMESTAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S.%f")
TRF_WORKING_DIR = f"trf_working_dir/{CURRENT_TIMESTAMP}"

TSV_OUTPUT_COLUMNS = [
    "Chrom",
    "Start0Based",
    "End1Based",
    "Locus",
    "LocusId",
    "INS_or_DEL",
    "Motif",
    "MotifInterruptionIndices",
    "CanonicalMotif",
    "MotifSize",
    "NumRepeatsInReference",
    "VcfPos",
    #"VcfRef",
    #"VcfAlt",
    "VcfGenotype",
    "SummaryString",
    "IsFoundInReference",
    "IsPureRepeat",
    "IsMultiallelic",
    "DetectionMode",
]

MAX_INDEL_SIZE = 100_000  # bp
TOO_MANY_Ns = "N"*100

FILTER_MORE_THAN_TWO_ALT_ALLELES = "more than two alt alleles"
FILTER_UNEXPECTED_GENOTYPE_FORMAT = "unexpected genotype format"
FILTER_ZERO_ALT_ALLELES = "variant has zero non-* alt alleles"
#FILTER_MULTIALLELIC_VARIANT_WITH_HOMOZYGOUS_GENOTYPE = "multiallelic variant with homozygous genotype"

FILTER_ALLELE_WITH_N_BASES = "contains Ns in the variant sequence"
FILTER_ALLELE_WITH_TOO_MANY_Ns_IN_FLANKS = "contains too many Ns in the flanks"
FILTER_ALLELE_TR_SPANS_TOO_MANY_BASES = "spans too many bases"
FILTER_ALLELE_SNV_OR_MNV = "SNV/MNV"
FILTER_ALLELE_MNV_INDEL = "complex multinucleotide indel"
FILTER_ALLELE_INDEL_WITHOUT_REPEATS = "INDEL without repeats"
FILTER_ALLELE_TOO_BIG = f"INDEL > {MAX_INDEL_SIZE}bp"
FILTER_TR_ALLELE_NOT_ENOUGH_REPEATS = "contains only %d full repeats"
FILTER_TR_ALLELE_DOESNT_SPAN_ENOUGH_BASE_PAIRS = "spans < %d bp"

#FILTER_TR_ALLELE_PARTIAL_REPEAT = "ends in partial repeat"

FILTER_TR_ALLELE_REPEAT_UNIT_TOO_SHORT = "repeat unit < %d bp"
FILTER_TR_ALLELE_REPEAT_UNIT_TOO_LONG = "repeat unit > %d bp"

FILTER_VARIANT_WITH_TR_ALLELES_WITH_DIFFERENT_MOTIFS = "tandem repeat alleles with different motifs"
FILTER_VARIANT_WITH_TR_ALLELES_WITH_DIFFERENT_INTERRUPTION_PATTERNS = "tandem repeat alleles with different interruption patterns"
#FILTER_VARIANT_WITH_TR_ALLELES_WITH_DIFFERENT_COORDS = "TR alleles with different coords"
FILTER_DUPLICATE_TR_LOCUS = "locus already detected due to a previous variant in the VCF"
FILTER_TR_LOCUS_THAT_IS_BETTER_REPRESENTED_BY_AN_OVERLAPPING_TR = "locus overlaps and can be better represented by another TR locus"
FILTER_TR_LOCUS_THAT_HAS_OVERLAPPING_VARIANTS = "locus overlaps more than one variant"

WILL_RUN_TRF_ON_THIS_ALLELE_IN_2ND_PASS = "will run TRF on this allele during the 2nd pass"

def parse_args():
    """Parse command-line arguments."""

    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument("-R", "--reference-fasta-path", help="Reference genome fasta path.", required=True)
    p.add_argument("--dont-allow-interruptions", action="store_true", help="Only detect perfect repeats. This implicitly "
                   "also enables --dont-run-trf since detection of pure repeats does not require running TandemRepeatFinder (TRF).")
    p.add_argument("--dont-run-trf", action="store_true", help="Don't use TandemRepeatFinder (TRF) to help detect imperfect TRs. Instead, only "
                   "use the simpler algorithm that allows one position within the repeat unit to vary across repeats.")
    p.add_argument("--trf-executable-path", help="Path to the TandemRepeatFinder (TRF) executable. This is "
                   "required unless --dont-run-trf is specified.")
    p.add_argument("-t", "--trf-threads", default=12, type=int, help="Number of TandemRepeatFinder (TRF) instances to run in parallel.")
    p.add_argument("--min-indel-size-to-run-trf", default=7, type=int, help="Only run TandemRepeatFinder (TRF) "
                    "on insertions and deletions that are at least this many base pairs.")

    p.add_argument("--min-tandem-repeat-length", type=int, default=9, help="Only detect tandem repeat variants that are at least this long (in base pairs). "
                   "This threshold will be applied to the total repeat sequence including any repeats in the flanking sequence to the left "
                   "and right of the variant in addition to the inserted or deleted bases")
    p.add_argument("--min-repeats", type=int, default=3, help="Only detect tandem repeat loci that consist of at least this many repeats. "
                   "This threshold will be applied to the total repeat sequence including any repeats in the flanking sequence to "
                   "the left and right of the variant in addition to the inserted or deleted bases")
    p.add_argument("--min-repeat-unit-length", type=int, default=1, help="Minimum repeat unit length in base pairs.")
    p.add_argument("--max-repeat-unit-length", type=int, default=10**9, help="Max repeat unit length in base pairs.")
    p.add_argument("--show-progress-bar", help="Show a progress bar in the terminal when processing variants.",
                   action="store_true")
    p.add_argument("-v", "--verbose", help="Print detailed logs.", action="store_true")
    p.add_argument( "--debug", help="Print any debugging info and don't delete intermediate files.", action="store_true")

    p.add_argument("--offset", default=0, type=int, help="Skip the first N variants in the VCF file. This is useful for testing ")
    p.add_argument("-n", type=int, help="Only process N rows from the VCF (after applying --offset). Useful for testing.")

    p.add_argument("-o", "--output-prefix", help="Output file prefix. If not specified, it will be computed based on "
                   "the input vcf filename")
    
    p.add_argument("-info", "--copy-info-field-keys-to-tsv", help="Optionally specify the name of an INFO fields from the "
                   "input VCF that should be copied to a separate column in the output TSV file. This option can be " 
                   "specified more than once.", action="append")

    p.add_argument("--write-vcf-file", help="Output a VCF file with all variants that were found to be TRs.", action="store_true")
    p.add_argument("--write-vcf-with-filtered-out-variants", help="Output a VCF file with variants where one allele was found to be an TR, "
                   "but that were still filtered out for reasons such as being multiallelic and having alleles with different motifs, "
                   "or because one allele was an TR while the other was an SNV. These types of variants are filtered out to reduce complexity "
                   "in downstream analyses.",
                   action="store_true")
    p.add_argument("--write-fasta", help="Output a FASTA file containing all TR alleles", action="store_true")

    p.add_argument("-L", "--interval", help="Only process variants in this genomic interval (format: chrom:start-end)",
                   action="append")
    
    p.add_argument("input_vcf_path", help="Input single-sample VCF file path. This script was designed and tested on VCFs produced by DipCall, "
                   "but should work with any single-sample VCF.")

    args = p.parse_args()

    if not args.dont_run_trf and not args.trf_executable_path:
        p.error(f"Must specify --trf-executable-path or --dont-run-trf")

    if args.dont_allow_interruptions:
        # drop some output columns that are only relevant for interrupted repeats
        for column in "NumPureRepeats", "PureRepeatSize (bp)", "FractionPureRepeats", "MotifInterruptionIndices":
            if column in TSV_OUTPUT_COLUMNS:
                TSV_OUTPUT_COLUMNS.remove(column)

    if args.copy_info_field_keys_to_tsv:
        TSV_OUTPUT_COLUMNS.extend(args.copy_info_field_keys_to_tsv)

    if args.write_vcf_with_filtered_out_variants and not args.write_vcf_file:
        print("NOTE: enabling --write-vcf because it's necessary for --write-vcf-with-filtered-out-variants")
        args.write_vcf_file = True

    args.write_bed_file = True

    return args


def get_trf_working_dir(args):
    input_vcf_prefix = re.sub(".vcf(.gz)?$", "", os.path.basename(args.input_vcf_path))
    return os.path.join(TRF_WORKING_DIR, f"{input_vcf_prefix}.h{abs(hash(args.input_vcf_path) % 10**9)}")


def tandem_repeat_allele_failed_filters(args, repeat_unit, total_repeats, counters):
    """Check if the given tandem repeat allele (represented by its repeat_unit and total_repeats) passes the filters
    specified in the command-line arguments.

    Args:
        args (argparse.Namespace): command-line arguments parsed by parse_args()
        repeat_unit (str): The repeat unit of the tandem repeat allele.
        total_repeats (int): The total number of repeats in the allele, including repeats in the flanking sequences.
        counters (dict): Dictionary of counters to collect summary stats about the number of TR variants found, etc.

    Return:
        str: A string describing the reason why the allele failed filters, or None if the allele passed all filters.
    """
    if total_repeats == 1:
        # no repeat unit found in this allele
        counters[f"allele filter: {FILTER_ALLELE_INDEL_WITHOUT_REPEATS}"] += 1
        return FILTER_ALLELE_INDEL_WITHOUT_REPEATS
    elif total_repeats < args.min_repeats:
        counters[f"allele filter: allele includes only {int(total_repeats)} full repeat(s)"] += 1
        return FILTER_TR_ALLELE_NOT_ENOUGH_REPEATS % int(total_repeats)
    elif total_repeats * len(repeat_unit) < args.min_tandem_repeat_length:
        # it has more than one repeat, so a repeat unit was found, but it was less than the minimum threshold
        counters[f"allele filter: allele sequence spans < {args.min_tandem_repeat_length}bp"] += 1
        return FILTER_TR_ALLELE_DOESNT_SPAN_ENOUGH_BASE_PAIRS % args.min_tandem_repeat_length
    elif len(repeat_unit) < args.min_repeat_unit_length:
        counters[f"allele filter: repeat unit is shorter than {args.min_repeat_unit_length}bp"] += 1
        return FILTER_TR_ALLELE_REPEAT_UNIT_TOO_SHORT % args.min_repeat_unit_length
    elif len(repeat_unit) > args.max_repeat_unit_length:
        counters[f"allele filter: repeat unit is longer than {args.max_repeat_unit_length}bp"] += 1
        return FILTER_TR_ALLELE_REPEAT_UNIT_TOO_LONG % args.max_repeat_unit_length

    return None  # did not fail filters


def check_if_allele_is_tandem_repeat(
    vcf_chrom,
    vcf_pos,
    vcf_ref,
    alt_allele,
    variant_bases,
    left_flanking_reference_sequence,
    right_flanking_reference_sequence,
    left_flank_end,
    right_flank_start_1based,
    args,
    counters,
    only_generate_trf_fasta,
    variants_to_process_using_trf,
    variants_per_trf_fasta,
    detection_mode=DETECTION_MODE_PURE_REPEATS,
    rerun_trf_on_longer_flanking_sequences=False,
):
    """Determine if the given allele is a tandem repeat expansion or contraction or neither.

    Args:
        vcf_chrom (str): chromosome
        vcf_pos (int): 1-based position
        vcf_ref (str): ref sequence
        alt_allele (str): alt allele sequence
        variant_bases (str): bases that were inserted or deleted in the variant allele
        left_flanking_reference_sequence (str): reference sequence to the left of the variant allele
        right_flanking_reference_sequence (str): reference sequence to the right of the variant allele
        left_flank_end (int): 1-based end position of the left flanking sequence
        right_flank_start_1based (int): 1-based start position of the right flanking sequence
        args (argparse.Namespace): command-line arguments parsed by parse_args()
        counters (dict): Dictionary of counters to collect summary stats about the number of TR variants found, etc.
        only_generate_trf_fasta (bool): if True, only generate the TRF fasta file for this variant and return None.
        variants_to_process_using_trf (dict): maps variant IDs that should be processed using TRF to their sequence
            number in the TRF fasta file.
        variants_per_trf_fasta (dict): maps TRF fasta input file path to the total number of variants in that fasta file.
        detection_mode (str): Should be either DETECTION_MODE_PURE_REPEATS, DETECTION_MODE_ALLOW_INTERRUPTIONS,
            DETECTION_MODE_TRF, or DETECTION_MODE_TRF_RERUN:

            DETECTION_MODE_PURE_REPEATS will only detect perfect repeats.
            DETECTION_MODE_ALLOW_INTERRUPTIONS will detect repeats where one position can vary across repeats (similar to
                known disease-associated loci like RUNX2).
            DETECTION_MODE_TRF will use TRF to detect more complicated imperfect repeats such as those often seen in VNTRs
        rerun_trf_on_longer_flanking_sequences (bool): if True, use a special fasta filename to process a single variant
    Return:
        2-tuple (dict, str):
            dict: if the allele represents a tandem repeat, this will be a JSON record describing the tandem repeat, 
                otherwise it will be None. The JSON record will contain the following keys:
            - "Chrom": chromosome name
            - "DetectionMode": either DETECTION_MODE_PURE_REPEATS or DETECTION_MODE_ALLOW_INTERRUPTIONS
            - "RepeatUnit": the repeat unit of the tandem repeat allele
            - "MotifInterruptionIndices": indices of interruptions in the repeat unit within the variant bases, or None if no interruptions were found
            - "IsPureRepeat": True if the allele is a pure repeat (no interruptions), False otherwise
            - "Start0Based": 0-based start position of the tandem repeat allele in the reference genome
            - "End1Based": 1-based end position of the tandem repeat allele in the reference genome
            - "NumRepeatsRef": total number of repeats in the reference sequence
            - "NumRepeatsAlt": total number of repeats in the alt allele, including flanking sequences
            - "AlleleRepeatSequence": the sequence of the tandem repeat allele
            
            str: if the allele is not a tandem repeat, this will be a string describing the reason why the allele failed tandem repeat filters, or otherwise None if it passed all filters.
    """

    result = {
        "Chrom": vcf_chrom,
        "DetectionMode": detection_mode,
    }

    if detection_mode == DETECTION_MODE_PURE_REPEATS:
        # check whether this variant allele + flanking sequences represent a pure tandem repeat expansion or contraction
        (
            repeat_unit,
            num_total_repeats_in_variant_bases,
            _,
        ) = find_repeat_unit_without_allowing_interruptions(variant_bases)

        num_total_repeats_left_flank = extend_repeat_into_sequence_without_allowing_interruptions(
            repeat_unit[::-1],
            left_flanking_reference_sequence[::-1])
        num_total_repeats_right_flank = extend_repeat_into_sequence_without_allowing_interruptions(
            repeat_unit,
            right_flanking_reference_sequence)

        result["RepeatUnit"] = repeat_unit
        result["MotifInterruptionIndices"] = None
        result["IsPureRepeat"] = True

        tandem_repeat_bases_in_left_flank = num_total_repeats_left_flank * len(repeat_unit)
        tandem_repeat_bases_in_right_flank = num_total_repeats_right_flank * len(repeat_unit)

    elif detection_mode == DETECTION_MODE_ALLOW_INTERRUPTIONS:
        (
            repeat_unit,
            num_pure_repeats_in_variant_bases,
            num_total_repeats_in_variant_bases,
            repeat_unit_interruption_index,
            _
        ) = find_repeat_unit_allowing_interruptions(variant_bases, allow_partial_repeats=False)

        reversed_repeat_unit_interruption_index = None
        if repeat_unit_interruption_index is not None:
            reversed_repeat_unit_interruption_index = (len(repeat_unit) - 1 - repeat_unit_interruption_index)

        num_pure_repeats_left_flank, num_total_repeats_left_flank, reversed_repeat_unit_interruption_index = extend_repeat_into_sequence_allowing_interruptions(
            repeat_unit[::-1],
            left_flanking_reference_sequence[::-1],
            repeat_unit_interruption_index=reversed_repeat_unit_interruption_index)

        if reversed_repeat_unit_interruption_index is not None:
            # reverse the repeat_unit_interruption_index
            repeat_unit_interruption_index = len(repeat_unit) - 1 - reversed_repeat_unit_interruption_index

        num_pure_repeats_right_flank, num_total_repeats_right_flank, repeat_unit_interruption_index = extend_repeat_into_sequence_allowing_interruptions(
            repeat_unit,
            right_flanking_reference_sequence,
            repeat_unit_interruption_index=repeat_unit_interruption_index)

        result["RepeatUnit"] = repeat_unit
        result["MotifInterruptionIndices"] = repeat_unit_interruption_index
        result["IsPureRepeat"] = False

        tandem_repeat_bases_in_left_flank = num_total_repeats_left_flank * len(repeat_unit)
        tandem_repeat_bases_in_right_flank = num_total_repeats_right_flank * len(repeat_unit)

    elif detection_mode == DETECTION_MODE_TRF:

        variant_id = f"{vcf_chrom}-{vcf_pos}-"
        variant_id += vcf_ref if len(vcf_ref) < 24 else (vcf_ref[0] + ".." + str(len(vcf_ref)) + ".." + vcf_ref[-1])
        variant_id += "-"
        variant_id += alt_allele if len(alt_allele) < 24 else (alt_allele[0] + ".." + str(len(alt_allele)) + ".." + alt_allele[-1])
        variant_id += "-h" + str(abs(hash(f"{vcf_chrom}-{vcf_pos}-{vcf_ref}-{alt_allele}")))  # add a hash to make the ID unique

        # determine which thread will process (or has processed) this variant
        if rerun_trf_on_longer_flanking_sequences:
            trf_thread_index = 0
            trf_fasta_filename = (
                f"{variant_id}_"
                f"and_left_{len(left_flanking_reference_sequence)}bp_"
                f"and_right_{len(right_flanking_reference_sequence)}bp.fa")
        else:
            if variant_id not in variants_to_process_using_trf:
                trf_thread_index = len(variants_to_process_using_trf) % args.trf_threads
            else:
                trf_thread_index, _ = variants_to_process_using_trf[variant_id]
            trf_fasta_filename = f"thread{trf_thread_index}.fa"

        trf_fasta_path = os.path.join(get_trf_working_dir(args), trf_fasta_filename)

        if only_generate_trf_fasta:
            if variant_id in variants_to_process_using_trf:
                raise ValueError(f"variant ID '{variant_id}' already exists in variants_to_process_using_trf")

            with open(trf_fasta_path, "at") as f:
                # reverse the left flanking sequence + variant bases so that TRF starts detecting repeats from the
                # end of the variant sequence rather than some random point in the flanking region. This will
                # make it so the repeat unit is detected in the correct orientation (after being reversed back).
                # For example, if the variant was T > TCAGCAGCAGCAG , we want TRF to detect the motif as CAG rather
                # than AGC of GCA.

                left_flank_and_variant_bases_reversed = f"{left_flanking_reference_sequence}{variant_bases}"[::-1]
                f.write(f">{variant_id}$left\n")
                f.write(f"{left_flank_and_variant_bases_reversed}\n")

                variant_bases_and_right_flank = f"{variant_bases}{right_flanking_reference_sequence}"
                f.write(f">{variant_id}$right\n")
                f.write(f"{variant_bases_and_right_flank}\n")

            trf_fasta_sequence_number = variants_per_trf_fasta[trf_fasta_path]
            variants_to_process_using_trf[variant_id] = (trf_thread_index, trf_fasta_sequence_number)
            variants_per_trf_fasta[trf_fasta_path] += 1

            return None, WILL_RUN_TRF_ON_THIS_ALLELE_IN_2ND_PASS

        # this is the 2nd pass where we need to parse the TRF output
        if variant_id not in variants_to_process_using_trf:
            raise ValueError(f"variant ID '{variant_id}' not found in variants_to_process_using_trf")


        # check that the results overlap with the variant bases
        trf_runner = TRFRunner(args.trf_executable_path, html_mode=True,
                       min_motif_size=args.min_repeat_unit_length,
                       max_motif_size=args.max_repeat_unit_length)

        if rerun_trf_on_longer_flanking_sequences:
            trf_runner.run_trf_on_fasta_file(trf_fasta_path)

        _, trf_fasta_sequence_number = variants_to_process_using_trf[variant_id]

        motif_size_to_matching_trf_results = collections.defaultdict(dict)

        for left_or_right in "left", "right":
            trf_results = trf_runner.parse_html_results(
                trf_fasta_path,
                sequence_number=2 * trf_fasta_sequence_number + (0 if left_or_right == "left" else 1) + 1,
                total_sequences=2 * variants_per_trf_fasta[trf_fasta_path])

            if not args.debug:
                # delete the TRF output html files to save space
                files_to_delete = set()
                for trf_result in trf_results:
                    for trf_output_file_name in "html_file_path", "html_file2_path", "html_summary_file_path":
                        files_to_delete.add(trf_result[trf_output_file_name])
                for file_to_delete in files_to_delete:
                    if os.path.isfile(file_to_delete):
                        os.remove(file_to_delete)

            for trf_result in trf_results:
                if trf_result["sequence_name"] != f"{variant_id}${left_or_right}":
                    raise ValueError(f"TRF result sequence name '{trf_result['sequence_name']}' does not match the expected variant ID '{variant_id}${left_or_right}'")

                fuzz = int(math.log10(trf_result["repeat_unit_length"]) * 2) # for large motifs, allow 1 or 2 bases of fuzziness in repeat boundaries
                if trf_result["start_0based"] > fuzz or trf_result["end_1based"] < len(variant_bases) - fuzz:
                    # repeat must start at or very close to the start of the variant bases and end at or after the end of the variant bases
                    continue

                # check that repeat boundary matches the variant breakpoint
                tandem_repeat_bases = trf_result["start_0based"]
                repeat_count = 0
                while repeat_count < len(trf_result["repeats"]) and tandem_repeat_bases < len(variant_bases):
                    current_repeat = trf_result["repeats"][repeat_count].replace("-", "")
                    tandem_repeat_bases += len(current_repeat)
                    repeat_count += 1

                if tandem_repeat_bases + fuzz < len(variant_bases) or tandem_repeat_bases - fuzz > len(variant_bases):
                    continue

                trf_result["num_repeats_in_variant"] = repeat_count
                motif_size_to_matching_trf_results[trf_result["repeat_unit_length"]][left_or_right] = trf_result

        if len(motif_size_to_matching_trf_results) == 0:
            return None, FILTER_ALLELE_INDEL_WITHOUT_REPEATS

        # get the entry with the smallest motif size
        best_trf_results = None
        for motif_size, matching_trf_results in sorted(motif_size_to_matching_trf_results.items()):
            best_motif_size = motif_size
            best_trf_results = matching_trf_results
            break

        # if the smallest motif size is 2 or less, then pick the entry with the largest alignment score instead
        if best_motif_size <= 2 and len(motif_size_to_matching_trf_results) > 1:
            # get the entry that has the max alignment score
            max_alignment_score = 0
            for motif_size, matching_trf_results in motif_size_to_matching_trf_results.items():
                alignment_score = 0
                for _trf_result in matching_trf_results.values():
                    alignment_score += _trf_result["alignment_score"]
                if alignment_score > max_alignment_score:
                    max_alignment_score = alignment_score
                    best_trf_results = matching_trf_results

        num_total_repeats_in_variant_bases = 0
        num_total_repeats_left_flank = tandem_repeat_bases_in_left_flank = 0
        num_total_repeats_right_flank = tandem_repeat_bases_in_right_flank = 0


        if best_trf_results.get("left"):
            tandem_repeat_bases_in_left_flank = max(0, best_trf_results["left"]["end_1based"] - len(variant_bases))
            num_total_repeats_left_flank = tandem_repeat_bases_in_left_flank / best_trf_results["left"]["repeat_unit_length"]
            num_total_repeats_in_variant_bases = best_trf_results["left"]["num_repeats_in_variant"]
        elif best_trf_results.get("right"):
            tandem_repeat_bases_in_right_flank = max(0, best_trf_results["right"]["end_1based"] - len(variant_bases))
            num_total_repeats_right_flank = tandem_repeat_bases_in_right_flank / best_trf_results["right"]["repeat_unit_length"]
            num_total_repeats_in_variant_bases = best_trf_results["right"]["num_repeats_in_variant"]
        else:
            return None, FILTER_ALLELE_INDEL_WITHOUT_REPEATS

        if best_trf_results.get("right") is not None:
            result["RepeatUnit"] = repeat_unit = best_trf_results["right"]["repeat_unit"]  # use values from "right" because it's not reversed
            result["MotifInterruptionIndices"] = ",".join(map(str, sorted([position_in_motif for position_in_motif, count in best_trf_results["right"]["motif_positions_with_interruptions"].items()])))
        elif best_trf_results.get("left") is not None:
            result["RepeatUnit"] = repeat_unit = best_trf_results["left"]["repeat_unit"][::-1]
            result["MotifInterruptionIndices"] = ",".join(map(str, sorted([len(repeat_unit) - position_in_motif + 1 for position_in_motif, count in best_trf_results["left"]["motif_positions_with_interruptions"].items()])))
        else:
            return None, FILTER_ALLELE_INDEL_WITHOUT_REPEATS

        result["IsPureRepeat"] = False
    else:
        raise ValueError(f"Invalid detection_mode: '{detection_mode}'. It must be 'pure_repeats' or 'allow_interruptions'.")


    total_repeats = num_total_repeats_left_flank + num_total_repeats_in_variant_bases + num_total_repeats_right_flank

    tandem_repeat_allele_failed_filters_reason = tandem_repeat_allele_failed_filters(args, repeat_unit, total_repeats, counters=counters)
    if tandem_repeat_allele_failed_filters_reason is not None:
        return None, tandem_repeat_allele_failed_filters_reason

    # the allele passed TR filters
    result["Start0Based"] = left_flank_end - tandem_repeat_bases_in_left_flank
    result["End1Based"] = right_flank_start_1based + tandem_repeat_bases_in_right_flank

    result["NumRepeatsRef"] = num_total_repeats_left_flank + num_total_repeats_right_flank
    result["NumRepeatsAlt"] = result["NumRepeatsRef"]
    if len(alt_allele) < len(vcf_ref):
        # add deleted repeats to the count of repeats in the reference
        result["NumRepeatsRef"] += num_total_repeats_in_variant_bases
    elif len(alt_allele) > len(vcf_ref):
        # add inserted repeats to the count of repeats in the variant bases
        result["NumRepeatsAlt"] += num_total_repeats_in_variant_bases

    result["AlleleRepeatSequence"] = ""
    if tandem_repeat_bases_in_left_flank > 0:
        result["AlleleRepeatSequence"] += left_flanking_reference_sequence[-tandem_repeat_bases_in_left_flank:]
    if len(alt_allele) > len(vcf_ref):
        result["AlleleRepeatSequence"] += variant_bases
    if tandem_repeat_bases_in_right_flank > 0:
        result["AlleleRepeatSequence"] += right_flanking_reference_sequence[:tandem_repeat_bases_in_right_flank]

    return result, None


def process_vcf_allele(
    fasta_obj,
    vcf_chrom,
    vcf_pos,
    vcf_ref,
    alt_allele,
    args,
    counters,
    only_generate_trf_fasta,
    variants_to_process_using_trf,
    variants_per_trf_fasta,
):
    """Process a single VCF allele and determine if it represents a tandem repeat expansion or contraction.

    Args:
        fasta_obj (pyfaidx.Fasta): Fasta object containing the reference genome.
        vcf_chrom (str): Chromosome name from the VCF file.
        vcf_pos (int): 1-based position of the variant in the VCF file.
        vcf_ref (str): Reference sequence at the variant position.
        alt_allele (str): Alternate allele sequence from the VCF file.
        args (argparse.Namespace): Command-line arguments parsed by parse_args().
        counters (dict, optional): Dictionary of counters to collect summary stats about the number of TR variants found, etc.
        only_generate_trf_fasta (bool): instead of running TRF, just generate a fasta file with the nucleotide sequence to run TRF on later.
        variants_to_process_using_trf (dict): maps variant IDs that should be processed using TRF to their sequence number in the TRF fasta file.
        variants_per_trf_fasta (dict): maps TRF fasta input file path to the total number of variants in that fasta file.

    Return:
        2-tuple (dict, str):
            dict: if the allele represents a tandem repeat, this will be a JSON record describing the tandem repeat, otherwise it will be None.
            str: if the allele is not a tandem repeat, this will be a string describing the reason why the allele failed tandem repeat filters, or otherwise None if it passed all filters.
    """

    # check if it's an indel variant and compute variant_bases
    if len(vcf_ref) == len(alt_allele):
        counters[f"allele filter: {'SNV' if len(alt_allele) == 1 else 'MNV'}"] += 1
        return None, FILTER_ALLELE_SNV_OR_MNV

    elif len(vcf_ref) < len(alt_allele):
        counters[f"allele counts: INS alleles"] += 1
        ins_or_del = "INS"
        if alt_allele.startswith(vcf_ref):
            variant_bases = alt_allele[len(vcf_ref):]
        else:
            counters[f"allele filter: complex MNV deletion/insertion"] += 1
            return None, FILTER_ALLELE_MNV_INDEL

    elif len(alt_allele) < len(vcf_ref):
        counters[f"allele counts: DEL alleles"] += 1
        ins_or_del = "DEL"
        if vcf_ref.startswith(alt_allele):
            variant_bases = vcf_ref[len(alt_allele):]
        else:
            counters[f"allele filter: complex MNV insertion/deletion"] += 1
            return None, FILTER_ALLELE_MNV_INDEL


    if len(variant_bases) > MAX_INDEL_SIZE:
        # this is a very large indel, so we don't want to process it
        counters[f"allele filter: {FILTER_ALLELE_TOO_BIG}"] += 1
        return None, FILTER_ALLELE_TOO_BIG

    # retrieve flanking sequences to the left and right of the variant

    # A VCF line like this:
    #      chr11	86734686	CACAATATACATATATATTATATATATAT	C,CACAATATACATATATATTATATATATATACAATATACATATATATTATATATATAT    GT   1/2
    # represents two alleles at pos = 86734686:
    #    allele 1 is a deletion where:
    #       left_flank = chr11_sequence[pos - num_flanking_bases: pos]
    #       variant_bases = "ACAATATACATATATATTATATATATAT"
    #       right_flank = chr11_sequence[pos: pos + len(variant_bases) + num_flanking_bases]
    #    allele 2 is an insertion:
    #       left_flank = chr11_sequence[pos - num_flanking_bases: pos + len(vcf_ref) - 1]
    #       variant_bases = "ACAATATACATATATATTATATATATAT"
    #       right_flank = chr11_sequence[pos + len(vcf_ref) - 1: pos + len(vcf_ref) - 1 + num_flanking_bases]

    chrom_size = len(fasta_obj[vcf_chrom])

    detection_modes_to_try = [DETECTION_MODE_PURE_REPEATS]
    flanking_sequence_size_multipliers = [3, 30, 100, 300, 1000, 3000]
    if not args.dont_allow_interruptions:
        detection_modes_to_try.append(DETECTION_MODE_ALLOW_INTERRUPTIONS)
        if not args.dont_run_trf and len(variant_bases) >= args.min_indel_size_to_run_trf:
            detection_modes_to_try.append(DETECTION_MODE_TRF)

    rerun_trf_on_longer_flanking_sequences = False
    current_variants_to_process_using_trf = variants_to_process_using_trf
    current_variants_per_trf_fasta = variants_per_trf_fasta
    for flanking_sequence_size_multiplier in flanking_sequence_size_multipliers:
        # start with relatively small flanking sequence sizes and increase them if it turns out that the entire flank
        # is covered by a tandem repeat
        num_flanking_bases = flanking_sequence_size_multiplier * max(len(variant_bases), 100)

        if len(vcf_ref) < len(alt_allele):
            left_flank_end = vcf_pos + len(vcf_ref) - 1
            right_flank_start_1based = left_flank_end
        elif len(alt_allele) < len(vcf_ref):
            # For deletions, the right-flanking sequence starts to the right of the variant bases
            left_flank_end = vcf_pos + len(alt_allele) - 1
            right_flank_start_1based = vcf_pos + len(vcf_ref) - 1

        left_flank_start_1based = max(left_flank_end - num_flanking_bases, 0)
        right_flank_end = min(right_flank_start_1based + num_flanking_bases, chrom_size)

        left_flanking_reference_sequence = str(fasta_obj[vcf_chrom][left_flank_start_1based : left_flank_end]).upper()
        right_flanking_reference_sequence = str(fasta_obj[vcf_chrom][right_flank_start_1based : right_flank_end]).upper()

        if TOO_MANY_Ns in left_flanking_reference_sequence or TOO_MANY_Ns in right_flanking_reference_sequence:
            counters[f"allele filter: {FILTER_ALLELE_WITH_TOO_MANY_Ns_IN_FLANKS}"] += 1
            return None, FILTER_ALLELE_WITH_TOO_MANY_Ns_IN_FLANKS

        for detection_mode_i, detection_mode in enumerate(detection_modes_to_try):
            if rerun_trf_on_longer_flanking_sequences and detection_mode == DETECTION_MODE_TRF:
                # This is a special, relatively rare case where we need to rerun TRF on the same variant, but passing in
                # longer flanking sequences. This special handling is triggered if, after running TRF in bulk mode, we
                # find that the detected repeats entirely cover (and probably exceeded) the left and/or right flanking
                # sequences.
                # This causes us to fall back on running TRF on one variant at a time (for just this variant).
                # Since this would be too slow to do for all variants by default, so we start with bulk processing, and
                # then fall back on this special-case handling where we output a new fasta file that contains just
                # the current variant with the longer flanks, and then immediately invoke TRF on this single-variant
                # fasta in the second call to check_if_allele_is_tandem_repeat(..) immediately below this if statement.
                check_if_allele_is_tandem_repeat(
                    vcf_chrom,
                    vcf_pos,
                    vcf_ref,
                    alt_allele,
                    variant_bases,
                    left_flanking_reference_sequence,
                    right_flanking_reference_sequence,
                    left_flank_end,
                    right_flank_start_1based,
                    args,
                    counters=None,
                    only_generate_trf_fasta=True,
                    variants_to_process_using_trf=current_variants_to_process_using_trf,
                    variants_per_trf_fasta=current_variants_per_trf_fasta,
                    detection_mode=DETECTION_MODE_TRF,
                    rerun_trf_on_longer_flanking_sequences=rerun_trf_on_longer_flanking_sequences)

            tandem_repeat_allele, tandem_repeat_allele_failed_filters_reason = check_if_allele_is_tandem_repeat(
                vcf_chrom,
                vcf_pos,
                vcf_ref,
                alt_allele,
                variant_bases,
                left_flanking_reference_sequence,
                right_flanking_reference_sequence,
                left_flank_end,
                right_flank_start_1based,
                args,
                counters=counters,
                only_generate_trf_fasta=only_generate_trf_fasta,
                variants_to_process_using_trf=current_variants_to_process_using_trf,
                variants_per_trf_fasta=current_variants_per_trf_fasta,
                detection_mode=detection_mode,
                rerun_trf_on_longer_flanking_sequences=rerun_trf_on_longer_flanking_sequences)

            if tandem_repeat_allele_failed_filters_reason is None:
                if (tandem_repeat_allele["Start0Based"] - len(tandem_repeat_allele["RepeatUnit"]) <= left_flank_start_1based) or (
                    tandem_repeat_allele["End1Based"] + len(tandem_repeat_allele["RepeatUnit"]) >= right_flank_end):
                    # The flanking sequence was not long enough to contain the entire repeat so we need to extend the flanks and try again
                    #print(f"Extending flanking sequence beyond {num_flanking_bases}bp for {tandem_repeat_allele['Chrom']}:{tandem_repeat_allele['Start0Based']}-{tandem_repeat_allele['End1Based']} ")

                    if detection_mode == DETECTION_MODE_TRF:
                        rerun_trf_on_longer_flanking_sequences = True
                        current_variants_to_process_using_trf = {}
                        current_variants_per_trf_fasta = collections.defaultdict(int)

                    print(f"NOTE: With detection_mode == {detection_mode}, the TR variant {vcf_chrom}:{vcf_pos} {vcf_ref} > {alt_allele} "
                          f"spanned the entire left and/or right flanking sequences (variant size: {len(variant_bases):,d}bp, "
                          f"flank size: {num_flanking_bases:,d}bp). Skipping...")

                    break

                # Found this allele to be a tandem repeat expansion/contraction
                tandem_repeat_allele["INS_or_DEL"] = ins_or_del
                return tandem_repeat_allele, None

            elif detection_mode_i == len(detection_modes_to_try) - 1:
                # this is the last detection mode, so the allele was not found to be a tandem repeat expansion/contraction
                return None, tandem_repeat_allele_failed_filters_reason
        else:
            raise RuntimeError("State error: should always break out of the loop before it completes")

    counters[f"allele filter: {FILTER_ALLELE_TR_SPANS_TOO_MANY_BASES}"] += 1
    return None, FILTER_ALLELE_TR_SPANS_TOO_MANY_BASES


def postprocess_multiallelic_tandem_repeats(tandem_repeat_alleles, counters):

    assert len(tandem_repeat_alleles) > 1, f"called with only {len(tandem_repeat_alleles)} allele"
    assert tandem_repeat_alleles[0]["Chrom"] == tandem_repeat_alleles[1]["Chrom"]

    # make sure the RepeatUnit is consistent across alleles
    if tandem_repeat_alleles[0]["MotifSize"] != tandem_repeat_alleles[1]["MotifSize"] or (
        tandem_repeat_alleles[0]["MotifSize"] <= 6
        and tandem_repeat_alleles[0]["CanonicalMotif"] != tandem_repeat_alleles[1]["CanonicalMotif"]
    ):
        return tandem_repeat_alleles, FILTER_VARIANT_WITH_TR_ALLELES_WITH_DIFFERENT_MOTIFS

    if tandem_repeat_alleles[0]["MotifInterruptionIndices"] is None:
        tandem_repeat_alleles[0]["MotifInterruptionIndices"] = tandem_repeat_alleles[1]["MotifInterruptionIndices"]
    elif tandem_repeat_alleles[1]["MotifInterruptionIndices"] is None:
        tandem_repeat_alleles[1]["MotifInterruptionIndices"] = tandem_repeat_alleles[0]["MotifInterruptionIndices"]

    if tandem_repeat_alleles[0]["DetectionMode"] != DETECTION_MODE_TRF and tandem_repeat_alleles[1]["DetectionMode"] != DETECTION_MODE_TRF:
        if tandem_repeat_alleles[0]["MotifInterruptionIndices"] != tandem_repeat_alleles[1]["MotifInterruptionIndices"]:
            return tandem_repeat_alleles, FILTER_VARIANT_WITH_TR_ALLELES_WITH_DIFFERENT_INTERRUPTION_PATTERNS

    # make sure the "IsPureRepeat" and "MotifInterruptionIndices" fields are consistent across alleles
    if tandem_repeat_alleles[0]["IsPureRepeat"] != tandem_repeat_alleles[1]["IsPureRepeat"]:
        tandem_repeat_alleles[0]["IsPureRepeat"] = tandem_repeat_alleles[1]["IsPureRepeat"] = False

    # if the locus definitions are not the same, use the larger one
    if tandem_repeat_alleles[0]["Start0Based"] != tandem_repeat_alleles[1]["Start0Based"] or tandem_repeat_alleles[0]["End1Based"] != tandem_repeat_alleles[1]["End1Based"]:
        copy_from, copy_to = (0, 1) if tandem_repeat_alleles[0]["NumRepeatsRef"] > tandem_repeat_alleles[1]["NumRepeatsRef"] else (1, 0)
        tandem_repeat_alleles[copy_to]["Start0Based"] = tandem_repeat_alleles[copy_from]["Start0Based"]
        tandem_repeat_alleles[copy_to]["End1Based"] = tandem_repeat_alleles[copy_from]["End1Based"]
        # increase the number of repeats in the alt allele to account for the added repeats in the reference
        tandem_repeat_alleles[copy_to]["NumRepeatsAlt"] += tandem_repeat_alleles[copy_from]["NumRepeatsRef"] - tandem_repeat_alleles[copy_to]["NumRepeatsRef"]
        tandem_repeat_alleles[copy_to]["NumRepeatsRef"] = tandem_repeat_alleles[copy_from]["NumRepeatsRef"]
        tandem_repeat_alleles[copy_to]["DetectionMode"] = tandem_repeat_alleles[copy_from]["DetectionMode"]
        counters["TR variant counts: multi-allelic variant: adjusted locus coordinates"] += 1

    return tandem_repeat_alleles, None


def is_found_in_reference(tandem_repeat_allele):
    """Return True if the allele has more than zero repeats in the reference genome, False otherwise."""

    return tandem_repeat_allele["End1Based"] - tandem_repeat_allele["Start0Based"] >= len(tandem_repeat_allele["RepeatUnit"])


def get_num_repeats_in_allele(tandem_repeat_alleles, genotype_index):
    """Looks up the number of TR repeats found in a given VCF record's allele(s) based on the genotype.

    Args:
        tandem_repeat_alleles (list): list containing the TR allele spec(s) returned by check_if_allele_is_tandem_repeat for each
            of the allele(s) in this variant.
        genotype_index (int): which genotype to look up.

    Return
        int: The number of repeats in this allele
    """
    if genotype_index == 0:
        return tandem_repeat_alleles[0]["NumRepeatsRef"]

    return tandem_repeat_alleles[genotype_index - 1]["NumRepeatsAlt"]


def compute_summary_string(tandem_repeat_alleles):
    """Returns a short easy-to-read string summarizing the TR locus.

    Args:
        tandem_repeat_alleles (list): list of 1 or more allele spec dictionaries

    Return:
        str: short summary of the variant
    """

    repeat_unit = tandem_repeat_alleles[0]["RepeatUnit"]
    summary_string = f"{len(repeat_unit)}bp:"
    if len(repeat_unit) > 30:
        summary_string += f"{repeat_unit[:30]}...:"
    else:
        summary_string += f"{repeat_unit}:"
        
    ins_or_del = [a["INS_or_DEL"] for a in tandem_repeat_alleles]
    summary_string += ",".join(ins_or_del)

    if tandem_repeat_alleles[0]["IsPureRepeat"]:
        summary_string += ":pure"
    else:
        summary_string += ":not-pure"

    return ":".join(ins_or_del), summary_string


def process_vcf_line(
    vcf_line_i,
    vcf_fields,
    fasta_obj,
    vcf_writer,
    tsv_writer,
    bed_writer,
    fasta_writer,
    args,
    counters,
    variant_intervals,
    only_generate_trf_fasta,
    variants_to_process_using_trf,
    variants_per_trf_fasta,
    tr_locus_ids,
):
    """Utility method for processing a single line from the input vcf

    Args:
        vcf_line_i (int): line number of the VCF file
        vcf_fields (list): list of fields from the VCF line
        fasta_obj (pyfaidx.Fasta): Fasta object for the reference genome
        vcf_writer (open file): open file for writing VCF output
        tsv_writer (open file): open file for writing locus TSV
        bed_writer (open file): open file for writing BED format TR loci
        fasta_writer (open file): open file for writing TR alleles in FASTA format
        args (argparse.Namespace): parsed command-line arguments
        counters (dict): dictionary of counters
        variant_intervals (dict): maps chrom to list of intervaltree.Intervals representing vcf variants for subsequent overlap detection
        only_generate_trf_fasta (bool): instead of running TRF, just generate a fasta file with the nucleotide sequence to run TRF on later.
        variants_to_process_using_trf (dict): maps variant IDs that should be processed using TRF to their sequence number in the TRF fasta file.
        variants_per_trf_fasta (dict): maps TRF fasta input file path to the total number of variants in that fasta file.
        tr_locus_ids (set): locus ids of TR loci identified so far

    Return:
        str: a string explaining the reason this variant is not a tandem repeat, or None if the variant is a tandem repeat
    """

    # parse the ALT allele(s)
    vcf_chrom = vcf_fields[0]
    vcf_chrom_without_chr_prefix = vcf_chrom.replace("chr", "")
    vcf_pos = int(vcf_fields[1])
    vcf_ref = vcf_fields[3].upper()
    vcf_alt = vcf_fields[4].upper()
    alt_alleles = vcf_alt.split(",")

    if vcf_chrom not in fasta_obj:
        raise ValueError(f"Chromosome '{vcf_chrom}' not found in the reference fasta")

    counters["variant counts: TOTAL variants"] += 1
    counters["allele counts: TOTAL alleles"] += len(alt_alleles)

    # check for N's in the ref or alt sequences
    if "N" in vcf_ref or "N" in vcf_alt:
        counters[f"allele filter: {FILTER_ALLELE_WITH_N_BASES}"] += 1
        return FILTER_ALLELE_WITH_N_BASES

    if not vcf_alt:
        raise ValueError(f"No ALT allele found in line {vcf_fields}")

    variant_id = f"{vcf_chrom_without_chr_prefix}-{vcf_pos}-{vcf_ref}-{vcf_alt}"

    if len(alt_alleles) > 2:
        if args.verbose:
            print(f"WARNING: vcf row #{vcf_line_i:,d}: {variant_id}: multi-allelic variant has "
                  f"{len(alt_alleles)} alt alleles. This script doesn't support more than 2 alt alleles. Skipping...")
        counters[f"WARNING: multi-allelic variant with {len(alt_alleles)} alt alleles"] += 1
        return FILTER_MORE_THAN_TWO_ALT_ALLELES

    # Handle '*' alleles
    star_allele_index = None
    if "*" in alt_alleles:
        # if this variant has 1 regular allele and 1 "*" allele (which represents an overlapping deletion), discard the
        # "*" allele and recode the genotype as haploid
        #counters["variant counts: removed * allele and converted to homozygous genotype"] += 1
        star_allele_index = alt_alleles.index("*") + 1
        alt_alleles = [a for a in alt_alleles if a != "*"]

    if len(alt_alleles) == 0:
        return FILTER_ZERO_ALT_ALLELES

    # Check if the allele(s) are tandem repeat expansions or contractions
    tandem_repeat_alleles = []
    for alt_allele in alt_alleles:
        tandem_repeat_allele, filter_reason = process_vcf_allele(
            fasta_obj,
            vcf_chrom,
            vcf_pos,
            vcf_ref,
            alt_allele,
            args,
            counters,
            only_generate_trf_fasta,
            variants_to_process_using_trf,
            variants_per_trf_fasta,
        )

        if filter_reason is not None and filter_reason != WILL_RUN_TRF_ON_THIS_ALLELE_IN_2ND_PASS:
            return filter_reason

        # update counters
        if tandem_repeat_allele is not None:
            counters[f"TR allele counts: TOTAL"] += 1
            counters[f"TR allele counts: {'INS' if len(vcf_ref) < len(alt_allele) else 'DEL'}"] += 1
            allele_repeat_unit = tandem_repeat_allele["RepeatUnit"]
            allele_detection_mode = tandem_repeat_allele["DetectionMode"]
            counters[f"TR allele motif size: {len(allele_repeat_unit) if len(allele_repeat_unit) < 9 else '9+'} bp"] += 1
            counters[f"TR allele detected by: {allele_detection_mode}"] += 1

        tandem_repeat_alleles.append(tandem_repeat_allele)

    if only_generate_trf_fasta:
        return WILL_RUN_TRF_ON_THIS_ALLELE_IN_2ND_PASS

    # Compute extra fields
    for tandem_repeat_allele in tandem_repeat_alleles:
        tandem_repeat_allele["MotifSize"] = len(tandem_repeat_allele["RepeatUnit"])
        tandem_repeat_allele["CanonicalMotif"] = compute_canonical_motif(tandem_repeat_allele["RepeatUnit"], include_reverse_complement=True)

    # Post-process multiallelic tandem repeat variants
    if len(tandem_repeat_alleles) > 1:
        tandem_repeat_alleles, filter_reason = postprocess_multiallelic_tandem_repeats(tandem_repeat_alleles, counters)

        if filter_reason is not None:
            return filter_reason

    # Parse the genotype
    vcf_genotype_format = vcf_fields[8].split(":")
    if vcf_genotype_format[0] != "GT":
        if args.verbose:
            print(f"WARNING: vcf row #{vcf_line_i:,d}: Unexpected genotype field in row #{vcf_line_i:,d}: "
                  f"{vcf_fields[8]}  (variant: {variant_id})")
        counters[f"WARNING: {FILTER_UNEXPECTED_GENOTYPE_FORMAT}"] += 1
        return FILTER_UNEXPECTED_GENOTYPE_FORMAT

    vcf_genotype = next(iter(vcf_fields[9].split(":")))

    counters["TR variant counts: TOTAL"] += 1

    # Get repeat unit, start, end coords for the variant
    repeat_unit = tandem_repeat_alleles[0]["RepeatUnit"]
    canonical_motif = tandem_repeat_alleles[0]["CanonicalMotif"]
    start_0based = tandem_repeat_alleles[0]["Start0Based"]
    end_1based = tandem_repeat_alleles[0]["End1Based"]
    num_repeats_in_reference = tandem_repeat_alleles[0]["NumRepeatsRef"]
    is_pure_repeat = tandem_repeat_alleles[0]["IsPureRepeat"]
    repeat_unit_interruption_indices = tandem_repeat_alleles[0]["MotifInterruptionIndices"]

    variant_detection_mode = tandem_repeat_alleles[0]["DetectionMode"]
    if len(tandem_repeat_alleles) > 1 and DETECTION_MODE_ORDER.index(variant_detection_mode) < DETECTION_MODE_ORDER.index(tandem_repeat_alleles[1]["DetectionMode"]):
        # if the detection mode of the first allele is lower than the second, use the second's detection mode
        variant_detection_mode = tandem_repeat_alleles[1]["DetectionMode"]

    locus_id = f"{vcf_chrom_without_chr_prefix}-{start_0based}-{end_1based}-{repeat_unit}"

    if locus_id in tr_locus_ids:
        return FILTER_DUPLICATE_TR_LOCUS

    tr_locus_ids.add(locus_id)

    # Generate output records
    variant_interval = intervaltree.Interval(start_0based, max(start_0based + 1, end_1based), data={
        "LocusId": locus_id,
        "IsPureRepeat": is_pure_repeat,
        "RepeatUnit": repeat_unit,
        "CanonicalMotif": canonical_motif,
    })

    variant_intervals[vcf_chrom].append(variant_interval)

    variant_ins_or_del, summary_string = compute_summary_string(tandem_repeat_alleles)

    # Generate the output records
    tsv_record = {
        "LocusId": locus_id,
        "Locus": f"{vcf_chrom_without_chr_prefix}:{start_0based}-{end_1based}",
        "Motif": repeat_unit,
        "CanonicalMotif": canonical_motif,
        "NumRepeatsInReference": num_repeats_in_reference,
        "IsPureRepeat": is_pure_repeat,
        "DetectionMode": variant_detection_mode,
    }

    if not is_pure_repeat:
        tsv_record["MotifInterruptionIndices"] = repeat_unit_interruption_indices if repeat_unit_interruption_indices is not None else ""

    if args.copy_info_field_keys_to_tsv:
        if vcf_fields[7] and vcf_fields[7] != ".":
            for info_key_value in vcf_fields[7].split(";"):
                info_field_tokens = info_key_value.split("=")
                info_field_key = info_field_tokens[0]
                if not any(info_field_key.lower() != k.lower() for k in args.copy_info_field_keys_to_tsv):
                    continue
                if len(info_field_tokens) > 1:
                    tsv_record[info_field_key] = info_field_tokens[1]
                else:
                    tsv_record[info_field_key] = True

    if vcf_writer is not None:
        # write results to a VCF file
        vcf_fields[2] = summary_string
        vcf_fields[4] = ",".join(alt_alleles)
        vcf_fields[7] = ";".join([f"{key}={value}" for key, value in tsv_record.items() if value is not None and value != ""])
        vcf_fields[9] = vcf_genotype
        vcf_writer.write("\t".join(vcf_fields) + "\n")

    # write results to TSVs
    if len(alt_alleles) > 1:
        counters[f"TR variant counts: multi-allelic"] += 1

    if is_pure_repeat:
        counters[f"TR variant counts: pure repeats"] += 1
    else:
        counters[f"TR variant counts: interrupted repeats"] += 1

    tsv_record.update({
        "Chrom": vcf_chrom,
        "Start0Based": start_0based,
        "End1Based": end_1based,
        "MotifSize": len(repeat_unit),
        "VcfPos": vcf_pos,
        #"VcfRef": vcf_ref,
        #"VcfAlt": ",".join(alt_alleles),
        "SummaryString": summary_string,
        "VcfGenotype": vcf_genotype,
        "INS_or_DEL": variant_ins_or_del,

        "IsMultiallelic": len(alt_alleles) > 1,
        "IsFoundInReference": any(is_found_in_reference(a) for a in tandem_repeat_alleles),
    })

    if tsv_writer is not None:
        tsv_writer.write("\t".join([str(tsv_record.get(c, "")) for c in TSV_OUTPUT_COLUMNS]) + "\n")

    if bed_writer is not None:
        bed_writer.write("\t".join(map(str, [vcf_chrom, start_0based, end_1based, repeat_unit, "."])) + "\n")

    assert len(alt_alleles) == len(tandem_repeat_alleles), f"{len(alt_alleles)} alt alleles, but {len(tandem_repeat_alleles)} tandem repeat alleles"

    if fasta_writer is not None:
        for i, (alt_allele, tandem_repeat_allele) in enumerate(zip(alt_alleles, tandem_repeat_alleles)):
            fasta_writer.write(f">{vcf_chrom}-{start_0based}-{end_1based}-{repeat_unit}-allele{i+1}: {summary_string}\n")
            fasta_writer.write(f"{tandem_repeat_allele['AlleleRepeatSequence']}\n")


def process_tandem_repeat_loci_that_have_overlapping_variants(
        file_path, locus_ids_with_overlapping_variants, filtered_out_variants_vcf_writer):
    """Process tandem repeat loci that turned out to overlap more than one variant in the VCF.
    Their true genotype is ambiguous and cannot be easily determined here, so just remove them from the output file.

    Args:
        file_path (str): Path to VCF or TSV file to rewrite in order to remove loci with the given locus ids which overlap multiple variants.
        locus_ids_with_overlapping_variants (dict): Locus ids that, in previous steps, were found to overlap multiple
            variants in the input vcf, mapped to a description of the specific reason to filter out that particular locus.
        filtered_out_variants_vcf_writer (VcfWriter): VCF writer for filtered out variants
    """

    if not any(file_path.endswith(suffix) for suffix in (".bed", ".vcf", ".tsv")):
        raise ValueError(f"Unexpected file type: {file_path}")

    filtered_count = total = 0
    with open(file_path, "rt") as f, open(f"{file_path}.temp", "wt") as fo:
        header = None
        for line in f:
            if file_path.endswith(".vcf") and line.startswith("#"):
                fo.write(line)
                continue
            if file_path.endswith(".tsv") and header is None:
                header = line.strip().split("\t")
                if "LocusId" not in header:
                    raise ValueError(f"LocusId column not found in {file_path} header: {header}")
                fo.write(line)
                continue

            total += 1
            fields = line.strip().split("\t")
            if file_path.endswith(".tsv"):
                row = dict(zip(header, fields))
                locus_id = row["LocusId"]
            elif file_path.endswith(".vcf"):
                info = dict([info_record.split("=") for info_record in fields[7].split(";") if "=" in info_record])
                if "LocusId" not in info:
                    raise ValueError(f"Unexpected info field format in {file_path} line: {fields}")
                locus_id = info["LocusId"]

            elif file_path.endswith(".bed"):
                chrom = fields[0].replace("chr", "")
                motif = fields[3]
                locus_id = f"{chrom}-{fields[1]}-{fields[2]}-{motif}"

            if locus_id in locus_ids_with_overlapping_variants:
                filtered_count += 1
                if file_path.endswith(".vcf") and filtered_out_variants_vcf_writer is not None:
                    # update the FILTER field
                    n_alleles = len(fields[4].split(","))
                    filter_reason = locus_ids_with_overlapping_variants[locus_id]
                    fields[6] = ";".join([filter_reason]*n_alleles)
                    filtered_out_variants_vcf_writer.write("\t".join(fields) + "\n")
            else:
                fo.write("\t".join(fields) + "\n")

    os.rename(f"{file_path}.temp", file_path)
    print(f"Removed {filtered_count:,d} out of {total:,d} rows from {file_path} because they were better represented by another overlapping TR variant")


def print_stats(counters):
    """Print out all the counters"""

    key_prefixes = set()
    for key, _ in counters.items():
        tokens = key.split(":")
        key_prefixes.add(f"{tokens[0]}:")

    for print_totals_only in True, False:
        for key_prefix in sorted(key_prefixes):
            if print_totals_only ^ (key_prefix in ("variant counts:", "allele counts:")):
                continue
            current_counter = [(key, count) for key, count in counters.items() if key.startswith(key_prefix)]
            current_counter = sorted(current_counter, key=lambda x: (-x[1], x[0]))
            print("--------------")
            for key, value in current_counter:
                if key_prefix.startswith("TR"):
                    total_key = "TR variant counts: TOTAL" if "variant" in key_prefix else "TR allele counts: TOTAL"
                else:
                    total_key = "variant counts: TOTAL variants" if "variant" in key_prefix else "allele counts: TOTAL alleles"

                total = counters[total_key]
                percent = f"{100*value / total:5.1f}%" if total > 0 else ""

                if print_totals_only:
                    print(f"{value:10,d}  {key}")
                else:
                    print(f"{value:10,d} out of {total:10,d} ({percent}) {key}")


def run(command):
    """Run a shell command"""

    os.system(command)


def main(args, only_generate_trf_fasta=False, variants_to_process_using_trf=None, variants_per_trf_fasta=None):
    """Main function to process the input VCF and identify tandem repeat variants.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.
        only_generate_trf_fasta (bool): If True, only generate the TRF fasta file without running TRF.
        variants_to_process_using_trf (dict): maps variant IDs that should be processed using TRF to their sequence number in the TRF fasta file.
        variants_per_trf_fasta (dict): maps TRF fasta input file path to the total number of variants in that fasta file.
    """
    counters = collections.defaultdict(int)

    if args.debug: print("-"*150)
    variant_intervals = None
    if not only_generate_trf_fasta:
        variant_intervals = collections.defaultdict(list)  # maps each chromosome to a list of intervaltree.Intervals representing all variants in the input VCF

    input_files_to_close = []

    # open the reference genome fasta
    fasta_obj = pyfaidx.Fasta(args.reference_fasta_path, one_based_attributes=False, as_raw=True)
    input_files_to_close.append(fasta_obj)

    # open output files for writing
    vcf_writer = None
    if args.write_vcf_file and not only_generate_trf_fasta:
        vcf_writer = open(f"{args.output_prefix}.vcf", "wt")

    filtered_out_variants_vcf_writer = None
    if args.write_vcf_with_filtered_out_variants and not only_generate_trf_fasta:
        filtered_out_variants_vcf_writer = open(f"{args.output_prefix}.filtered_out_rows.vcf", "wt")

    tsv_writer = None
    if not only_generate_trf_fasta:
        tsv_writer = open(f"{args.output_prefix}.loci.tsv", "wt")
        tsv_writer.write("\t".join(TSV_OUTPUT_COLUMNS) + "\n")

    bed_writer = None
    if args.write_bed_file and not only_generate_trf_fasta:
        bed_writer = open(f"{args.output_prefix}.loci.bed", "wt")

    fasta_writer = None
    if args.write_fasta and not only_generate_trf_fasta:
        fasta_writer = gzip.open(f"{args.output_prefix}.fasta.gz", "wt")

    # open the input VCF
    if args.interval:
        print(f"Fetching interval(s):", ", ".join(args.interval))
        tabix_file = pysam.TabixFile(args.input_vcf_path)
        vcf_iterator = (line for interval in args.interval for line in tabix_file.fetch(interval))
        input_files_to_close.append(tabix_file)
    else:
        vcf_iterator = open_file(args.input_vcf_path, is_text_file=True)
        input_files_to_close.append(vcf_iterator)

    if args.show_progress_bar:
        vcf_iterator = tqdm.tqdm(vcf_iterator, unit=" variants", unit_scale=True)

    # iterate over all VCF rows
    tr_locus_ids = set()  # keep a running set of detected TR loci to avoid duplicates
    vcf_line_i = 0
    for line in vcf_iterator:
        if line.startswith("##"):
            # copy the VCF header to the output
            if vcf_writer is not None:
                vcf_writer.write(line)
            if filtered_out_variants_vcf_writer is not None:
                filtered_out_variants_vcf_writer.write(line)

            continue

        vcf_fields = line.strip().split("\t")

        if line.startswith("#"):
            # process the last line of the VCF header
            sample_ids = vcf_fields[9:]
            if len(sample_ids) != 1:
                raise ValueError(f"{args.input_vcf_path} contains {len(sample_ids)} samples, but this script only "
                                 f"supports single-sample VCFs.")

            if vcf_writer is not None:
                vcf_writer.write(line)
            if filtered_out_variants_vcf_writer is not None:
                filtered_out_variants_vcf_writer.write(line)

            continue

        if vcf_line_i < args.offset:
            vcf_line_i += 1
            continue
        if args.n is not None and vcf_line_i >= args.offset + args.n:
            break

        vcf_line_i += 1

        filter_string = process_vcf_line(
            vcf_line_i, vcf_fields, fasta_obj, vcf_writer, 
            tsv_writer, bed_writer, fasta_writer, args, counters,
            variant_intervals, only_generate_trf_fasta, variants_to_process_using_trf, variants_per_trf_fasta, tr_locus_ids)

        if filter_string and filtered_out_variants_vcf_writer is not None:
            # if this variant is an indel that was filtered out, write it to the filtered out variants VCF

            vcf_fields[6] = filter_string

            if filter_string not in (
                FILTER_ALLELE_SNV_OR_MNV,
                ";".join([FILTER_ALLELE_SNV_OR_MNV]*2),
            ):
                filtered_out_variants_vcf_writer.write("\t".join(vcf_fields) + "\n")

    for file_obj in input_files_to_close:
        try:
            file_obj.close()
        except Exception as e:
            print(f"WARNING: unable to close {file_obj.name}: {e}")

    if only_generate_trf_fasta:
        return

    handle_overlapping_loci(variant_intervals, tsv_writer, vcf_writer, bed_writer, filtered_out_variants_vcf_writer, args, counters)

    # close and post-process all output files
    for writer in filter(None, [tsv_writer, vcf_writer, bed_writer, filtered_out_variants_vcf_writer]):
        writer.close()   # close writers (it's ok to close a file that's already been closed)

        if writer.name.endswith(".tsv"):
            run(f"bgzip -f {writer.name}")
        elif writer.name.endswith(".vcf"):
            if writer is filtered_out_variants_vcf_writer:
                # sort the vcf
                run(f"cat {writer.name} | grep ^# > {writer.name}.sorted.vcf")
                run(f"cat {writer.name} | grep -v ^# | sort -k1,1V -k2,2n >> {writer.name}.sorted.vcf")
                os.rename(f"{writer.name}.sorted.vcf", writer.name)
            run(f"bgzip -f {writer.name}")
            run(f"tabix -f {writer.name}.gz")
        elif writer.name.endswith(".bed"):
            run(f"bedtools sort -i {bed_writer.name} | bgzip > {bed_writer.name}.gz")
            run(f"tabix -f {bed_writer.name}.gz")
            os.remove(bed_writer.name)

        print(f"Finished writing to {writer.name}.gz")

    if not args.dont_run_trf and not args.verbose:
        # delete the TRF fasta file directory
        trf_working_dir = get_trf_working_dir(args)
        shutil.rmtree(trf_working_dir)
        print(f"Deleted {trf_working_dir} directory")

    print_stats(counters)

def handle_overlapping_loci(variant_intervals, tsv_writer, vcf_writer, bed_writer, filtered_out_variants_vcf_writer, args, counters):

    # prepare to filter out detected TR loci that overlap one or more other variants in the VCF (ie. other TRs, non-TR indels, or SNVs).
    print(f"Starting overlap detection")
    locus_ids_with_overlapping_variants = {}
    chrom_intervals_iter = variant_intervals.items()
    if args.show_progress_bar:
        chrom_intervals_iter = tqdm.tqdm(chrom_intervals_iter, unit=" chromosomes", unit_scale=True)

    for chrom, intervals in chrom_intervals_iter:
        # create an interval tree for this chromosome
        interval_tree = intervaltree.IntervalTree(intervals)
        for interval in intervals:
            locus_id = interval.data["LocusId"]
            if args.debug: print(f"#### Checking overlap with {locus_id}")
            repeat_unit = interval.data["RepeatUnit"]
            overlapping_intervals = interval_tree.overlap(interval.begin - len(repeat_unit), interval.end + len(repeat_unit))

            assert len(overlapping_intervals) > 0, f"ERROR: Expected the locus interval to overlap with itself: {interval.data}"

            if len(overlapping_intervals) == 1:
                continue

            for overlapping_interval in overlapping_intervals:
                if overlapping_interval.data["LocusId"] == locus_id:
                    # don't compare to self
                    continue

                # Filter out TR that overlaps another TR locus that has the same canonical motif but wider locus boundaries.
                # This can happen when the locus was identified by detecting pure repeats, while the other (wider) locus was
                # detected via interrupted repeats or by using TRF
                overlapping_locus_is_wider = (interval.end - interval.begin) <= (overlapping_interval.end - overlapping_interval.begin)
                overlapping_locus_already_filtered_out = overlapping_interval.data["LocusId"] in locus_ids_with_overlapping_variants
                loci_have_similar_motif = interval.data["CanonicalMotif"] == overlapping_interval.data["CanonicalMotif"] or (
                        len(interval.data["CanonicalMotif"]) > 6
                        and len(interval.data["CanonicalMotif"]) == len(overlapping_interval.data["CanonicalMotif"])) or (
                                                  len(interval.data["CanonicalMotif"]) > 12
                                                  and len(interval.data["CanonicalMotif"]) - len(overlapping_interval.data["CanonicalMotif"]) <= 2
                                          )

                if args.debug: print(f"#### Checking overlap between {locus_id} and", overlapping_interval.data["LocusId"], ":",
                                     "loci_have_similar_motif =", loci_have_similar_motif,
                                     "overlapping_locus_is_wider =", overlapping_locus_is_wider,
                                     "overlapping_locus_already_filtered_out = ", overlapping_locus_already_filtered_out,

                                     )
                if  loci_have_similar_motif and overlapping_locus_is_wider and not overlapping_locus_already_filtered_out:
                    counters[f"variant filter: {FILTER_TR_LOCUS_THAT_HAS_OVERLAPPING_VARIANTS}"] += 1
                    locus_ids_with_overlapping_variants[locus_id] = FILTER_TR_LOCUS_THAT_HAS_OVERLAPPING_VARIANTS
                    break


    # remove locus_ids_with_overlapping_variants from all output files
    if len(locus_ids_with_overlapping_variants) > 0:
        for writer in filter(None, [tsv_writer, vcf_writer, bed_writer]):
            writer.close()
            #print(f"Filtering out {len(locus_ids_with_overlapping_variants):,d} loci from {writer.name} because they overlap "
            #      f"more than one variant in the input VCF")

            process_tandem_repeat_loci_that_have_overlapping_variants(
                writer.name, locus_ids_with_overlapping_variants, filtered_out_variants_vcf_writer)


def run_trf(args, trf_fasta_path):
    trf_runner = TRFRunner(args.trf_executable_path, html_mode=True,
                           min_motif_size=args.min_repeat_unit_length,
                           max_motif_size=args.max_repeat_unit_length)

    trf_runner.run_trf_on_fasta_file(os.path.basename(trf_fasta_path))


if __name__ == "__main__":
    args = parse_args()

    if not args.output_prefix:
        args.output_prefix = re.sub(".vcf(.gz)?", "", os.path.basename(args.input_vcf_path)) + ".TRs"

    variants_to_process_using_trf = None
    variants_per_trf_fasta = None
    if not args.dont_run_trf:
        start_time = datetime.datetime.now()

        variants_to_process_using_trf = {}
        variants_per_trf_fasta = collections.defaultdict(int)

        trf_working_dir = get_trf_working_dir(args)
        if not os.path.isdir(trf_working_dir):
            os.makedirs(trf_working_dir)

        main(args, only_generate_trf_fasta=True, variants_to_process_using_trf=variants_to_process_using_trf,
             variants_per_trf_fasta=variants_per_trf_fasta)

        # start multithreaded TRF processing
        original_working_dir = os.getcwd()
        os.chdir(trf_working_dir)

        n_threads = min(args.trf_threads, len(variants_to_process_using_trf))
        print(f"Launching {n_threads} TRF instance(s) to process {len(variants_to_process_using_trf):,d} insertion & deletion alleles")
        threads = []
        for thread_i in range(0, n_threads):
            trf_fasta_path = f"thread{thread_i}.fa"

            thread = threading.Thread(target=run_trf, args=(args, trf_fasta_path))
            thread.start()
            threads.append(thread)

        for thread in threads:
            thread.join()

        os.chdir(original_working_dir)

        print(f"Done running TRF. Running time: {datetime.datetime.now() - start_time}. Parsing results...")

    main(args, only_generate_trf_fasta=False, variants_to_process_using_trf=variants_to_process_using_trf,
         variants_per_trf_fasta=variants_per_trf_fasta)



