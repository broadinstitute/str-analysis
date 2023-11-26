#!/usr/bin/env python3

"""
This script takes a single-sample .vcf and filters it to insertions and deletions where either the REF or ALT allele
is a short tandem repeat (STR).
"""

import argparse
import collections
import gzip
import os
from pprint import pformat
import pyfaidx
import pysam
import re
import tqdm

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.find_repeat_unit import find_repeat_unit_allowing_interruptions
from str_analysis.utils.find_repeat_unit import find_repeat_unit_without_allowing_interruptions
from str_analysis.utils.find_repeat_unit import extend_repeat_into_sequence_allowing_interruptions
from str_analysis.utils.find_repeat_unit import extend_repeat_into_sequence_without_allowing_interruptions


COMMON_TSV_OUTPUT_COLUMNS = [
    "Chrom",
    "Start1Based",
    "End1Based",
    "Locus",
    "LocusId",
    "INS_or_DEL",
    "Motif",
    "MotifInterruptionIndex",
    "CanonicalMotif",
    "MotifSize",
    "NumRepeatsInReference",
    "VcfPos",
    "VcfRef",
    "VcfAlt",
    "VcfGenotype",
    "SummaryString",
    "IsFoundInReference",
    "IsPureRepeat",
    "IsMultiallelic",
]

VARIANT_TSV_OUTPUT_COLUMNS = COMMON_TSV_OUTPUT_COLUMNS + [
    "HET_or_HOM_or_HEMI_or_MULTI",
    "NumRepeatsShortAllele",
    "NumRepeatsLongAllele",
    "RepeatSizeShortAllele (bp)",
    "RepeatSizeLongAllele (bp)",
]

ALLELE_TSV_OUTPUT_COLUMNS = COMMON_TSV_OUTPUT_COLUMNS + [
    "NumRepeats",
    "RepeatSize (bp)",
    "NumPureRepeats",
    "PureRepeatSize (bp)",
    "FractionPureRepeats",
]

FILTER_MORE_THAN_TWO_ALT_ALLELES = "more than two alt alleles"
FILTER_UNEXPECTED_GENOTYPE_FORMAT = "unexpected genotype format"
FILTER_ZERO_ALT_ALLELES = "variant has zero non-* alt alleles"
FILTER_MULTIALLELIC_VARIANT_WITH_HOMOZYGOUS_GENOTYPE = "multiallelic variant with homozygous genotype"

FILTER_ALLELE_WITH_N_BASES = "contains N"
FILTER_ALLELE_SNV_OR_MNV = "SNV/MNV"
FILTER_ALLELE_MNV_INDEL = "complex multinucleotide insertion + deletion"
FILTER_ALLELE_NON_STR_INDEL = "INDEL without repeats"
FILTER_STR_ALLELE_NOT_ENOUGH_REPEATS = "is only %d repeats"
FILTER_STR_ALLELE_DOESNT_SPAN_ENOUGH_BASE_PAIRS = "spans < %d bp"
FILTER_STR_ALLELE_PARTIAL_REPEAT = "ends in partial repeat"

FILTER_STR_ALLELE_REPEAT_UNIT_TOO_SHORT = "repeat unit < %d bp"
FILTER_STR_ALLELE_REPEAT_UNIT_TOO_LONG = "repeat unit > %d bp"

FILTER_VARIANT_WITH_STR_ALLELES_WITH_DIFFERENT_MOTIFS = "STR alleles with different motifs"
FILTER_VARIANT_WITH_STR_ALLELES_WITH_DIFFERENT_INTERRUPTION_PATTERNS = "STR alleles with different interruption patterns"
FILTER_VARIANT_WITH_STR_ALLELES_WITH_DIFFERENT_COORDS = "STR alleles with different coords"
FILTER_STR_LOCUS_WITH_MULTIPLE_STR_VARIANTS = "locus overlaps more than one STR variant"


def parse_args():
    """Parse command-line arguments."""

    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument("-R", "--reference-fasta-path", help="Reference genome fasta path.", required=True)
    p.add_argument("--allow-interruptions", help="Whether to allow interruptions in the repeat sequence. There are 3 "
                   "options: 1) 'no' disallows interruptions and only looks for pure repeats. "
                   "2) 'only-if-pure-repeats-not-found' checks for interrupted repeats only for variants that don't "
                   "pass filters as pure repeats. 3) 'always' checks all variants for interrupted repeats, extending "
                   "their locus start and end coordinates to include interrupted repeats in the reference even when "
                   "the variant sequence contains only pure repeats.",
                   choices=["no", "only-if-pure-repeats-not-found", "always"], required=True)
    p.add_argument("--min-str-length", type=int, default=9, help="Minimum STR length in base pairs. This threshold "
                   "applies to the total repeat length comprising any repeats in the flanking sequence to the left "
                   "and right of the variant + in the inserted or deleted bases themselves")
    p.add_argument("--min-str-repeats", type=int, default=3, help="Minimum STR size in number of repeats. This "
                   "threshold applies to the total repeat length comprising any repeats in the flanking sequence to "
                   "the left and right of the variant + in the inserted or deleted bases themselves")
    p.add_argument("--min-repeat-unit-length", type=int, default=1, help="Minimum repeat unit length in base pairs.")
    p.add_argument("--max-repeat-unit-length", type=int, default=10**9, help="Max repeat unit length in base pairs.")
    p.add_argument("--show-progress-bar", help="Show a progress bar in the terminal when processing variants.",
                   action="store_true")
    p.add_argument("-v", "--verbose", help="Print detailed logs.", action="store_true")
    p.add_argument("-n", type=int, help="Only process the first N rows of the VCF. Useful for testing.")

    p.add_argument("-o", "--output-prefix", help="Output file prefix. If not specified, it will be computed based on "
                   "the input vcf filename")
    p.add_argument("-ik", "--copy-info-field-keys-to-tsv", help="Copy the values of these INFO fields from the input VCF to "
                                                                "the output TSV files.", action="append")
    p.add_argument("--write-bed-file", help="Whether to output a .bed file containing the STR variants. This requires "
                   "bedtools, bgzip and tabix tools to be available in the shell environment.",
                   action="store_true")

    p.add_argument("--write-vcf-with-filtered-out-variants", help="Output VCF files with variants that were filtered "
                   "out even though at least one allele was found to be an STR for reasons such as being multiallelic "
                   "and having alleles with different motifs, or one allele being an STR while the other is an SNV. "
                   "These alleles are filtered out to reduce complexity for downstream analyses.",
                   action="store_true")
    p.add_argument("-L", "--interval", help="Only process variants in this genomic interval (format: chrN:start-end)",
                   action="append")
    p.add_argument("input_vcf_path")

    args = p.parse_args()

    if not args.allow_interruptions != "no":
        # drop some output columns
        for header in VARIANT_TSV_OUTPUT_COLUMNS, ALLELE_TSV_OUTPUT_COLUMNS:
            for column in "NumPureRepeats", "PureRepeatSize (bp)", "FractionPureRepeats", "MotifInterruptionIndex":
                if column in header:
                    header.remove(column)

    if args.copy_info_field_keys_to_tsv:
        VARIANT_TSV_OUTPUT_COLUMNS.extend(args.copy_info_field_keys_to_tsv)
        ALLELE_TSV_OUTPUT_COLUMNS.extend(args.copy_info_field_keys_to_tsv)

    return args


def get_flanking_reference_sequences(fasta_obj, chrom, pos, ref, alt, num_flanking_bases=None):
    """Takes a single insertion or deletion and returns flanking reference sequence around the variant.

    Args:
        fasta_obj: pyfaidx wrapper of the reference genome fasta
        chrom (str): variant chromosome
        pos (str): variant position
        ref (str): variant ref allele
        alt (str): variant alt allele
        num_flanking_bases (int): (optional) num flanking bases to return. If not specified, it will be set to 5x the
            allele length.

    Return:
        3-tuple: left_flanking_reference_sequence,
            variant_bases (the bases that were inserted or deleted by this variant),
            right_flanking_reference_sequence
    """
    variant_bases = compute_indel_variant_bases(ref, alt)
    if variant_bases is None:
        raise ValueError(f"Invalid variant: {chrom}-{pos}-{ref}-{alt} is not an expansion or contraction.")

    if num_flanking_bases is None:
        num_flanking_bases = 5 * len(variant_bases)

    if chrom not in fasta_obj:
        raise ValueError(f"Chromosome '{chrom}' not found in reference fasta")

    left_flank_start_1based = max(0, pos - num_flanking_bases)
    left_flank_end = pos
    right_flank_start_1based = pos
    right_flank_end = min(pos + num_flanking_bases, len(fasta_obj[chrom]))
    if len(ref) > len(alt):
        # For deletions, the right-flanking sequence starts to the right of the variant bases
        right_flank_start_1based += len(variant_bases)
        right_flank_end = min(right_flank_end + len(variant_bases), len(fasta_obj[chrom]))

    left_flanking_reference_sequence = str(fasta_obj[chrom][left_flank_start_1based : left_flank_end]).upper()
    right_flanking_reference_sequence = str(fasta_obj[chrom][right_flank_start_1based : right_flank_end]).upper()

    return left_flanking_reference_sequence, variant_bases, right_flanking_reference_sequence


def determine_reason_indel_allele_failed_str_filter(
        variant_bases, repeat_unit, min_str_length, min_repeat_unit_length, max_repeat_unit_length,
        num_total_repeats_in_str, allow_interruptions, counters):
    """Utility function for refining the reason why a given indel allele did not pass the STR filter.

    Return:
        str: A string describing the reason why the allele was filtered out. If the allele passed the filter, returns None.
    """
    if len(repeat_unit) < min_repeat_unit_length:
        counters[f"allele filter: repeat unit is shorter than {min_repeat_unit_length}bp"] += 1
        return FILTER_STR_ALLELE_REPEAT_UNIT_TOO_SHORT % min_repeat_unit_length
    if len(repeat_unit) > max_repeat_unit_length:
        counters[f"allele filter: repeat unit is longer than {max_repeat_unit_length}bp"] += 1
        return FILTER_STR_ALLELE_REPEAT_UNIT_TOO_LONG % max_repeat_unit_length

    if num_total_repeats_in_str > 1:
        # it has more than one repeat, so a repeat unit was found, but it was less than the minimum threshold
        if num_total_repeats_in_str * len(repeat_unit) < min_str_length:
            counters[f"allele filter: allele sequence spans < {min_str_length}bp"] += 1
            return FILTER_STR_ALLELE_DOESNT_SPAN_ENOUGH_BASE_PAIRS % min_str_length
        else:
            counters[f"allele filter: allele consists of only {num_total_repeats_in_str} repeats"] += 1
            return FILTER_STR_ALLELE_NOT_ENOUGH_REPEATS % num_total_repeats_in_str

    # check whether it consists of some repeats if it's allowed to end in a partial repeat
    if allow_interruptions:
        repeat_unit_allowing_partial, _, num_total_repeats_allowing_partial, _, has_partial_repeats = find_repeat_unit_allowing_interruptions(
            variant_bases,
            allow_partial_repeats=True,
        )
    else:
        repeat_unit_allowing_partial, num_total_repeats_allowing_partial, has_partial_repeats = find_repeat_unit_without_allowing_interruptions(
            variant_bases,
            allow_partial_repeats=True,
        )

    if num_total_repeats_allowing_partial > 1:
        # found that it consists of some repeats if it's allowed to end in a partial repeat
        if not has_partial_repeats:
            # sanity check
            raise ValueError(f"Unexpected return value: has_partial_repeats is False for {variant_bases}, "
                             f"allow_interruptions={allow_interruptions}")
        counters[f"allele filter: {FILTER_STR_ALLELE_PARTIAL_REPEAT}"] += 1
        return FILTER_STR_ALLELE_PARTIAL_REPEAT

    # no repeat unit found in this allele
    counters[f"allele filter: {FILTER_ALLELE_NON_STR_INDEL}"] += 1
    return FILTER_ALLELE_NON_STR_INDEL


def compute_indel_variant_bases(ref, alt):
    if len(ref) < len(alt):
        if alt.startswith(ref):
            return alt[len(ref):]
        else:
            return None
    elif len(alt) < len(ref):
        if ref.startswith(alt):
            return ref[len(alt):]
        else:
            return None
    else:
        raise ValueError(f"Invalid ref, alt allele pair: {ref} > {alt}")


def check_if_allele_is_str(
    fasta_obj, chrom, pos, ref, alt,
    min_str_repeats, min_str_length, min_repeat_unit_length, max_repeat_unit_length,
    allow_interruptions=False,
    counters=None,
):
    """Determine if the given allele (represented by chrom/pos/ref/alt) is an STR expansion or contraction or neither.

    Args:
        fasta_obj (object): pyfasta object for accessing reference sequence.
        chrom (str): chromosome
        pos (int): 1-based position
        ref (str): VCF ref sequence
        alt (str): VCF alt sequence
        min_str_repeats (int): The min number of repeats that must be found in the variant sequence + flanking
            reference sequences for the variant to be considered an STR.
        min_str_length (int): The repeats in the variant sequence + the flanking reference sequence must cover at least
            this many base-pairs for the variant to be considered an STR.
        min_repeat_unit_length (int): If the repeat unit is shorter than this, the variant will be filtered out.
        max_repeat_unit_length (int): If the repeat unit is longer than this, the variant will be filtered out.
        allow_interruptions (bool): Whether to allow interruptions in the repeat sequence
        counters (dict): Dictionary of counters to collect summary stats about the number of STR variants found, etc.

    Return:
        dict: Fields describing the results of checking whether this variant is an STR
    """
    null_result = {
        "RepeatUnit": None,
        "FilterReason": None,
    }

    if "N" in ref.upper() or "N" in alt.upper():
        counters[f"allele filter: N bases"] += 1
        null_result["FilterReason"] = FILTER_ALLELE_WITH_N_BASES
        return null_result

    if set(alt.upper()) - {"A", "C", "G", "T"}:
        counters[f"allele filter: non-ACGT bases"] += 1

    if len(ref) == len(alt):
        counters[f"allele filter: "+("SNV" if len(alt) == 1 else "MNV")] += 1
        null_result["FilterReason"] = FILTER_ALLELE_SNV_OR_MNV
        return null_result

    variant_bases = compute_indel_variant_bases(ref, alt)
    if variant_bases is None:
        counters[f"allele filter: complex MNV indel"] += 1
        null_result["FilterReason"] = FILTER_ALLELE_MNV_INDEL
        return null_result

    if len(ref) == len(alt):
        raise ValueError(f"Unexpected ref, alt pair: len('{ref}') == len('{alt}')")
    elif len(ref) < len(alt):
        ins_or_del = "INS"
    else:
        ins_or_del = "DEL"
    counters[f"allele counts: {ins_or_del} alleles"] += 1

    left_flanking_reference_sequence, variant_bases2, right_flanking_reference_sequence = get_flanking_reference_sequences(
        fasta_obj, chrom, pos, ref, alt, num_flanking_bases=2000)

    assert variant_bases == variant_bases2

    if allow_interruptions:
        # check for interrupted repeats
        (
            repeat_unit,
            num_pure_repeats_within_variant_bases,
            num_total_repeats_within_variant_bases,
            repeat_unit_interruption_index,
            has_partial_repeats
        ) = find_repeat_unit_allowing_interruptions(variant_bases)

        reversed_repeat_unit_interruption_index = (len(repeat_unit) - 1 - repeat_unit_interruption_index) if repeat_unit_interruption_index is not None else None

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

    else:
        # check for pure repeats
        (
            repeat_unit,
            num_pure_repeats_within_variant_bases,
            has_partial_repeats,
        ) = find_repeat_unit_without_allowing_interruptions(variant_bases)

        repeat_unit_interruption_index = None
        num_total_repeats_within_variant_bases = num_pure_repeats_within_variant_bases

        num_pure_repeats_left_flank = num_total_repeats_left_flank = extend_repeat_into_sequence_without_allowing_interruptions(
            repeat_unit[::-1],
            left_flanking_reference_sequence[::-1])

        num_pure_repeats_right_flank = num_total_repeats_right_flank = extend_repeat_into_sequence_without_allowing_interruptions(
            repeat_unit,
            right_flanking_reference_sequence)

    # even though the VCF position is 1-based, it represents the location of the base preceding the variant bases, so
    # add 1 to get the 1-based position of the true base
    start_1based = pos + 1 - num_total_repeats_left_flank * len(repeat_unit)
    end_1based = pos + (len(variant_bases) if len(ref) > len(alt) else 0) + num_total_repeats_right_flank * len(repeat_unit)

    pure_start_1based = pos + 1 - num_pure_repeats_left_flank * len(repeat_unit)
    pure_end_1based = pos + (len(variant_bases) if len(ref) > len(alt) else 0) + num_pure_repeats_right_flank * len(repeat_unit)

    num_pure_repeats_ref = num_pure_repeats_alt = num_pure_repeats_left_flank + num_pure_repeats_right_flank
    num_total_repeats_ref = num_total_repeats_alt = num_total_repeats_left_flank + num_total_repeats_right_flank
    if len(ref) < len(alt):
        # insertion
        num_pure_repeats_alt += num_pure_repeats_within_variant_bases
        num_total_repeats_alt += num_total_repeats_within_variant_bases
    else:
        # deletion
        num_pure_repeats_ref += num_pure_repeats_within_variant_bases
        num_total_repeats_ref += num_total_repeats_within_variant_bases

    num_pure_repeats_in_variant_plus_flanks = (
        num_pure_repeats_left_flank + num_pure_repeats_within_variant_bases + num_pure_repeats_right_flank)
    num_total_repeats_in_variant_plus_flanks = (
        num_total_repeats_left_flank + num_total_repeats_within_variant_bases + num_total_repeats_right_flank)

    if num_total_repeats_in_variant_plus_flanks < min_str_repeats \
            or num_total_repeats_in_variant_plus_flanks * len(repeat_unit) < min_str_length \
            or len(repeat_unit) < min_repeat_unit_length \
            or len(repeat_unit) > max_repeat_unit_length:
        # this allele didn't pass filters. Determine the detailed reason it's being filtered out.
        null_result["FilterReason"] = determine_reason_indel_allele_failed_str_filter(
            variant_bases, repeat_unit, min_str_length, min_repeat_unit_length, max_repeat_unit_length,
            num_total_repeats_in_variant_plus_flanks, allow_interruptions, counters)
        return null_result

    result = {
        "Chrom": chrom,
        "Pos": pos,
        "Ref": ref,
        "Alt": alt,
        "Start1Based": start_1based,
        "End1Based": end_1based,
        "RepeatUnit": repeat_unit,
        "NumRepeatsRef": num_total_repeats_ref,
        "NumRepeatsAlt": num_total_repeats_alt,
        "NumRepeatsLeftFlank": num_total_repeats_left_flank,
        "NumRepeatsRightFlank": num_total_repeats_right_flank,
        "NumRepeatsInVariant": num_total_repeats_within_variant_bases,
        "IsPureRepeat": num_pure_repeats_in_variant_plus_flanks == num_total_repeats_in_variant_plus_flanks,
        "FractionPureRepeatsRef": num_pure_repeats_ref / num_total_repeats_ref if num_total_repeats_ref > 0 else None,
        "FractionPureRepeatsAlt": num_pure_repeats_alt / num_total_repeats_alt if num_total_repeats_alt > 0 else None,
        "MotifInterruptionIndex": repeat_unit_interruption_index,
        "PureStart1Based": pure_start_1based,
        "PureEnd1Based": pure_end_1based,
        "NumPureRepeatsRef": num_pure_repeats_ref,
        "NumPureRepeatsAlt": num_pure_repeats_alt,
        "NumPureRepeatsLeftFlank": num_pure_repeats_left_flank,
        "NumPureRepeatsRightFlank": num_pure_repeats_right_flank,
        "NumPureRepeatsInVariant": num_pure_repeats_within_variant_bases,
        "FilterReason": None,
    }

    if result["IsPureRepeat"] and result["MotifInterruptionIndex"] is not None:
        raise ValueError(f"Allele has IsPureRepeat set to True, but MotifInterruptionIndex is not None:" + pformat(result))

    # update counters
    if counters:
        if num_total_repeats_left_flank > 0 and num_total_repeats_right_flank > 0:
            left_or_right = 'both left and right'
        elif num_total_repeats_left_flank > 0:
            left_or_right = 'left'
        elif num_total_repeats_right_flank > 0:
            left_or_right = 'right'
        else:
            left_or_right = 'no'

        if len(variant_bases) < 500:
            num_base_pairs_within_variant_bases = f"{25*int(len(variant_bases)/25)}-{25*(1+int(len(variant_bases)/25))}bp"
        else:
            num_base_pairs_within_variant_bases = "500+bp"

        counters[f"STR allele counts: TOTAL"] += 1
        counters[f"STR allele counts: {ins_or_del}"] += 1
        counters[f"STR allele delta: {num_total_repeats_within_variant_bases if num_total_repeats_within_variant_bases < 9 else '9+'} repeats"] += 1
        counters[f"STR allele motif size: {len(repeat_unit) if len(repeat_unit) < 9 else '9+'} bp"] += 1
        counters[f"STR allele size: {num_base_pairs_within_variant_bases}"] += 1
        counters[f"STR allele reference repeats: with {left_or_right} matching ref. repeat"] += 1

    return result


def check_if_variant_is_str(
        fasta_obj,
        vcf_chrom,
        vcf_pos,
        vcf_ref,
        alt_alleles,
        min_str_repeats,
        min_str_length,
        min_repeat_unit_length,
        max_repeat_unit_length,
        counters,
        interruptions="no",
        vcf_line_i=None,
        verbose=False,
):

    # first check whether the alleles are both pure STRs
    if interruptions != "always":
        current_counters = collections.defaultdict(int)
        str_specs_with_pure_repeats = [check_if_allele_is_str(
            fasta_obj, vcf_chrom, vcf_pos, vcf_ref, alt_allele,
            min_str_repeats=min_str_repeats, min_str_length=min_str_length,
            min_repeat_unit_length=min_repeat_unit_length, max_repeat_unit_length=max_repeat_unit_length,
            allow_interruptions=False,
            counters=current_counters,
        ) for alt_allele in alt_alleles]

        str_specs_with_pure_repeats, variant_filter_string = postprocess_str_variant(
            vcf_line_i, str_specs_with_pure_repeats, allow_interruptions=False, counters=current_counters, verbose=verbose)

        if interruptions == "no" or not variant_filter_string:
            # found that both alleles are pure STRs, so no need to check for repeats with interruptions
            # if interruptions are not allowed, then return the pure STR results
            for counter_key, count in current_counters.items():
                counters[counter_key] += count

            return str_specs_with_pure_repeats, variant_filter_string

    # check again, now allowing interruptions
    current_counters = collections.defaultdict(int)
    str_specs_allowing_interruptions = [check_if_allele_is_str(
        fasta_obj, vcf_chrom, vcf_pos, vcf_ref, alt_allele,
        min_str_repeats=min_str_repeats, min_str_length=min_str_length,
        min_repeat_unit_length=min_repeat_unit_length, max_repeat_unit_length=max_repeat_unit_length,
        allow_interruptions=True,
        counters=current_counters,
    ) for alt_allele in alt_alleles]

    str_specs_allowing_interruptions, variant_filter_string = postprocess_str_variant(
        vcf_line_i, str_specs_allowing_interruptions, allow_interruptions=True, counters=current_counters, verbose=verbose)

    for counter_key, count in current_counters.items():
        counters[counter_key] += count

    return str_specs_allowing_interruptions, variant_filter_string


def postprocess_str_variant(vcf_line_i, str_allele_specs, allow_interruptions=False, counters=None, verbose=False):
    """Check whether the STR alleles are valid, or whether they should be filtered out.

    Return:
        str: the reason a variant should be filtered out, or None if the variant should not be filtered out
    """
    if len(str_allele_specs) == 1:
        # handle mono-allelic (rather than multi-allelic) variants
        if str_allele_specs[0]["FilterReason"] is not None:
            counters[f"variant filter: " + str_allele_specs[0]["FilterReason"]] += 1
            
        return str_allele_specs, str_allele_specs[0]["FilterReason"]

    # Filter out multi-allelics where one allele is an STR and the other is not
    if any(allele_spec["RepeatUnit"] is None for allele_spec in str_allele_specs):
        if str_allele_specs[0]["FilterReason"] is None or str_allele_specs[1]["FilterReason"] is None:
            counters[f"variant filter: multi-allelic: one STR, one non-STR allele"] += 1
        elif str_allele_specs[0]["FilterReason"] == str_allele_specs[1]["FilterReason"]:
            counters[f"variant filter: multi-allelic: " + str_allele_specs[0]["FilterReason"]] += 1
        else:
            counters[f"variant filter: multi-allelic: both alleles fail STR filters for different reasons"] += 1

        return str_allele_specs, ";".join((allele_spec["FilterReason"] or f"STR allele") for allele_spec in str_allele_specs)

    if str_allele_specs[0]["IsPureRepeat"] != str_allele_specs[1]["IsPureRepeat"]:
        str_allele_specs[0]["IsPureRepeat"] = str_allele_specs[1]["IsPureRepeat"] = False

        if str_allele_specs[0]["MotifInterruptionIndex"] is None:
            str_allele_specs[0]["MotifInterruptionIndex"] = str_allele_specs[1]["MotifInterruptionIndex"]
        elif str_allele_specs[1]["MotifInterruptionIndex"] is None:
            str_allele_specs[1]["MotifInterruptionIndex"] = str_allele_specs[0]["MotifInterruptionIndex"]

    # Filter out multi-allelics where the STR alleles have different repeat units or coordinates
    for attribute, filter_reason, error_when_different in [
        ("Chrom", "", True),
        ("IsPureRepeat", "", True),
        ("RepeatUnit", FILTER_VARIANT_WITH_STR_ALLELES_WITH_DIFFERENT_MOTIFS, False),
        ("MotifInterruptionIndex", FILTER_VARIANT_WITH_STR_ALLELES_WITH_DIFFERENT_INTERRUPTION_PATTERNS, True if not allow_interruptions else False),
        ("Start1Based", FILTER_VARIANT_WITH_STR_ALLELES_WITH_DIFFERENT_COORDS, False),
        ("End1Based", FILTER_VARIANT_WITH_STR_ALLELES_WITH_DIFFERENT_COORDS, False),
    ]:
        if str_allele_specs[0][attribute] == str_allele_specs[1][attribute]:
            continue

        chrom = str_allele_specs[0]['Chrom'].replace('chr', '')

        variant_id = f"{chrom}-{str_allele_specs[0]['Start1Based']}-{str_allele_specs[0]['End1Based']}-"
        variant_id += f"{str_allele_specs[0]['Ref']}-"
        variant_id += f"{str_allele_specs[0]['Alt']} (RU:{str_allele_specs[0]['RepeatUnit']}"
        variant_id += f":{str_allele_specs[0]['MotifInterruptionIndex']})" if allow_interruptions else ")"
        variant_id += f", {str_allele_specs[1]['Alt']} (RU:{str_allele_specs[1]['RepeatUnit']}"
        variant_id += f":{str_allele_specs[1]['MotifInterruptionIndex']})" if allow_interruptions else ")"

        if error_when_different:
            raise ValueError(f"STR alleles have different {attribute} values: "
                             f"{str_allele_specs[0][attribute]} vs {str_allele_specs[1][attribute]}")

        if verbose:
            print(f"WARNING: vcf row #{vcf_line_i:,d}: {variant_id}: {'interrupted' if allow_interruptions else 'pure'} "
                  f"multi-allelic STR has alleles with different {attribute}:",
                  str_allele_specs[0][attribute], " vs ", str_allele_specs[1][attribute], "  Skipping...")

        counters[f"variant filter: multi-allelic: has {filter_reason}"] += 1

        filter_string = ";".join([filter_reason] * len(str_allele_specs))  # add filter reason to both alleles
        return str_allele_specs, filter_string

    return str_allele_specs, None


def is_found_in_reference(alt_STR_allele):
    """Return True if the allele is a deletion, since it by definition means some of the repeats were present in the
    reference genome. Otherwise, return True if the allele had STRs in the left or right flanking sequence.
    """
    is_deletion = len(alt_STR_allele["Ref"]) > len(alt_STR_allele["Alt"])
    if is_deletion:
        return True

    num_repeats_in_left_or_right_flank = sum([alt_STR_allele[f"NumRepeats{s}Flank"] for s in ("Left", "Right")])
    if num_repeats_in_left_or_right_flank > 0:
        return True

    return False


def get_num_repeats_in_allele(str_allele_specs, genotype_index):
    """Looks up the number of STR repeats found in a given VCF record's allele(s) based on the genotype.

    Args:
        str_allele_specs (list): list containing the STR allele spec(s) returned by check_if_allele_is_str for each
            of the allele(s) in this variant.
        genotype_index (int): which genotype to look up.

    Return
        int: The number of repeats in this allele
    """
    if not isinstance(str_allele_specs, list):
        raise ValueError(f"str_allele_specs argument is not of type list: {str_allele_specs}")
    if not isinstance(genotype_index, int):
        raise ValueError(f"genotype_index argument is not of type int: {genotype_index}")

    if genotype_index == 0:
        return str_allele_specs[0]["NumRepeatsRef"]

    if genotype_index - 1 >= len(str_allele_specs):
        raise ValueError(f"genotype_index {genotype_index} is larger than the str_allele_specs list: {str_allele_specs}")

    return str_allele_specs[genotype_index - 1]["NumRepeatsAlt"]


def compute_variant_summary_string(str_allele_specs, het_or_hom_or_multi):
    """Returns a short easy-to-read string summarizing the STR variant.

    Args:
        str_allele_specs (list): list of 1 or more allele spec dictionaries
        het_or_hom_or_multi (str): describes the genotype as "HET" or "HOM" or "MULTI" meaning multi-allelic
    Return:
        str: short summary of the variant
    """

    repeat_unit = str_allele_specs[0]["RepeatUnit"]
    summary_string = "{}bp:{}:".format(len(repeat_unit), repeat_unit)

    ins_or_del = []
    for i in range(0, len(str_allele_specs)):
        if len(str_allele_specs[i]["Ref"]) > len(str_allele_specs[i]["Alt"]):
            ins_or_del.append("DEL")
        elif len(str_allele_specs[i]["Ref"]) < len(str_allele_specs[i]["Alt"]):
            ins_or_del.append("INS")
        elif len(str_allele_specs[i]["Ref"]) == len(str_allele_specs[i]["Alt"]):
            raise ValueError(f"Unexpected ref, alt pair: len('{str_allele_specs[i]['Ref']}') == len('{str_allele_specs[i]['Alt']}')")

    summary_string += ",".join(ins_or_del) + ":"
    summary_string += str(str_allele_specs[i]["NumRepeatsRef"]) + "=>" + ",".join([
        str(str_allele_specs[i]["NumRepeatsAlt"]) for i in range(0, len(str_allele_specs))
    ])
    summary_string += ":" + het_or_hom_or_multi

    if str_allele_specs[0]["IsPureRepeat"]:
        summary_string += ":pure"
    else:
        summary_string += ":not-pure"

    return ":".join(ins_or_del), summary_string


def process_vcf_line(
        vcf_line_i,
        vcf_fields,
        fasta_obj,
        vcf_writer,
        variants_tsv_writer,
        alleles_tsv_writer,
        bed_writer,
        args,
        counters,
        indels_per_locus_counter,
):
    """Utility method for processing a single line from the input vcf

    Args:
        vcf_line_i (int): line number of the VCF file
        vcf_fields (list): list of fields from the VCF line
        fasta_obj (pyfaidx.Fasta): Fasta object for the reference genome
        vcf_writer (vcf.Writer): VCF writer object
        variants_tsv_writer (csv.DictWriter): TSV writer object for the variants file
        alleles_tsv_writer (csv.DictWriter): TSV writer object for the alleles file
        bed_writer (csv.DictWriter): BED writer object
        args (argparse.Namespace): parsed command-line arguments
        counters (dict): dictionary of counters
        indels_per_locus_counter (collections.Counter): counter for the number of indels per locus

    Return:
        str: a string explaining the reason this variant is not a tandem repeat, or None if the variant is a tandem repeat
    """

    # parse the ALT allele(s)
    vcf_chrom = vcf_fields[0]
    vcf_chrom_without_chr_prefix = vcf_chrom.replace("chr", "")
    vcf_pos = int(vcf_fields[1])
    vcf_ref = vcf_fields[3].upper()
    vcf_alt = vcf_fields[4].upper()
    if not vcf_alt:
        raise ValueError(f"No ALT allele found in line {vcf_fields}")

    variant_id = f"{vcf_chrom_without_chr_prefix}-{vcf_pos}-{vcf_ref}-{vcf_alt}"

    if args.interval:
        print(f"Processing variant {variant_id}...")

    alt_alleles = vcf_alt.split(",")

    counters["variant counts: TOTAL variants"] += 1
    counters["allele counts: TOTAL alleles"] += len(alt_alleles)

    if len(alt_alleles) > 2:
        if args.verbose:
            print(f"WARNING: vcf row #{vcf_line_i:,d}: {variant_id}: multi-allelic variant has "
                  f"{len(alt_alleles)} alt alleles. This script doesn't support more than 2 alt alleles. Skipping...")
        counters[f"WARNING: multi-allelic variant with {len(alt_alleles)} alt alleles"] += 1
        return FILTER_MORE_THAN_TWO_ALT_ALLELES

    # parse the genotype
    vcf_genotype_format = vcf_fields[8].split(":")
    if vcf_genotype_format[0] != "GT":
        if args.verbose:
            print(f"WARNING: vcf row #{vcf_line_i:,d}: Unexpected genotype field in row #{vcf_line_i:,d}: "
                  f"{vcf_fields[8]}  (variant: {variant_id})")
        counters[f"WARNING: {FILTER_UNEXPECTED_GENOTYPE_FORMAT}"] += 1
        return FILTER_UNEXPECTED_GENOTYPE_FORMAT

    vcf_genotype = next(iter(vcf_fields[9].split(":")))
    vcf_genotype_indices = re.split(r"[/|\\]", vcf_genotype)
    if len(vcf_genotype_indices) not in (1, 2) or any(not gi.isdigit() for gi in vcf_genotype_indices if gi != "."):
        if args.verbose:
            print(f"WARNING: vcf row #{vcf_line_i:,d}:  Unexpected genotype GT format in row #{vcf_line_i:,d}: "
                  f"{vcf_genotype}  (variant: {variant_id})")
        counters[f"WARNING: {FILTER_UNEXPECTED_GENOTYPE_FORMAT}"] += 1
        return FILTER_UNEXPECTED_GENOTYPE_FORMAT

    vcf_genotype_indices = [int(gi) for gi in vcf_genotype_indices if gi != "."]
    vcf_genotype_separator = "|" if "|" in vcf_genotype else ("\\" if "\\" in vcf_genotype else "/")
    vcf_genotype = vcf_genotype_separator.join(map(str, vcf_genotype_indices))

    # check that genotype indices correctly correspond to alt allele(s)
    for genotype_index in vcf_genotype_indices:
        try:
            genotype_index = int(genotype_index)
        except ValueError:
            raise ValueError(f"Unexpected genotype field has non-numeric genotype indices {vcf_genotype_indices}")
        if genotype_index > len(alt_alleles):
            raise ValueError(f"Unexpected genotype field has a genotype index that is larger than the number of "
                             f"alt alleles {alt_alleles}: {vcf_genotype_indices}")

    # handle '*' alleles
    if "*" in alt_alleles:
        #counters["variant counts: removed * allele and converted to homozygous genotype"] += 1
        # if this variant has 1 regular allele and 1 "*" allele (which represents an overlapping deletion), discard the
        # "*" allele and recode the genotype as haploid
        star_allele_index = alt_alleles.index("*") + 1
        alt_alleles = [a for a in alt_alleles if a != "*"]
        vcf_genotype_indices = [gi for gi in vcf_genotype_indices if gi != star_allele_index]
        vcf_genotype = vcf_genotype_separator.join(map(str, vcf_genotype_indices))

    if len(alt_alleles) == 0 or len(vcf_genotype_indices) == 0:
        return FILTER_ZERO_ALT_ALLELES

    # confirm that there are no more "*" alleles
    if "*" in alt_alleles:
        raise ValueError(f"Unexpected '*' allele in row #{vcf_line_i:,d}: {variant_id}")

    is_homozygous = len(vcf_genotype_indices) == 2 and vcf_genotype_indices[0] == vcf_genotype_indices[1]
    is_hemizygous = len(vcf_genotype_indices) == 1
    het_or_hom_or_hemi_or_multi = "HOM" if is_homozygous else ("MULTI" if len(alt_alleles) > 1 else ("HEMI" if is_hemizygous else "HET"))

    # since these variants are from a single sample VCF, a homozygous genotype and multiple alleles are
    # an error of some sort.
    if len(alt_alleles) > 1 and len(set(vcf_genotype_indices)) == 1:
        return FILTER_MULTIALLELIC_VARIANT_WITH_HOMOZYGOUS_GENOTYPE

    str_allele_specs, filter_string = check_if_variant_is_str(
        fasta_obj,
        vcf_chrom,
        vcf_pos,
        vcf_ref,
        alt_alleles,
        args.min_str_repeats,
        args.min_str_length,
        args.min_repeat_unit_length,
        args.max_repeat_unit_length,
        counters=counters,
        interruptions=args.allow_interruptions,
        vcf_line_i=vcf_line_i,
        verbose=args.verbose,
    )

    if filter_string:
        return filter_string

    variant_ins_or_del, variant_summary_string = compute_variant_summary_string(
        str_allele_specs, het_or_hom_or_hemi_or_multi)

    # get repeat unit, start, end coords for the variant
    repeat_unit = str_allele_specs[0]["RepeatUnit"]
    start_1based = str_allele_specs[0]["Start1Based"]
    end_1based = str_allele_specs[0]["End1Based"]
    is_pure_repeat = str_allele_specs[0]["IsPureRepeat"]
    repeat_unit_interruption_index = str_allele_specs[0]["MotifInterruptionIndex"]

    # get allele sizes based on the variant's genotype in the VCF
    try:
        num_repeats_in_allele1 = get_num_repeats_in_allele(str_allele_specs, vcf_genotype_indices[0])
        if len(vcf_genotype_indices) == 2:
            num_repeats_in_allele2 = get_num_repeats_in_allele(str_allele_specs, vcf_genotype_indices[1])
            num_repeats_short_allele = min(num_repeats_in_allele1, num_repeats_in_allele2)
            num_repeats_long_allele  = max(num_repeats_in_allele1, num_repeats_in_allele2)
        else:
            num_repeats_short_allele = num_repeats_in_allele1
            num_repeats_long_allele  = ""
    except Exception as e:
        raise ValueError(f"{e} {vcf_ref} {alt_alleles} {vcf_genotype}")


    counters["STR variant counts: TOTAL"] += 1

    locus_id = f"{vcf_chrom_without_chr_prefix}-{start_1based - 1}-{end_1based}-{repeat_unit}"
    indels_per_locus_counter[locus_id] += 1

    num_repeats_in_reference = (end_1based - start_1based + 1)/len(repeat_unit)

    # write results to a VCF file
    vcf_fields[2] = variant_summary_string
    vcf_fields[4] = ",".join(alt_alleles)

    info_field_dict = {}
    if vcf_fields[7] and vcf_fields[7] != ".":
        for info_key_value in vcf_fields[7].split(";"):
            info_field_tokens = info_key_value.split("=")
            if len(info_field_tokens) > 1:
                info_field_dict[info_field_tokens[0]] = info_field_tokens[1]
            else:
                info_field_dict[info_field_tokens[0]] = True

    info_field_dict.update({
        f"LocusId": f"{locus_id}",
        f"Locus": f"{vcf_chrom_without_chr_prefix}:{start_1based}-{end_1based}",
        f"Motif": f"{repeat_unit}",
        f"NumRepeatsShortAllele": f"{num_repeats_short_allele}",
        f"NumRepeatsLongAllele": f"{num_repeats_long_allele}",
        f"NumRepeatsInReference": f"{num_repeats_in_reference}",
        f"IsPureRepeat": f"{is_pure_repeat}",
    })

    if not is_pure_repeat:
        info_field_dict["MotifInterruptionIndex"] = str(repeat_unit_interruption_index)

    vcf_fields[7] = ";".join([f"{key}={value}" for key, value in info_field_dict.items()])
    vcf_fields[9] = vcf_genotype
    vcf_writer.write("\t".join(vcf_fields) + "\n")

    # remove "NumRepeatsShortAllele" and "NumRepeatsLongAllele" keys from info_fields_dict since it doesn't make sense
    # to add them to the alleles tsv
    del info_field_dict["NumRepeatsShortAllele"]
    del info_field_dict["NumRepeatsLongAllele"]

    # write results to TSVs
    if len(alt_alleles) > 1:
        counters[f"STR variant counts: multi-allelic"] += 1

    tsv_record = info_field_dict
    tsv_record.update({
        "Chrom": vcf_chrom,
        "Start1Based": start_1based,
        "End1Based": end_1based,
        "Locus": f"{vcf_chrom_without_chr_prefix}:{start_1based}-{end_1based}",
        "LocusId": locus_id,
        "NumRepeatsInReference": num_repeats_in_reference,
        "Motif": repeat_unit,
        "CanonicalMotif": compute_canonical_motif(repeat_unit, include_reverse_complement=True),
        "MotifSize": len(repeat_unit),
        "VcfPos": vcf_pos,
        "VcfRef": vcf_ref,
        "VcfGenotype": vcf_genotype,
        "IsPureRepeat": is_pure_repeat,
        "MotifInterruptionIndex": repeat_unit_interruption_index if repeat_unit_interruption_index is not None else "",
        "IsMultiallelic": len(alt_alleles) > 1,
        "IsFoundInReference": any(is_found_in_reference(spec) for spec in str_allele_specs),
    })

    if num_repeats_short_allele < 0 or (num_repeats_long_allele and num_repeats_long_allele < 0):
        raise ValueError(f"Short or long allele size is < 0: "
                         f"{num_repeats_short_allele}, {num_repeats_long_allele}  {pformat(tsv_record)}")

    if is_pure_repeat:
        counters[f"STR variant counts: pure repeats"] += 1
    else:
        counters[f"STR variant counts: interrupted repeats"] += 1

    variant_tsv_record = dict(tsv_record)
    variant_tsv_record.update({
        "VcfAlt": ",".join(alt_alleles),
        "INS_or_DEL": variant_ins_or_del,
        "HET_or_HOM_or_HEMI_or_MULTI": het_or_hom_or_hemi_or_multi,
        "SummaryString": variant_summary_string,
        "NumRepeatsShortAllele": num_repeats_short_allele,
        "NumRepeatsLongAllele": num_repeats_long_allele,
        "RepeatSizeShortAllele (bp)": num_repeats_short_allele * len(repeat_unit),
        "RepeatSizeLongAllele (bp)": num_repeats_long_allele * len(repeat_unit) if num_repeats_long_allele != "" else "",
    })
    variants_tsv_writer.write("\t".join([str(variant_tsv_record.get(c, "")) for c in VARIANT_TSV_OUTPUT_COLUMNS]) + "\n")

    if bed_writer is not None:
        bed_writer.write("\t".join(map(str, [vcf_chrom, start_1based - 1, end_1based, variant_summary_string, "."])) + "\n")

    for alt_allele, alt_STR_allele_spec in zip(alt_alleles, str_allele_specs):
        allele_tsv_record = dict(tsv_record)
        ins_or_del, summary_string = compute_variant_summary_string([alt_STR_allele_spec], het_or_hom_or_hemi_or_multi)
        allele_tsv_record.update({
            "VcfAlt": alt_allele,
            "INS_or_DEL": ins_or_del,
            "SummaryString": summary_string,
            "NumRepeats": alt_STR_allele_spec["NumRepeatsAlt"],
            "RepeatSize (bp)": alt_STR_allele_spec["NumRepeatsAlt"] * len(repeat_unit),
            "NumPureRepeats": alt_STR_allele_spec["NumPureRepeatsAlt"],
            "PureRepeatSize (bp)": alt_STR_allele_spec["NumPureRepeatsAlt"] * len(repeat_unit),
            "FractionPureRepeats": ("%0.3f" % alt_STR_allele_spec["FractionPureRepeatsAlt"])
                                   if alt_STR_allele_spec["FractionPureRepeatsAlt"] is not None else "",
        })

        alleles_tsv_writer.write("\t".join([str(allele_tsv_record.get(c, "")) for c in ALLELE_TSV_OUTPUT_COLUMNS]) + "\n")


def process_STR_loci_with_multiple_indels(file_path, locus_ids_with_multiple_indels, filtered_out_indels_vcf_writer):
    """Go back and filter out STR loci that turned out to overlap more than one indel since their genotype is unclear.

    Args:
        file_path (str): Path to VCF or TSV file to rewrite in order to remove loci that have this multiple-indels
        locus_ids_with_multiple_indels (set): Set of locus IDs that, in previous steps, were found to overlap multiple indels
        filtered_out_indels_vcf_writer (VcfWriter): VCF writer for filtered out indels
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
                summary_string_fields = fields[3].split(":")
                if len(summary_string_fields) < 2:
                    raise ValueError(f"Unexpected summary string format in {file_path} line: {fields}")
                motif = summary_string_fields[1]
                locus_id = f"{chrom}-{fields[1]}-{fields[2]}-{motif}"

            if locus_id in locus_ids_with_multiple_indels:
                filtered_count += 1
                if file_path.endswith(".vcf") and filtered_out_indels_vcf_writer:
                    # update the FILTER field
                    n_alleles = len(fields[4].split(","))
                    fields[6] = ";".join([FILTER_STR_LOCUS_WITH_MULTIPLE_STR_VARIANTS]*n_alleles)
                    filtered_out_indels_vcf_writer.write("\t".join(fields) + "\n")
            else:
                fo.write("\t".join(fields) + "\n")

    os.rename(f"{file_path}.temp", file_path)
    print(f"Discarded {filtered_count} out of {total} rows from {file_path} as loci that contained more than one STR variant")


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
                if key_prefix.startswith("STR"):
                    total_key = "STR variant counts: TOTAL" if "variant" in key_prefix else "STR allele counts: TOTAL"
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


def main():
    args = parse_args()

    if not args.output_prefix:
        args.output_prefix = re.sub(".vcf(.gz)?", "", os.path.basename(args.input_vcf_path)) + ".STRs"

    counters = collections.defaultdict(int)
    indels_per_locus_counter = collections.defaultdict(int)

    # open the reference genome fasta
    fasta_obj = pyfaidx.Fasta(args.reference_fasta_path, one_based_attributes=False, as_raw=True)

    # open output files for writing
    fopen = gzip.open if args.input_vcf_path.endswith("gz") else open
    vcf_writer = open(f"{args.output_prefix}.vcf", "wt")
    filtered_out_indels_vcf_writer = open(f"{args.output_prefix}.filtered_out_indels.vcf", "wt") if args.write_vcf_with_filtered_out_variants else None

    variants_tsv_writer = open(f"{args.output_prefix}.variants.tsv", "wt")
    variants_tsv_writer.write("\t".join(VARIANT_TSV_OUTPUT_COLUMNS) + "\n")
    alleles_tsv_writer = open(f"{args.output_prefix}.alleles.tsv", "wt")
    alleles_tsv_writer.write("\t".join(ALLELE_TSV_OUTPUT_COLUMNS) + "\n")
    bed_writer = open(f"{args.output_prefix}.variants.bed", "wt") if args.write_bed_file else None

    # open the input VCF
    if args.interval:
        tabix_file = pysam.TabixFile(args.input_vcf_path)
        print(f"Fetching interval(s):", ", ".join(args.interval))
        vcf_iterator = (line for interval in args.interval for line in tabix_file.fetch(interval))
    else:
        vcf_iterator = fopen(args.input_vcf_path, "rt")

    if args.show_progress_bar:
        vcf_iterator = tqdm.tqdm(vcf_iterator, unit=" rows")

    # iterate over all VCF rows
    vcf_line_i = 0
    for line in vcf_iterator:
        if line.startswith("##"):
            # copy the VCF header to the output
            vcf_writer.write(line)
            if args.write_vcf_with_filtered_out_variants:
                filtered_out_indels_vcf_writer.write(line)

            continue

        vcf_fields = line.strip().split("\t")

        if line.startswith("#"):
            # process the last line of the VCF header
            sample_ids = vcf_fields[9:]
            if len(sample_ids) != 1:
                raise ValueError(f"{args.input_vcf_path} contains {len(sample_ids)} samples, but this script only "
                                 f"supports single-sample VCFs.")

            vcf_writer.write(line)
            if args.write_vcf_with_filtered_out_variants:
                filtered_out_indels_vcf_writer.write(line)
            continue

        if args.n is not None and vcf_line_i >= args.n:
            break

        vcf_line_i += 1

        filter_string = process_vcf_line(vcf_line_i, vcf_fields, fasta_obj, vcf_writer, variants_tsv_writer,
                                         alleles_tsv_writer, bed_writer, args, counters, indels_per_locus_counter)

        if filter_string and args.write_vcf_with_filtered_out_variants:
            # if this variant was filtered out, record it in output VCFs

            vcf_fields[6] = filter_string

            if filter_string not in (
                FILTER_MORE_THAN_TWO_ALT_ALLELES,
                FILTER_UNEXPECTED_GENOTYPE_FORMAT,
                FILTER_ALLELE_SNV_OR_MNV,
                ";".join([FILTER_ALLELE_SNV_OR_MNV]*2),
            ):
                filtered_out_indels_vcf_writer.write("\t".join(vcf_fields) + "\n")

    locus_ids_with_multiple_indels = set()
    for locus_id, count in [item for item in indels_per_locus_counter.items() if item[1] > 1]:
        if args.verbose:
            print(f"WARNING: {locus_id} locus contained {count} different STR INDEL variants. Skipping...")
        counters[f"variant filter: {FILTER_STR_LOCUS_WITH_MULTIPLE_STR_VARIANTS}"] += 1
        locus_ids_with_multiple_indels.add(locus_id)

    # close and post-process all output files
    for writer, discard_loci_with_multiple_indels in [
        (vcf_writer, True),
        (variants_tsv_writer, True),
        (alleles_tsv_writer, True),
        (bed_writer, True),
        (filtered_out_indels_vcf_writer, False),
    ]:
        if writer is None:
            continue
        writer.close()
        if len(locus_ids_with_multiple_indels) > 0:
            print(f"Filtering {len(locus_ids_with_multiple_indels)} loci from {writer.name} because they overlap "
                  f"more than one STR INDEL variant")
            if discard_loci_with_multiple_indels:
                process_STR_loci_with_multiple_indels(
                    writer.name, locus_ids_with_multiple_indels, filtered_out_indels_vcf_writer)

        if writer.name.endswith(".tsv"):
            run(f"bgzip -f {writer.name}")
        elif writer.name.endswith(".vcf"):
            if not discard_loci_with_multiple_indels:
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

    print_stats(counters)


if __name__ == "__main__":
    main()

