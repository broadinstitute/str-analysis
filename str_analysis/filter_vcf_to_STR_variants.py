"""
This script takes a .vcf and filters it to insertions and deletions  where either the REF or ALT allele
is a short tandem repeat (STR).
"""
import argparse
import collections
import gzip
import os
from pprint import pformat, pprint
import pyfaidx
import re
import tqdm

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.find_repeat_unit import find_repeat_unit, extend_repeat_into_sequence

MIN_SEQUENCE_LENGTH_FOR_TRYING_REVERSE_SEQUENCE = 9
MIN_SEQUENCE_LENGTH_FOR_RUNNING_TRF = 12

COMMON_TSV_OUTPUT_COLUMNS = [
    "Chrom",
    "Start1Based",
    "End1Based",
    "Locus",
    "LocusId",
    "INS_or_DEL",
    "HET_or_HOM",
    "Motif",
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
    "RepeatUnitInterruptionIndex",
]


def parse_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument("-R", "--reference-fasta-path", help="Reference genome fasta path.", required=True)
    p.add_argument("--min-str-length", type=int, default=9, help="Minimum STR length in base pairs. This threshold "
                   "applies to the total repeat length comprising any repeats in the flanking sequence to the left "
                   "and right of the variant + in the inserted or deleted bases themselves")
    p.add_argument("--min-str-repeats", type=int, default=3, help="Minimum STR size in number of repeats. This "
                   "threshold applies to the total repeat length comprising any repeats in the flanking sequence to "
                   "the left and right of the variant + in the inserted or deleted bases themselves")

    p.add_argument("--allow-interruptions", help="Whether to allow interruptions in the repeat sequence.", action="store_true")
    p.add_argument("--discard-loci-with-multiple-indels", help="Sometimes a VCF will have multiple INDEL variants within the same STR locus. "
                   "This occurs more often within interrupted repeats since left alignment of INDELs doesn't work as well there."
                   "If this flag is used, STR loci that contain more than one INDEL variant will be filtered out.", action="store_true")

    p.add_argument("--show-progress-bar", help="Show a progress bar in the terminal when processing variants.", action="store_true")
    p.add_argument("-n", type=int, help="Only process the first N rows of the VCF. Useful for testing.")
    p.add_argument("-o", "--output-prefix", help="Output vcf prefix. If not specified, it will be computed based on "
                                                 "the input vcf filename")
    p.add_argument("--write-bed-file", help="Whether to output a .bed file containing the STR variants. This requires "
                   "bedtools, bgzip and tabix tools to be available in the shell environment.", action="store_true")
    p.add_argument("-L", "--interval", help="Only process variants in this genomic interval (format: chrN:start-end)", action="append")
    p.add_argument("input_vcf_path")

    args = p.parse_args()

    if not args.allow_interruptions:
        # drop some output columns
        for header in VARIANT_TSV_OUTPUT_COLUMNS, ALLELE_TSV_OUTPUT_COLUMNS:
            for column in "NumPureRepeats", "PureRepeatSize (bp)", "FractionPureRepeats", "RepeatUnitInterruptionIndex":
                if column in header:
                    header.remove(column)

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
    if len(ref) < len(alt):
        variant_bases = alt[len(ref):]
    elif len(ref) > len(alt):
        variant_bases = ref[len(alt):]
    else:
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


def check_if_allele_is_str(
    fasta_obj, chrom, pos, ref, alt,
    min_str_repeats, min_str_length,
    counters,
    allow_interruptions=False,
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
        counters (dict): Dictionary of counters to collect summary stats about the number of STR variants found, etc.
        allow_interruptions (bool): Whether to allow interruptions in the repeat sequence

    Return:
        dict: Fields describing the results of checking whether this variant is an STR
    """
    null_result = {
        "RepeatUnit": None,
    }

    counter_key_suffix = " (interruptions allowed)" if allow_interruptions else ""

    if len(ref) == len(alt):
        counters[f"non-STR allele: SNV{counter_key_suffix}" if len(ref) == 1 else "non-STR allele: MNV"] += 1
        return null_result

    ins_or_del = "INS" if len(ref) < len(alt) else "DEL"
    counters[f"allele counts: {ins_or_del} alleles{counter_key_suffix}"] += 1

    variant_bases = alt[len(ref):] if len(ref) < len(alt) else ref[len(alt):]
    if len(variant_bases) == 1:
        counters[f"non-STR allele: 1bp {ins_or_del}{counter_key_suffix}"] += 1
        return null_result

    repeat_unit, num_pure_repeats_within_variant_bases, num_total_repeats_within_variant_bases, repeat_unit_interruption_index = find_repeat_unit(
        variant_bases,
        allow_interruptions=allow_interruptions,
        allow_partial_repeats=False,
    )

    if len(repeat_unit) == 1:
        counters[f"non-STR allele: homopolymer{counter_key_suffix}"] += 1
        return null_result

    left_flanking_reference_sequence, variant_bases2, right_flanking_reference_sequence = get_flanking_reference_sequences(
        fasta_obj, chrom, pos, ref, alt,
        num_flanking_bases=len(repeat_unit) * 1000)

    if variant_bases2 != variant_bases:
        raise ValueError(f"variant_bases: {variant_bases} != {variant_bases2} variant bases returned by "
                         f"get_flanking_reference_sequences for variant {chrom}-{pos}-{ref}-{alt} ")

    reversed_repeat_unit_interruption_index = None
    if repeat_unit_interruption_index is not None:
        reversed_repeat_unit_interruption_index = len(repeat_unit) - 1 - repeat_unit_interruption_index
    num_pure_repeats_left_flank, num_total_repeats_left_flank, reversed_repeat_unit_interruption_index = extend_repeat_into_sequence(
        repeat_unit[::-1],
        left_flanking_reference_sequence[::-1],
        allow_interruptions=allow_interruptions,
        repeat_unit_interruption_index=reversed_repeat_unit_interruption_index)
    if reversed_repeat_unit_interruption_index is not None:
        # reverse the repeat_unit_interruption_index
        repeat_unit_interruption_index = len(repeat_unit) - 1 - reversed_repeat_unit_interruption_index

    num_pure_repeats_right_flank, num_total_repeats_right_flank, repeat_unit_interruption_index = extend_repeat_into_sequence(
        repeat_unit,
        right_flanking_reference_sequence,
        allow_interruptions=allow_interruptions,
        repeat_unit_interruption_index=repeat_unit_interruption_index)

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

    num_pure_repeats_in_str = num_pure_repeats_left_flank + num_pure_repeats_within_variant_bases + num_pure_repeats_right_flank
    num_total_repeats_in_str = max(num_total_repeats_ref, num_total_repeats_alt)

    if not allow_interruptions and num_pure_repeats_in_str != num_total_repeats_in_str:
        # sanity check
        raise Exception("Repeat has interruptions even though allow_interruptions is False")

    if num_total_repeats_in_str < min_str_repeats or num_total_repeats_in_str * len(repeat_unit) < min_str_length:
        counters[f"non-STR allele: allele repeat sequence is too short{counter_key_suffix}"] += 1
        return null_result

    is_pure_repeat = num_pure_repeats_in_str == num_total_repeats_in_str
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
        "IsPureRepeat": is_pure_repeat,
        "FractionPureRepeats": num_pure_repeats_in_str / num_total_repeats_in_str,
        "RepeatUnitInterruptionIndex": repeat_unit_interruption_index,
        "PureStart1Based": pure_start_1based,
        "PureEnd1Based": pure_end_1based,
        "NumPureRepeatsRef": num_pure_repeats_ref,
        "NumPureRepeatsAlt": num_pure_repeats_alt,
        "NumPureRepeatsInStr": num_pure_repeats_in_str,
        "NumPureRepeatsLeftFlank": num_pure_repeats_left_flank,
        "NumPureRepeatsRightFlank": num_pure_repeats_right_flank,
        "NumPureRepeatsInVariant": num_pure_repeats_within_variant_bases,
        "NumRepeatsInStr": num_total_repeats_in_str,
    }

    # update counters
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

    counters[f"STR allele counts{counter_key_suffix}: TOTAL"] += 1
    counters[f"STR allele counts{counter_key_suffix}: {ins_or_del}"] += 1
    counters[f"STR allele delta{counter_key_suffix}: {num_total_repeats_within_variant_bases if num_total_repeats_within_variant_bases < 9 else '9+'} repeats"] += 1
    counters[f"STR allele motif size{counter_key_suffix}: {len(repeat_unit) if len(repeat_unit) < 9 else '9+'} bp"] += 1
    counters[f"STR allele size{counter_key_suffix}: {num_base_pairs_within_variant_bases}"] += 1
    counters[f"STR allele reference repeats{counter_key_suffix}: with {left_or_right} matching ref. repeat"] += 1
    return result


def found_in_reference(alt_STR_allele):
    # if the allele is a deletion, it automatically means some of the repeats are found in the reference
    # otherwise, return True if the allele had STRs in the left or right flanking sequence
    is_deletion = len(alt_STR_allele["Ref"]) > len(alt_STR_allele["Alt"])
    if is_deletion:
        return True

    num_repeats_in_left_or_right_flank = sum([alt_STR_allele[f"NumRepeats{s}Flank"] for s in ("Left", "Right")])
    if num_repeats_in_left_or_right_flank > 0:
        return True

    return False


def get_num_repeats_in_allele(alt_STR_allele_specs, genotype_index):
    """Looks up the number of STR repeats found in a given VCF record's allele(s) based on the genotype.

    Args:
        alt_STR_allele_specs (list): list containing the STR allele spec(s) returned by check_if_allele_is_str for each
            of the allele(s) in this variant.
        genotype_index (int): which genotype to look up.

    Return
        int: The number of repeats in this allele
    """
    if not isinstance(alt_STR_allele_specs, list):
        raise ValueError(f"alt_STR_allele_specs argument is not of type list: {alt_STR_allele_specs}")
    if not isinstance(genotype_index, int):
        raise ValueError(f"genotype_index argument is not of type int: {genotype_index}")

    if genotype_index == 0:
        return alt_STR_allele_specs[0]["NumRepeatsRef"]

    if genotype_index - 1 >= len(alt_STR_allele_specs):
        raise ValueError(f"genotype_index argument is larger than the alt_STR_allele_specs list size: "
                         f"{genotype_index} vs {alt_STR_allele_specs}")

    return alt_STR_allele_specs[genotype_index - 1]["NumRepeatsAlt"]


def compute_variant_summary_string(alt_STR_allele_specs, het_or_hom):
    """Returns a short easy-to-read string summarizing the STR variant.

    Args:
        alt_STR_allele_specs (list): list of 1 or more allele spec dictionaries
        het_or_hom (str): describes the genotype as "HET" or "HOM"
    Return:
        str: short summary of the variant
    """

    repeat_unit = alt_STR_allele_specs[0]["RepeatUnit"]
    summary_string = "{}bp:{}:".format(len(repeat_unit), repeat_unit)

    ins_or_del = []
    for i in range(0, len(alt_STR_allele_specs)):
        if len(alt_STR_allele_specs[i]["Ref"]) > len(alt_STR_allele_specs[i]["Alt"]):
            ins_or_del.append("DEL")
        elif len(alt_STR_allele_specs[i]["Ref"]) < len(alt_STR_allele_specs[i]["Alt"]):
            ins_or_del.append("INS")
        elif len(alt_STR_allele_specs[i]["Ref"]) == len(alt_STR_allele_specs[i]["Alt"]):
            ins_or_del.append("REF")

    summary_string += ",".join(ins_or_del) + ":"
    summary_string += str(alt_STR_allele_specs[i]["NumRepeatsRef"]) + "=>" + ",".join([
        str(alt_STR_allele_specs[i]["NumRepeatsAlt"]) for i in range(0, len(alt_STR_allele_specs))
    ])
    summary_string += ":" + het_or_hom

    if not alt_STR_allele_specs[0]["IsPureRepeat"]:
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
    """Utility method for processing a single line from the input vcf"""

    # parse the ALT allele(s)
    vcf_chrom = vcf_fields[0]
    vcf_chrom_without_chr_prefix = vcf_chrom.replace("chr", "")
    vcf_pos = int(vcf_fields[1])
    vcf_ref = vcf_fields[3]
    vcf_alt = vcf_fields[4]
    if not vcf_alt:
        print(f"No ALT allele found in line {vcf_fields}")
        counters[f"ERROR: no alt allele in VCF row"] += 1
        return

    variant_id = f"{vcf_chrom_without_chr_prefix}-{vcf_pos}-{vcf_ref}-{vcf_alt}"

    alt_alleles = vcf_alt.split(",")
    if len(alt_alleles) > 2:
        print(f"WARNING: Multi-allelic variant with {len(alt_alleles)} alt alleles found: row #{vcf_line_i}: {variant_id}. "
              "This script doesn't support more than 2 alt alleles. Skipping...")
        counters[f"WARNING: multi-allelic variant with {len(alt_alleles)} alt alleles"] += 1
        return

    # parse the genotype
    vcf_genotype_format = vcf_fields[8].split(":")
    if vcf_genotype_format[0] != "GT":
        raise ValueError(f"Unexpected genotype field in row #{vcf_line_i}: {vcf_fields[8]}  (variant: {variant_id})")
    vcf_genotype = next(iter(vcf_fields[9].split(":")))
    vcf_genotype_indices = re.split(r"[/|\\]", vcf_genotype)
    if len(vcf_genotype_indices) != 2 or any(not gi.isdigit() for gi in vcf_genotype_indices):
        raise ValueError(f"Unexpected genotype GT format in row #{vcf_line_i}: {vcf_genotype}  (variant: {variant_id})")

    vcf_genotype_indices = [int(gi) for gi in vcf_genotype_indices]

    # check that genotype indices correctly correspond to alt allele(s)
    for genotype_index in vcf_genotype_indices:
        try:
            genotype_index = int(genotype_index)
        except ValueError:
            raise ValueError(f"Unexpected genotype field has non-numeric genotype indices {vcf_genotype_indices}")
        if genotype_index > len(alt_alleles):
            raise ValueError(f"Unexpected genotype field has a genotype index that is larger than the number of "
                             f"alt alleles {alt_alleles}: {vcf_genotype_indices}")

    # if this variant has 1 regular allele and 1 "*" allele (which represents an overlapping deletion), discard the
    # "*" allele and recode the genotype as homozygous
    if len(alt_alleles) == 2 and alt_alleles[1] == "*":
        #counters["variant counts: removed * allele and converted to homozygous genotype"] += 1
        alt_alleles = [alt_alleles[0]]
        vcf_genotype_indices = [1,1]
        vcf_genotype_separator = "|" if "|" in vcf_genotype else ("\\" if "\\" in vcf_genotype else "/")
        vcf_genotype = vcf_genotype_separator.join(map(str, vcf_genotype_indices))

    # confirm that there are no more "*" alleles
    if "*" in alt_alleles:
        raise ValueError(f"Unexpected '*' allele in row #{vcf_line_i}: {variant_id}")

    is_homozygous = vcf_genotype_indices[0] == vcf_genotype_indices[1]
    het_or_hom = "HOM" if is_homozygous else "HET"

    counters["variant counts: TOTAL variants"] += 1

    for allow_interruptions in ([False, True] if args.allow_interruptions else [False]):
        # check whether the allele(s) are STRs
        alt_STR_allele_specs = []
        for current_alt in alt_alleles:
            if set(current_alt) - set("ACGTN"):
                raise ValueError(f"Invalid bases in ALT allele: {current_alt}  in row #{vcf_line_i}: {variant_id}")

            counters["allele counts: TOTAL alleles"] += 1

            str_spec = check_if_allele_is_str(
                fasta_obj, vcf_chrom, vcf_pos, vcf_ref, current_alt,
                min_str_repeats=args.min_str_repeats,
                min_str_length=args.min_str_length,
                counters=counters,
                allow_interruptions=allow_interruptions,
            )

            if str_spec["RepeatUnit"] is not None:
                alt_STR_allele_specs.append(str_spec)
            else:
                # append None to indicate this is not an STR allele
                alt_STR_allele_specs.append(None)

        if any(allele_spec is not None for allele_spec in alt_STR_allele_specs):
            # found at least one STR allele
            break

    if all(allele_spec is None for allele_spec in alt_STR_allele_specs):
        counters[f"skipped variant: no repeat units found in variant"] += 1
        return

    if len(alt_STR_allele_specs) > 1 and any(allele_spec is None for allele_spec in alt_STR_allele_specs):
        # multiallelics with one non-STR allele (such as an SNV) and one STR allele are difficult to process downstream.
        # There are very few of these, so just discard them.
        counters[f"skipped variant: mixed STR/non-STR multi-allelic variant"] += 1
        return

    if len(alt_STR_allele_specs) > 1 and is_homozygous:
        # since these variants are from a single diploid sample, a homozygous genotype and multiple alleles are
        # an error of some sort.
        raise ValueError(f"Multi-allelic variant is homozygous in row #{vcf_line_i}: {vcf_genotype} ({vcf_ref} {alt_alleles})")

    is_multiallelic = len(alt_STR_allele_specs) > 1
    variant_ins_or_del, variant_summary_string = compute_variant_summary_string(
        alt_STR_allele_specs, "MULTI" if is_multiallelic else het_or_hom)

    # check fields that should be the same for both STR alleles in a multi-allelic STR variant
    if len(alt_STR_allele_specs) > 1:
        if alt_STR_allele_specs[0]["Chrom"] != alt_STR_allele_specs[1]["Chrom"]:
            # sanity check
            raise Exception("Alleles in a multiallelic STR have different chromosomes")
        for attribute in "RepeatUnit", "Start1Based", "End1Based", "IsPureRepeat":
            if alt_STR_allele_specs[0][attribute] != alt_STR_allele_specs[1][attribute]:
                print(f"WARNING: Multi-allelic STR {variant_ins_or_del} {variant_summary_string} has alleles with "
                      f"different {attribute}:",
                      alt_STR_allele_specs[0][attribute], " vs ",
                      alt_STR_allele_specs[1][attribute], "  Skipping...")
                counters[f"skipped variant: multi-allelic with different {attribute}"] += 1
                return

    # get repeat unit, start, end coords
    repeat_unit = alt_STR_allele_specs[0]["RepeatUnit"]
    start_1based = alt_STR_allele_specs[0]["Start1Based"]
    end_1based = alt_STR_allele_specs[0]["End1Based"]
    is_pure_repeat = alt_STR_allele_specs[0]["IsPureRepeat"]

    # get allele sizes
    try:
        num_repeats_in_allele1 = get_num_repeats_in_allele(alt_STR_allele_specs, vcf_genotype_indices[0])
        num_repeats_in_allele2 = get_num_repeats_in_allele(alt_STR_allele_specs, vcf_genotype_indices[1])
    except Exception as e:
        raise ValueError(f"{e} {vcf_ref} {alt_alleles} {vcf_genotype}")

    variant_short_allele_size = min(num_repeats_in_allele1, num_repeats_in_allele2)
    variant_long_allele_size  = max(num_repeats_in_allele1, num_repeats_in_allele2)

    counters["STR variant counts: TOTAL"] += 1

    locus_id = f"{vcf_chrom_without_chr_prefix}-{start_1based - 1}-{end_1based}-{repeat_unit}"
    indels_per_locus_counter[locus_id] += 1

    # output the vcf record
    INFO_field = ";".join([
        f"Motif={repeat_unit}",
        f"NumRepeats1={num_repeats_in_allele1}",
        f"NumRepeats2={num_repeats_in_allele2}",
        f"LocusId={locus_id}",
    ])
    vcf_fields[2] = variant_summary_string
    vcf_fields[4] = ",".join(alt_alleles)
    vcf_fields[7] = INFO_field
    vcf_fields[9] = vcf_genotype
    vcf_writer.write("\t".join(vcf_fields) + "\n")

    # write results to TSVs
    if len(alt_alleles) > 1:
        counters[f"STR variant counts: multi-allelic"] += 1

    tsv_record = {
        "Chrom": vcf_chrom,
        "Start1Based": start_1based,
        "End1Based": end_1based,
        "Locus": f"{vcf_chrom_without_chr_prefix}:{start_1based}-{end_1based}",
        "LocusId": locus_id,
        "NumRepeatsInReference": (end_1based - start_1based + 1)/len(repeat_unit),
        "Motif": repeat_unit,
        "CanonicalMotif": compute_canonical_motif(repeat_unit, include_reverse_complement=True),
        "MotifSize": len(repeat_unit),
        "VcfPos": vcf_pos,
        "VcfRef": vcf_ref,
        "VcfGenotype": vcf_genotype,
        "HET_or_HOM": het_or_hom,
        "IsPureRepeat": is_pure_repeat,
        "IsMultiallelic": len(alt_alleles) > 1,
        "IsFoundInReference": any(found_in_reference(spec) for spec in alt_STR_allele_specs),
    }

    if variant_short_allele_size < 0 or variant_long_allele_size < 0:
        raise ValueError(f"Short or long allele size is < 0: "
                         f"{variant_short_allele_size}, {variant_long_allele_size}  {pformat(tsv_record)}")

    counters[f"STR variant counts: pure repeats"] += 1 if is_pure_repeat else 0

    variant_tsv_record = dict(tsv_record)
    variant_tsv_record.update({
        "VcfAlt": ",".join(alt_alleles),
        "INS_or_DEL": variant_ins_or_del,
        "SummaryString": variant_summary_string,
        "NumRepeatsShortAllele": variant_short_allele_size,
        "NumRepeatsLongAllele": variant_long_allele_size,
        "RepeatSizeShortAllele (bp)": variant_short_allele_size * len(repeat_unit),
        "RepeatSizeLongAllele (bp)": variant_long_allele_size * len(repeat_unit),
    })
    variants_tsv_writer.write("\t".join([str(variant_tsv_record[c]) for c in VARIANT_TSV_OUTPUT_COLUMNS]) + "\n")

    if bed_writer is not None:
        bed_writer.write("\t".join(map(str, [vcf_chrom, start_1based - 1, end_1based, variant_summary_string, "."])) + "\n")

    for alt_allele, alt_STR_allele_spec in zip(alt_alleles, alt_STR_allele_specs):
        allele_tsv_record = dict(tsv_record)
        ins_or_del, summary_string = compute_variant_summary_string([alt_STR_allele_spec], het_or_hom)
        allele_tsv_record.update({
            "VcfAlt": alt_allele,
            "INS_or_DEL": ins_or_del,
            "SummaryString": summary_string,
            "NumRepeats": alt_STR_allele_spec["NumRepeatsAlt"],
            "RepeatSize (bp)": alt_STR_allele_spec["NumRepeatsAlt"] * len(repeat_unit),
            "NumPureRepeats": alt_STR_allele_spec["NumPureRepeatsAlt"],
            "PureRepeatSize (bp)": alt_STR_allele_spec["NumPureRepeatsAlt"] * len(repeat_unit),
            "FractionPureRepeats": "%0.3f" % alt_STR_allele_spec["FractionPureRepeats"],
            "RepeatUnitInterruptionIndex": alt_STR_allele_spec["RepeatUnitInterruptionIndex"],
        })

        alleles_tsv_writer.write("\t".join([str(allele_tsv_record[c]) for c in ALLELE_TSV_OUTPUT_COLUMNS]) + "\n")


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

    variants_tsv_writer = open(f"{args.output_prefix}.variants.tsv", "wt")
    variants_tsv_writer.write("\t".join(VARIANT_TSV_OUTPUT_COLUMNS) + "\n")
    alleles_tsv_writer = open(f"{args.output_prefix}.alleles.tsv", "wt")
    alleles_tsv_writer.write("\t".join(ALLELE_TSV_OUTPUT_COLUMNS) + "\n")
    bed_writer = open(f"{args.output_prefix}.variants.bed", "wt") if args.write_bed_file else None

    # open the input VCF
    vcf_iterator = fopen(args.input_vcf_path, "rt")
    if args.show_progress_bar:
        vcf_iterator = tqdm.tqdm(vcf_iterator, unit=" rows")

    # iterate over all VCF rows
    vcf_line_i = 0
    for line in vcf_iterator:
        if line.startswith("##"):
            # copy the VCF header to the output
            vcf_writer.write(line)
            continue

        vcf_fields = line.strip().split("\t")

        if line.startswith("#"):
            # process the last line of the VCF header
            vcf_writer.write(line)
            sample_ids = vcf_fields[9:]
            if len(sample_ids) != 1:
                raise ValueError(f"{args.input_vcf_path} contains {len(sample_ids)} samples, but this script only supports single-sample VCFs.")
            continue

        if args.n is not None and vcf_line_i >= args.n:
            break
        vcf_line_i += 1

        process_vcf_line(vcf_line_i, vcf_fields, fasta_obj,
            vcf_writer, variants_tsv_writer, alleles_tsv_writer, bed_writer, args, counters, indels_per_locus_counter)

    locus_ids_with_multiple_indels = set()
    for locus_id, count in indels_per_locus_counter.items():
        if count > 1:
            print(f"WARNING: {locus_id} locus contained {count} different INDEL variants in the input vcf")
            locus_ids_with_multiple_indels.add(locus_id)

    for writer in vcf_writer, variants_tsv_writer, alleles_tsv_writer, bed_writer:
        if writer is None:
            continue
        writer.close()
        if args.discard_loci_with_multiple_indels and len(locus_ids_with_multiple_indels) > 0:
            print(f"Filtering {len(locus_ids_with_multiple_indels)} loci from {writer.name} because they have "
                  f"more than one STR variant")
            discard_loci_with_multiple_indels(writer.name, locus_ids_with_multiple_indels)

        print("Finished writing to", writer.name)

    print_stats(counters)

    # final post-processing steps
    os.system(f"bgzip -f {variants_tsv_writer.name}")
    os.system(f"bgzip -f {alleles_tsv_writer.name}")

    os.system(f"bgzip -f {args.output_prefix}.vcf")
    os.system(f"tabix -f {args.output_prefix}.vcf.gz")

    if args.write_bed_file:
        os.system(f"bedtools sort -i {bed_writer.name} | bgzip > {bed_writer.name}.gz")
        os.system(f"tabix -f {bed_writer.name}.gz")


def discard_loci_with_multiple_indels(file_path, locus_ids_with_multiple_indels):
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
                continue

            total += 1
            fields = line.strip().split("\t")
            if file_path.endswith(".tsv"):
                row = dict(zip(header, fields))
                locus_id = row["LocusId"]
            elif file_path.endswith(".vcf"):
                info = dict([info_record.split("=") for info_record in fields[7].split(";")])
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
                continue

            fo.write(line)

    os.rename(f"{file_path}.temp", file_path)

    print(f"Discarded {filtered_count} out of {total} rows from {file_path} as loci that contained more than one STR variant")


def print_stats(counters):
    """Print out all the counters"""

    key_prefixes = set()
    for key, _ in counters.items():
        tokens = key.split(":")
        key_prefixes.add(f"{tokens[0]}:")

    for key_prefix in sorted(key_prefixes):
        if key_prefix in ("variant counts:", "allele counts:"):
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

            print(f"{value:10,d} out of {total:10,d} ({percent}) {key}")


if __name__ == "__main__":
    main()

