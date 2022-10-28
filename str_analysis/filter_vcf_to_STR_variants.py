"""
This script takes a .vcf and filters it to insertions and deletions  where either the REF or ALT allele
is a short tandem repeat (STR).
"""

import argparse
import collections
import logging
import os
from pprint import pformat
import pyfaidx
import re
import tqdm
import vcf

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.find_repeat_unit import find_repeat_unit_using_trf, find_repeat_unit

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')

logger = logging.getLogger()

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
]


def parse_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument("-R", "--reference-fasta-path", help="Reference genome fasta path.", required=True)
    p.add_argument("--min-fraction-of-variant-covered-by-repeat", type=float, default=0.90, help="Fraction of inserted "
                   "or deleted bases that must be covered by repeats. Setting this to 1 will require perfect repeats.")
    p.add_argument("--min-str-length", type=int, default=9, help="Minimum STR length in base pairs. This threshold "
                   "applies to the total repeat length comprising any repeats in the flanking sequence to the left "
                   "and right of the variant + in the inserted or deleted bases themselves")
    p.add_argument("--min-str-repeats", type=int, default=3, help="Minimum STR size in number of repeats. This "
                   "threshold applies to the total repeat length comprising any repeats in the flanking sequence to "
                   "the left and right of the variant + in the inserted or deleted bases themselves")
    p.add_argument("--use-trf", action="store_true", help="Use TandemRepeatFinder to look for repeats in the variant sequence")
    p.add_argument("--trf-path", help="TandemRepeatFinder executable path", default="~/bin/trf409.macosx")
    p.add_argument("--show-progress-bar", help="Show a progress bar in the terminal when processing variants.", action="store_true")
    p.add_argument("-n", type=int, help="Only process the first N rows of the VCF. Useful for testing.")
    p.add_argument("-o", "--output-prefix", help="Output vcf prefix. If not specified, it will be computed based on "
                                                 "the input vcf filename")
    p.add_argument("--write-bed-file", help="Whether to output a .bed file containing the STR variants. This requires "
                   "bgzip and tabix tools to be available in the shell environment.", action="store_true")
    p.add_argument("-L", "--interval", help="Only process variants in this genomic interval (format: chrN:start-end)", action="append")
    p.add_argument("input_vcf_path")

    return p.parse_args()


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

    left_flanking_reference_sequence = fasta_obj[chrom][left_flank_start_1based : left_flank_end]
    right_flanking_reference_sequence = fasta_obj[chrom][right_flank_start_1based : right_flank_end]

    return left_flanking_reference_sequence, variant_bases, right_flanking_reference_sequence


def extend_repeat_into_sequence(
    repeat_unit,
    sequence,
    min_fraction_covered_by_repeat=1,
    max_number_of_mismatching_repeats_in_a_row=3,
):
    """This method walks along the given sequence from left to right, one repeat unit length at a time
    (defined by the given repeat unit), and checks for the longest stretch of repeats where at least
    min_fraction_covered_by_repeat fraction of the repeats exactly matches the given repeat_unit.

    Args:
        repeat_unit (str): For example "CAG"
        sequence (str): a longer sequence that may contain repeats of the given repeat unit starting at the left end,
            before switching to random other sequence or repeats of a different repeat unit.
            For example: "CAGCAGCAGCTAGTGCAGTGACAGT"
        min_fraction_covered_by_repeat (float): Look for repeats that correspond to at least.
        max_number_of_mismatching_repeats_in_a_row (int): If this many repeats in a row mismatch the repeat unit,
            stop looking further to the right.

    Returns the number of repeats, as well as whether it's a perfect repeat.
    """
    i = 0
    num_repeats = 0
    num_matching_repeats = 0
    num_mismatching_repeats_in_a_row = 0
    largest_number_of_repeats_above_threshold = 0
    is_perfect_repeat = True
    while num_mismatching_repeats_in_a_row < max_number_of_mismatching_repeats_in_a_row and i <= len(sequence) - len(repeat_unit):
        current_repeat = sequence[i:i+len(repeat_unit)]
        num_repeats += 1
        if current_repeat == repeat_unit:
            num_matching_repeats += 1
            if num_matching_repeats / num_repeats >= min_fraction_covered_by_repeat:
                largest_number_of_repeats_above_threshold = num_repeats
                if num_mismatching_repeats_in_a_row > 0:
                    is_perfect_repeat = False
            num_mismatching_repeats_in_a_row = 0
        else:
            num_mismatching_repeats_in_a_row += 1

        i += len(repeat_unit)

    return largest_number_of_repeats_above_threshold, is_perfect_repeat


def check_if_allele_is_str(
    fasta_obj, chrom, pos, ref, alt,
    min_str_repeats, min_str_length,
    min_fraction_of_variant_covered_by_repeat,
    counters,
    use_trf=False,
    trf_path=None,
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
        min_fraction_of_variant_covered_by_repeat (float): To allow for imperfect repeats, not all repeats need to match
            the given repeat unit. This is the minimum fraction of repeats that must exactly match the repeat unit.
        counters (dict): Dictionary of counters to collect summary stats about the number of STR variants found, etc.
        use_trf (bool): Whether to run TandemRepeatFinder on the variant sequence to check for repeats.
        trf_path (str): Path of TandemRepeatFinder executable.

    Return:
        dict: Fields describing the results of checking whether this variant is an STR
    """
    null_result = {
        "RepeatUnit": None,
    }

    if len(ref) == len(alt):
        counters["non-STR allele: SNV" if len(ref) == 1 else "non-STR allele: MNV"] += 1
        return null_result

    ins_or_del = "INS" if len(ref) < len(alt) else "DEL"
    counters[f"allele counts: {ins_or_del} alleles"] += 1

    variant_bases = alt[len(ref):] if len(ref) < len(alt) else ref[len(alt):]
    if len(variant_bases) == 1:
        counters[f"non-STR allele: 1bp {ins_or_del}"] += 1
        return null_result

    found_by = "NoRepeatFoundInVariantBases"
    repeat_unit, num_repeats_within_variant_bases = find_repeat_unit(
        variant_bases,
        min_fraction_covered_by_repeat=min_fraction_of_variant_covered_by_repeat,
        allow_indel_interruptions=False,
    )
    if num_repeats_within_variant_bases > 1:
        found_by = "BruteForce"

    if num_repeats_within_variant_bases <= 1 and min_fraction_of_variant_covered_by_repeat < 1 and len(variant_bases) >= MIN_SEQUENCE_LENGTH_FOR_TRYING_REVERSE_SEQUENCE:
        # see if we can find a repeat by reversing the sequence
        repeat_unit, num_repeats_within_variant_bases = find_repeat_unit(
            variant_bases[::-1],
            min_fraction_covered_by_repeat=min_fraction_of_variant_covered_by_repeat,
            allow_indel_interruptions=False,
        )
        repeat_unit = repeat_unit[::-1]
        if num_repeats_within_variant_bases > 1:
            found_by = "BruteForceReverseSequence"

    if use_trf and num_repeats_within_variant_bases <= 1 and len(variant_bases) >= MIN_SEQUENCE_LENGTH_FOR_RUNNING_TRF:
        # run TRF to find repeats
        repeat_unit, num_repeats_within_variant_bases = find_repeat_unit_using_trf(
            variant_bases,
            min_fraction_covered_by_repeat=min_fraction_of_variant_covered_by_repeat,
            trf_path=trf_path,
        )
        if num_repeats_within_variant_bases > 1:
            found_by = "TRF"

    if len(repeat_unit) == 1:
        counters["non-STR allele: homopolymer"] += 1
        return null_result

    left_flanking_reference_sequence, variant_bases2, right_flanking_reference_sequence = get_flanking_reference_sequences(
        fasta_obj, chrom, pos, ref, alt,
        num_flanking_bases=len(repeat_unit) * 1000)

    if variant_bases2 != variant_bases:
        raise ValueError(f"variant_bases: {variant_bases} != {variant_bases2} variant bases returned by "
                         f"get_flanking_reference_sequences for variant {chrom}-{pos}-{ref}-{alt} ")

    num_repeats_left_flank, is_perfect_repeat_in_left_flank = extend_repeat_into_sequence(
        repeat_unit[::-1],
        left_flanking_reference_sequence[::-1],
        min_fraction_covered_by_repeat=min_fraction_of_variant_covered_by_repeat)
    num_repeats_right_flank, is_perfect_repeat_in_right_flank = extend_repeat_into_sequence(
        repeat_unit,
        right_flanking_reference_sequence,
        min_fraction_covered_by_repeat=min_fraction_of_variant_covered_by_repeat)

    # even though the VCF position is 1-based, it represents the location of the base preceding the variant bases, so
    # add 1 to get the 1-based position of the true base
    start_1based = pos + 1 - num_repeats_left_flank * len(repeat_unit)
    end_1based = pos + (len(variant_bases) if len(ref) > len(alt) else 0) + num_repeats_right_flank * len(repeat_unit)

    is_perfect_repeat = variant_bases == (repeat_unit * num_repeats_within_variant_bases)
    if not is_perfect_repeat:
        counters[f"STR allele purity: interrupted repeats"] += 1
        if num_repeats_left_flank == 0 and num_repeats_right_flank == 0:
            # only accept imperfect repeats in the variant when there's also support for the same repeat unit in the
            # flanking reference sequence
            return null_result

    num_repeats_ref = num_repeats_alt = num_repeats_right_flank + num_repeats_left_flank
    if len(ref) < len(alt):
        # insertion
        num_repeats_alt += num_repeats_within_variant_bases
    else:
        # deletion
        num_repeats_ref += num_repeats_within_variant_bases

    num_repeats_in_str = max(num_repeats_ref, num_repeats_alt)
    if num_repeats_in_str < min_str_repeats or num_repeats_in_str * len(repeat_unit) < min_str_length:
        counters["non-STR allele: allele repeat sequence is too short"] += 1
        return null_result

    result = {
        "Chrom": chrom,
        "Pos": pos,
        "Ref": ref,
        "Alt": alt,
        "Start1Based": start_1based,
        "End1Based": end_1based,
        "RepeatUnit": repeat_unit,
        "NumRepeatsRef": num_repeats_ref,
        "NumRepeatsAlt": num_repeats_alt,
        "NumRepeatsLeftFlank": num_repeats_left_flank,
        "NumRepeatsRightFlank": num_repeats_right_flank,
        "NumRepeatsInVariant": num_repeats_within_variant_bases,
        "IsPureRepeat": "Yes" if is_perfect_repeat else "No",
    }

    # pprint(result)

    # update counters
    if num_repeats_left_flank > 0 and num_repeats_right_flank > 0:
        left_or_right = 'both left and right'
    elif num_repeats_left_flank > 0:
        left_or_right = 'left'
    elif num_repeats_right_flank > 0:
        left_or_right = 'right'
    else:
        left_or_right = 'no'

    if len(variant_bases) < 500:
        num_base_pairs_within_variant_bases = f"{25*int(len(variant_bases)/25)}-{25*(1+int(len(variant_bases)/25))}bp"
    else:
        num_base_pairs_within_variant_bases = "500+bp"

    counters[f"STR allele counts: TOTAL"] += 1
    counters[f"STR allele counts: {ins_or_del}"] += 1
    if found_by:
        counters[f"STR allele source: found by {found_by}"] += 1
    counters[f"STR allele delta: {num_repeats_within_variant_bases if num_repeats_within_variant_bases < 9 else '9+'} repeats"] += 1
    counters[f"STR allele motif size: {len(repeat_unit) if len(repeat_unit) < 9 else '9+'} bp"] += 1
    counters[f"STR allele size: {num_base_pairs_within_variant_bases}"] += 1
    counters[f"STR allele reference repeats: with {left_or_right} matching ref. repeat"] += 1
    counters[f"STR allele purity: perfect repeats"] += 1 if is_perfect_repeat else 0

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


def compute_summary_string(alt_STR_allele_specs):
    """Returns a short easy-to-read string summarizing the STR variant"""

    repeat_unit = alt_STR_allele_specs[0]["RepeatUnit"]
    summary_string = "RU{}:{}:".format(len(repeat_unit), repeat_unit)

    ins_or_del = []
    for i in range(0, len(alt_STR_allele_specs)):
        if len(alt_STR_allele_specs[i]["Ref"]) > len(alt_STR_allele_specs[i]["Alt"]):
            ins_or_del.append("DEL")
        elif len(alt_STR_allele_specs[i]["Ref"]) < len(alt_STR_allele_specs[i]["Alt"]):
            ins_or_del.append("INS")

    summary_string += ",".join(ins_or_del) + ":"
    summary_string += str(alt_STR_allele_specs[i]["NumRepeatsRef"]) + "=>" + ",".join([
        str(alt_STR_allele_specs[i]["NumRepeatsAlt"]) for i in range(0, len(alt_STR_allele_specs))
    ])

    return ":".join(ins_or_del), summary_string


def main():
    args = parse_args()

    if not args.output_prefix:
        args.output_prefix = re.sub(".vcf(.gz)?", "", os.path.basename(args.input_vcf_path)) + ".STRs"

    # create the input VCF reader and open output files for writing
    vcf_reader = vcf.VCFReader(filename=args.input_vcf_path, encoding="UTF-8")
    vcf_writer = vcf.VCFWriter(open(f"{args.output_prefix}.vcf", "w"), vcf_reader)

    variants_bed_records = set()
    variants_tsv_writer = open(f"{args.output_prefix}.variants.tsv", "wt")
    variants_tsv_writer.write("\t".join(VARIANT_TSV_OUTPUT_COLUMNS) + "\n")
    alleles_tsv_writer = open(f"{args.output_prefix}.alleles.tsv", "wt")
    alleles_tsv_writer.write("\t".join(ALLELE_TSV_OUTPUT_COLUMNS) + "\n")

    if args.interval:
        vcf_iterator = []
        for interval in args.interval:
            vcf_iterator.extend(vcf_reader.fetch(*parse_interval(interval)))
    else:
        vcf_iterator = vcf_reader

    if len(vcf_reader.samples) != 1:
        raise ValueError(f"{args.input_vcf_path} contains {len(vcf_reader.samples)} samples. This script only supports single-sample VCFs.")

    if args.show_progress_bar:
        vcf_iterator = tqdm.tqdm(vcf_iterator, unit=" rows")

    # open the reference genome fasta
    fasta_obj = pyfaidx.Fasta(args.reference_fasta_path, one_based_attributes=False, as_raw=True)

    # iterate over all VCF rows
    counters = collections.defaultdict(int)
    for row_i, row in enumerate(vcf_iterator):

        if args.n is not None and row_i >= args.n:
            break

        # parse the ALT allele(s)
        if row.ALT is None or row.ALT[0] is None:
            print(f"No ALT allele found in line {pformat(row.to_dict())}")
            counters[f"ERROR: no alt allele in VCF row"] += 1
            continue

        ref = str(row.REF)
        alt_alleles = [str(a) for a in row.ALT]
        variant_id = f"{row.CHROM}-{row.POS}-{ref}-" + ",".join(alt_alleles)

        if len(alt_alleles) > 2:
            print(f"WARNING: Multi-allelic variant with {len(row.ALT)} alt alleles found: row #{row_i + 1}: {variant_id}. " 
                  "This script doesn't support more than 2 alt alleles. Skipping...")
            counters[f"WARNING: multi-allelic variant with {len(row.ALT)} alt alleles"] += 1
            continue

        # parse the genotype
        vcf_genotype = row.genotype(row.samples[0].sample)
        vcf_genotype = vcf_genotype.data[0]
        vcf_genotype_indices = re.split(r"[/|\\]", vcf_genotype)
        if len(vcf_genotype_indices) != 2 or any(not gi.isdigit() for gi in vcf_genotype_indices):
            raise ValueError(f"Unexpected genotype format in row #{row_i + 1}: {variant_id} {vcf_genotype}")

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
            alt_alleles = [alt_alleles[0]]
            vcf_genotype_indices = [1,1]
            vcf_genotype_separator = "|" if "|" in vcf_genotype else ("\\" if "\\" in vcf_genotype else "/")
            vcf_genotype = vcf_genotype_separator.join(map(str, vcf_genotype_indices))

        # confirm that there are no more "*" alleles
        if "*" in alt_alleles:
            raise ValueError(f"Unexpected '*' allele in row #{row_i + 1}: {variant_id}")

        is_homozygous = vcf_genotype_indices[0] == vcf_genotype_indices[1]
        het_or_hom = "HOM" if is_homozygous else "HET"

        counters["variant counts: TOTAL variants"] += 1
        # check whether the allele(s) are STRs
        alt_STR_allele_specs = []
        for current_alt in alt_alleles:
            if set(current_alt) - {"A", "C", "G", "T", "N"}:
                raise ValueError(f"Invalid bases in ALT allele: {current_alt}  in row #{row_i + 1}: {variant_id}")

            counters["allele counts: TOTAL alleles"] += 1

            str_spec = check_if_allele_is_str(
                fasta_obj, row.CHROM, row.POS, ref, current_alt,
                min_str_repeats=args.min_str_repeats,
                min_str_length=args.min_str_length,
                min_fraction_of_variant_covered_by_repeat=args.min_fraction_of_variant_covered_by_repeat,
                counters=counters,
                use_trf=args.use_trf,
                trf_path=args.trf_path,
            )

            if str_spec["RepeatUnit"] is not None:
                alt_STR_allele_specs.append(str_spec)
            else:
                alt_STR_allele_specs.append(None)

        if all(allele_spec is None for allele_spec in alt_STR_allele_specs):
            counters[f"skipped variant: no repeat units found in variant"] += 1
            continue

        if any(allele_spec is None for allele_spec in alt_STR_allele_specs):
            counters[f"skipped variant: mixed STR/non-STR multi-allelic variant"] += 1
            continue

        if len(alt_STR_allele_specs) > 1 and is_homozygous:
            raise ValueError(f"Multi-allelic variant is homozygous in row #{row_i + 1}: {vcf_genotype} ({ref} {alt_alleles})")

        # check fields that should be the same for both STR alleles in a multi-allelic STR variant
        if len(alt_STR_allele_specs) > 1:
            skip_this_variant = False
            for field in "Chrom", "RepeatUnit":
                if alt_STR_allele_specs[0][field] != alt_STR_allele_specs[1][field]:
                    print(f"WARNING: Multi-allelic STR {variant_ins_or_del} {variant_summary_string} has "
                          f"alleles with different " + re.sub(r'(?<!^)(?=[A-Z])', '_', field).lower() + "s:",
                          alt_STR_allele_specs[0][field], ",", alt_STR_allele_specs[1][field], ". Skipping...")
                    counters[f"skipped variant: multi-allelic with different {field}s"] += 1
                    skip_this_variant = True
                    break

            if skip_this_variant:
                continue

        repeat_unit = alt_STR_allele_specs[0]["RepeatUnit"]
        variant_locus_start_1based = min(spec["Start1Based"] for spec in alt_STR_allele_specs)
        variant_locus_end_1based = max(spec["End1Based"] for spec in alt_STR_allele_specs)

        variant_ins_or_del, variant_summary_string = compute_summary_string(alt_STR_allele_specs)

        try:
            num_repeats_in_allele1 = get_num_repeats_in_allele(alt_STR_allele_specs, vcf_genotype_indices[0])
            num_repeats_in_allele2 = get_num_repeats_in_allele(alt_STR_allele_specs, vcf_genotype_indices[1])
        except Exception as e:
            raise ValueError(f"{e} {ref} {alt_alleles} {vcf_genotype}")

        variant_short_allele_size = min(num_repeats_in_allele1, num_repeats_in_allele2)
        variant_long_allele_size  = max(num_repeats_in_allele1, num_repeats_in_allele2)

        counters["STR variant counts: TOTAL"] += 1

        # add info about the repeat to both the ID and the INFO field for convenience & IGV
        new_INFO = dict(row.INFO)
        new_INFO.update({
            "Motif": repeat_unit,
            "NumRepeats1": num_repeats_in_allele1,
            "NumRepeats2": num_repeats_in_allele2,
        })
        record = vcf.model._Record(
            row.CHROM,
            row.POS,
            variant_summary_string,
            row.REF,
            row.ALT,
            row.QUAL,
            row.FILTER,
            new_INFO,
            row.FORMAT,
            sample_indexes={c.sample: c for c in row.samples},
            samples=row.samples)

        vcf_writer.write_record(record)

        # write results to TSVs
        if len(alt_alleles) > 1:
            counters[f"STR variant counts: multi-allelic"] += 1

        tsv_record = {
            "Chrom": row.CHROM,
            "Motif": repeat_unit,
            "CanonicalMotif": compute_canonical_motif(repeat_unit, include_reverse_complement=True),
            "MotifSize": len(repeat_unit),
            "VcfPos": row.POS,
            "VcfRef": ref,
            "VcfGenotype": vcf_genotype,
            "HET_or_HOM": het_or_hom,
            "IsMultiallelic": "Yes" if len(alt_alleles) > 1 else "No",
            "IsFoundInReference": "Yes" if any(found_in_reference(spec) for spec in alt_STR_allele_specs) else "No",
        }

        if variant_short_allele_size < 0 or variant_long_allele_size < 0:
            raise ValueError(f"Short or long allele size is < 0: "
                             f"{variant_short_allele_size}, {variant_long_allele_size}  {pformat(tsv_record)}")

        variant_tsv_record = dict(tsv_record)
        variant_tsv_record.update({
            "VcfAlt": ",".join(alt_alleles),
            "INS_or_DEL": variant_ins_or_del,
            "SummaryString": variant_summary_string,
            "IsPureRepeat": all(spec["IsPureRepeat"] for spec in alt_STR_allele_specs),
            "NumRepeatsShortAllele": variant_short_allele_size,
            "NumRepeatsLongAllele": variant_long_allele_size,
            "RepeatSizeShortAllele (bp)": variant_short_allele_size * len(repeat_unit),
            "RepeatSizeLongAllele (bp)": variant_long_allele_size * len(repeat_unit),
            "Start1Based": variant_locus_start_1based,
            "End1Based": variant_locus_end_1based,
            "Locus": f"{row.CHROM}:{variant_locus_start_1based}-{variant_locus_end_1based}",
            "LocusId": f"{row.CHROM}-{variant_locus_start_1based - 1}-{variant_locus_end_1based}-{repeat_unit}",
            "NumRepeatsInReference": (variant_locus_end_1based - variant_locus_start_1based + 1)/len(repeat_unit),
        })
        variants_tsv_writer.write("\t".join([str(variant_tsv_record[c]) for c in VARIANT_TSV_OUTPUT_COLUMNS]) + "\n")
        variants_bed_records.add(
            (row.CHROM, variant_locus_start_1based - 1, variant_locus_end_1based, variant_summary_string))

        for alt_allele, alt_STR_allele_spec in zip(alt_alleles, alt_STR_allele_specs):
            allele_tsv_record = dict(tsv_record)
            ins_or_del, summary_string = compute_summary_string([alt_STR_allele_spec])
            allele_locus_start_1based = alt_STR_allele_spec['Start1Based']
            allele_locus_end_1based = alt_STR_allele_spec['End1Based']
            allele_tsv_record.update({
                "VcfAlt": alt_allele,
                "INS_or_DEL": ins_or_del,
                "SummaryString": summary_string,
                "IsPureRepeat": alt_STR_allele_spec["IsPureRepeat"],
                "NumRepeats": alt_STR_allele_spec["NumRepeatsAlt"],
                "RepeatSize (bp)": alt_STR_allele_spec["NumRepeatsAlt"] * len(repeat_unit),
                "Start1Based": allele_locus_start_1based,
                "End1Based": allele_locus_end_1based,
                "Locus": f"{row.CHROM}:{allele_locus_start_1based}-{allele_locus_end_1based}",
                "LocusId": f"{row.CHROM}-{allele_locus_start_1based - 1}-{allele_locus_end_1based}-{repeat_unit}",
                "NumRepeatsInReference": (allele_locus_end_1based - allele_locus_start_1based + 1)/len(repeat_unit),
            })
            alleles_tsv_writer.write("\t".join([str(allele_tsv_record[c]) for c in ALLELE_TSV_OUTPUT_COLUMNS]) + "\n")

    # print stats
    key_prefixes = set()
    for key, _ in counters.items():
        tokens = key.split(":")
        key_prefixes.add(f"{tokens[0]}:")

    for key_prefix in sorted(key_prefixes):
        if key_prefix in ("variant counts:", "allele counts:"):
            continue
        current_counter = [(key, count) for key, count in counters.items() if key.startswith(key_prefix)]
        current_counter = sorted(current_counter, key=lambda x: (-x[1], x[0]))
        logger.info("--------------")
        for key, value in current_counter:
            if key_prefix.startswith("STR"):
                total_key = "STR variant counts: TOTAL" if "variant" in key_prefix else "STR allele counts: TOTAL"
            else:
                total_key = "variant counts: TOTAL variants" if "variant" in key_prefix else "allele counts: TOTAL alleles"

            total = counters[total_key]
            percent = f"{100*value / total:5.1f}%" if total > 0 else ""

            logger.info(f"{value:10,d} out of {total:10,d} ({percent}) {key}")

    vcf_writer.close()
    print("Finished writing to ", variants_tsv_writer.name)
    variants_tsv_writer.close()
    print("Finished writing to ", alleles_tsv_writer.name)
    alleles_tsv_writer.close()

    os.system(f"bgzip -f {variants_tsv_writer.name}")
    os.system(f"bgzip -f {alleles_tsv_writer.name}")

    os.system(f"bgzip -f {args.output_prefix}.vcf")
    os.system(f"tabix -f {args.output_prefix}.vcf.gz")

    if args.write_bed_file:
        output_bed_path = f"{args.output_prefix}.variants.bed"
        with open(output_bed_path, "wt") as f:
            for chrom, start_0based, end_1based, name in sorted(variants_bed_records):
                f.write("\t".join(map(str, [chrom, start_0based, end_1based, name, "."])) + "\n")
        os.system(f"bgzip -f {output_bed_path}")
        os.system(f"tabix -f {output_bed_path}.gz")
        print(f"Wrote {len(variants_bed_records)} to {output_bed_path}.gz")


if __name__ == "__main__":
    main()

