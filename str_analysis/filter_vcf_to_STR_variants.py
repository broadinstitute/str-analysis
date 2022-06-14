"""
This script takes a .vcf and filters it to variants where either the REF or ALT allele contains a tandem repeat.
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

from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.find_repeat_unit import find_repeat_unit_using_trf, find_repeat_unit

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')

logger = logging.getLogger()

MIN_SEQUENCE_LENGTH_FOR_TRYING_REVERSE_SEQUENCE = 9
MIN_SEQUENCE_LENGTH_FOR_RUNNING_TRF = 12

OUTPUT_TSV_COLUMNS = [
    "Chrom",
    "Start1Based",
    "End1Based",
    "RepeatUnit",
    "Locus",
    "LocusId",
    "InsDel",
    "HETvsHOM",
    "IsFoundInReference",
    "Motif",
    "MotifSize",
    "NumRepeatsShortAllele",
    "NumRepeatsLongAllele",
    "RepeatSizeShortAllele (bp)",
    "RepeatSizeLongAllele (bp)",
    "VcfPos",
    "VcfRef",
    "VcfAlt",
    "VcfGenotype",
    "SummaryString",
    "FoundBy",
    "IsPerfectRepeat",
    "IsMultiAllelic",
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
    p.add_argument("-n", type=int, help="Only process the first N rows of the VCF. Useful for testing.")
    p.add_argument("-o", "--output-prefix", help="Output vcf prefix. If not specified, it will be computed based on "
                                                 "the input vcf filename")
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


def check_if_variant_is_str(
    fasta_obj, chrom, pos, ref, alt,
    min_str_repeats, min_str_length,
    min_fraction_of_variant_covered_by_repeat,
    counters,
    use_trf=False,
    trf_path=None,
):
    """Determine if the given chrom/pos/ref/alt variant represents an STR expansion or contraction or neither.

    Args:

    """
    null_result = {
        "RepeatUnit": None,
    }

    counters["TOTAL"] += 1
    if len(ref) == len(alt):
        counters["skipped: SNV" if len(ref) == 1 else "skipped: MNV"] += 1
        return null_result

    ins_or_del = "INS" if len(ref) < len(alt) else "DEL"
    counters[ins_or_del] += 1

    variant_bases = alt[len(ref):] if len(ref) < len(alt) else ref[len(alt):]
    if len(variant_bases) == 1:
        counters["skipped: 1bp INDEL"] += 1
        return null_result

    found_by = None
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
        counters["skipped: homopolymer"] += 1
        return null_result

    left_flanking_reference_sequence, variant_bases2, right_flanking_reference_sequence = get_flanking_reference_sequences(
        fasta_obj, chrom, pos, ref, alt,
        num_flanking_bases=len(repeat_unit) * 1000)

    if variant_bases2 != variant_bases:
        raise ValueError(f"variant_bases: {variant_bases} != {variant_bases2} variant bases returned by "
                         f"get_flanking_reference_sequences for variant {chrom}-{pos}-{ref}-{alt} ")

    #num_repeats_left_flank, is_perfect_repeat_in_left_flank = extend_repeat_into_sequence(
    #   repeat_unit[::-1], left_flanking_reference_sequence[::-1], min_fraction_covered_by_repeat=min_fraction_of_variant_covered_by_repeat)
    #num_repeats_right_flank, is_perfect_repeat_in_right_flank = extend_repeat_into_sequence(
    #   repeat_unit, right_flanking_reference_sequence, min_fraction_covered_by_repeat=min_fraction_of_variant_covered_by_repeat)

    # check the left flanking sequence for additional repeats matching the repeat unit found in the variant
    i = len(left_flanking_reference_sequence)
    num_repeats_left_flank = 0
    while left_flanking_reference_sequence[i-len(repeat_unit):i] == repeat_unit:
        i -= len(repeat_unit)
        num_repeats_left_flank += 1

    # check the right flanking sequence for additional repeats matching the repeat unit found in the variant
    i = 0
    num_repeats_right_flank = 0
    while right_flanking_reference_sequence[i:i+len(repeat_unit)] == repeat_unit:
        i += len(repeat_unit)
        num_repeats_right_flank += 1

    # even though the VCF position is 1-based, it represents the location of the base preceding the variant bases, so
    # add 1 to get the 1-based position of the true base
    start_1based = (pos + 1) - num_repeats_left_flank * len(repeat_unit)
    end_1based = (pos + 1) + (len(variant_bases) if len(ref) > len(alt) else 0) + num_repeats_right_flank * len(repeat_unit)

    is_perfect_repeat = variant_bases == (repeat_unit * num_repeats_within_variant_bases)
    if not is_perfect_repeat:
        counters[f"STR imperfect repeats"] += 1
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
        counters["skipped: STR too short"] += 1
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
        "FoundBy": found_by,
        "IsPerfectRepeat": is_perfect_repeat,
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

    counters[f"STR TOTAL"] += 1
    counters[f"STR {ins_or_del}"] += 1
    if found_by:
        counters[f"STR found by {found_by}"] += 1
    counters[f"STR delta {num_repeats_within_variant_bases if num_repeats_within_variant_bases < 9 else '9+'} repeats"] += 1
    counters[f"STR motif size {len(repeat_unit) if len(repeat_unit) < 9 else '9+'} bp"] += 1
    counters[f"STR size {num_base_pairs_within_variant_bases}"] += 1
    counters[f"STR with {left_or_right} matching ref. repeat"] += 1
    counters[f"STR perfect repeats"] += 1 if is_perfect_repeat else 0

    return result


def compute_summary_string(alt_STR_alleles):
    summary_string = "RU{}:{}:".format(len(alt_STR_alleles[0]["RepeatUnit"]), alt_STR_alleles[0]["RepeatUnit"])

    ins_or_del = []
    for i in range(0, len(alt_STR_alleles)):
        if len(alt_STR_alleles[i]["Ref"]) > len(alt_STR_alleles[i]["Alt"]):
            ins_or_del.append("DEL")
        elif len(alt_STR_alleles[i]["Ref"]) < len(alt_STR_alleles[i]["Alt"]):
            ins_or_del.append("INS")

    summary_string += ",".join(ins_or_del) + ":"
    summary_string += str(alt_STR_alleles[i]["NumRepeatsRef"]) + "=>" + ",".join([
        str(alt_STR_alleles[i]["NumRepeatsAlt"]) for i in range(0, len(alt_STR_alleles))
    ])

    return ":".join(ins_or_del), summary_string


def main():
    args = parse_args()

    if not args.output_prefix:
        args.output_prefix = re.sub(".vcf(.gz)?", "", os.path.basename(args.input_vcf_path)) + ".STRs"

    vcf_reader = vcf.VCFReader(filename=args.input_vcf_path, encoding="UTF-8")
    vcf_writer = vcf.VCFWriter(open(f"{args.output_prefix}.vcf", "w"), vcf_reader)

    tsv_writer = open(f"{args.output_prefix}.tsv", "wt")
    tsv_writer.write("\t".join(OUTPUT_TSV_COLUMNS) + "\n")

    counters = collections.defaultdict(int)
    fasta_obj = pyfaidx.Fasta(args.reference_fasta_path, one_based_attributes=False, as_raw=True)
    if args.interval:
        vcf_iterator = []
        for interval in args.interval:
            vcf_iterator.extend(vcf_reader.fetch(*parse_interval(interval)))
    else:
        vcf_iterator = vcf_reader

    if len(vcf_reader.samples) != 1:
        raise ValueError(f"{args.input_vcf_path} contains {len(vcf_reader.samples)} samples. "
                         f"This script only supports single-sample VCFs.")

    for row_i, row in tqdm.tqdm(enumerate(vcf_iterator), unit=" rows"):

        if args.n is not None and row_i >= args.n:
            break

        if row.ALT is None or row.ALT[0] is None:
            print(f"No ALT allele found in line {pformat(row.to_dict())}")
            counters[f"ERROR: no alt allele in VCF row"] += 1
            continue

        if len(row.ALT) > 2:
            print(f"WARNING: Multi-allelic variant with {len(row.ALT)} alt alleles found. This script supports no more than 2. Skipping...")
            counters[f"multi-allelic with {len(row.ALT)} alt alleles"] += 1
            continue

        chrom = row.CHROM
        ref = str(row.REF)
        alt_alleles = list(map(str, row.ALT))

        vcf_genotype = row.genotype(row.samples[0].sample)
        vcf_genotype = vcf_genotype.data[0]
        vcf_genotype_fields = re.split(r"[/|\\]", vcf_genotype)
        if len(vcf_genotype_fields) != 2 or any(not s.isdigit() for s in vcf_genotype_fields):
            raise ValueError(f"Unexpected genotype format in row #{row_i + 1}: {vcf_genotype}")

        # if this variant has a * allele (which represents an overlapping deletion) discard it and recode the genotype to homozygous
        if len(alt_alleles) == 2 and alt_alleles[1] == "*":
            alt_alleles = [alt_alleles[0]]
            vcf_genotype_fields[1] = vcf_genotype_fields[0]
            vcf_genotype_separator = "|" if "|" in vcf_genotype else ("\\" if "\\" in vcf_genotype else "/")
            vcf_genotype = vcf_genotype_separator.join(vcf_genotype_fields)

        is_homozygous = vcf_genotype_fields[0] == vcf_genotype_fields[1]
        het_or_hom = "HOM" if is_homozygous else "HET"

        alt_STR_alleles = []
        found_non_STR_allele = False
        for current_alt in alt_alleles:
            current_alt = str(current_alt)
            if set(current_alt) - {"A", "C", "G", "T", "N"}:
                raise ValueError(f"Invalid bases in ALT allele: {current_alt}")

            str_spec = check_if_variant_is_str(
                fasta_obj, chrom, row.POS, ref, current_alt,
                min_str_repeats=args.min_str_repeats,
                min_str_length=args.min_str_length,
                min_fraction_of_variant_covered_by_repeat=args.min_fraction_of_variant_covered_by_repeat,
                counters=counters,
                use_trf=args.use_trf,
                trf_path=args.trf_path,
            )

            if str_spec["RepeatUnit"] is not None:
                alt_STR_alleles.append(str_spec)
            else:
                found_non_STR_allele = True

        if found_non_STR_allele:
            if len(alt_STR_alleles) == 2 and len(alt_STR_alleles) > 0:
                counters[f"skipped: multi-allelic with 1 STR allele and 1 non-STR allele"] += 1
            continue

        # check fields that should be the same for both STR alleles in a multi-allelic variant
        if len(alt_STR_alleles) == 2:
            skip_this_variant = False
            for field in "RepeatUnit", "Start1Based", "End1Based":
                if alt_STR_alleles[0][field] != alt_STR_alleles[1][field]:
                    #print(f"WARNING: Multi-allelic variant STR has alleles with different {field} values:",
                    #      alt_STR_alleles[0][field],  alt_STR_alleles[1][field], ". Skipping...")
                    counters[f"skipped: multi-allelic with different {field}"] += 1
                    skip_this_variant = True
            if skip_this_variant:
                continue

        repeat_unit = alt_STR_alleles[0]["RepeatUnit"]
        locus_start_1based = alt_STR_alleles[0]["Start1Based"]
        locus_end_1based = alt_STR_alleles[0]["End1Based"]

        if len(alt_STR_alleles) == 2:
            if is_homozygous:
                raise ValueError(f"Multi-allelic variant is homozygous in row #{row_i + 1}: {vcf_genotype} ({ref} {alt_alleles})")

            short_allele_size = min(alt_STR_alleles[0]["NumRepeatsAlt"], alt_STR_alleles[1]["NumRepeatsAlt"])
            long_allele_size  = max(alt_STR_alleles[0]["NumRepeatsAlt"], alt_STR_alleles[1]["NumRepeatsAlt"])
        else:
            if is_homozygous:
                short_allele_size = long_allele_size = alt_STR_alleles[0]["NumRepeatsAlt"]
            else:
                short_allele_size = min(alt_STR_alleles[0]["NumRepeatsRef"], alt_STR_alleles[0]["NumRepeatsAlt"])
                long_allele_size  = max(alt_STR_alleles[0]["NumRepeatsRef"], alt_STR_alleles[0]["NumRepeatsAlt"])

        # compute misc. other summary stats
        if any(len(alt_STR_alleles[i]["Ref"]) > len(alt_STR_alleles[i]["Alt"]) for i in range(0, len(alt_STR_alleles))):
            # if either allele is a deletion, it automatically means some of the repeats are found in the reference
            is_found_in_reference = True
        else:
            # check if either allele had STRs in the left or right flanking sequence
            is_found_in_reference = sum([
                alt_STR_alleles[i][f"NumRepeats{s}Flank"]
                for i in range(0, len(alt_STR_alleles))
                for s in ("Left", "Right")
            ]) > 0

        ins_or_del, summary_string = compute_summary_string(alt_STR_alleles)
        found_by = ",".join(sorted(set(map(str, [alt_STR_alleles[i]["FoundBy"] for i in range(0, len(alt_STR_alleles))]))))
        is_perfect_repeat = all(alt_STR_alleles[i]["IsPerfectRepeat"] for i in range(0, len(alt_STR_alleles)))
        is_multiallelic = len(alt_STR_alleles) > 1

        if is_multiallelic:
            counters[f"STR multiallelic"] += 1

        # add info about the repeat to both the ID and the INFO field for convenience & IGV
        new_INFO = dict(row.INFO)
        new_INFO.update({"RU": repeat_unit, "NumRepeats1": short_allele_size, "NumRepeats2": long_allele_size})
        record = vcf.model._Record(
            row.CHROM,
            row.POS,
            summary_string,
            row.REF,
            row.ALT,
            row.QUAL,
            row.FILTER,
            new_INFO,
            row.FORMAT,
            sample_indexes={c.sample: c for c in row.samples},
            samples=row.samples)

        vcf_writer.write_record(record)

        tsv_record = dict(new_INFO)
        #short_allele_size, long_allele_size = compute_short_and_long_allele_size(
        #    locus_start_1based, locus_end_1based, ref, alt, repeat_unit_length=len(repeat_unit), genotype_num_alt=num_alt,
        #)

        tsv_record = {
            "Chrom": chrom,
            "Start1Based": locus_start_1based,
            "End1Based": locus_end_1based,
            "RepeatUnit": repeat_unit,
            "Locus": f"{chrom}:{locus_start_1based}-{locus_end_1based}",
            "LocusId": f"{chrom}:{locus_start_1based - 1}-{locus_end_1based}-{repeat_unit}",
            "InsDel": ins_or_del,
            "HETvsHOM": het_or_hom,
            "IsFoundInReference": "Reference" if is_found_in_reference else "DeNovo",
            "Motif": repeat_unit,
            "MotifSize": len(repeat_unit),
            "NumRepeatsShortAllele": short_allele_size,
            "NumRepeatsLongAllele": long_allele_size,
            "RepeatSizeShortAllele (bp)": short_allele_size * len(repeat_unit),
            "RepeatSizeLongAllele (bp)": long_allele_size * len(repeat_unit),
            "VcfPos": row.POS,
            "VcfRef": ref,
            "VcfAlt": ",".join(alt_alleles),
            "VcfGenotype": vcf_genotype,
            "SummaryString": summary_string,
            "FoundBy": found_by,
            "IsPerfectRepeat": is_perfect_repeat,
            "IsMultiAllelic": is_multiallelic,
        }

        if short_allele_size < 0 or long_allele_size < 0:
            raise ValueError(f"Short or long allele size is < 0: {short_allele_size}, {long_allele_size}  {pformat(tsv_record)}")

        tsv_writer.write("\t".join([str(tsv_record[c]) for c in OUTPUT_TSV_COLUMNS]) + "\n")

    counter_sort_order = lambda x: (x[0].startswith("STR"), len(x[0]), x[0])
    str_keys = False
    for key, value in sorted(counters.items(), key=counter_sort_order):
        if key.startswith("STR") and not str_keys:
            str_keys = True
            logger.info("--------------")
        try:
            percent = 100*value
            percent /= counters["STR TOTAL"] if key.startswith("STR") else counters["TOTAL"]
            percent = f"{percent:5.1f}%"
        except:
            percent = ""

        logger.info(f"{value:10d}  ({percent}) {key}")

    vcf_writer.close()
    tsv_writer.close()


if __name__ == "__main__":
    main()

