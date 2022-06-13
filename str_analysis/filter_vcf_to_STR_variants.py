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


MIN_SEQUENCE_LENGTH_FOR_RUNNING_TRF = 12


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
    null_result=(None, None, None, None, 0, 0, 0, 0, None, None)
):
    """Determine if the given chrom/pos/ref/alt variant represents an STR expansion or contraction or neither.

    Args:

    """

    counters["TOTAL"] += 1
    if len(ref) == len(alt):
        counters["skipped: SNV" if len(ref) == 1 else "skipped: MNV"] += 1
        return null_result

    ins_or_del = "INS" if len(ref) < len(alt) else "DEL"
    counters[ins_or_del] += 1

    variant_bases = alt[1:] if len(ref) < len(alt) else ref[1:]
    if len(variant_bases) == 1:
        counters["skipped: 1bp INDEL"] += 1
        return null_result

    repeat_unit, num_repeats_within_variant_bases = find_repeat_unit(
        variant_bases,
        min_fraction_covered_by_repeat=min_fraction_of_variant_covered_by_repeat,
        allow_indel_interruptions=False,
    )

    found_by_TRF = False
    if use_trf and num_repeats_within_variant_bases <= 1 and len(variant_bases) >= MIN_SEQUENCE_LENGTH_FOR_RUNNING_TRF:
        repeat_unit, num_repeats_within_variant_bases = find_repeat_unit_using_trf(
            variant_bases,
            min_fraction_covered_by_repeat=min_fraction_of_variant_covered_by_repeat,
            trf_path=trf_path,
        )
        found_by_TRF = True

    if len(repeat_unit) == 1:
        counters["skipped: homopolymer"] += 1
        return null_result

    left_flanking_reference_sequence, variant_bases2, right_flanking_reference_sequence = get_flanking_reference_sequences(
        fasta_obj, chrom, pos, ref, alt,
        num_flanking_bases=len(repeat_unit) * 1000)

    assert variant_bases2 == variant_bases

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

    start_0based = pos - num_repeats_left_flank * len(repeat_unit)
    start_1based = start_0based + 1
    end_1based = pos + (len(variant_bases) if len(ref) > len(alt) else 0) + num_repeats_right_flank * len(repeat_unit)

    is_perfect_repeat = variant_bases == repeat_unit * num_repeats_within_variant_bases
    if not is_perfect_repeat:
        counters[f"STR imperfect repeats"] += 1
        #print("TRF:" if found_by_TRF else "---", variant_bases, "is",
        #       f"{num_repeats_within_variant_bases} x {repeat_unit} ... "
        #       f"flanks: [{num_repeats_left_flank}, {num_repeats_right_flank}]")
        if num_repeats_left_flank == 0 and num_repeats_right_flank == 0:
            return null_result

    num_repeats_ref = num_repeats_alt = num_repeats_right_flank + num_repeats_left_flank
    if len(ref) < len(alt):
        # insertion
        num_repeats_alt += num_repeats_within_variant_bases
    else:
        # deletion
        num_repeats_ref += num_repeats_within_variant_bases


    #pprint({
    #    "ref": ref,
    #    "alt": alt,
    #    "repeat_unit": repeat_unit,
    #    "left_flank": left_flanking_reference_sequence[-1*abs(num_repeats_left_flank)*len(repeat_unit):] if num_repeats_left_flank else "",
    #    "right_flank": right_flanking_reference_sequence[:num_repeats_right_flank*len(repeat_unit)],
    #    "num_repeats_within_variant_bases": num_repeats_within_variant_bases,
    #    "num_repeats_ref": num_repeats_ref,
    #    "num_repeats_alt": num_repeats_alt,
    #    "num_repeats_left_flank": num_repeats_left_flank,
    #    "num_repeats_right_flank": num_repeats_right_flank,
    #    "found_by_TRF": found_by_TRF,
    #    "is_perfect_repeat": is_perfect_repeat,
    #})

    num_repeats_in_str = max(num_repeats_ref, num_repeats_alt)
    if num_repeats_in_str < min_str_repeats or num_repeats_in_str * len(repeat_unit) < min_str_length:
        counters["skipped: STR too short"] += 1
        return null_result

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
    counters[f"STR found by TRF"] += 1 if found_by_TRF else 0
    counters[f"STR delta {num_repeats_within_variant_bases if num_repeats_within_variant_bases < 9 else '9+'} repeats"] += 1
    counters[f"STR motif size {len(repeat_unit) if len(repeat_unit) < 9 else '9+'} bp"] += 1
    counters[f"STR size {num_base_pairs_within_variant_bases}"] += 1
    counters[f"STR with {left_or_right} matching ref. repeat"] += 1
    counters[f"STR perfect repeats"] += 1 if is_perfect_repeat else 0

    return chrom, start_1based, end_1based, repeat_unit, num_repeats_ref, num_repeats_alt, \
           num_repeats_left_flank, num_repeats_right_flank, found_by_TRF, is_perfect_repeat


def parse_num_alt_from_genotype(genotype):
    GT = [int(gt) for gt in re.split(r"[/|\\]", genotype) if gt != "."]
    if len(GT) > 0 and max(GT) > 1:
        raise ValueError(f"Found multi-allelic genotype with GT elements > 1: {genotype}")

    num_alt = sum(GT)

    return num_alt


def compute_short_and_long_allele_size(start_1based, end_1based, ref, alt, repeat_unit_length, genotype_num_alt):
    start_0based = start_1based - 1
    locus_size_in_base_pairs = (end_1based - start_0based)
    if locus_size_in_base_pairs % repeat_unit_length != 0:
        raise ValueError(f"Unexpected locus size (bp): {locus_size_in_base_pairs} for repeat unit length: {repeat_unit_length}")

    if len(ref) < len(alt):
        # insertion
        alt_allele_size_in_base_pairs = (end_1based - start_0based) + (len(alt) - 1)
        long_allele_size = int(alt_allele_size_in_base_pairs/repeat_unit_length)
        if genotype_num_alt == 2:
            short_allele_size = long_allele_size
        else:
            short_allele_size = int(locus_size_in_base_pairs/repeat_unit_length)

    else:
        # deletion
        alt_allele_size_in_base_pairs = (end_1based - start_0based) - (len(ref) - 1)
        short_allele_size = int(alt_allele_size_in_base_pairs/repeat_unit_length)
        if genotype_num_alt == 2:
            long_allele_size = short_allele_size
        else:
            long_allele_size = int(locus_size_in_base_pairs/repeat_unit_length)

    if alt_allele_size_in_base_pairs % repeat_unit_length != 0:
        raise ValueError(f"Unexpected alt allele size (bp): {alt_allele_size_in_base_pairs} for repeat unit length: {repeat_unit_length}")

    return short_allele_size, long_allele_size


def main():
    args = parse_args()

    if not args.output_prefix:
        args.output_prefix = re.sub(".vcf(.gz)?", "", os.path.basename(args.input_vcf_path)) + ".STRs"

    vcf_reader = vcf.VCFReader(filename=args.input_vcf_path, encoding="UTF-8")
    vcf_writer = vcf.VCFWriter(open(f"{args.output_prefix}.vcf", "w"), vcf_reader)
    tsv_writer = open(f"{args.output_prefix}.tsv", "wt")
    tsv_columns = [
        "Chrom", "Pos", "Start1Based", "End1Based", "Locus",
        "ExpansionContraction", "HETvsHOM", "IsFoundInReference",
        "Motif", "MotifSize",
        "NumRepeatsShortAllele",  "NumRepeatsLongAllele", "NumRepeatsInVariant",
        "RepeatSizeShortAllele (bp)", "RepeatSizeLongAllele (bp)", "RepeatSizeVariantBases (bp)",
        "Genotype", "Ref", "Alt", "SummaryField",
        "NumRepeatsLeftFlank", "NumRepeatsRightFlank",
        "FoundByTRF",
    ]

    tsv_writer.write("\t".join(tsv_columns) + "\n")
    counters = collections.defaultdict(int)
    fasta_obj = pyfaidx.Fasta(args.reference_fasta_path, one_based_attributes=False, as_raw=True)
    if args.interval:
        vcf_iterator = []
        for interval in args.interval:
            vcf_iterator.extend(vcf_reader.fetch(*parse_interval(interval)))
    else:
        vcf_iterator = vcf_reader

    for i, row in tqdm.tqdm(enumerate(vcf_iterator), unit=" rows"):
        if args.n is not None and i >= args.n:
            break

        if len(row.ALT) > 1:
            raise ValueError("Multi-allelic variant found. VCF should be split and normalized first.")

        if row.ALT is None or row.ALT[0] is None:
            print(f"No ALT allele found in line {pformat(row.to_dict())}")
            continue

        chrom, pos, ref, alt = row.CHROM, row.POS, str(row.REF), str(row.ALT[0])

        (
            locus_chrom,
            locus_start_1based,
            locus_end_1based,
            repeat_unit,
            num_repeats_ref,
            num_repeats_alt,
            num_repeats_left_flank,
            num_repeats_right_flank,
            found_by_TRF,
            is_perfect_repeat,
        ) = check_if_variant_is_str(
            fasta_obj, chrom, pos, ref, alt,
            min_str_repeats=args.min_str_repeats,
            min_str_length=args.min_str_length,
            min_fraction_of_variant_covered_by_repeat=args.min_fraction_of_variant_covered_by_repeat,
            counters=counters,
            use_trf=args.use_trf,
            trf_path=args.trf_path,
        )

        if repeat_unit is None:
            continue

        # add info about the repeat to both the ID and the INFO field for convenience & IGV
        id_field = "RU{}:{}:".format(len(repeat_unit), repeat_unit)
        if len(ref) > len(alt):
            id_field += "DEL:"
        elif len(ref) < len(alt):
            id_field += "INS:"
            
        id_field += "{}=>{}".format(num_repeats_ref, num_repeats_alt)

        new_INFO = dict(row.INFO)
        new_INFO.update({
            "RU": repeat_unit,
            "NumRepeatsRef": num_repeats_ref,
            "NumRepeatsAlt": num_repeats_alt,
            "NumRepeatsLeftFlank": num_repeats_left_flank,
            "NumRepeatsRightFlank": num_repeats_right_flank,
        })

        record = vcf.model._Record(
            chrom,
            pos,
            id_field,
            row.REF,
            row.ALT,
            row.QUAL,
            row.FILTER,
            new_INFO,
            row.FORMAT,
            sample_indexes={c.sample: c for c in row.samples},
            samples=row.samples)

        vcf_writer.write_record(record)

        genotype = row.genotype(row.samples[0].sample)
        genotype = genotype.data[0]
        num_alt = parse_num_alt_from_genotype(genotype)
        if num_alt != 1 and num_alt != 2:
            raise ValueError(f"Unexpected num_alt value: {num_alt}")

        tsv_record = dict(new_INFO)
        short_allele_size, long_allele_size = compute_short_and_long_allele_size(
            locus_start_1based, locus_end_1based, ref, alt, repeat_unit_length=len(repeat_unit), genotype_num_alt=num_alt,
        )

        if len(row.samples) > 1:
            raise ValueError(f"The input vcf contains more than 1 sample: {len(row.samples)}")

        tsv_record.update({
            "Chrom": chrom,
            "Pos": pos,
            "Start1Based": locus_start_1based,
            "End1Based": locus_end_1based,
            "Locus": f"{chrom}:{locus_start_1based}-{locus_end_1based}",
            "ExpansionContraction": "Expansion" if len(ref) < len(alt) else "Contraction",
            "HETvsHOM": "HOM" if num_alt == 2 else ("HET" if num_alt == 1 else num_alt),
            "IsFoundInReference": "Reference" if num_repeats_left_flank + num_repeats_right_flank > 0 else "DeNovo",
            "Motif": repeat_unit,
            "MotifSize": len(repeat_unit),
            "NumRepeatsInVariant": abs(len(ref) - len(alt))/len(repeat_unit),
            "NumRepeatsShortAllele": short_allele_size,
            "NumRepeatsLongAllele": long_allele_size,
            "RepeatSizeShortAllele (bp)": short_allele_size * len(repeat_unit),
            "RepeatSizeLongAllele (bp)": long_allele_size * len(repeat_unit),
            "RepeatSizeVariantBases (bp)": abs(len(ref) - len(alt)),
            "Genotype": genotype,
            "Ref": row.REF,
            "Alt": row.ALT[0],
            "SummaryField": id_field,
            "FoundByTRF": found_by_TRF,
            "IsPerfectRepeat": is_perfect_repeat,
        })

        if short_allele_size < 0 or long_allele_size < 0:
            raise ValueError(f"Short or long allele size is < 0: {short_allele_size}, {long_allele_size}  {pformat(tsv_record)}")

        tsv_writer.write("\t".join([str(tsv_record[c]) for c in tsv_columns]) + "\n")

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

