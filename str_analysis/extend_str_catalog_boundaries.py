"""Extend tandem repeat locus boundaries outward into their flanking sequence using a gap-purity heuristic.
The EP400 locus served as a motivating example for this (described in https://gnomad.broadinstitute.org/news/2026-07-str-data-update/).


Algorithm pseudo-code, applied to each side of the existing locus definition:

    boundary = the locus's current start (going left) or end (going right)
    loop:
		# 1. scan the flank for additional copies of the repeat motif
        scan outward from boundary in steps of one motif length
            stop the scan once an exact copy of the motif is found  -> anchor
            give up if you see more than MAX_REPEATS_IN_GAP motif copies that don't match the motif
        if no anchor was found:
            stop, this side is done
            
        # 2. extend to include all exact copies of the motif that follow the anchor
        run_end = anchor
        while the copy starting at run_end is exact:
            run_end += motif length
            
        # 3. check if the candidate extension sequence passes a purity threshold
        purity = fraction of bases in the [boundary, run_end) interval that match a perfect repeat sequence of the motif
        if purity >= MIN_PURITY_OF_NEW_SEQUENCE:  # default min purity = 0.9 (matching ExpansionHunter's internal threshold for detecting repetitive reads)
            boundary = run_end          # accept, then try to extend further from here
        else:
            stop, this side is done     # rejection is final for this side

Sequence is compared against the motif tiled in the locus's existing frame, so a repeat in the flank
that is out of phase with the boundary will not match. 

Output:

By default an extended definition replaces original locus definition. Specifying
--keep-original-definitions-of-extended-loci causes the original locus definition to be included 
in the output alongside the extended definition, making the output a superset of the input.
Any loci that didn't get extended are included in the output regardless.
Also, extended definitions are deduplicated in both modes - eg. if neighboring loci in a
catalog get extended to produce identical locus definitions. In that case, the neighboring loci collapse 
to the single extended definition they share. An extended definition is also dropped when 
an identical definition is already in the output as an input locus.


Usage:

    python3 -m str_analysis.extend_str_catalog_boundaries -R hg38.fa catalog.bed.gz

"""

import argparse
import collections
import os
import re

import pysam
from tqdm import tqdm

from str_analysis.utils.eh_catalog_utils import get_variant_catalog_iterator, \
    parse_motifs_from_locus_structure, convert_json_records_to_bed_format_tuples
from str_analysis.utils.fasta_utils import normalize_chrom_using_pysam_fasta
from str_analysis.utils.gap_purity_extension import compute_extension, \
    MAX_REPEATS_IN_GAP, MIN_PURITY_OF_NEW_SEQUENCE
from str_analysis.utils.misc_utils import parse_interval


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extend tandem repeat locus boundaries into their flanking sequence using the "
                    "gap-purity rule.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-R", "--reference-fasta", required=True, help="Reference fasta file path.")
    parser.add_argument("-o", "--output-path", help="Output BED path. Defaults to the input path with "
                        "'.extended' added before the file extension. Always bgzipped and tabix-indexed.")
    parser.add_argument("--max-repeats-in-gap", type=int, default=MAX_REPEATS_IN_GAP,
                        help="How many consecutive interrupted motif copies may separate the current "
                             "boundary from the next exact copy of the motif. Raising this lets the "
                             "rule reach across longer interruptions.")
    parser.add_argument("--min-purity-of-new-sequence", type=float, default=MIN_PURITY_OF_NEW_SEQUENCE,
                        help="The fraction of bases in the newly added sequence, interrupted copies "
                             "included, that must match a perfect tiling of the motif. Lowering this "
                             "accepts extensions that are less clean.")
    parser.add_argument("-k", "--keep-original-definitions-of-extended-loci", action="store_true",
                        help="If specified, the original locus definitions will be included in the "
                             "output alongside the extended locus definitions. If not, any locus "
                             "definitions that are extended will be replaced with the new extended "
                             "definitions. Regardless, if two adjacent loci in the input catalog yield "
                             "the same extended locus definition (same start, end, and motif), the "
                             "output will be deduplicated and only one copy of this extended "
                             "definition will be included.")
    parser.add_argument("-n", "--process-n-loci", type=int, help="Only process the first N loci.")
    parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar.")
    parser.add_argument("--verbose", action="store_true", help="Print every accept/reject decision the "
                        "rule makes, with the purity behind it. Intended for looking at a handful of "
                        "loci, not a whole catalog.")
    parser.add_argument("catalog_bed", help="A catalog of tandem repeats in BED format.")
    args = parser.parse_args()

    if not 0 < args.min_purity_of_new_sequence <= 1:
        parser.error(f"--min-purity-of-new-sequence must be > 0 and <= 1, got "
                     f"{args.min_purity_of_new_sequence}")
    if args.max_repeats_in_gap < 0:
        parser.error(f"--max-repeats-in-gap cannot be negative, got {args.max_repeats_in_gap}")
    if not os.path.isfile(args.reference_fasta):
        parser.error(f"Reference fasta not found: {args.reference_fasta}")
    if ".json" in os.path.basename(args.catalog_bed):
        parser.error(f"This script takes a BED catalog, not JSON: {args.catalog_bed}")

    if not args.output_path:
        args.output_path = compute_output_path(args.catalog_bed)
    else:
        # The output is bgzipped and tabix-indexed, and bgzip always writes <prefix>.gz, so a path
        # without that suffix, or one ending in .bgz, would name a file that never appears. Fixing it
        # here keeps every later message pointing at the file that is actually produced.
        args.output_path = re.sub(r"\.bgz$", ".gz", args.output_path)
        if not args.output_path.endswith(".gz"):
            args.output_path = f"{args.output_path}.gz"

    return args


def compute_output_path(catalog_bed):
    """<input prefix>.extended.bed.gz"""
    filename_prefix = re.sub(r"\.bed(\.b?gz)?$", "", os.path.basename(catalog_bed))
    return f"{filename_prefix}.extended.bed.gz"


def catalog_record_key(record):
    """What makes two locus definitions the same definition: the interval and the motif.

    Deliberately ignores LocusId, since two records that describe the same interval with the same motif
    are duplicates however they are named.
    """
    return record["ReferenceRegion"], record["LocusStructure"]


def compute_extended_locus_id(record, taken_locus_ids):
    """Name an extended definition after the interval it now covers, matching the chrom-start-end-motif
    convention this repo already uses for BED-derived locus ids.
    """
    chrom, start_0based, end_1based = parse_interval(record["ReferenceRegion"])
    motif = parse_motifs_from_locus_structure(record["LocusStructure"])[0]
    locus_id = f"{chrom.replace('chr', '')}-{start_0based}-{end_1based}-{motif}"
    # A collision here means some other locus already carries this name while describing something
    # else, since an identical definition would have been dropped as a duplicate before we got here.
    return locus_id if locus_id not in taken_locus_ids else f"{locus_id}_extended"


def compute_span(record):
    """Base pairs the record's interval covers."""
    _, start_0based, end_1based = parse_interval(record["ReferenceRegion"])
    return end_1based - start_0based


def extend_catalog_record(record, fasta_obj, args, counters):
    """Return an extended copy of the record, or None if the rule declines to move either boundary."""
    motifs = parse_motifs_from_locus_structure(record["LocusStructure"])
    if len(motifs) != 1:
        # A TRGT-format BED row can define several motifs sharing one interval. Where each motif starts
        # inside that interval is not recorded, and the motif phase of the right boundary depends on
        # exactly that, so the locus cannot be extended correctly, nor written back out as one BED row.
        # Stop now rather than after processing the whole catalog.
        raise ValueError(f"This script handles one motif per locus, but {record['LocusId']} has "
                         f"{len(motifs)} ({record['LocusStructure']}). Split it into one row per motif.")

    chrom, start_0based, end_1based = parse_interval(record["ReferenceRegion"])
    normalized_chrom = normalize_chrom_using_pysam_fasta(fasta_obj, chrom)

    if args.verbose:
        print(f"{record['LocusId']}: {record['ReferenceRegion']}")
    left_extension, left_hit_limit = compute_extension(
        fasta_obj, normalized_chrom, start_0based, end_1based, motifs[0], "left",
        args.max_repeats_in_gap, args.min_purity_of_new_sequence, verbose=args.verbose)
    right_extension, right_hit_limit = compute_extension(
        fasta_obj, normalized_chrom, start_0based, end_1based, motifs[0], "right",
        args.max_repeats_in_gap, args.min_purity_of_new_sequence, verbose=args.verbose)

    if left_hit_limit or right_hit_limit:
        counters["loci whose extension was still growing when the search window stopped widening"] += 1

    if left_extension == 0 and right_extension == 0:
        return None

    extended_record = dict(record)
    extended_record["ReferenceRegion"] = (f"{chrom}:{start_0based - left_extension}"
                                          f"-{end_1based + right_extension}")
    return extended_record


def build_output_records(records_and_extensions, keep_original_definitions, counters):
    """Assemble the output catalog from each input record paired with its extended form or None.

    Args:
        records_and_extensions (list): (input record, extended record or None) pairs, in input order.
        keep_original_definitions (bool): keep the original definition of a locus that was extended,
            in addition to its extended definition, rather than replacing it.
        counters (collections.Counter): updated with what was added and skipped.

    Return:
        list: the output catalog records.
    """
    # Which definitions survive as themselves. An extended definition that matches one of these would
    # be a duplicate of a locus already in the output, so it is not added. This has to be computed over
    # the whole catalog up front, because the locus it collides with may come later in the file.
    kept_original_keys = {catalog_record_key(record)
                          for record, extended in records_and_extensions
                          if extended is None or keep_original_definitions}
    taken_locus_ids = {record["LocusId"] for record, extended in records_and_extensions
                       if extended is None or keep_original_definitions}

    output_records = []
    added_extended_keys = set()
    for record, extended in records_and_extensions:
        if extended is None or keep_original_definitions:
            output_records.append(record)

        if extended is None:
            continue

        extended_key = catalog_record_key(extended)
        if extended_key in kept_original_keys:
            counters["extended definitions not added because that definition was already in the catalog"] += 1
            continue
        if extended_key in added_extended_keys:
            counters["extended definitions not added because another locus extended to the same interval"] += 1
            continue

        added_extended_keys.add(extended_key)
        counters["extended definitions added to the output"] += 1
        counters["base pairs added"] += compute_span(extended) - compute_span(record)
        if keep_original_definitions:
            # The original keeps its own id and stays in the output, so the extended definition needs a
            # different one. When the original is being replaced its id is free, and reusing it is what
            # preserves meaningful names that a coordinate string would throw away.
            extended = dict(extended)
            extended["LocusId"] = compute_extended_locus_id(extended, taken_locus_ids)
            taken_locus_ids.add(extended["LocusId"])
        output_records.append(extended)

    return output_records


def write_output_catalog(output_records, output_path):
    """Write the BED rows, then bgzip and tabix them."""
    uncompressed_path = re.sub(r"\.b?gz$", "", output_path)
    print(f"Writing {output_path}")
    with open(uncompressed_path, "wt") as output_bed:
        for bed_record in sorted(convert_json_records_to_bed_format_tuples(output_records)):
            output_bed.write("\t".join(map(str, bed_record)) + "\n")
    os.system(f"bgzip -f {uncompressed_path}")
    os.system(f"tabix -f -p bed {uncompressed_path}.gz")


def main():
    args = parse_args()
    fasta_obj = pysam.FastaFile(args.reference_fasta)

    print(f"Parsing {args.catalog_bed}")
    catalog_iterator = get_variant_catalog_iterator(args.catalog_bed)
    if args.show_progress_bar:
        catalog_iterator = tqdm(catalog_iterator, unit=" loci", unit_scale=True)

    counters = collections.Counter()
    records_and_extensions = []
    for i, record in enumerate(catalog_iterator):
        if args.process_n_loci is not None and i >= args.process_n_loci:
            break
        counters["total"] += 1
        extended_record = extend_catalog_record(record, fasta_obj, args, counters)
        if extended_record is not None:
            counters["extended"] += 1
            if args.verbose:
                print(f"  -> {extended_record['ReferenceRegion']}")
        records_and_extensions.append((record, extended_record))

    output_records = build_output_records(
        records_and_extensions, args.keep_original_definitions_of_extended_loci, counters)

    locus_ids = [record["LocusId"] for record in output_records]
    duplicates = [locus_id for locus_id, count in collections.Counter(locus_ids).items() if count > 1]
    if duplicates:
        # These can only have come in from the catalog, since compute_extended_locus_id gives every id
        # it generates a suffix rather than let it collide. Say so and carry on: refusing to write the
        # output would throw away the whole run over a pre-existing problem in the input.
        print(f"WARNING: {len(duplicates):,d} LocusId(s) appear more than once in the output, carried "
              f"over from duplicates in the input catalog, starting with {duplicates[:3]}")

    write_output_catalog(output_records, args.output_path)

    print(f"The rule extended {counters['extended']:,d} out of {counters['total']:,d} loci "
          f"({counters['extended'] / max(counters['total'], 1):.1%}). "
          f"{counters['extended definitions added to the output']:,d} extended definitions reached the "
          f"output, adding {counters['base pairs added']:,d}bp")
    print(f"Wrote {len(output_records):,d} locus definitions to {args.output_path}")
    already_reported = ("total", "extended", "extended definitions added to the output",
                        "base pairs added")
    for label, count in counters.items():
        if label not in already_reported:
            print(f"{count:,d} {label}")


if __name__ == "__main__":
    main()
