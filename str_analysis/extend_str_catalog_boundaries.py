"""Extend tandem repeat locus boundaries outward into their flanking sequence using the gap-purity rule.

WHAT THIS DOES
Catalog boundaries routinely stop short of where the repeat actually ends, because the tools that
generated them stop at the first imperfect copy. This script walks outward from each boundary and
pulls in flanking sequence that is still recognizably the same repeat, even when a copy or two along
the way is interrupted by a substitution.

WHY IT MATTERS
A boundary that stops early leaves real repeat sequence in the flank, where a genotyper treats it as
unique anchoring sequence rather than as part of the repeat. Extending the locus to cover it measurably
improves how often ExpansionHunter's genotypes match assembly-derived truth genotypes.

WHY IT IS CONSERVATIVE, AND WHERE IT STOPS
The rule only ever grows a locus, never shrinks or shifts it, and it moves in whole motif copies, so a
locus's count of leftover partial-copy bases is exactly what it was before. That matters because
ExpansionHunter models a locus as a whole number of copies, so changing the leftover count moves its
quality score for reasons unrelated to whether the boundary is better. The purity test below is what
keeps it from running away: a single mismatched base has to be paid for by roughly ten bases of clean
repeat, so an isolated interrupted copy with nothing good after it is refused.

ALGORITHM, PER SIDE OF EACH LOCUS, INDEPENDENTLY

    boundary = the locus's current start (going left) or end (going right)
    loop:
        # 1. look outward for the next exact copy of the motif, in the locus's own frame,
        #    tolerating a short run of interrupted copies on the way there
        scan outward from boundary in steps of one motif length
            stop the scan once an exact copy of the motif is found  -> anchor
            give up if more than MAX_REPEATS_IN_GAP copies pass without an exact one
        if no anchor was found:
            stop, this side is done
        # 2. take the whole unbroken run of exact copies that follows the anchor
        run_end = anchor
        while the copy starting at run_end is exact:
            run_end += motif length
        # 3. the candidate addition is everything from the boundary through that run,
        #    the imperfect gap included, and it is judged as a whole
        purity = fraction of bases in [boundary, run_end) matching a perfect tiling of the motif
        if purity >= MIN_PURITY_OF_NEW_SEQUENCE:
            boundary = run_end          # accept, then try to extend further from here
        else:
            stop, this side is done     # rejection is final for this side

Sequence is compared against the motif tiled in the locus's existing frame, so a repeat in the flank
that is out of phase with the boundary will not match. IUPAC ambiguity codes in the motif (RFC1's
AARRG, RUNX2's GCN) match any base they represent. One input path does not reach them: the shared BED
reader filters PLAIN BED rows to ACGTN motifs, so a plain BED row whose motif uses any of R, Y, S, W,
K, M, B, D, H or V is dropped with a warning before this script sees it. N survives that filter, and
JSON catalogs and TRGT-format BED rows are not filtered at all, so any code reaches the rule through
those.

WHAT THE OUTPUT CONTAINS
By default an extended definition replaces the locus it came from, so the output has one definition per
input locus. With --keep-original-definitions-of-extended-loci the original is kept as well and the
extended definition is added alongside it, making the output a superset of the input, which is useful
when you want to compare the two or not commit to the extension yet. A locus the rule declines to
extend is passed through unchanged either way.

Extended definitions are deduplicated on their interval and motif in both modes. Neighbouring loci in a
catalog frequently describe overlapping pieces of the same underlying repeat, and extending them can
land them all on one interval, so without this the output would carry several identical copies of it.
An extended definition is also dropped when an identical definition is already in the output as an
input locus.

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
from str_analysis.utils.export_json import export_json
from str_analysis.utils.fasta_utils import normalize_chrom_using_pysam_fasta
from str_analysis.utils.gap_purity_extension import compute_extension, \
    MAX_REPEATS_IN_GAP, MIN_PURITY_OF_NEW_SEQUENCE
from str_analysis.utils.misc_utils import parse_interval

# Annotations that annotate_and_filter_str_catalog.py computes from a locus's interval. Once the
# interval moves they describe the wrong sequence, and none of them can be recomputed here: the
# mappability ones need a bigWig track and the purity ones need a different pass over the reference.
# Dropping them leaves the field absent, which is honest, where keeping them would be quietly wrong.
# Re-running annotate_and_filter_str_catalog.py on the output restores all of them correctly.
ANNOTATIONS_INVALIDATED_BY_EXTENSION = (
    "NumRepeatsInReference", "ReferenceRepeatPurity", "NsInFlanks",
    "LeftFlankMappability", "RightFlankMappability", "FlanksAndLocusMappability",
    "HighestPurityMotif", "HighestPurityMotifPurity", "HighestPurityMotifQuality",
    "TandemRepeatFinderMotif", "TandemRepeatFinderMotifQuality",
    # set only when the boundaries match a known disease locus exactly, which extending them undoes
    "KnownDiseaseAssociatedLocus",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extend tandem repeat locus boundaries into their flanking sequence using the "
                    "gap-purity rule.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-R", "--reference-fasta", required=True, help="Reference fasta file path.")
    parser.add_argument("-o", "--output-path", help="Output path. Defaults to the input path with "
                        "'.extended' added before the file extension.")
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
    parser.add_argument("catalog_json_or_bed", help="A catalog of tandem repeats in ExpansionHunter "
                        "JSON format or BED format.")
    args = parser.parse_args()

    if not 0 < args.min_purity_of_new_sequence <= 1:
        parser.error(f"--min-purity-of-new-sequence must be > 0 and <= 1, got "
                     f"{args.min_purity_of_new_sequence}")
    if args.max_repeats_in_gap < 0:
        parser.error(f"--max-repeats-in-gap cannot be negative, got {args.max_repeats_in_gap}")
    if not os.path.isfile(args.reference_fasta):
        parser.error(f"Reference fasta not found: {args.reference_fasta}")

    if not args.output_path:
        args.output_path = compute_output_path(args.catalog_json_or_bed)
    elif ".json" not in args.output_path:
        # BED output is bgzipped and tabix-indexed, and bgzip always writes <prefix>.gz, so a path
        # without that suffix, or one ending in .bgz, would name a file that never appears. Fixing it
        # here keeps every later message pointing at the file that is actually produced.
        args.output_path = re.sub(r"\.bgz$", ".gz", args.output_path)
        if not args.output_path.endswith(".gz"):
            args.output_path = f"{args.output_path}.gz"

    return args


def compute_output_path(catalog_json_or_bed):
    """<input prefix>.extended.<input extension>.gz, keeping the input's json or bed format."""
    filename = os.path.basename(catalog_json_or_bed)
    filename_prefix = re.sub(r"\.(json|bed)(\.b?gz)?$", "", filename)
    extension = "json" if ".json" in filename else "bed"
    return f"{filename_prefix}.extended.{extension}.gz"


def get_reference_regions(record):
    """The record's ReferenceRegion as a list, whether it holds one region or several."""
    reference_regions = record["ReferenceRegion"]
    return reference_regions if isinstance(reference_regions, list) else [reference_regions]


def compute_span(record):
    """Total base pairs the record's reference regions cover."""
    return sum(end - start for _, start, end in map(parse_interval, get_reference_regions(record)))


def catalog_record_key(record):
    """What makes two locus definitions the same definition: their intervals and their motifs.

    Deliberately ignores LocusId, since two records that describe the same interval with the same motif
    are duplicates however they are named.
    """
    return (tuple(get_reference_regions(record)),
            tuple(parse_motifs_from_locus_structure(record["LocusStructure"])))


def rename_variant_ids(variant_ids, old_locus_id, new_locus_id):
    """Point a record's VariantId at its new LocusId, keeping any per-region suffixes intact.

    Returns None when no entry mentions the old LocusId, since then substitution cannot tell the two
    records apart and leaving the value alone would put the same VariantId on both. The caller drops
    the field in that case, which is safe because VariantId is optional and consumers fall back to
    LocusId.
    """
    is_list = isinstance(variant_ids, list)
    renamed = [variant_id.replace(old_locus_id, new_locus_id)
               for variant_id in (variant_ids if is_list else [variant_ids])]
    if renamed == (variant_ids if is_list else [variant_ids]):
        return None
    return renamed if is_list else renamed[0]


def compute_extended_locus_id(record, taken_locus_ids):
    """Name an extended definition after the interval it now covers, matching the chrom-start-end-motif
    convention this repo already uses for BED-derived locus ids.
    """
    reference_regions = get_reference_regions(record)
    chrom, start_0based, _ = parse_interval(reference_regions[0])
    _, _, end_1based = parse_interval(reference_regions[-1])
    motifs = parse_motifs_from_locus_structure(record["LocusStructure"])
    locus_id = f"{chrom.replace('chr', '')}-{start_0based}-{end_1based}-{'-'.join(motifs)}"
    # A collision here means some other locus already carries this name while describing something
    # else, since an identical definition would have been dropped as a duplicate before we got here.
    return locus_id if locus_id not in taken_locus_ids else f"{locus_id}_extended"


def extend_catalog_record(record, fasta_obj, args, counters):
    """Return an extended copy of the record, or None if the rule declines to move either boundary.

    A record with several adjacent ReferenceRegions is extended only at its two outer boundaries: the
    interior boundaries are where one repeat hands off to the next, so moving them would overlap the
    neighbouring region and no longer describe the same locus.
    """
    if record.get("Diseases"):
        # Thresholds like NormalMax and PathogenicMin count repeats from this exact boundary, so
        # widening it would shift every genotype called here and quietly invalidate them.
        counters["disease-associated loci left unextended, since their repeat-count thresholds are "
                 "tied to the current boundary"] += 1
        return None

    reference_regions = get_reference_regions(record)
    is_single_region = not isinstance(record["ReferenceRegion"], list)

    motifs = parse_motifs_from_locus_structure(record["LocusStructure"])
    if len(motifs) != len(reference_regions):
        # A TRGT-format BED row parses into several motifs sharing one spanning ReferenceRegion. Where
        # each motif starts inside that span is not recorded, and the motif phase of the right boundary
        # depends on exactly that, so there is no way to extend the locus correctly. Pass it through
        # unextended rather than aborting the run or guessing.
        counters["loci passed through unextended because their motifs and reference regions do not "
                 "correspond one-to-one (convert them to a JSON catalog with one region per motif)"] += 1
        return None

    intervals = [parse_interval(reference_region) for reference_region in reference_regions]
    chrom = intervals[0][0]
    normalized_chrom = normalize_chrom_using_pysam_fasta(fasta_obj, chrom)

    if args.verbose:
        print(f"{record['LocusId']}: {record['ReferenceRegion']}")
    left_extension, left_hit_limit = compute_extension(
        fasta_obj, normalized_chrom, intervals[0][1], intervals[0][2], motifs[0], "left",
        args.max_repeats_in_gap, args.min_purity_of_new_sequence, verbose=args.verbose)
    right_extension, right_hit_limit = compute_extension(
        fasta_obj, normalized_chrom, intervals[-1][1], intervals[-1][2], motifs[-1], "right",
        args.max_repeats_in_gap, args.min_purity_of_new_sequence, verbose=args.verbose)

    if left_hit_limit or right_hit_limit:
        counters["loci whose extension was still growing when the search window stopped widening"] += 1

    if left_extension == 0 and right_extension == 0:
        return None

    intervals[0] = (intervals[0][0], intervals[0][1] - left_extension, intervals[0][2])
    intervals[-1] = (intervals[-1][0], intervals[-1][1], intervals[-1][2] + right_extension)

    extended_regions = [f"{region_chrom}:{start}-{end}" for region_chrom, start, end in intervals]
    extended_record = dict(record)
    extended_record["ReferenceRegion"] = extended_regions[0] if is_single_region else extended_regions
    # Catalogs annotated by annotate_and_filter_str_catalog.py also carry MainReferenceRegion. Left
    # alone it would still name the pre-extension interval, so the record's two coordinate fields would
    # disagree. Only the outer regions moved, so only those two need remapping.
    main_reference_region = record.get("MainReferenceRegion")
    if main_reference_region == reference_regions[0]:
        extended_record["MainReferenceRegion"] = extended_regions[0]
    elif main_reference_region == reference_regions[-1]:
        extended_record["MainReferenceRegion"] = extended_regions[-1]

    for annotation in ANNOTATIONS_INVALIDATED_BY_EXTENSION:
        if extended_record.pop(annotation, None) is not None:
            counters[f"extended loci whose stale '{annotation}' annotation was dropped"] += 1

    return extended_record


def build_output_records(records_and_extensions, keep_original_definitions, counters):
    """Assemble the output catalog from each input record paired with its extended form or None.

    Args:
        records_and_extensions (list): (input record, extended record or None) pairs, in input order.
        keep_original_definitions (bool): keep the original definition of a locus that was extended,
            in addition to its extended definition, rather than replacing it.
        counters (collections.Counter): updated with what was kept, added and skipped.

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
    emitted_keys = set()
    for record, extended in records_and_extensions:
        if extended is None or keep_original_definitions:
            output_records.append(record)
            emitted_keys.add(catalog_record_key(record))

        if extended is None:
            continue

        extended_key = catalog_record_key(extended)
        if extended_key in kept_original_keys or extended_key in emitted_keys:
            if extended_key in kept_original_keys:
                counters["extended definitions not added because that definition was already in the catalog"] += 1
            else:
                counters["extended definitions not added because another locus extended to the same interval"] += 1
            if not keep_original_definitions:
                # Replacing this locus with a definition that is already in the output would delete the
                # locus outright, so fall back to its own original definition. That fallback is subject
                # to the same dedup as everything else: if the original is itself already in the
                # output, adding it would create the duplicate this function exists to prevent.
                original_key = catalog_record_key(record)
                if original_key in kept_original_keys or original_key in emitted_keys:
                    counters["loci contributing no definition of their own, since both their original "
                             "and extended definitions were already in the output"] += 1
                else:
                    output_records.append(record)
                    emitted_keys.add(original_key)
            continue

        emitted_keys.add(extended_key)
        counters["extended definitions added to the output"] += 1
        counters["base pairs added"] += compute_span(extended) - compute_span(record)
        if keep_original_definitions:
            # The original keeps its own id and stays in the output, so the extended definition needs a
            # different one. When the original is being replaced its id is free, and reusing it is what
            # preserves meaningful names like HTT or ATXN7 that a coordinate string would throw away.
            extended = dict(extended)
            extended["LocusId"] = compute_extended_locus_id(extended, taken_locus_ids)
            taken_locus_ids.add(extended["LocusId"])
            # VariantId is derived from LocusId, so leaving it alone would put the same VariantId on
            # two different locus definitions, which downstream tools use as a unique key.
            if "VariantId" in extended:
                renamed_variant_ids = rename_variant_ids(
                    extended["VariantId"], record["LocusId"], extended["LocusId"])
                if renamed_variant_ids is None:
                    del extended["VariantId"]
                    counters["extended definitions whose VariantId was dropped, since it does not "
                             "derive from the LocusId and would have collided with the original"] += 1
                else:
                    extended["VariantId"] = renamed_variant_ids
        output_records.append(extended)

    return output_records


def write_output_catalog(output_records, output_path):
    """Write JSON directly, or write BED plain and then bgzip and tabix it."""
    if ".json" in output_path:
        export_json(output_records, output_path)
        return

    # convert_json_records_to_bed_format_tuples needs one region per motif, since a BED row carries a
    # single interval and a single motif. A TRGT-format BED input parses into records that violate that,
    # and they were already passed through unextended for the same reason, so say so here with the way
    # out rather than letting the shared converter raise from three frames down.
    unrepresentable = sum(1 for record in output_records
                          if len(parse_motifs_from_locus_structure(record["LocusStructure"])) !=
                          len(get_reference_regions(record)))
    if unrepresentable:
        raise ValueError(f"{unrepresentable:,d} of the {len(output_records):,d} locus definitions have "
                         f"several motifs sharing one reference region, which BED format cannot "
                         f"represent. Re-run with a .json output path to keep them.")

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

    print(f"Parsing {args.catalog_json_or_bed}")
    catalog_iterator = get_variant_catalog_iterator(args.catalog_json_or_bed)
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
