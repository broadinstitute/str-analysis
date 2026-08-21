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

A locus is only put through the rule at all when its existing definition is itself pure enough:
at least MIN_REFERENCE_PURITY of the bases it already covers must match a perfect tiling of its own
motif. Where a locus is a poor match for the motif it is annotated with, its boundaries are not
something this rule can place sensibly, so it is left exactly as it is.

Polishing a locus that was extended:

The loop above moves a boundary by whole motif copies only, so it can stop with a run of bases that
do continue the repeat, but are fewer than one copy, still outside the locus. Every locus the rule
extended is therefore polished afterwards:

    1. each boundary is pushed outward, one base at a time, for as long as the flanking bases keep
       matching the motif's tiling, which adds up to len(motif) - 1 bases per side
    2. the motif is re-phased to the copy that now starts at the new left boundary, so the locus
       still satisfies the convention that its first base is the motif's first base. A CAG locus
       whose left boundary picks up one base becomes a GCA locus, for instance.

Loci the rule did not extend are left exactly as they were, motif included.

Output:

By default an extended definition replaces original locus definition. Specifying
--keep-original-definitions-of-extended-loci causes the original locus definition to be included
in the output alongside the extended definition, making the output a superset of the input.
Any loci that didn't get extended are included in the output regardless.

An extended definition is dropped before any of that when the input catalog already described the
same repeat with more: a locus that shares its canonical motif, overlaps it by any amount, and covers
more base pairs than it does. The locus it grew from is left as it was.

Extended definitions are deduplicated against each other in both modes: among the extended
definitions that overlap and share a canonical motif, only the longest is kept, taking them
longest-first so that a definition is dropped only when it overlaps one that was actually kept.
Definitions of equal length are ranked by coordinate, motif and locus id, so one of those can lose
to a definition that is no longer than it is.
Neighboring loci that grow onto the same repeat therefore collapse, though a run of definitions
that overlap only through each other collapses to as many as are needed to cover it, not to one. An
extended definition is also dropped when an identical definition is already in the output as an
input locus. When an extended definition is dropped and the original it came from is not being kept,
the original is put back unless what was kept for its repeat already describes every base it covered,
so the output never covers less sequence than the input and never repeats a definition it already has.


Usage:

    python3 -m str_analysis.extend_str_catalog_boundaries -R hg38.fa catalog.bed.gz

"""

import argparse
import bisect
import collections
import itertools
import os
import re

import pysam
from tqdm import tqdm

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.eh_catalog_utils import get_variant_catalog_iterator, \
    parse_motifs_from_locus_structure, convert_json_records_to_bed_format_tuples
from str_analysis.utils.fasta_utils import normalize_chrom_using_pysam_fasta
from str_analysis.utils.gap_purity_extension import compute_extension, compute_reference_purity, \
    extend_to_capture_partial_copy, MAX_REPEATS_IN_GAP, MIN_PURITY_OF_NEW_SEQUENCE, \
    MIN_REFERENCE_PURITY
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
    parser.add_argument("--min-reference-purity", type=float, default=MIN_REFERENCE_PURITY,
                        help="Only consider extending a locus whose existing definition is at least "
                             "this pure, meaning this fraction of the bases it already covers match a "
                             "perfect tiling of its own motif. A locus that is a poor match for its "
                             "own motif is not one whose boundaries this rule can place sensibly. "
                             "Set to 0 to attempt every locus.")
    parser.add_argument("-k", "--keep-original-definitions-of-extended-loci", action="store_true",
                        help="If specified, the original locus definitions will be included in the "
                             "output alongside the extended locus definitions. If not, any locus "
                             "definitions that are extended will be replaced with the new extended "
                             "definitions. Regardless, the extended definitions are deduplicated "
                             "against each other: taken longest-first, one is dropped when it overlaps "
                             "a definition already kept that shares its chromosome and canonical "
                             "motif, so a chain of overlapping definitions can still leave more than "
                             "one survivor. When an extended definition is dropped and the original is "
                             "not being kept, the original is put back unless what was kept already "
                             "covers it.")
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
    if not 0 <= args.min_reference_purity <= 1:
        parser.error(f"--min-reference-purity must be between 0 and 1, got "
                     f"{args.min_reference_purity}")
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


def rotate_motif(motif, bases_added_on_the_left):
    """Re-phase the motif onto a left boundary that moved out by the given number of bases.

    The locus's first base is its motif's first base, so moving the left boundary out by a whole
    number of copies leaves the motif alone, while moving it out by part of a copy starts the locus
    part-way through the motif and needs the motif rotated to match.
    """
    rotation = (-bases_added_on_the_left) % len(motif)
    return motif[rotation:] + motif[:rotation]


def polish_extended_locus(fasta_obj, chrom, start_0based, end_1based, motif, core_length):
    """Grow an extended locus over the partial motif copies its boundaries stopped short of.

    Args:
        fasta_obj (pysam.FastaFile): the reference genome.
        chrom (str): chromosome name, already matching the reference's naming.
        start_0based (int): the extended locus start.
        end_1based (int): the extended locus end.
        motif (str): the locus motif, in the phase it had before the extension.
        core_length (int): length of the locus before it was extended, which fixes the motif phase.
            The gap-purity rule moves each boundary by whole copies, so this still gives the right
            phase at the extended boundaries.

    Return:
        tuple: (start_0based, end_1based, motif) of the polished locus.
    """
    left_flank = fasta_obj.fetch(chrom, max(0, start_0based - len(motif) + 1), start_0based)[::-1]
    right_flank = fasta_obj.fetch(
        chrom, end_1based, min(fasta_obj.get_reference_length(chrom), end_1based + len(motif) - 1))
    bases_added_on_the_left = extend_to_capture_partial_copy(left_flank, motif, core_length, "left")
    bases_added_on_the_right = extend_to_capture_partial_copy(right_flank, motif, core_length, "right")

    return (start_0based - bases_added_on_the_left,
            end_1based + bases_added_on_the_right,
            rotate_motif(motif, bases_added_on_the_left))


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

    reference_purity = compute_reference_purity(
        fasta_obj.fetch(normalized_chrom, start_0based, end_1based), motifs[0])
    if reference_purity < args.min_reference_purity:
        counters["loci not considered because the existing definition is not pure enough"] += 1
        if args.verbose:
            print(f"    reference purity={reference_purity:.3f} < {args.min_reference_purity} "
                  f"-- leave this locus alone")
        return None

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

    polished_start, polished_end, polished_motif = polish_extended_locus(
        fasta_obj, normalized_chrom, start_0based - left_extension, end_1based + right_extension,
        motifs[0], end_1based - start_0based)
    if polished_motif != motifs[0]:
        counters["extended definitions whose motif was re-phased onto the new left boundary"] += 1

    extended_record = dict(record)
    extended_record["ReferenceRegion"] = f"{chrom}:{polished_start}-{polished_end}"
    extended_record["LocusStructure"] = f"({polished_motif})*"
    return extended_record


def select_longest_non_overlapping(extended_records):
    """Choose which extended definitions to keep when several grow onto the same repeat.

    Definitions are considered longest-first, and one is kept whenever it doesn't overlap a
    definition that was already kept and shares its chromosome and canonical motif. Going
    longest-first means a definition is only ever dropped for one that ranked ahead of it and really
    does overlap it, so a chain of overlaps doesn't collapse further than it has to: of [100,200),
    [150,260) and [250,500), the middle one is dropped and both ends survive. Definitions of the same
    length rank by the tie-breakers below, so one of those can lose to an equally long definition.

    Args:
        extended_records (list): the extended definitions to choose between.

    Return:
        tuple: (the indices into extended_records of the definitions to keep,
                (chrom, canonical motif) -> the (starts, ends) of what was kept for it, sorted and
                disjoint)
    """
    kept_intervals_per_group = collections.defaultdict(lambda: ([], []))
    kept_indices = set()

    definitions = [parse_interval(record["ReferenceRegion"])
                   + (parse_motifs_from_locus_structure(record["LocusStructure"])[0],
                      record["LocusId"], index)
                   for index, record in enumerate(extended_records)]
    # Length, then coordinate, motif and locus id, so that equally long definitions are resolved the
    # same way every run rather than in whatever order the catalog happened to list them. The locus id
    # only settles definitions that are identical, where the choice is between two records that
    # describe the same thing, so it decides which record is returned rather than anything about the
    # catalog that gets written.
    definitions.sort(key=lambda d: (-(d[2] - d[1]), d[0], d[1], d[3], d[4]))

    for chrom, start_0based, end_1based, motif, _, index in definitions:
        starts, ends = kept_intervals_per_group[chrom, compute_canonical_motif(motif)]
        # What's been kept for a group never overlaps, so it stays sorted and disjoint. Everything in
        # it is also at least as long as this definition, so it cannot hold an interval that fits
        # inside this one, and the only two that can overlap it are the interval starting at or
        # before it and the one after that.
        position = bisect.bisect_right(starts, start_0based) - 1
        if any(0 <= i < len(starts) and starts[i] < end_1based and ends[i] > start_0based
               for i in (position, position + 1)):
            continue
        starts.insert(position + 1, start_0based)
        ends.insert(position + 1, end_1based)
        kept_indices.add(index)

    return kept_indices, kept_intervals_per_group


def is_covered_by_kept_definitions(record, kept_intervals_per_group):
    """Whether the definitions that were kept reach over the whole of the record's interval.

    What was kept for a group is disjoint but can touch, and two definitions that touch describe the
    sequence between them just as fully as one would, so this walks the run of touching intervals
    rather than looking for a single interval that contains the record.
    """
    chrom, start_0based, end_1based = parse_interval(record["ReferenceRegion"])
    motif = parse_motifs_from_locus_structure(record["LocusStructure"])[0]
    starts, ends = kept_intervals_per_group.get((chrom, compute_canonical_motif(motif)), ([], []))

    position = bisect.bisect_right(starts, start_0based) - 1
    if position < 0 or ends[position] <= start_0based:
        return False
    while (ends[position] < end_1based and position + 1 < len(starts)
           and starts[position + 1] == ends[position]):
        position += 1
    return ends[position] >= end_1based


def discard_extensions_outdone_by_the_input(records_and_extensions, counters):
    """Drop an extended definition where the input catalog already described its repeat with more.

    An extension is only worth having when it says more about the repeat than the catalog already did.
    Where an input locus shares the extension's canonical motif, overlaps it, and covers more base
    pairs than it does, the catalog was already the better description, so the extension is dropped
    and the locus it came from is left exactly as it was. A locus's own original can never trigger
    this, since an extension always covers strictly more than what it grew from.

    Any overlap counts, however small, and an input locus gets its say whether or not it is pure
    enough for the rule to have extended it. Both are deliberate: a locus is being read here as a
    statement about where a repeat is, which does not depend on the rule being able to place its
    boundaries. On the v2 catalog this drops 54,195 extensions, 95.7% of them to an input locus that
    covers the whole extension.

    Args:
        records_and_extensions (list): (input record, extended record or None) pairs, in input order.
        counters (collections.Counter): updated with how many extensions were dropped.

    Return:
        list: the same pairs, with a dropped extension replaced by None.
    """
    # Indexing the extensions and reading the catalog past them, rather than the other way around:
    # there are far fewer extensions than input loci, and only the extensions need to be looked up.
    canonical_motifs = {}

    def canonical_motif_of(record):
        motif = parse_motifs_from_locus_structure(record["LocusStructure"])[0]
        if motif not in canonical_motifs:
            canonical_motifs[motif] = compute_canonical_motif(motif)
        return canonical_motifs[motif]

    extensions_per_group = collections.defaultdict(list)
    for index, (_, extended) in enumerate(records_and_extensions):
        if extended is not None:
            chrom, start_0based, end_1based = parse_interval(extended["ReferenceRegion"])
            extensions_per_group[chrom, canonical_motif_of(extended)].append(
                (start_0based, end_1based, index))

    index_per_group = {}
    for group, extensions in extensions_per_group.items():
        extensions.sort()
        starts, ends, indices = (list(column) for column in zip(*extensions))
        # The highest end among everything up to each position, which is what lets a backwards scan
        # stop as soon as nothing further back can still reach the locus being checked.
        index_per_group[group] = (starts, ends, indices,
                                  list(itertools.accumulate(ends, max)))

    outdone = set()
    for record, _ in records_and_extensions:
        chrom, start_0based, end_1based = parse_interval(record["ReferenceRegion"])
        if (chrom, canonical_motif_of(record)) not in index_per_group:
            continue
        starts, ends, indices, highest_end = index_per_group[chrom, canonical_motif_of(record)]
        span = end_1based - start_0based
        position = bisect.bisect_left(starts, end_1based) - 1
        while position >= 0 and highest_end[position] > start_0based:
            if ends[position] > start_0based and ends[position] - starts[position] < span:
                outdone.add(indices[position])
            position -= 1

    counters["extended definitions dropped because a longer locus in the input already covers "
             "that repeat"] += len(outdone)
    return [(record, None if index in outdone else extended)
            for index, (record, extended) in enumerate(records_and_extensions)]


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

    # Deciding this needs every extended definition in hand, since the one that displaces a given
    # definition can come from anywhere in the catalog.
    candidates = [extended for _, extended in records_and_extensions
                  if extended is not None and catalog_record_key(extended) not in kept_original_keys]
    kept_indices, kept_intervals_per_group = select_longest_non_overlapping(candidates)
    # Keyed by identity, not equality: two loci can extend to the very same definition, and which
    # copy of it survives is exactly what has been decided here.
    kept_candidates = {id(candidates[index]) for index in kept_indices}

    output_records = []
    for record, extended in records_and_extensions:
        if extended is None or keep_original_definitions:
            output_records.append(record)

        if extended is None:
            continue

        if catalog_record_key(extended) in kept_original_keys:
            counters["extended definitions not added because that definition was already in the catalog"] += 1
            continue
        if id(extended) not in kept_candidates:
            counters["extended definitions dropped in favor of an overlapping extended definition that was kept"] += 1
            if not keep_original_definitions and not is_covered_by_kept_definitions(
                    record, kept_intervals_per_group):
                # The definitions that displaced this one only have to overlap it, so they don't
                # always reach over everything this locus covered. Put the locus back when they
                # don't, rather than let the catalog lose sequence it used to describe. Asking the
                # whole group, rather than only what displaced this one definition, is what keeps a
                # locus that some other kept definition already covers from coming back as a
                # duplicate.
                output_records.append(record)
            continue

        counters["extended definitions added to the output"] += 1
        if keep_original_definitions:
            # The original keeps its own id and stays in the output, so the extended definition needs a
            # different one. When the original is being replaced its id is free and the extended
            # definition inherits it. Ids only matter to callers working with these records: the BED
            # this script writes has no name column, so none of them reach the output file.
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
                print(f"  -> {extended_record['ReferenceRegion']} "
                      f"{extended_record['LocusStructure']}")
        records_and_extensions.append((record, extended_record))

    records_and_extensions = discard_extensions_outdone_by_the_input(records_and_extensions, counters)
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

    # Taken from the finished catalog rather than accumulated per extension, since a locus can also
    # leave the output when its extension is dropped for an overlapping one, and a figure that only
    # ever added would report growth the catalog didn't get.
    base_pairs_added = (sum(compute_span(record) for record in output_records)
                        - sum(compute_span(record) for record, _ in records_and_extensions))
    print(f"The rule extended {counters['extended']:,d} out of {counters['total']:,d} loci "
          f"({counters['extended'] / max(counters['total'], 1):.1%}). "
          f"{counters['extended definitions added to the output']:,d} extended definitions reached the "
          f"output, whose locus definitions now span {base_pairs_added:,d}bp more in total than the "
          f"input's did")
    print(f"Wrote {len(output_records):,d} locus definitions to {args.output_path}")
    already_reported = ("total", "extended", "extended definitions added to the output")
    for label, count in counters.items():
        if label not in already_reported:
            print(f"{count:,d} {label}")


if __name__ == "__main__":
    main()
