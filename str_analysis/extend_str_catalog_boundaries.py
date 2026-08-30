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

A locus whose motif contains IUPAC codes is never put through the rule. Such a motif stands for a set
of sequences rather than one, so it matches flanking sequence that a concrete motif would not, and the
rule would read that as the repeat continuing.

A locus is only put through the rule at all when its existing definition is itself pure enough:
at least MIN_REFERENCE_PURITY of the bases it already covers must match a perfect tiling of its own
motif. Where a locus is a poor match for the motif it is annotated with, its boundaries are not
something this rule can place sensibly, so it is left exactly as it is.

Polishing a locus that was extended:

The loop above moves a boundary by whole motif copies only, so it can stop with a run of bases that
do continue the repeat, but are fewer than one copy, still outside the locus. Every locus the rule
extended therefore has its right boundary pushed outward afterwards, one base at a time, for as long
as the flanking bases keep matching the motif's tiling, which adds up to at most len(motif) - 1 bases.

The left boundary is never moved by part of a copy, so an extended locus is still read in the same
frame as the locus it grew from. The locus's first base is its motif's first base, and taking a
partial copy on the left would leave the locus starting part-way through the motif, with the motif
having to be rotated to match: a CAG locus that picked up one base on the left would become a GCA
locus. At EP400, for instance, chr12:132062548-132062611 (CAG)* extends to
chr12:132062524-132062611 (CAG)*, stopping short of the G at 132062523 that a rotated motif could
have taken.

Loci the rule did not extend are left exactly as they were.

Output:

There are three output modes. By default an extended definition replaces the original locus
definition. Specifying --keep-original-definitions-of-extended-loci causes the original locus
definition to be included in the output alongside the extended definition, making the output a
superset of the input. In both of those the output is a catalog in its own right, and any loci that
didn't get extended are included regardless.

--only-output-extended-loci is the third, and writes only the definitions the rule produced, leaving
out every locus the input catalog held. What it writes is therefore not a catalog on its own: it is
the set of definitions this rule adds to the one it was given, meant to be appended to that catalog
as a separate source of loci. The two paragraphs below about putting a displaced locus back do not
apply to it, since the locus in question never left the catalog this output gets appended to.

An extended definition is dropped before any of that when the input catalog already described the
same repeat with more: a locus that shares its canonical motif, overlaps it by any amount, and covers
more base pairs than it does. The locus it grew from is left as it was.

Extended definitions are deduplicated against each other in every mode: among the extended
definitions that overlap and share a canonical motif, only the longest is kept, taking them
longest-first so that a definition is dropped only when it overlaps one that was actually kept.
Definitions of equal length are ranked by coordinate, motif and locus id, so one of those can lose
to a definition that is no longer than it is.
Neighboring loci that grow onto the same repeat therefore collapse, though a run of definitions
that overlap only through each other collapses to as many as are needed to cover it, not to one. An
extended definition is also dropped when an identical definition is already in the output as an
input locus. In the two whole-catalog modes, when an extended definition is dropped and the original
it came from is not being kept, the original is put back unless what was kept for its repeat already
describes every base it covered, so the output never covers less sequence than the input and never
repeats a definition it already has.


The input catalog can be in BED or JSON format, and the output is written in the same format as the
input. Every locus must be defined by one interval and one motif: a JSON record that gives a locus
several adjacent intervals, or a BED row that lists several motifs for one interval, stops the run,
since which interval each boundary belongs to and the motif phase at each one would have to be
tracked separately. Split such a catalog first with
str_analysis.split_adjacent_loci_in_expansion_hunter_catalog.

BED carries only the interval and the motif, so a BED run drops every other field, while a JSON
run passes all of them through untouched. ReferenceRegion and LocusId are the only fields this tool
rewrites, so fields computed from the original boundaries, such as NumRepeatsInReference or
ReferenceRepeatPurity, survive a JSON run while describing the locus as it was before it grew, and
need recomputing afterwards.


Usage:

    python3 -m str_analysis.extend_str_catalog_boundaries -R hg38.fa catalog.bed.gz
    python3 -m str_analysis.extend_str_catalog_boundaries -R hg38.fa catalog.json.gz
    python3 -m str_analysis.extend_str_catalog_boundaries -R hg38.fa --only-output-extended-loci catalog.bed.gz

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
from str_analysis.utils.export_json import export_json
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
    parser.add_argument("-o", "--output-path", help="Output path. Defaults to the input path with "
                        "'.extended' added before the file extension, in the same format as the input. "
                        "A BED output is always bgzipped and tabix-indexed, and a JSON output gzipped.")
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
    parser.add_argument("--only-output-extended-loci", action="store_true",
                        help="Write only the definitions the rule produced, leaving out every locus "
                             "from the input catalog. The output is then not a catalog in its own "
                             "right: it is the set of definitions this rule adds to the one it was "
                             "given, meant to be appended to that catalog as a separate source of "
                             "loci. Extended definitions are still dropped when the input already "
                             "describes the same repeat, so nothing here duplicates the input.")
    parser.add_argument("-n", "--process-n-loci", type=int, help="Only process the first N loci.")
    parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar.")
    parser.add_argument("--verbose", action="store_true", help="Print every accept/reject decision the "
                        "rule makes, with the purity behind it. Intended for looking at a handful of "
                        "loci, not a whole catalog.")
    parser.add_argument("catalog_json_or_bed", help="A catalog of tandem repeats in BED or JSON "
                        "format, with one interval and one motif per locus.")
    args = parser.parse_args()

    if not 0 < args.min_purity_of_new_sequence <= 1:
        parser.error(f"--min-purity-of-new-sequence must be > 0 and <= 1, got "
                     f"{args.min_purity_of_new_sequence}")
    if not 0 <= args.min_reference_purity <= 1:
        parser.error(f"--min-reference-purity must be between 0 and 1, got "
                     f"{args.min_reference_purity}")
    if args.only_output_extended_loci and args.keep_original_definitions_of_extended_loci:
        parser.error("--only-output-extended-loci writes none of the input's own definitions, so it "
                     "cannot be combined with --keep-original-definitions-of-extended-loci, which "
                     "asks for them to be kept.")
    if args.max_repeats_in_gap < 0:
        parser.error(f"--max-repeats-in-gap cannot be negative, got {args.max_repeats_in_gap}")
    if not os.path.isfile(args.reference_fasta):
        parser.error(f"Reference fasta not found: {args.reference_fasta}")
    if not args.output_path:
        args.output_path = compute_output_path(args.catalog_json_or_bed)
    else:
        # The output is compressed, and bgzip and gzip both write <prefix>.gz, so a path without that
        # suffix, or one ending in .bgz, would name a file that never appears. Fixing it here keeps
        # every later message pointing at the file that is actually produced.
        args.output_path = re.sub(r"\.bgz$", ".gz", args.output_path)
        if not args.output_path.endswith(".gz"):
            args.output_path = f"{args.output_path}.gz"
        if is_json_path(args.output_path) != is_json_path(args.catalog_json_or_bed):
            # The format written is taken from the output path, and BED keeps only the interval and
            # the motif, so a JSON catalog sent to a path that doesn't say .json would quietly lose
            # every other field.
            expected_format = "JSON" if is_json_path(args.catalog_json_or_bed) else "BED"
            parser.error(f"--output-path must name a {expected_format} file, to match the input "
                         f"catalog: got {args.output_path}")

    return args


def compute_output_path(catalog_json_or_bed):
    """<input prefix>.extended.json.gz or <input prefix>.extended.bed.gz, matching the input format."""
    filename = os.path.basename(catalog_json_or_bed)
    if is_json_path(filename):
        return re.sub(r"\.json(\.b?gz)?$", "", filename) + ".extended.json.gz"

    return re.sub(r"\.bed(\.b?gz)?$", "", filename) + ".extended.bed.gz"


def is_json_path(path):
    """Whether the path names a JSON catalog rather than a BED one. Matches how
    get_variant_catalog_iterator decides which parser to use, so input and output stay consistent.
    """
    return ".json" in path


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


def polish_right_boundary(fasta_obj, chrom, start_0based, end_1based, motif):
    """Grow an extended locus over the partial motif copy its right boundary stopped short of.

    Args:
        fasta_obj (pysam.FastaFile): the reference genome.
        chrom (str): chromosome name, already matching the reference's naming.
        start_0based (int): the extended locus start.
        end_1based (int): the extended locus end.
        motif (str): the locus motif.

    Return:
        int: the polished end_1based.
    """
    right_flank = fasta_obj.fetch(
        chrom, end_1based, min(fasta_obj.get_reference_length(chrom), end_1based + len(motif) - 1))

    # The locus starts on the motif's first base, both before and after the extension, so its own
    # length is the phase the flank continues from.
    return end_1based + extend_to_capture_partial_copy(right_flank, motif, end_1based - start_0based)


def extend_catalog_record(record, fasta_obj, args, counters):
    """Return an extended copy of the record, or None if the rule declines to move either boundary."""
    if isinstance(record["ReferenceRegion"], list):
        # A JSON catalog can define one locus as several adjacent intervals. Which interval each
        # boundary belongs to, and the motif phase at each one, would have to be tracked separately,
        # so stop now rather than mangle the definition.
        raise ValueError(f"This script handles one interval per locus, but {record['LocusId']} has "
                         f"{len(record['ReferenceRegion'])}. Split it into one record per interval "
                         f"with str_analysis.split_adjacent_loci_in_expansion_hunter_catalog.")

    motifs = parse_motifs_from_locus_structure(record["LocusStructure"])
    if len(motifs) != 1:
        # A TRGT-format BED row can define several motifs sharing one interval. Where each motif starts
        # inside that interval is not recorded, and the motif phase of the right boundary depends on
        # exactly that, so the locus cannot be extended correctly, nor written back out as one BED row.
        # Stop now rather than after processing the whole catalog.
        raise ValueError(f"This script handles one motif per locus, but {record['LocusId']} has "
                         f"{len(motifs)} ({record['LocusStructure']}). Split it into one row per motif.")

    if set(motifs[0]) - set("ACGT"):
        # An IUPAC code stands for a set of bases, so it matches more of the flanking sequence than a
        # concrete motif would, and the rule would read that as the repeat continuing. GCN, for
        # instance, matches any GC followed by anything at all. Where the motif is written this
        # loosely its boundaries are not something this rule can place, so the locus is left alone.
        counters["loci not considered because the motif contains IUPAC codes"] += 1
        if args.verbose:
            print(f"{record['LocusId']}: motif {motifs[0]} contains IUPAC codes "
                  f"-- leave this locus alone")
        return None

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

    extended_start = start_0based - left_extension
    extended_end = polish_right_boundary(fasta_obj, normalized_chrom, extended_start,
                                         end_1based + right_extension, motifs[0])

    extended_record = dict(record)
    extended_record["ReferenceRegion"] = f"{chrom}:{extended_start}-{extended_end}"
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


def build_output_records(records_and_extensions, keep_original_definitions, counters,
                         only_extended_definitions=False):
    """Assemble the output catalog from each input record paired with its extended form or None.

    Args:
        records_and_extensions (list): (input record, extended record or None) pairs, in input order.
        keep_original_definitions (bool): keep the original definition of a locus that was extended,
            in addition to its extended definition, rather than replacing it.
        counters (collections.Counter): updated with what was added and skipped.
        only_extended_definitions (bool): write only the definitions the rule produced, leaving out
            every input locus. The result is not a catalog in its own right, it is the set of
            definitions this rule adds to the one it was given.

    Return:
        list: the output catalog records.
    """
    # Which of the input's own definitions are still around for an extended definition to collide
    # with. An extended definition that matches one of these would be a duplicate, so it is not added.
    # This has to be computed over the whole catalog up front, because the locus it collides with may
    # come later in the file. With only_extended_definitions none of them are written here, but the
    # output is meant to be appended to the catalog it came from, where all of them remain.
    surviving_records = [record for record, extended in records_and_extensions
                         if only_extended_definitions or extended is None or keep_original_definitions]
    kept_original_keys = {catalog_record_key(record) for record in surviving_records}
    taken_locus_ids = {record["LocusId"] for record in surviving_records}

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
        if not only_extended_definitions and (extended is None or keep_original_definitions):
            output_records.append(record)

        if extended is None:
            continue

        if catalog_record_key(extended) in kept_original_keys:
            counters["extended definitions not added because that definition was already in the catalog"] += 1
            continue
        if id(extended) not in kept_candidates:
            counters["extended definitions dropped in favor of an overlapping extended definition that was kept"] += 1
            if not only_extended_definitions and not keep_original_definitions \
                    and not is_covered_by_kept_definitions(record, kept_intervals_per_group):
                # The definitions that displaced this one only have to overlap it, so they don't
                # always reach over everything this locus covered. Put the locus back when they
                # don't, rather than let the catalog lose sequence it used to describe. Asking the
                # whole group, rather than only what displaced this one definition, is what keeps a
                # locus that some other kept definition already covers from coming back as a
                # duplicate.
                output_records.append(record)
            continue

        counters["extended definitions added to the output"] += 1
        # An extended definition is always named after the interval it now covers, whatever the
        # original was called. The original id described the locus before it grew, and the original
        # may still be in the output holding that id, so the extended definition takes a
        # chrom-start-end-motif name of its own.
        extended = dict(extended)
        extended["LocusId"] = compute_extended_locus_id(extended, taken_locus_ids)
        taken_locus_ids.add(extended["LocusId"])
        output_records.append(extended)

    return output_records


def sort_catalog_records(output_records):
    """Order the catalog the way its consumers require: each chromosome contiguous, keeping the order
    the chromosomes arrived in, and non-decreasing coordinates within each one.

    Sorting is not cosmetic here. An extension can move a locus's left boundary, so an extended
    definition can start before the locus it grew from, and appending it after that locus leaves the
    catalog locally out of order. str_analysis.utils.eh_catalog_utils.group_overlapping_loci, which
    the annotation steps downstream run on this output, raises on the first pair of records whose
    coordinates decrease.
    """
    chrom_order = {}
    for record in output_records:
        chrom, _, _ = parse_interval(record["ReferenceRegion"])
        if chrom not in chrom_order:
            chrom_order[chrom] = len(chrom_order)

    def sort_key(record):
        chrom, start_0based, end_1based = parse_interval(record["ReferenceRegion"])
        return (chrom_order[chrom], start_0based, end_1based, record["LocusStructure"],
                record["LocusId"])

    return sorted(output_records, key=sort_key)


def write_output_catalog(output_records, output_path):
    """Write the catalog, as gzipped JSON that keeps every field, or as bgzipped and tabixed BED rows
    that keep only the interval and the motif.
    """
    if is_json_path(output_path):
        export_json(sort_catalog_records(output_records), output_path)
        return

    uncompressed_path = re.sub(r"\.b?gz$", "", output_path)
    print(f"Writing {output_path}")
    with open(uncompressed_path, "wt") as output_bed:
        for bed_record in sorted(convert_json_records_to_bed_format_tuples(output_records)):
            output_bed.write("\t".join(map(str, bed_record)) + "\n")
    if os.system(f"bgzip -f {uncompressed_path}") != 0:
        raise RuntimeError(f"bgzip failed for {uncompressed_path}")
    if os.system(f"tabix -f -p bed {uncompressed_path}.gz") != 0:
        raise RuntimeError(f"tabix failed for {uncompressed_path}.gz")


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
                print(f"  -> {extended_record['ReferenceRegion']} "
                      f"{extended_record['LocusStructure']}")
        records_and_extensions.append((record, extended_record))

    records_and_extensions = discard_extensions_outdone_by_the_input(records_and_extensions, counters)
    output_records = build_output_records(
        records_and_extensions, args.keep_original_definitions_of_extended_loci, counters,
        only_extended_definitions=args.only_output_extended_loci)

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
          f"{counters['extended definitions added to the output']:,d} extended definitions reached "
          f"the output")
    if not args.only_output_extended_loci:
        # Taken from the finished catalog rather than accumulated per extension, since a locus can
        # also leave the output when its extension is dropped for an overlapping one, and a figure
        # that only ever added would report growth the catalog didn't get. It says nothing when the
        # input's own definitions were left out of the output.
        base_pairs_added = (sum(compute_span(record) for record in output_records)
                            - sum(compute_span(record) for record, _ in records_and_extensions))
        print(f"Locus definitions in the output span {base_pairs_added:,d}bp more in total than the "
              f"input's did")
    print(f"Wrote {len(output_records):,d} locus definitions to {args.output_path}")
    already_reported = ("total", "extended", "extended definitions added to the output")
    for label, count in counters.items():
        if label not in already_reported:
            print(f"{count:,d} {label}")


if __name__ == "__main__":
    main()
