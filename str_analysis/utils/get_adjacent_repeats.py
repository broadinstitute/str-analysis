import argparse
import intervaltree

from str_analysis.utils.misc_utils import parse_interval

# allow up to 20bp of reference sequence between adjacent repeats
MAX_DISTANCE_BETWEEN_REPEATS = 6

# look for adjacent repeats that are no more than 1000bp to the left or right of the main repeat.
MAX_TOTAL_ADJACENT_REGION_SIZE = 1000

# TandemRepeatFinder can output repeat stretch that overlap by several base pairs - for example:
# ATXN8OS repeats db records are chr13|70139352|70139385|TAC|..|ATXN8OS and chr13|70139384|70139429|CTG|..|ATXN8OS
# Allow overlap for this overlap.
MAX_OVERLAP_BETWEEN_ADJACENT_REPEATS = 3


def get_repeat_unit_from_fasta(chrom, start_1based, end_1based, repeat_unit_length, pysam_fasta_file):
    """Get the reference sequence from the genome fasta file at the given coordinates."""

    ref_fasta_sequence = pysam_fasta_file.fetch(chrom, start_1based - 1, end_1based)  # fetch uses 0-based coords
    repeat_unit = ref_fasta_sequence[:repeat_unit_length]  # get repeat unit from reference that matches these coords
    repeat_unit = repeat_unit.upper()  # convert to upper-case since the reference contains some lower-case regions

    return repeat_unit


def compute_locus_id(chrom, start_1based, end_1based, repeat_unit):
    return f"{chrom}:{start_1based}-{end_1based}-{repeat_unit}"


def get_adjacent_repeats(locus_interval_0based, repeat_unit, pysam_fasta_file, interval_tree_0based,
    max_distance_between_adjacent_repeats=MAX_DISTANCE_BETWEEN_REPEATS,
    max_total_adjacent_region_size=MAX_TOTAL_ADJACENT_REGION_SIZE,
    max_overlap_between_adjacent_repeats=MAX_OVERLAP_BETWEEN_ADJACENT_REPEATS,
    max_adjacent_repeats=None,
):
    """Find adjacent repeats on the left and right of the given locus interval.

    Args:
        locus_interval_0based (str): locus interval in 0-based format, e.g. "chr1:100-200"
        repeat_unit (str): repeat unit
        pysam_fasta_file (pysam.FastaFile): pysam FastaFile object
        interval_tree (intervaltree.IntervalTree): an IntervalTree containing 0-based intervals that represent the
            reference start and end coordinates of all TR loci found on the same chromosome as the given record.
            Adjacent repeats will be queried from this data structure.
        max_distance_between_adjacent_repeats (int): maximum distance between adjacent repeats
        max_total_adjacent_region_size (int): maximum total size of the region to search for adjacent repeats
        max_overlap_between_adjacent_repeats (int): maximum overlap between adjacent repeats
        max_adjacent_repeats (int): maximum number of adjacent repeats to return

    Returns:
         3-tuple: (adjacent_repeats_on_left, locus_structure, adjacent_repeats_on_right)
    """
    chrom, start_0based, end_1based = parse_interval(locus_interval_0based)

    repeat_units_already_added = set()

    repeat_unit = get_repeat_unit_from_fasta(chrom, start_0based + 1, end_1based, len(repeat_unit), pysam_fasta_file)
    locus_structure = f"({repeat_unit})*"

    repeat_units_already_added.add(repeat_unit)

    overlapping_intervals_left = interval_tree_0based.overlap(
        intervaltree.Interval(
            start_0based - max_total_adjacent_region_size,
            start_0based + max_overlap_between_adjacent_repeats))
    overlapping_intervals_left = sorted(overlapping_intervals_left, key=lambda i: i.end, reverse=True)

    adjacent_repeats_on_left = []
    current_left_coord_1based = start_0based + 1
    for adj_repeat_0based in overlapping_intervals_left:
        if current_left_coord_1based - adj_repeat_0based.begin < 2*len(adj_repeat_0based.data):
            # skip the overlapping interval if it's (almost) the same as the previous interval to the right of it
            continue

        adj_repeat = argparse.Namespace(
            start_1based=adj_repeat_0based.begin + 1,
            end_1based=adj_repeat_0based.end,
            repeat_unit=adj_repeat_0based.data)

        if adj_repeat.end_1based + max_distance_between_adjacent_repeats + 1 < current_left_coord_1based:
            break

        if adj_repeat.end_1based >= current_left_coord_1based:
            adj_repeat.end_1based = current_left_coord_1based - 1

        # If the adjacent repeat length is not an exact multiple of the repeat unit length, trim extra bases
        # since ExpansionHunter works slightly better with repeats that are exact multiples of the repeat unit
        adj_repeat.start_1based += (adj_repeat.end_1based - adj_repeat.start_1based + 1) % len(adj_repeat.repeat_unit)

        if max_adjacent_repeats and len(adjacent_repeats_on_left) + 1 >= max_adjacent_repeats:
            break

        # record the ref spacer (if one exists) between the next repeat on the left, and the current repeat
        # (which, on the 1st iteration, is the main repeat)
        ref_spacer = ""
        if adj_repeat.end_1based < current_left_coord_1based - 1:
            ref_spacer = pysam_fasta_file.fetch(chrom, adj_repeat.end_1based, current_left_coord_1based - 1)
            ref_spacer = ref_spacer.upper() # convert to upper-case since the reference contains some lower-case regions

        # record the next repeat to the left of the current repeat (which, on the 1st iteration, is the main repeat)
        repeat_unit = get_repeat_unit_from_fasta(
            chrom, adj_repeat.start_1based, adj_repeat.end_1based, len(adj_repeat.repeat_unit), pysam_fasta_file)
        if repeat_unit in repeat_units_already_added:
            break  # ExpansionHunter can't handle the same repeat unit being specified more than once
        repeat_units_already_added.add(repeat_unit)

        adj_repeat_left_string = compute_locus_id(
            chrom, adj_repeat.start_1based - 1, adj_repeat.end_1based, repeat_unit)
        adjacent_repeats_on_left.append(adj_repeat_left_string)

        # move current_left_coord to the start_1based position of this next repeat on the left.
        current_left_coord_1based = adj_repeat.start_1based
        locus_structure = f"({repeat_unit})*{ref_spacer}{locus_structure}"

    adjacent_repeats_on_left.reverse()  # Reverse adjacent-left repeats to put them in genomic coordinate order

    overlapping_intervals_right = interval_tree_0based.overlap(
        intervaltree.Interval(
            end_1based - max_overlap_between_adjacent_repeats,
            end_1based + max_total_adjacent_region_size))

    # sort overlapping intervals?
    adjacent_repeats_on_right = []
    current_right_coord_1based = end_1based
    for adj_repeat_0based in sorted(overlapping_intervals_right, key=lambda i: i.begin):
        if adj_repeat_0based.end - current_right_coord_1based < 2*len(adj_repeat_0based.data):
            # skip the overlapping interval if it's (almost) the same as the previous interval to the left of it
            continue

        adj_repeat = argparse.Namespace(
            start_1based=adj_repeat_0based.begin + 1,
            end_1based=adj_repeat_0based.end,
            repeat_unit=adj_repeat_0based.data)

        if adj_repeat.start_1based - max_distance_between_adjacent_repeats - 1 > current_right_coord_1based:
            break

        if adj_repeat.start_1based <= current_right_coord_1based:
            adj_repeat.start_1based = current_right_coord_1based + 1

        # If the adjacent repeat length is not an exact multiple of the repeat unit length, trim extra bases
        # since ExpansionHunter works slightly better with repeats that are exact multiples of the repeat unit
        adj_repeat.end_1based -= (adj_repeat.end_1based - adj_repeat.start_1based + 1) % len(adj_repeat.repeat_unit)

        if max_adjacent_repeats and len(adjacent_repeats_on_right) + 1 >= max_adjacent_repeats:
            break

        # record the ref spacer (if one exists) between the next repeat on the right, and the current repeat
        # (which, on the 1st iteration, is the main repeat)
        ref_spacer = ""
        if adj_repeat.start_1based > current_right_coord_1based + 1:
            ref_spacer = pysam_fasta_file.fetch(chrom, current_right_coord_1based, adj_repeat.start_1based - 1)
            ref_spacer = ref_spacer.upper() # convert to upper-case since the reference contains some lower-case regions

        # record the next repeat to the right of the current repeat (which, on the 1st iteration, is the main repeat)
        repeat_unit = get_repeat_unit_from_fasta(
            chrom, adj_repeat.start_1based, adj_repeat.end_1based, len(adj_repeat.repeat_unit), pysam_fasta_file)
        if repeat_unit in repeat_units_already_added:
            break # ExpansionHunter can't handle the same repeat unit being specified more than once
        repeat_units_already_added.add(repeat_unit)

        adj_repeat_right_string = compute_locus_id(
            chrom, adj_repeat.start_1based - 1, adj_repeat.end_1based, repeat_unit)
        adjacent_repeats_on_right.append(adj_repeat_right_string)

        current_right_coord_1based = adj_repeat.end_1based
        locus_structure = f"{locus_structure}{ref_spacer}({repeat_unit})*"

    return adjacent_repeats_on_left, locus_structure, adjacent_repeats_on_right

