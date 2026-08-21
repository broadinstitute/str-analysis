"""The gap-purity rule for extending a tandem repeat locus boundary outward into its flanking sequence.

The rule tolerates a short run of imperfect motif copies in the flank rather than requiring an
immediate exact match, as long as the whole newly added stretch of sequence stays highly pure. See
extend_str_catalog_boundaries.py for the pseudocode and the reasoning behind it.

Convention used throughout: the annotated repeat region is assumed to be phase-aligned such that
core_seq[0] == motif[0]. This matches how ExpansionHunter anchors reference-repeat purity in
RepeatPurity.cpp, which performs no rotation search. The motif phase at any position relative to the
core's start is therefore just `relative_position % motif_length`, and Python's modulo already
returns a non-negative result for the negative relative positions that left-flank bases need.
"""

import numpy as np

MAX_REPEATS_IN_GAP = 2
MIN_PURITY_OF_NEW_SEQUENCE = 0.9
MIN_REFERENCE_PURITY = 0.9

# How much flank to fetch on the first attempt, and how many times to widen it. Each attempt triples
# the window, so the last is 81x the first. The first window is motif-dependent (see
# compute_extension) and at least INITIAL_FLANK_WINDOW, which for a motif of 100bp or less gives
# 300bp, 900bp, 2700bp, 8100bp, 24300bp.
INITIAL_FLANK_WINDOW = 300
MAX_RESCAN_ATTEMPTS = 5

# IUPAC ambiguity codes: the set of literal bases each code represents. A motif position using one of
# these (RFC1's "AARRG" where R = A-or-G, RUNX2's "GCN" where N = any base) matches any base in its
# set rather than only its own letter, since otherwise a real, biologically correct repeat copy is
# counted as mismatching at every degenerate position. Note this deliberately departs from
# ExpansionHunter's own computeRepeatSequencePurity, which compares characters literally and so never
# matches an 'N' in the motif against a real base.
IUPAC_CODE_TO_BASES = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


def _build_iupac_match_table():
    """256x256 boolean lookup: table[sequence_base][motif_base] is True when the sequence base is a
    literal base consistent with the motif base's IUPAC code.
    """
    table = np.zeros((256, 256), dtype=bool)
    for motif_code, represented_bases in IUPAC_CODE_TO_BASES.items():
        for sequence_base in represented_bases:
            table[ord(sequence_base), ord(motif_code)] = True
    return table


_IUPAC_MATCH_TABLE = _build_iupac_match_table()


def _flank_phases(core_length, motif_length, side, max_offset):
    """Motif phase (index into the motif) for each of the first `max_offset` bases of the given flank.

    side='right': flank[k] sits at relative position core_length + k.
    side='left': flank[k], where k=0 is the base immediately left of the core, sits at -(k + 1).
    """
    offsets = np.arange(max_offset)
    if side == "right":
        relative_positions = core_length + offsets
    elif side == "left":
        relative_positions = -(offsets + 1)
    else:
        raise ValueError(f"side must be 'left' or 'right', got {side!r}")
    return relative_positions % motif_length


def cumulative_mismatched_bases(flank_seq, motif, core_length, side):
    """Hamming distance from the region [boundary, boundary + offset) to a perfect tiled repeat of the
    motif, for offset in 1..len(flank_seq).

    Args:
        flank_seq (str): flank sequence in forward genomic order for 'right', and reversed
            (boundary-adjacent base first) for 'left'.
        motif (str): the locus motif, which may contain IUPAC ambiguity codes.
        core_length (int): length of the annotated repeat region, which fixes the motif phase.
        side (str): 'left' or 'right'.

    Return:
        np.ndarray: index i holds the mismatch count for offset i+1.
    """
    max_offset = len(flank_seq)
    phases = _flank_phases(core_length, len(motif), side, max_offset)
    motif_array = np.frombuffer(motif.upper().encode(), dtype=np.uint8)
    flank_array = np.frombuffer(flank_seq.upper().encode(), dtype=np.uint8)
    is_match = _IUPAC_MATCH_TABLE[flank_array, motif_array[phases]]
    return np.arange(1, max_offset + 1) - np.cumsum(is_match)


def compute_reference_purity(sequence, motif):
    """The fraction of the sequence's bases that match a perfect tiling of the motif from its start.

    The same measure the gap-purity rule applies to sequence it is considering adding, turned on the
    locus's own existing definition.

    Args:
        sequence (str): the locus's reference sequence.
        motif (str): the locus motif, which may contain IUPAC ambiguity codes.

    Return:
        float: between 0 and 1, and 0 for an empty sequence.
    """
    if not sequence:
        return 0.0
    # A repeat region tiled from its own first base has phase k % motif_length at offset k, which is
    # what cumulative_mismatched_bases works out for a right flank following a zero-length region.
    return 1 - int(cumulative_mismatched_bases(sequence, motif, 0, "right")[-1]) / len(sequence)


def extend_by_gap_purity(flank_seq, motif, core_length, side,
                         max_repeats_in_gap=MAX_REPEATS_IN_GAP,
                         min_purity_of_new_sequence=MIN_PURITY_OF_NEW_SEQUENCE,
                         verbose=False):
    """Return how far this one flank should be extended, in base pairs.

    Args:
        flank_seq (str): flank sequence in forward genomic order for 'right', and reversed
            (boundary-adjacent base first) for 'left'.
        motif (str): the locus motif.
        core_length (int): length of the annotated repeat region.
        side (str): 'left' or 'right'.
        max_repeats_in_gap (int): how many consecutive non-exact motif copies may separate the current
            boundary from the next exact copy.
        min_purity_of_new_sequence (float): the fraction of bases in the newly added sequence that must
            match a perfect tiling of the motif for the extension to be accepted.
        verbose (bool): print each accept/reject decision and the purity behind it.

    Return:
        int: the accepted extension length in base pairs, which is 0 when nothing is accepted.
    """
    motif = motif.upper()
    flank_seq = flank_seq.upper()
    motif_length = len(motif)
    max_offset = len(flank_seq)

    mismatched = cumulative_mismatched_bases(flank_seq, motif, core_length, side)

    def mismatches_in(start, end):
        if start >= end:
            return 0
        return int(mismatched[end - 1]) - (int(mismatched[start - 1]) if start > 0 else 0)

    def is_exact_copy(copy_start):
        copy_end = copy_start + motif_length
        return copy_end <= max_offset and mismatches_in(copy_start, copy_end) == 0

    accepted = 0
    while True:
        position = accepted
        copies_in_gap = 0
        anchor_start = None
        while copies_in_gap <= max_repeats_in_gap:
            if position + motif_length > max_offset:
                break
            if is_exact_copy(position):
                anchor_start = position
                break
            position += motif_length
            copies_in_gap += 1

        if anchor_start is None:
            if verbose:
                print(f"    no exact copy within {max_repeats_in_gap} imperfect one(s) -- stop at {accepted}")
            return accepted

        run_end = anchor_start
        while is_exact_copy(run_end):
            run_end += motif_length

        added_bases = run_end - accepted
        purity = 1 - mismatches_in(accepted, run_end) / added_bases
        if purity < min_purity_of_new_sequence:
            if verbose:
                print(f"    reject: newly-added={flank_seq[accepted:run_end]!r} ({added_bases}bp), "
                      f"purity={purity:.3f} < {min_purity_of_new_sequence} -- stop at {accepted}")
            return accepted

        if verbose:
            print(f"    accept: newly-added={flank_seq[accepted:run_end]!r} ({added_bases}bp = "
                  f"{copies_in_gap} imperfect copy/copies + "
                  f"{(run_end - anchor_start) // motif_length} exact), "
                  f"purity={purity:.3f} >= {min_purity_of_new_sequence} -- extend to {run_end}")
        accepted = run_end


def extend_to_capture_partial_copy(flank_seq, motif, core_length, side):
    """Return how many of the flank's boundary-adjacent bases continue the motif's tiling exactly.

    The gap-purity rule only ever moves a boundary by whole motif copies, so it can leave a run of
    bases that do continue the repeat, but are fewer than one copy, sitting just outside the locus.
    The EP400 locus is an example: its boundary lands one base short of the 'G' that starts the next
    copy. This picks those bases up.

    Args:
        flank_seq (str): flank sequence starting at the boundary, in forward genomic order for
            'right', and reversed (boundary-adjacent base first) for 'left'.
        motif (str): the locus motif, which may contain IUPAC ambiguity codes.
        core_length (int): length of the repeat region whose phase the flank continues.
        side (str): 'left' or 'right'.

    Return:
        int: how many bases to add, between 0 and len(motif) - 1.
    """
    # A whole copy's worth of matching bases is a full exact copy, which the gap-purity rule would
    # have taken already, so stopping one base short of that is not a limit this can run into.
    max_offset = min(len(flank_seq), len(motif) - 1)
    if max_offset <= 0:
        return 0

    motif = motif.upper()
    # An IUPAC code matches whatever base happens to sit under it. That is what keeps a real copy of
    # an RFC1 or RUNX2-style motif from being scored as a mismatch at every degenerate position, but
    # here it would pull a base into the locus on no evidence at all, so the walk stops before one.
    phases = _flank_phases(core_length, len(motif), side, max_offset)
    for offset, phase in enumerate(phases):
        if motif[phase] not in "ACGT":
            max_offset = offset
            break
    if max_offset == 0:
        return 0

    mismatched = cumulative_mismatched_bases(flank_seq[:max_offset], motif, core_length, side)
    # mismatched is cumulative, so its first non-zero entry is the first base that doesn't match, and
    # its index is the number of matching bases before it.
    first_mismatch = np.flatnonzero(mismatched)
    return int(first_mismatch[0]) if len(first_mismatch) > 0 else max_offset


def compute_extension(fasta_obj, chrom, start_0based, end_1based, motif, side,
                      max_repeats_in_gap=MAX_REPEATS_IN_GAP,
                      min_purity_of_new_sequence=MIN_PURITY_OF_NEW_SEQUENCE,
                      verbose=False):
    """Extend one side of a locus, widening the fetched flank until the answer stops changing.

    extend_by_gap_purity only ever sees the flank it was handed, so a window that is too small can cut
    short either the search for the next exact copy or the run of exact copies after it, in either case
    returning a value that looks like a genuine rule-determined stop but is not. The only way to tell
    the two apart from outside that function is to check whether a bigger window changes the answer.

    Args:
        fasta_obj (pysam.FastaFile): the reference genome.
        chrom (str): chromosome name, already matching the reference's naming.
        start_0based (int): the locus start.
        end_1based (int): the locus end.
        motif (str): the locus motif.
        side (str): 'left' or 'right'.
        max_repeats_in_gap (int): see extend_by_gap_purity.
        min_purity_of_new_sequence (float): see extend_by_gap_purity.
        verbose (bool): print each accept/reject decision and the purity behind it. Only the final,
            widest window's decisions are printed, since the earlier attempts are provisional.

    Return:
        tuple: (extension in base pairs, whether the widest window was still changing the answer).
    """
    core_length = end_1based - start_0based
    chrom_size = fasta_obj.get_reference_length(chrom)
    # The first window covers the search's full reach, (max_repeats_in_gap + 1) motif copies. A fixed
    # 300bp start is smaller than that for motifs over 100bp, and two windows too small to finish the
    # search can agree on a value neither of them was in a position to determine.
    window = max(INITIAL_FLANK_WINDOW, (max_repeats_in_gap + 1) * len(motif))
    previous_extension = None
    for _ in range(MAX_RESCAN_ATTEMPTS):
        if side == "left":
            flank = fasta_obj.fetch(chrom, max(0, start_0based - window), start_0based)[::-1]
        else:
            flank = fasta_obj.fetch(chrom, end_1based, min(chrom_size, end_1based + window))
        extension = extend_by_gap_purity(flank, motif, core_length, side,
                                         max_repeats_in_gap, min_purity_of_new_sequence)
        if extension == previous_extension:
            if verbose:
                # Re-run the converged window so the tracing shows the decisions behind the answer that
                # was actually returned, rather than those of every discarded attempt along the way.
                extend_by_gap_purity(flank, motif, core_length, side, max_repeats_in_gap,
                                     min_purity_of_new_sequence, verbose=True)
            return extension, False
        previous_extension = extension
        window *= 3
    return previous_extension, True
