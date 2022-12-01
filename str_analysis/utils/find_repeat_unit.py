import os
import re
from str_analysis.utils.trf_runner import TRFRunner

ALL_ALLOWED_BASES = {"A", "C", "G", "T"}
COMPILED_REGEX_CACHE = {}

DEFAULT_MIN_INTERRUPTED_REPEAT_UNIT_LENGTH = 3
DEFAULT_MAX_INTERRUPTED_REPEAT_UNIT_LENGTH = 24


def get_repeat_unit_regex_with_N_base(repeat_unit, i, allow_homopolymer=True):
    """Generate a regular expression (regex) that recognizes the motif, with the base at position i treated as an N (ie.
    any nucleotide is allowed at that position). If allow_homopolymer is False, any base is allowed at the N position
    as long as it doesn't turn the motif into a homopolymer. For example, the motif "AAT" with i=0 would yield the regex
    "[ACGT]AT", but "AAT" with i=2 and allow_homopolymer=False would yield the regex "AA[CGT]" (keep "A" out of the
    allowed nucleotides in the last position).

    Args:
        repeat_unit (str): A nucleotide sequence
        i (int): The position where to place the "N" base. Should be between 0 and len(motif).
        allow_homopolymer (bool): If False, base substitutions that convert the entire motif to a homopolymer will not
            be allowed by the returned regex (see example in the docs above).

    Return:
        str: regular expression
    """
    if i < 0 or i >= len(repeat_unit):
        raise ValueError(f"Argument i is out of bounds: {i} vs motif length: {len(repeat_unit)}")

    left_subsequence = repeat_unit[:i]
    right_subsequence = repeat_unit[i+1:]
    allowed_bases = "[ACGT]"

    if not allow_homopolymer:
        existing_bases = set(left_subsequence) | set(right_subsequence)
        if len(existing_bases) == 1:
            allowed_bases = "[" + "".join(ALL_ALLOWED_BASES - existing_bases) + "]"

    return left_subsequence + allowed_bases + right_subsequence


def find_repeat_unit(
    sequence,
    allow_interruptions=False,
    min_interrupted_repeat_unit_length=DEFAULT_MIN_INTERRUPTED_REPEAT_UNIT_LENGTH,
    max_interrupted_repeat_unit_length=DEFAULT_MAX_INTERRUPTED_REPEAT_UNIT_LENGTH,
    allow_partial_repeats=False,
):
    """Check whether the given nucleotide sequence is a continuous set of repeats made up of some smaller
    repeat unit (eg. CAGCAGCAGCAGCAG = 5xCAG). If no such repeat unit is found, return the input sequence itself.
    This method first checks whether the input sequence consists entirely of perfect repeats (eg. CAGCAGCAG).
    If not, the method checks for repeats where 1 position can vary (eg. CAGCAACAGCAA
    would be counted as a CAG repeat with 2 interruptions). Interruptions are allowed in motifs between length 3bp and
    max_interrupted_repeat_unit_length.

    Args:
        sequence (str): nucleotide sequence
        allow_interruptions (bool): Whether to allow interruptions in the repeat sequence.
        min_interrupted_repeat_unit_length (int): Interruptions will not be allowed in motifs shorter than this threshold.
            This avoids degenerate repeat sequences with messy interrupted  homopolymer or dinucleotide motifs.
        max_interrupted_repeat_unit_length (int): Interruptions will not be allowed in motifs longer than this threshold.
            Allowing for interruptions increases compute time proportional to motif length, so this threshold is
            necessary to prevent the method from running too slowly.
        allow_partial_repeats (bool): If False, partial repeats will not be allowed. For example, an input sequence of
            "CAGCAGCA" would not be considered as made up of "CAG" repeats since the last repeat is incomplete. If True,
            "CAGCAGCA" would be counted as consisting of 2 "CAG" repeats.
    Return:
         4-tuple (str, int, int, int): repeat unit, number of pure repeats, num total repeats,
            which position in the repeat unit varies across
            repeats or None if the input sequence was found to consist entirely of pure repeats. For example, a return
            value of (CAG, 7, 10, 0) means the input sequence consists of 7 pure CAGs and 3 CAGs with an interruption at
            position 0 of the repeat unit (whould could be "AAG", "GAG" or "TAG").
    """

    # find the smallest repeat unit that covers the entire motif
    repeat_unit_length = 1
    while repeat_unit_length <= len(sequence)/2:
        repeat_unit = sequence[:repeat_unit_length]
        num_repeats = sequence.count(repeat_unit)
        partial_repeat_length = len(sequence) % repeat_unit_length
        if num_repeats * repeat_unit_length == len(sequence) or (
                allow_partial_repeats
                and num_repeats * repeat_unit_length == len(sequence) - partial_repeat_length
                and sequence[num_repeats * repeat_unit_length:] == repeat_unit[:partial_repeat_length]
        ):
            return repeat_unit, num_repeats, num_repeats, None

        repeat_unit_length += 1

    # pure repeat not found, so check for repeats with interruptions
    if allow_interruptions:
        repeat_unit_length = min_interrupted_repeat_unit_length
        while repeat_unit_length <= len(sequence)/2:
            repeat_unit = sequence[:repeat_unit_length]

            repeat_unit_length += 1
            if repeat_unit_length > max_interrupted_repeat_unit_length:
                # don't check regexp for very large motifs since this takes too long
                break

            if len(set(repeat_unit)) == 1:
                # skip homopolymer repeat units
                continue

            num_total_repeats_expected = int(len(sequence) / len(repeat_unit))
            for i in range(len(repeat_unit)):
                regex = get_repeat_unit_regex_with_N_base(repeat_unit, i, allow_homopolymer=False)
                regex = "^(" + regex + "){%d}" % num_total_repeats_expected
                if allow_partial_repeats:
                    partial_repeat_length = len(sequence) % len(repeat_unit)
                    regex += repeat_unit[:partial_repeat_length]
                regex += "$"

                if regex not in COMPILED_REGEX_CACHE:
                    COMPILED_REGEX_CACHE[regex] = re.compile(regex)

                if COMPILED_REGEX_CACHE[regex].match(sequence):
                    num_pure_repeats = sequence.count(repeat_unit)
                    return repeat_unit, num_pure_repeats, num_total_repeats_expected, i


    # no repeat unit found, so return the input sequence itself as the repeat unit
    repeat_unit = sequence
    num_pure_repeats = num_total_repeats_expected = 1
    return repeat_unit, num_pure_repeats, num_total_repeats_expected, None


def extend_repeat_into_sequence(
        repeat_unit, 
        sequence,
        allow_interruptions=False,
        min_interrupted_repeat_unit_length=DEFAULT_MIN_INTERRUPTED_REPEAT_UNIT_LENGTH,
        max_interrupted_repeat_unit_length=DEFAULT_MAX_INTERRUPTED_REPEAT_UNIT_LENGTH,
        repeat_unit_interruption_index=None,
):
    """This method walks along the given sequence from left to right, one repeat unit length at a time
    (defined by the given repeat unit), and checks for the longest stretch of repeats that exactly
    matches the given repeat_unit.

    Args:
        repeat_unit (str): For example "CAG"
        sequence (str): a longer sequence that may contain repeats of the given repeat unit starting at the left end,
            before switching to random other sequence or repeats of a different repeat unit.
            For example: "CAGCAGCAGCTAGTGCAGTGACAGT"
        allow_interruptions (bool): Whether to allow interruptions in the repeat sequence.
        min_interrupted_repeat_unit_length (int): Interruptions will not be allowed in motifs shorter than this threshold.
            This avoids degenerate repeat sequences with messy interrupted  homopolymer or dinucleotide motifs.
        max_interrupted_repeat_unit_length (int): Interruptions will not be allowed in motifs longer than this threshold.
            Allowing for interruptions increases compute time proportional to motif length, so this threshold is
            necessary to avoid very slow compute times.
        repeat_unit_interruption_index (int): This should be an integer between 0 and len(repeat_unit) or None.
            If not None, this position in the motif will be treated as N (eg. allowed to be any nucleotide A,C,G, or T)
            If None, the sequence will be extended allowing for up to 1 base in the motif to be treated as N (eg.
            allowed to be any nucleotide A,C,G, or T). Even with interruptions though, the repeat sequence must end in
            at least one pure copy of the given repeat_unit (eg. the "CAT" interruption will be included if
            the extended sequence is "CAGCAGCATCAG" , but not if it's just "CAGCAGCAT")

    Return:
        3-tuple: (int, int, int) num_pure_repeats, num_total_repeats, repeat_unit_interruption_index
    """

    # compute the number of pure repeats
    i = 0
    num_pure_repeats = 0
    while i <= len(sequence) - len(repeat_unit):
        if sequence[i:i+len(repeat_unit)] != repeat_unit:
            break
        num_pure_repeats += 1
        i += len(repeat_unit)

    num_total_repeats = num_pure_repeats
    
    # extend allowing for interruptions
    if allow_interruptions:
        if repeat_unit_interruption_index is None:
            if len(repeat_unit) < min_interrupted_repeat_unit_length or len(repeat_unit) > max_interrupted_repeat_unit_length:
                return num_pure_repeats, num_pure_repeats, repeat_unit_interruption_index
    
            repeat_unit_interruption_indices_list = range(len(repeat_unit))
    
        else:
            repeat_unit_interruption_indices_list = [repeat_unit_interruption_index]
    
        for current_index in repeat_unit_interruption_indices_list:
            regex = get_repeat_unit_regex_with_N_base(repeat_unit, current_index, allow_homopolymer=False)
            regex = f"^((?:{regex})*{repeat_unit})(.*)"
            if regex not in COMPILED_REGEX_CACHE:
                COMPILED_REGEX_CACHE[regex] = re.compile(regex)
    
            match = COMPILED_REGEX_CACHE[regex].match(sequence)
            if match:
                matching_impure_sequence = match.group(1)
                current_num_total_repeats = int(len(matching_impure_sequence) / len(repeat_unit))
                if current_num_total_repeats > num_total_repeats:
                    num_total_repeats = current_num_total_repeats
                    if repeat_unit_interruption_index is None:
                        repeat_unit_interruption_index = current_index

                rest_of_sequence = match.group(2)
                if len(rest_of_sequence) == 0:
                    print(f"WARNING: extend_repeat_into_sequence({repeat_unit}) reached the end of the "
                          f"{len(sequence)}bp sequence: {num_total_repeats}x{repeat_unit} ==> {sequence}")

    return num_pure_repeats, num_total_repeats, repeat_unit_interruption_index


def find_repeat_unit_kmer(sequence, min_fraction_covered_by_repeat=1, allow_indel_interruptions=False):
    """Check whether the given nucleotide sequence is a contiguous set of repeats made up of some smaller
    repeat unit (eg. CAGCAGCAGCAGCAG = 5xCAG). If no such repeat unit is found, return the input sequence itself.
    Partial repeats are not recognized: an input sequence of "CAGCAGCA" would not be conrepeat_unitsidered as made up of "CAG"
    repeats since the last repeat is incomplete. Although this is very strict, it makes biological sense.

    Args:
        sequence (str): nucleotide sequence
        min_fraction_covered_by_repeat (float): Repeat units must cover at least this fraction of the input sequence.
            For example, if this parameter is set to 0.9, an input sequence consisting of 9 "CAG" repeats and 1 "CAT"
            would yield a return value of ("CAG", 10), but an input sequence of 8 "CAG"s and 1 "CAT" would have a
            return value of ("CAGCAGCAGCAGCAGCAGCAGCAGCAT", 1) since no smaller repeat unit satisfies the 0.9 threshold.
        allow_indel_interruptions (bool): Whether the repeat unit can be some
    Return:
         2-tuple (str, int): repeat unit, number of repeats
    """
    if len(sequence) <= 1:
        raise ValueError(f"len(s) must be > 1: {sequence}")

    result_repeat_unit = sequence
    result_num_repeats = 1
    max_fraction_covered_so_far = 0

    # find the smallest repeat unit that covers the largest fraction of the input sequence
    repeat_unit_length = 1
    while repeat_unit_length <= len(sequence)/2:
        if not allow_indel_interruptions and len(sequence) % repeat_unit_length != 0:
            repeat_unit_length += 1
            continue

        repeat_unit = sequence[:repeat_unit_length]
        num_repeats = sequence.count(repeat_unit)
        fraction_covered = num_repeats * repeat_unit_length / len(sequence)
        if fraction_covered >= min_fraction_covered_by_repeat and fraction_covered > max_fraction_covered_so_far:
            max_fraction_covered_so_far = fraction_covered
            result_repeat_unit = repeat_unit
            result_num_repeats = int(len(sequence) / len(repeat_unit))

        repeat_unit_length += 1

    return result_repeat_unit, result_num_repeats


def find_repeat_unit_using_trf(
        sequence, min_fraction_covered_by_repeat=0.9, trf_path="~/bin/trf409.macosx",
):
    """Run TandemRepeatFinder [Benson 1999] to find the repeat unit allowing for mismatches before falling back on the
    simpler find_repeat_unit(..) method.

    Args:
        sequence (str): The nucleotide sequence to search for repeats
        min_fraction_covered_by_repeat (float): A detected repeat is discarded unless it covers this fraction of s.
        trf_path (str): Path of TRF executable.

    Return:
        2-tuple: repeat unit, number of repeats found in the sequence. If no repeats are found, the input sequence is
        returned without modification, along with repeat count = 1.
    """
    if len(sequence) <= 1:
        raise ValueError(f"len(s) must be > 1: {sequence}")

    # if the simple approach didn't find repeats within the sequence, try using TRF to find imperfect repeats
    trf_runner = TRFRunner(
        os.path.expanduser(trf_path),
        match_score=2,
        mismatch_penalty=7,
        indel_penalty=10**6,
        minscore=len(sequence) * 2 * min_fraction_covered_by_repeat,
    )

    dat_records = list(trf_runner.run_TRF(sequence))

    for dat_record in dat_records:
        if len(sequence) % len(dat_record.repeat_unit) == 0:
            repeat_count = int(len(sequence) / len(dat_record.repeat_unit))
            repeat_unit = dat_record.repeat_unit
            assert repeat_count * len(repeat_unit) / len(sequence) >= min_fraction_covered_by_repeat

            return repeat_unit, repeat_count

    return sequence, 1

