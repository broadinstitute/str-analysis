import collections
import os
import re
from str_analysis.utils.trf_runner import TRFRunner

"""Interruptions will not be allowed in repeat units shorter than this threshold. Allowing for interruptions increases 
compute time proportional to repeat unit length, so this threshold is necessary to reduce overall compute time.
"""
MAX_INTERRUPTED_REPEAT_UNIT_LENGTH = 24

"""The maximum fraction of bases that can be interruptions before it's no longer considered a homopolymer repeat."""
MIN_PURE_REPEAT_FRACTION_FOR_HOMOPOLYMER_REPEATS = 9/10

"""The maximum fraction of repeats that can have interruptions before it's no longer considered a dinucleotide repeat."""
MIN_PURE_REPEAT_FRACTION_FOR_DINUCLEOTIDE_REPEATS = 4/5

"""Cache of compiled regular expressions for matching repeats."""
COMPILED_REGEX_CACHE = {}


def get_repeat_unit_regex_with_N_base(repeat_unit, i, allow_homopolymer=True):
    """Generate a regular expression (regex) that recognizes the repeat unit, with the base at position i set to N (ie.
    any nucleotide is allowed at that position). If allow_homopolymer is False, any base is allowed at the N position
    as long as it doesn't turn the repeat unit into a homopolymer. For example, the repeat unit "AAT" with i=0 would
    yield the regex "[ACGT]AT", but "AAT" with i=2 and allow_homopolymer=False would yield the regex "AA[CGT]" (keep "A"
    out of the allowed nucleotides in the last position).

    Args:
        repeat_unit (str): A nucleotide sequence
        i (int): The position where to place the "N" base. Should be between 0 and len(repeat_unit).
        allow_homopolymer (bool): If False, base substitutions that convert the entire repeat unit to a homopolymer
            will not be allowed by the returned regex (see example in the docs above).
    Return:
        str: regular expression
    """
    if i < 0 or i >= len(repeat_unit):
        raise ValueError(f"Argument i is out of bounds: {i} vs repeat unit length: {len(repeat_unit)}")

    if len(repeat_unit) == 1:
        raise ValueError("Repeat unit length must be > 1")

    left_subsequence = repeat_unit[:i]
    right_subsequence = repeat_unit[i+1:]
    allowed_bases = "[ACGT]"

    if not allow_homopolymer:
        existing_bases = set(left_subsequence) | set(right_subsequence)
        if len(existing_bases) == 1:
            allowed_bases = "[" + "".join([b for b in "ACGT" if b not in existing_bases]) + "]"

    repeat_unit_regex = left_subsequence + allowed_bases + right_subsequence

    return repeat_unit_regex


def get_most_common_repeat_unit(sequence, repeat_unit_length):
    """Find the most common repeat unit of the given length in the given sequence. For example, for
    "CAGCAGCAT" and repeat unit length 3, the most common repeat unit is "CAG" (2 repeats).

    Args:
        sequence (str): nucleotide sequence
        repeat_unit_length (int): length of the repeat unit to search for

    Return:
        str: the most common repeat unit of the given length in the given sequence
    """
    repeat_unit_counter = collections.defaultdict(int)
    most_common_repeat_unit = None
    most_common_repeat_unit_count = 0
    for i in range(0, len(sequence), repeat_unit_length):
        repeat_unit = sequence[i:i+repeat_unit_length]
        repeat_unit_counter[repeat_unit] += 1
        if repeat_unit_counter[repeat_unit] > most_common_repeat_unit_count:
            most_common_repeat_unit = repeat_unit
            most_common_repeat_unit_count = repeat_unit_counter[repeat_unit]

    return most_common_repeat_unit, most_common_repeat_unit_count


def find_repeat_unit_without_allowing_interruptions(sequence, allow_partial_repeats=False):
    """Check whether the given nucleotide sequence entirely consists of repeats of some smaller repeat unit
    (eg. CAGCAGCAGCAGCAG = 5xCAG). If no such repeat unit is found, return the input sequence itself.

    Args:
        sequence (str): nucleotide sequence
        allow_partial_repeats (bool): If False, partial repeats will not be allowed. For example, an input sequence of
            "CAGCAGCA" would not be counted as having repeat unit "CAG" repeats because the last repeat is incomplete.
            If True, "CAGCAGCA" would be counted as consisting of 2 "CAG" repeats.

    Return:
         3-tuple (str, int, bool):
            repeat unit,
            number of pure repeats,
            has_partial_repeats
    """

    # find the smallest repeat unit that covers the entire sequence
    repeat_unit_length = 1
    while repeat_unit_length <= len(sequence)/2:
        if not allow_partial_repeats and len(sequence) % repeat_unit_length != 0:
            repeat_unit_length += 1
            continue

        repeat_unit = sequence[:repeat_unit_length]
        num_repeats = sequence.count(repeat_unit)
        if num_repeats * repeat_unit_length == len(sequence):
            return repeat_unit, num_repeats, False

        if allow_partial_repeats:
            partial_repeat_length = len(sequence) % repeat_unit_length
            if num_repeats * repeat_unit_length == len(sequence) - partial_repeat_length \
                    and sequence[num_repeats * repeat_unit_length:] == repeat_unit[:partial_repeat_length]:
                return repeat_unit, num_repeats, True

        repeat_unit_length += 1

    # no repeat unit found, so return the input sequence itself as the repeat unit
    return sequence, 1, False


def find_repeat_unit_allowing_interruptions(sequence, allow_partial_repeats=False):
    """Check whether the given nucleotide sequence entirely consists of repeats of some smaller repeat unit
    (eg. CAGCAGCAGCAGCAG = 5xCAG), allowing 1 position to vary across repeats (eg. CAG.CAA.CAG.CAA would be counted as
    a CAG repeat with 2 interruptions). If no such repeat unit is found, return the input sequence itself.

    Args:
        sequence (str): nucleotide sequence
        allow_partial_repeats (bool): If False, partial repeats will not be allowed. For example, an input sequence of
            "CAGCAGCA" would not be counted as having repeat unit "CAG" repeats because the last repeat is incomplete.
            If True, "CAGCAGCA" would be counted as consisting of 2 "CAG" repeats.

    Return:
         5-tuple (str, int, int, int, bool):
            repeat unit,
            number of pure repeats,
            num total repeats,
            which position in the repeat unit varies across repeats or None if the input sequence was found to
                consist entirely of pure repeats. For example, a return value of (CAG, 7, 10, 0) means the input
                sequence consists of 7 pure CAGs and 3 CAGs with an interruption at position 0 of the repeat unit
                (could be "AAG", "GAG" or "TAG"),
            has_partial_repeats
    """

    # check for repeats with interruptions
    repeat_unit_length = 1
    while repeat_unit_length <= len(sequence)/2:
        if not allow_partial_repeats and len(sequence) % repeat_unit_length != 0:
            repeat_unit_length += 1
            continue

        repeat_unit, num_pure_repeats = get_most_common_repeat_unit(sequence, repeat_unit_length)
        num_total_repeats_expected = int(len(sequence) / repeat_unit_length)

        if repeat_unit_length == 1:
            # apply the MIN_PURE_REPEAT_FRACTION_FOR_HOMOPOLYMER_REPEATS threshold
            if num_pure_repeats < num_total_repeats_expected * MIN_PURE_REPEAT_FRACTION_FOR_HOMOPOLYMER_REPEATS:
                repeat_unit_length += 1
                continue

            # Don't allow homopolymer interruptions to be in consecutive bases. If 2 or more bases in a row deviate
            # from the homopolymer repeat, the sequence is no longer considered to be a homopolymer.
            homopolymer_regex_key = f"{repeat_unit}_homopolymer"
            if homopolymer_regex_key not in COMPILED_REGEX_CACHE:
                all_bases_except_homopolymer_base = "ACGT".replace(repeat_unit, "")
                regex = "["+all_bases_except_homopolymer_base+"]"
                regex = f"{regex}{regex}+"
                COMPILED_REGEX_CACHE[homopolymer_regex_key] = re.compile(regex)

            if COMPILED_REGEX_CACHE[homopolymer_regex_key].search(sequence):
                repeat_unit_length += 1
                continue

            repeat_unit_interruption_index = 0 if num_pure_repeats < num_total_repeats_expected else None
            return repeat_unit, num_pure_repeats, num_total_repeats_expected, repeat_unit_interruption_index, False

        if repeat_unit_length == 2:
            if num_pure_repeats < num_total_repeats_expected * MIN_PURE_REPEAT_FRACTION_FOR_DINUCLEOTIDE_REPEATS:
                # apply the MIN_PURE_REPEAT_FRACTION_FOR_DINUCLEOTIDE_REPEATS threshold
                repeat_unit_length += 1
                continue

        if repeat_unit_length <= MAX_INTERRUPTED_REPEAT_UNIT_LENGTH:
            # allow position i within the repeat unit to vary across repeats
            for i in range(repeat_unit_length):
                repeat_unit_regex = get_repeat_unit_regex_with_N_base(repeat_unit, i, allow_homopolymer=False)
                regex = "^(" + repeat_unit_regex + ")+"
                if allow_partial_repeats:
                    partial_repeat_length = len(sequence) % repeat_unit_length
                    regex += repeat_unit[:partial_repeat_length]
                regex += "$"

                if regex not in COMPILED_REGEX_CACHE:
                    COMPILED_REGEX_CACHE[regex] = re.compile(regex)

                if COMPILED_REGEX_CACHE[regex].match(sequence):
                    num_pure_repeats = count_pure_repeats(sequence, repeat_unit)
                    repeat_unit_interruption_index = i if num_pure_repeats < num_total_repeats_expected else None
                    has_partial_repeats = len(sequence) % repeat_unit_length > 0
                    return repeat_unit, num_pure_repeats, num_total_repeats_expected, repeat_unit_interruption_index, has_partial_repeats
        else:
            # check for pure repeats only
            repeat_unit = sequence[:repeat_unit_length]
            num_repeats = sequence.count(repeat_unit)
            if num_repeats * repeat_unit_length == len(sequence):
                return repeat_unit, num_repeats, num_repeats, None, False

            if allow_partial_repeats:
                partial_repeat_length = len(sequence) % repeat_unit_length
                if num_repeats * repeat_unit_length == len(sequence) - partial_repeat_length \
                        and sequence[num_repeats * repeat_unit_length:] == repeat_unit[:partial_repeat_length]:
                    return repeat_unit, num_repeats, num_repeats, None, True

        repeat_unit_length += 1

    # no repeat unit found, so return the input sequence itself as the repeat unit
    return sequence, 1, 1, None, False


def extend_repeat_into_sequence_without_allowing_interruptions(repeat_unit, sequence):
    """This method walks along the given sequence from left to right, one repeat unit length at a time, and returns
    the longest stretch of repeats that exactly matches the given repeat_unit.

    Args:
        repeat_unit (str): For example "CAG"
        sequence (str): a longer sequence that may contain repeats of the given repeat unit starting at the left end,
            before switching to random other sequence or repeats of a different repeat unit.
            For example: "CAGCAGCAGCTAGTGCAGTGACAGT"

    Return:
        int: The number of exact repeats of the given repeat unit found at the left end of the given sequence.
    """

    # compute the number of pure repeats
    i = 0
    num_pure_repeats = 0
    while i <= len(sequence) - len(repeat_unit):
        if sequence[i:i+len(repeat_unit)] != repeat_unit:
            break
        num_pure_repeats += 1
        i += len(repeat_unit)

    return num_pure_repeats


def extend_repeat_into_sequence_allowing_interruptions(
        repeat_unit, 
        sequence,
        repeat_unit_interruption_index=None,
):
    """This method walks along the given sequence from left to right, one repeat unit length at a time, and returns
    the longest stretch of repeats that matches the given repeat_unit, allowing for certain types of interruptions.

    Args:
        repeat_unit (str): For example "CAG"
        sequence (str): a longer sequence that may contain repeats of the given repeat unit starting at the left end,
            before switching to random other sequence or repeats of a different repeat unit.
            For example: "CAGCAGCAGCTAGTGCAGTGACAGT"
        repeat_unit_interruption_index (int): This should be an integer between 0 and len(repeat_unit) or None. If not
            None, this position in the repeat unit will be treated as N (eg. allowed to be any nucleotide A,C,G, or T)
            If None, the sequence will be extended allowing for up to 1 base in the repeat unit to be treated as N (eg.
            allowed to be any nucleotide A,C,G, or T). Even with interruptions though, the repeat sequence must end in
            at least one pure copy of the given repeat_unit (eg. the "CAT" interruption will be included if
            the extended sequence is "CAGCAGCATCAG" , but not if it's just "CAGCAGCAT")

    Return:
        3-tuple: (int, int, int) num_pure_repeats, num_total_repeats, repeat_unit_interruption_index
    """
    if len(repeat_unit) < 3 or len(repeat_unit) > MAX_INTERRUPTED_REPEAT_UNIT_LENGTH:
        # only allow interruptions for repeat units of length 3 up to MAX_INTERRUPTED_REPEAT_UNIT_LENGTH
        num_pure_repeats = extend_repeat_into_sequence_without_allowing_interruptions(repeat_unit, sequence)
        return num_pure_repeats, num_pure_repeats, repeat_unit_interruption_index

    # extend allowing for interruptions
    num_pure_repeats = num_total_repeats = 0

    if repeat_unit_interruption_index is None:
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
                num_pure_repeats = count_pure_repeats(matching_impure_sequence, repeat_unit)
                num_total_repeats = current_num_total_repeats
                if num_pure_repeats < num_total_repeats:
                    repeat_unit_interruption_index = current_index

            rest_of_sequence = match.group(2)
            if len(rest_of_sequence) == 0:
                print(f"WARNING: extend_repeat_into_sequence({repeat_unit}) reached the end of the "
                      f"{len(sequence)}bp sequence: {num_total_repeats}x{repeat_unit} ==> {sequence}. "
                      f"Consider increasing the length of the flanking sequencing.")

    return num_pure_repeats, num_total_repeats, repeat_unit_interruption_index


def count_pure_repeats(sequence, repeat_unit):
    """Count the number of pure repeats of the given repeat unit in the given sequence."""

    if not sequence or not repeat_unit or len(sequence) < len(repeat_unit):
        raise ValueError(f"Invalid sequence arg: '{sequence}' or repeat unit arg: '{repeat_unit}'")

    num_pure_repeats = 0
    for i in range(0, len(sequence) - len(repeat_unit) + 1, len(repeat_unit)):
        if sequence[i:i+len(repeat_unit)] == repeat_unit:
            num_pure_repeats += 1

    return num_pure_repeats


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

