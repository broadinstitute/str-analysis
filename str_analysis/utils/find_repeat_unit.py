import os
from str_analysis.utils.trf_runner import TRFRunner


# TODO add docs
def find_repeat_unit(s, min_fraction_covered_by_repeat=0.9, allow_indel_interruptions=False):
    if len(s) <= 1:
        raise ValueError(f"len(s) must be > 1: {s}")

    result_repeat_unit = s
    result_num_repeats = 1
    max_fraction_covered_so_far = 0

    # find the smallest repeat unit that covers the largest fraction of the input sequence
    repeat_unit_length = 1
    while repeat_unit_length <= len(s)/2:
        if not allow_indel_interruptions and len(s) % repeat_unit_length != 0:
            repeat_unit_length += 1
            continue

        repeat_unit = s[:repeat_unit_length]
        num_repeats = s.count(repeat_unit)
        fraction_covered = num_repeats * repeat_unit_length / len(s)
        if fraction_covered >= min_fraction_covered_by_repeat and fraction_covered > max_fraction_covered_so_far:
            max_fraction_covered_so_far = fraction_covered
            result_repeat_unit = repeat_unit
            result_num_repeats = int(len(s) / len(repeat_unit))

        repeat_unit_length += 1

    return result_repeat_unit, result_num_repeats


def find_repeat_unit_using_trf(
    s, min_fraction_covered_by_repeat=0.9, trf_path="~/bin/trf409.macosx",
):
    """Run TandemRepeatFinder [Benson 1999] to find the repeat unit allowing for mismatches before falling back on the
    simpler find_repeat_unit(..) method.

    Args:
        s (str): The nucleotide sequence to search for repeats
        min_fraction_covered_by_repeat (float): A detected repeat is discarded unless it covers this fraction of s.
        trf_path (str): Path of TRF executable.

    Return:
        2-tuple: repeat unit, number of repeats found in the sequence. If no repeats are found, the input sequence is
        returned without modification, along with repeat count = 1.
    """
    if len(s) <= 1:
        raise ValueError(f"len(s) must be > 1: {s}")

    # if the simple approach didn't find repeats within the sequence, try using TRF to find imperfect repeats
    trf_runner = TRFRunner(
        os.path.expanduser(trf_path),
        match_score=2,
        mismatch_penalty=7,
        indel_penalty=10**6,
        minscore=len(s)*2*min_fraction_covered_by_repeat,
    )

    dat_records = list(trf_runner.run_TRF(s))

    for dat_record in dat_records:
        if len(s) % len(dat_record.repeat_unit) == 0:
            repeat_count = int(len(s) / len(dat_record.repeat_unit))
            repeat_unit = dat_record.repeat_unit
            assert repeat_count * len(repeat_unit) / len(s) >= min_fraction_covered_by_repeat

            return repeat_unit, repeat_count

    return s, 1
