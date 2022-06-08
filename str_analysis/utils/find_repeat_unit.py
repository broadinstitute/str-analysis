import os
from utils.trf_runner import TRFRunner


def find_repeat_unit(s, min_fraction_covered_by_repeat=0.9):
    len_s = len(s)
    if len_s == 1:
        raise ValueError(f"len(s) must be > 1: {s}")

    repeat_unit_length = 1
    while repeat_unit_length <= len(s)/2:
        repeat_unit = s[:repeat_unit_length]
        num_repeats = s.count(repeat_unit)
        if num_repeats * repeat_unit_length / len_s >= min_fraction_covered_by_repeat:
            return repeat_unit, num_repeats
        repeat_unit_length += 1

    return s, 1


def find_repeat_unit_using_trf(
    s, min_fraction_covered_by_repeat=0.9, min_sequence_length_for_running_trf=12, trf_path="~/bin/trf409.macosx",
):
    """Run TandemRepeatFinder [Benson 1999] to find repeat unit allowing for mismatches before falling back on the
    simpler find_repeat_unit(..) method.

    Args:
        s (str): The nucleotide sequence to search for repeats
        min_fraction_covered_by_repeat (float): A detected repeat is discarded unless it covers this fraction of s.
        min_sequence_length_for_running_trf (str): TRF will only be run if s is at least this many nucleotides long.
        trf_path (str): Path of TRF executable.

    Return:
        3-tuple: repeat unit, number of repeats found in the sequence, whether the repeats were found by TRF or the
        simpler pure-repeat finding algorithm. If no repeats are found, the input sequence is returned with modification,
        along with repeat count = 1, and False,
    """
    len_s = len(s)
    if len_s == 1:
        raise ValueError(f"len(s) must be > 1: {s}")

    # try using TRF to find repeat sequence
    found_by_TRF = True
    if len_s >= min_sequence_length_for_running_trf:
        trf_runner = TRFRunner(os.path.expanduser(trf_path), mismatch_penalty=7, indel_penalty=100000, minscore=1)
        dat_records = list(trf_runner.run_TRF(s))
        for dat_record in dat_records:
            repeat_count = int(dat_record.repeat_count)
            if repeat_count * len(dat_record.repeat_unit) / len_s >= min_fraction_covered_by_repeat:
                return dat_record.repeat_unit, repeat_count, found_by_TRF

    found_by_TRF = False
    repeat_unit, repeat_count = find_repeat_unit(s, min_fraction_covered_by_repeat=min_fraction_covered_by_repeat)
    return repeat_unit, repeat_count, found_by_TRF
