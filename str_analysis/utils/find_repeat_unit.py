
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

