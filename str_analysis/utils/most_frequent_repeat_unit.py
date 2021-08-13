

def compute_most_frequent_repeat_unit(sequence, repeat_unit_size=5, min_occurrences=3, min_fraction_bases_covered=0.8):
    """Returns the most frequent repeat unit of the given size within the given sequence.

    Args:
        sequence (str): a sequence of dna bases
        repeat_unit_size (int): exact repeat unit size in bases
        min_occurrences (int): repeat_units that occur fewer than this many times will be ignored
        min_fraction_bases_covered (float): repeat units that cover less than this fraction the sequence will be ignored

    Returns:
        2-tuple: (repeat unit (str), number of times it occurs within sequence (int)) If no repeat unit satisfies the
            min_fraction_bases_covered constraint, the return value is (None, 0).
    """
    if len(sequence) < repeat_unit_size:
        raise ValueError(f"len(sequence) < repeat_unit_size: {len(sequence)} < {repeat_unit_size}")

    repeat_units = {}  # maps repeat_unit to the number of times it occurs in sequence
    for i in range(len(sequence) - repeat_unit_size + 1):
        current_repeat_unit = sequence[i:i+repeat_unit_size]
        if current_repeat_unit not in repeat_units:
            repeat_units[current_repeat_unit] = sequence.count(current_repeat_unit)

    most_frequent_repeat_unit = None
    most_frequent_repeat_unit_occurrences = 0
    for repeat_unit, count in repeat_units.items():
        if (
                count > most_frequent_repeat_unit_occurrences and
                count >= min_occurrences and
                float(count * repeat_unit_size)/len(sequence) >= min_fraction_bases_covered

        ):
            most_frequent_repeat_unit = repeat_unit
            most_frequent_repeat_unit_occurrences = count

    return most_frequent_repeat_unit, most_frequent_repeat_unit_occurrences
