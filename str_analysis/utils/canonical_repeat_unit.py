"""
Ported from ExpansionHunterDenovo code:

https://github.com/Illumina/ExpansionHunterDenovo/blob/master/source/reads/IrrFinder.cpp#L125
"""


COMPLEMENT = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
    'N': 'N',
}


def reverse_complement(dna):
    return "".join([COMPLEMENT[c] for c in dna[::-1]])


def _minimal_unit_under_shift(unit):
    minimal_unit = unit
    double_unit = unit + unit
    for i in range(len(unit)):
        current_unit = double_unit[i:i+len(unit)]
        if current_unit < minimal_unit:
            minimal_unit = current_unit
    return minimal_unit


def compute_canonical_repeat_unit(unit):
    minimal_unit = _minimal_unit_under_shift(unit)

    unit_rc = reverse_complement(unit)
    minimal_unit_rc = _minimal_unit_under_shift(unit_rc)

    if minimal_unit_rc < minimal_unit:
        return minimal_unit_rc

    return minimal_unit
