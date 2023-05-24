"""Ported from ExpansionHunterDenovo code:

https://github.com/Illumina/ExpansionHunterDenovo/blob/master/source/reads/IrrFinder.cpp#L125
"""


COMPLEMENT = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
    'N': 'N',
    'Y': 'R',   # source: https://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
    'R': 'Y',
    'S': 'S',
    'W': 'W',
    'M': 'K',
    'K': 'M',
    'B': 'V',
    'V': 'B',
    'D': 'H',
    'H': 'D',
}


def reverse_complement(dna):
    """Take a string representing a DNA sequence and return its reverse-complement"""
    return "".join([COMPLEMENT[c] for c in dna[::-1]])


def _alphabetically_first_motif_under_shift(motif):
    """Return the permutation of the given repeat motif that's first alphabetically.

    Args:
        motif (str): A repeat motif like "CAG".
    Return:
        str: The alphabetically first repeat motif.
    """
    minimal_motif = motif
    double_motif = motif + motif
    for i in range(len(motif)):
        current_motif = double_motif[i:i+len(motif)]
        if current_motif < minimal_motif:
            minimal_motif = current_motif
    return minimal_motif


def compute_canonical_motif(motif, include_reverse_complement=True):
    """Take an STR motif string like "GAA" and returns the "canonical" representation of the motif. This is the
    rearrangement of the motif bases (including reverse complements) that is alphabetically first. 
    
    For "GAA", this would return "AAG" since it's alphabetically first among "GAA", "AGA", "AAG", "TTC", "TCT", and "CTT".

    Args:
        motif (str): A repeat motif like "CAG".
        include_reverse_complement (bool): if False, will not consider the reverse-complement of the given motif.  
    Return:
        str: The alphabetically first repeat motif.
    """
    motif = motif.upper()
    minimal_motif = _alphabetically_first_motif_under_shift(motif)

    if include_reverse_complement:
        motif_rc = reverse_complement(motif)
        minimal_motif_rc = _alphabetically_first_motif_under_shift(motif_rc)
    
        if minimal_motif_rc < minimal_motif:
            return minimal_motif_rc

    return minimal_motif

