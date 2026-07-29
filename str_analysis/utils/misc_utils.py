import logging
import subprocess

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


def parse_interval(interval_string):
    """Parses interval string like "chr1:12345-54321" and returns 3-tuple (chrom, start, end)"""

    try:
        tokens = interval_string.split(":")
        chrom = ":".join(tokens[:-1])  # some super-contig names have : in them
        start, end = map(int, tokens[-1].split("-"))
    except Exception as e:
        raise ValueError(f"Unable to parse interval: '{interval_string}': {e}")

    return chrom, start, end


def run(command):
    """Run a shell command and return its output. Raises an exception if the command exits with a non-zero exit code"""

    logging.info(command)

    return subprocess.check_output(['/bin/bash', '-c', command]).decode('UTF-8')


