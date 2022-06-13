import logging
import subprocess

CHROMOSOME_ORDER = list(map(str, range(1,23))) + ["X", "Y", "M", "MT"]
CHROMOSOME_ORDER += [f"chr{s}" for s in CHROMOSOME_ORDER]


def intervals_in_genomic_sort_order(interval_strings):
    """Sorts a list of intervals by genomic coordinates. Takes a list of 'chr:start-end' interval strings like
    ['chr1:12345-54321', 'chr1:23456-65432', 'chr3:34567-76543'] and returns the same list, but sorted by genomic
    coordinates (based on the chromosome, start coordinate, end coordinate).
    """

    def sort_key(interval_string):
        chrom, positions = interval_string.split(":")
        if chrom in CHROMOSOME_ORDER:
            chrom_ordinal = CHROMOSOME_ORDER.index(chrom)
        else:
            chrom_ordinal = sum(ord(c)*10**(len(chrom) - i) for i, c in enumerate(chrom))

        start_pos, end_pos = positions.split("-")
        return chrom_ordinal, int(start_pos), int(end_pos)

    return sorted(interval_strings, key=sort_key)


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


def parse_key_value(key_value, delimiter="="):
    """Splits a string like "x=3" into a 2-tuple ("x", "3")

    Args:
        key_value (str): The string containing the key and value separated by the delimiter.
        delimiter (str): Separator between the key and value. Default is "=".

    Return:
        2-tuple: (key string, value string)
    """

    key_value_list = key_value.split(delimiter)
    if len(key_value_list) != 2 or not key_value_list[0] or not key_value_list[1]:
        raise ValueError(f"Invalid arg {key_value}")

    return tuple(key_value_list)


