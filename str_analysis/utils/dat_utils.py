"""
Utils for using the TandemRepeatsFinder [Benson 1999] to find repeats in nucleotide sequences.
"""

from collections import namedtuple
import gzip
import logging
import re


logger = logging.getLogger()

# Represents a single row of TRF output
DatRecord = namedtuple('DatRecord', [
    'chrom', 'start_0based', 'end_1based',
    'repeat_unit', 'repeat_count', 'repeat_unit_length',
    'percent_matches', 'percent_indels', 'alignment_score', 'entropy',
    'full_repeat_sequence', 'full_repeat_sequence_length',
    'left_flank', 'right_flank'])



def parse_dat_record(line, chromosome):
    """Parse a data line from Tandem Repeat Finder .dat output format. Similar to .fasta files, the
    chromosome name is on a separate header line, so must be passed in.
    """

    # TRF output data columns are:
    ### Indices of the repeat relative to the start of the sequence.
    ### Period size of the repeat.
    ### Number of copies aligned with the consensus pattern.
    ### Size of consensus pattern (may differ slightly from the period size).
    ### Percent of matches between adjacent copies overall.
    ### Percent of indels between adjacent copies overall.
    ### Alignment score.
    ### Percent composition for each of the four nucleotides.
    ### Entropy measure based on percent composition.

    fields = line.split()

    dat_record = {
        'chrom': chromosome,
        'start_0based': int(fields[0]) - 1,   # bed file start is 0-based, TRF output is 1-based
        'end_1based': int(fields[1]),
        'repeat_unit_length': int(fields[2]),
        'repeat_count': float(fields[3]),  # not always a whole number - because of inexact repeats and indels
        #'consensus_length': int(fields[4]),
        'percent_matches': int(fields[5]),
        'percent_indels': int(fields[6]),
        'alignment_score': int(fields[7]),
        'entropy': float(fields[12]),
        'repeat_unit': fields[13],
        'full_repeat_sequence': fields[14],
        'full_repeat_sequence_length': len(fields[14]),
        'left_flank': None,
        'right_flank': None,
    }

    if len(fields) > 16:
        dat_record['left_flank'] = fields[15]
        dat_record['right_flank'] = fields[16]

    return DatRecord(**dat_record)


def parse_dat_file(dat_file_path, dat_in_original_format=False, limit=None):
    """Parse a Tandem Repeat Finder output .dat file and yield DatRecord objects.

    Args:
        dat_file_path (str): .dat file
        dat_in_original_format (bool): set to True if the .dat file was generated without passing -ngs to Tandem Repeat Finder.
        limit (int): returns only the first N records. Useful for testing.
    Yield:
        DatRecord object for each row in the .dat file
    """
    chrom = None

    fopen = gzip.open if dat_file_path.endswith("gz") else open
    with fopen(dat_file_path, "rt") as input_dat_file:
        for i, line in enumerate(input_dat_file):
            line = line.strip()
            if not line:
                continue

            if limit is not None and i > limit:
                break

            # parse sequence header
            if line.startswith("@") or dat_in_original_format:
                if dat_in_original_format and any(line.startswith(header_line_prefix) for header_line_prefix in [
                    "Tandem", "Gary", "Program", "Boston", "Version", "Parameters"
                ]):
                    continue

                is_header_line = line.startswith("@") or line.startswith("Sequence:")
                if is_header_line:
                    if line.startswith("@"):
                        line = line[1:]
                    elif line.startswith("Sequence:"):
                        line = line[len("Sequence:"):]

                    if " chromosome " in line:
                        match = re.search("chromosome ([chrXYMT0-9]+)[, ]", line)
                        if not match:
                            raise ValueError(f"Couldn't parse chromosome from line: {line}")
                        chrom = match.group(1)
                    else:
                        chrom = line.strip()

                    continue

            # parse data row
            yield parse_dat_record(line, chrom)

