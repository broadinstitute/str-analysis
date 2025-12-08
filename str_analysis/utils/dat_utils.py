"""
Utils for using the TandemRepeatsFinder [Benson 1999] to find repeats in nucleotide sequences.
"""

import gzip
import logging


logger = logging.getLogger()

def parse_dat_record(line, sequence_name):
    """Parse a data line from Tandem Repeat Finder .dat output format. Similar to .fasta files, the
    chromosome name is on a separate header line, so must be passed in.
    """

    fields = line.split()

    dat_record = {
        'chrom': sequence_name,  # alias
        'sequence_name': sequence_name,
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
        'full_repeat_sequence_length': len(fields[14]),
        'full_repeat_sequence': None,  #fields[14],
        'left_flank': None,
        'right_flank': None,
    }

    #if len(fields) > 16:
    #    dat_record['left_flank'] = fields[15] if fields[15] != "." else None
    #    dat_record['right_flank'] = fields[16] if fields[16] != "." else None

    return dat_record


def parse_dat_file(dat_file_path, n=None):
    """Parse a Tandem Repeat Finder output .dat file and yield DatRecord objects.

    Args:
        dat_file_path (str): .dat file
        n (int): returns only the first N records. Useful for testing.
    Yield:
        DatRecord object for each row in the .dat file
    """

    fopen = gzip.open if dat_file_path.endswith("gz") else open
    with fopen(dat_file_path, "rt") as input_dat_file:
        current_sequence_name = None
        for i, line in enumerate(input_dat_file):
            line = line.strip()
            if n is not None and i > 2*n:
                break

            # parse sequence header
            if line.startswith("@"):
                current_sequence_name = line[1:]
                continue

            if current_sequence_name is None:
                raise ValueError(
                    f"{dat_file_path}: line #{i + 1} was not preceded by a line starting with '@' and "
                    f"the sequence name: {line}")

            yield parse_dat_record(line, current_sequence_name)
