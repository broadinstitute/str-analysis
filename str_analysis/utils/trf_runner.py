"""
Utils for using the TandemRepeatsFinder [Benson 1999] to find repeats in nucleotide sequences.
"""

import os
import subprocess
import tempfile
import unittest

from str_analysis.utils.dat_utils import parse_dat_file


class TRFRunner:
    """Class for piping nucleotide sequences to a background TandemRepeatFinder process."""

    def __init__(self,
                 trf_command_path,
                 match_score = 2,
                 mismatch_penalty = 7,
                 indel_penalty = 7,
                 pm = 80,
                 pi = 10,
                 minscore = 24,
                 ):

        self.trf_command_path = trf_command_path
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.indel_penalty = indel_penalty
        self.pm = pm
        self.pi = pi
        self.minscore = minscore


    def run_TRF(self, nucleotide_sequence, chromosome_name="chrN"):
        temp_fasta_path = _write_sequences_to_temp_fasta({chromosome_name : nucleotide_sequence})
        for dat_record in _run_tandem_repeats_finder(
                self.trf_command_path,
                temp_fasta_path,
                match_score = self.match_score,
                mismatch_penalty = self.mismatch_penalty,
                indel_penalty = self.indel_penalty,
                pm = self.pm,
                pi = self.pi,
                minscore = self.minscore,
        ):
            yield dat_record


def _write_sequences_to_temp_fasta(sequences):
    """Create a temporary .fasta file with the given nucleotide sequences.

    Args:
        sequences (dict): dictionary that maps chromosome names to nucleotide sequences

    Return:
        str: temp file path
    """
    temp_fasta_file = tempfile.NamedTemporaryFile("w+", suffix=".fasta", delete=False)
    for chromosome_name in sequences:
        temp_fasta_file.write(f">{chromosome_name}\n")
        temp_fasta_file.write(f"{sequences[chromosome_name]}\n")
    temp_fasta_file.close()

    return temp_fasta_file.name


def _run_tandem_repeats_finder(
        trf_path,
        input_fasta_path,
        match_score = 2,
        mismatch_penalty = 3,
        indel_penalty = 5,
        pm = 80,
        pi = 10,
        minscore = 8,
):
    """Runs Tandem Repeats Finder on the given input fasta, parses results and returns a generator of DatRecords

    Args:
        trf_path (str): TRF executable path
        input_fasta_path (str): Fasta file path to run on.
        match_score (int): see Tandem Repeats Finder docs
        mismatch_penalty (int): see Tandem Repeats Finder docs
        indel_penalty (int): see Tandem Repeats Finder docs
        pm (int): see Tandem Repeats Finder docs
        pi (int): see Tandem Repeats Finder docs
        minscore (int): see Tandem Repeats Finder docs

    Yield:
        DatRecord: each representing an output row from TandemRepeatFinder
    """

    output_dat_path = f"{input_fasta_path}.dat"

    command = f"{trf_path} {input_fasta_path} {match_score} {mismatch_penalty} {indel_penalty} {pm} {pi} {minscore} 2000 -h -l 6 -ngs > {output_dat_path}"

    #logger.info(f"Running: {command}")
    subprocess.check_output(command, shell=True)

    for dat_record in parse_dat_file(output_dat_path):
        yield dat_record

    os.remove(output_dat_path)



class Tests(unittest.TestCase):
    def test_trf_runner(self):

        trf_runner = TRFRunner(os.path.expanduser("~/bin/trf409.macosx"), mismatch_penalty=3, indel_penalty=5, minscore=8)

        repeat_sequence = "ACACACACAC"
        dat_records = list(trf_runner.run_TRF(repeat_sequence))
        self.assertEqual(len(dat_records), 1)
        self.assertEqual(dat_records[0].repeat_unit, "AC")
        self.assertEqual(dat_records[0].repeat_count, 5)

        repeat_sequence = "ACGAACGAACGAACGAATGA"
        dat_records = list(trf_runner.run_TRF(repeat_sequence))
        self.assertEqual(len(dat_records), 1)
        self.assertEqual(dat_records[0].repeat_unit, "ACGA")
        self.assertEqual(dat_records[0].repeat_count, 5)
