"""
Utils for using the TandemRepeatsFinder [Benson 1999] to find repeats in nucleotide sequences.
"""

import os
import subprocess
import tempfile
import unittest

from str_analysis.utils.dat_utils import parse_dat_file


class TRFRunner:
    """Class for running TandemRepeatFinder on a nucleotide sequence."""

    def __init__(self,
                 trf_command_path,
                 match_score = 2,
                 mismatch_penalty = 7,
                 indel_penalty = 7,
                 pm = 80,
                 pi = 10,
                 minscore = 24):

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
        mismatch_penalty = 7,
        indel_penalty = 7,
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

    # -l 6  = maximum TR length expected (in millions) (eg, -l 3 or -l=3 for 3 million). Human genome HG38 would need -l 6
    # -h = suppress html output
    # -ngs =  more compact .dat output on multisequence files, returns 0 on success. Output is printed to the screen.
    max_period = 2000  # maximum period size to report. Must be between 1 and 2000, inclusive
    command = f"{trf_path} {input_fasta_path} {match_score} {mismatch_penalty} {indel_penalty} {pm} {pi} {minscore} {max_period} "
    command += f"-h -l 6 -ngs "
    command += f"> {output_dat_path}"

    subprocess.check_output(command, shell=True)

    for dat_record in parse_dat_file(output_dat_path):
        yield dat_record

    os.remove(output_dat_path)



def parse_trf_html_output(input_seq, locus_name, trf_output_filename):
    """Parse the HTML output of TRF and return a list of motifs"""

    with open(trf_output_filename, "rt") as f:
        html_content = f.read()

    results = []
    for i, html_section in enumerate(html_content.split("<A NAME=")[1:]):
        #print(f"Parsing html section #{i+1}")

        result_row = {}

        # find the line that starts with "Indices:"
        html_section = html_section.split("Indices: ")[1]
        html_section = html_section.split("Statistics")[0].strip()

        lines = html_section.split("\n")

        # Parse indices and score using regex - example: "Indices: 1--5730  Score: 9516"
        indices_score_match = re.search(r"(\d+)--(\d+)\s+Score: (\d+)", lines[0])
        if not indices_score_match:
            raise ValueError(f"Could not parse indices and score from line: {lines[0]}")

        start_index = int(indices_score_match.group(1))
        end_index = int(indices_score_match.group(2))
        score = int(indices_score_match.group(3))
        indices = [start_index - 1, end_index]  # Adjust start index to be 0-based


        # Parse period, copynumber, and consensus size using regex - example: "Period size: 30  Copynumber: 191.0  Consensus size: 30"
        period_copynumber_consensus_size_match = re.search(r"Period size: (\d+)\s+Copynumber: (\d+\.?\d*)\s+Consensus size: (\d+)", lines[1])
        if not period_copynumber_consensus_size_match:
            raise ValueError(f"Could not parse period, copynumber, and consensus size from line: {lines[1]}")

        period = int(period_copynumber_consensus_size_match.group(1))
        copynumber = float(period_copynumber_consensus_size_match.group(2))
        consensus_size = int(period_copynumber_consensus_size_match.group(3))


        result_row[f"{locus_name}_trf_score"] = score
        #result_row[f"{locus_name}_trf_period"] = period
        result_row[f"{locus_name}_trf_copy_number"] = copynumber
        result_row[f"{locus_name}_trf_consensus_motif_size"] = consensus_size
        result_row[f"{locus_name}_trf_repeat_sequence_length"] = indices[1] - indices[0]
        result_row[f"{locus_name}_trf_seq_coverage"] = result_row[f"{locus_name}_trf_repeat_sequence_length"] / len(input_seq)
        result_row[f"{locus_name}_trf_uncovered_length"] = len(input_seq) - result_row[f"{locus_name}_trf_repeat_sequence_length"]

        if consensus_size >= 50:
            result_row[f"{locus_name}_trf_consensus_motif"] = None
            result_row[f"{locus_name}_trf_motifs"] = []
            result_row[f"{locus_name}_trf_motifs_and_start_indices"] = []
            continue

        """
                    *            *
        121 GACCCTGACCTTACTAGTTTACAACCACAC
          1 GACCCTGACCTGACTAGTTTACAATCACAC
        """
        lines = lines[3:]

        consensus_motif = None
        motifs = []
        motifs_and_start_indices = []
        block_counter = 0
        while i < len(lines):
            start_block_i = i
            while i < len(lines) and lines[i].strip():
                i += 1

            lines_in_block = lines[start_block_i:i]
            block_text = "\n".join(lines_in_block)
            i += 1
            if not block_text.strip() or len(lines_in_block) < 2:
                continue

            try:
                consensus_motif_line = lines_in_block[-1].strip()
                consensus_motif_match = re.search(r"(\d+)\s+([ACGTN-]+)", consensus_motif_line)
                if not consensus_motif_match:
                    raise ValueError(f"Could not parse consensus motif from line: {consensus_motif_line}")

                if consensus_motif is not None and consensus_motif != consensus_motif_match.group(2).replace("-", ""):
                    if i + 2 < len(lines):
                        print(f"WARNING: {trf_output_filename} consensus motif mismatch: {consensus_motif_match.group(2)} != {consensus_motif} at line: {consensus_motif_line}")
                    continue

                consensus_motif = consensus_motif_match.group(2).replace("-", "")

                motif_line = lines_in_block[-2].strip()
                #use regex to parse the motif line: 121 GACCCTGACCTTACTAGTTTACAACCACAC
                motif_match = re.search(r"(\d+)\s+([ACGTN-]+)", motif_line)
                if not motif_match:
                    raise ValueError(f"Could not parse motif from line: {motif_line}")
                motif_start_index = int(motif_match.group(1))
                motif_sequence = motif_match.group(2).replace("-", "")
                motifs_and_start_indices.append((motif_start_index, motif_sequence))
                motifs.append(motif_sequence)
            except Exception as e:
                print(f"Error parsing {trf_output_filename} at lines {start_block_i} to {i}: {block_text}\n{e}")
                continue

            block_counter += 1

        result_row[f"{locus_name}_trf_consensus_motif"] = consensus_motif
        result_row[f"{locus_name}_trf_motifs"] = motifs
        result_row[f"{locus_name}_trf_motifs_and_start_indices"] = motifs_and_start_indices

        results.append(result_row)

    result_with_max_score = max(results, key=lambda x: x[f"{locus_name}_trf_score"])
    max_score = result_with_max_score[f"{locus_name}_trf_score"]

    result_with_min_consensus_motif_size = min([t for t in results if t[f"{locus_name}_trf_score"] > 0.8 * max_score], key=lambda x: x[f"{locus_name}_trf_consensus_motif_size"])

    return result_with_min_consensus_motif_size


