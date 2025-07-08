"""
Utils for using the TandemRepeatsFinder [Benson 1999] to find repeats in nucleotide sequences.
"""

import os
import re
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

        self.trf_command_path = os.path.expanduser(trf_command_path)
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.indel_penalty = indel_penalty
        self.pm = pm
        self.pi = pi
        self.minscore = minscore


    def run_TRF(self, nucleotide_sequence, chromosome_name="chrN", min_motif_size=None, max_motif_size=None):
        """Run TRF on the given nucleotide sequence and return a generator of DatRecords"""
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
                output_format="dat",
                min_motif_size=min_motif_size,
                max_motif_size=max_motif_size,
        ):
            yield dat_record

    def run_TRF_and_get_motif_composition(self, nucleotide_sequence, chromosome_name="chrN", min_motif_size=None, max_motif_size=None):
        """Run TRF on the given nucleotide sequence and return a generator of DatRecords with motif composition"""
        temp_fasta_path = _write_sequences_to_temp_fasta({chromosome_name : nucleotide_sequence})
        for record in _run_tandem_repeats_finder(
                self.trf_command_path,
                temp_fasta_path,
                match_score = self.match_score,
                mismatch_penalty = self.mismatch_penalty,
                indel_penalty = self.indel_penalty,
                pm = self.pm,
                pi = self.pi,
                minscore = self.minscore,
                output_format="html",
                min_motif_size=min_motif_size,
                max_motif_size=max_motif_size,
        ):
            yield record


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
        output_format="dat",
        output_filename_prefix=None,
        min_motif_size=None,
        max_motif_size=None,
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
        output_format (str): this should be "dat" or "html"
        output_filename_prefix (str): if specified, the output will be saved to a file with this prefix.
        min_motif_size (int): if specified, only return loci with this motif size or larger.
        max_motif_size (int): if specified, only return loci with this motif size or smaller.

    Yield:
        if output_format was = "dat", yield DatRecords each of which represents an output row from TandemRepeatFinder
        if output_format was = "html", yield records representing 
    """
    # -l 6  = maximum TR length expected (in millions) (eg, -l 3 or -l=3 for 3 million). Human genome HG38 would need -l 6
    # -h = suppress html output
    # -ngs =  more compact .dat output on multisequence files, returns 0 on success. Output is printed to the screen.
    max_period = 2000  # maximum period size to report. Must be between 1 and 2000, inclusive
    command = f"{trf_path} {input_fasta_path} {match_score} {mismatch_penalty} {indel_penalty} {pm} {pi} {minscore} {max_period} -l 6 "

    if output_format == "dat":
        output_dat_path = f"{input_fasta_path}.dat"
        command += f"-h -ngs > {output_dat_path}"

        subprocess.check_output(command, shell=True)

        for dat_record in parse_dat_file(output_dat_path):
            if min_motif_size is not None and dat_record.repeat_unit_length < min_motif_size:
                continue
            if max_motif_size is not None and dat_record.repeat_unit_length > max_motif_size:
                continue

            yield dat_record

        os.remove(output_dat_path)

    elif output_format == "html":

        subprocess.run(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)

        # read html output
        trf_output_filename_prefix = f"{os.path.basename(input_fasta_path)}.{match_score}.{mismatch_penalty}.{indel_penalty}.{pm}.{pi}.{minscore}.{max_period}"
        i = 1
        html_file_path =  f"{trf_output_filename_prefix}.{i}.txt.html"
        while os.path.isfile(html_file_path):

            for record in _parse_trf_html_output(html_file_path, min_motif_size=min_motif_size, max_motif_size=max_motif_size):
                yield record

            i += 1
            html_file_path =  f"{trf_output_filename_prefix}.{i}.txt.html"

    else:
        raise ValueError(f"Invalid output_format: '{output_format}'. It must be either 'dat' or 'html'.")


def _parse_trf_html_output(html_file_path, min_motif_size=None, max_motif_size=None):
    """Parse the HTML output of TRF and return a list of motifs"""

    with open(html_file_path, "rt") as f:
        html_content = f.read()

    results = []
    for i, html_section in enumerate(html_content.split("<A NAME=")[1:]):
        #print(f"Parsing html section #{i+1}")


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

        if min_motif_size is not None and consensus_size < min_motif_size:
            continue
        if max_motif_size is not None and consensus_size > max_motif_size:
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
                        print(f"WARNING: {html_file_path} consensus motif mismatch: {consensus_motif_match.group(2)} != {consensus_motif} at line: {consensus_motif_line}")
                    continue

                consensus_motif = consensus_motif_match.group(2).replace("-", "")

                motif_line = lines_in_block[-2].strip()
                #use regex to parse the motif line: 121 GACCCTGACCTTACTAGTTTACAACCACAC
                motif_match = re.search(r"(\d+)\s+([-ACGTN\s]+)", motif_line)
                if not motif_match:
                    raise ValueError(f"Could not parse motif from line: {motif_line}")
                motif_start_index = int(motif_match.group(1)) - 1
                motif_sequence = motif_match.group(2).replace("-", "")
                for motif in motif_sequence.split(" "):
                    motifs_and_start_indices.append((motif_start_index, motif))
                    motifs.append(motif)
                    motif_start_index += len(motif)
            except Exception as e:
                print(f"Error parsing {html_file_path} at lines {start_block_i} to {i}: {block_text}\n{e}")
                continue

            block_counter += 1

        results.append({
            "motif": consensus_motif,
            "motif_size": len(consensus_motif),
            #"motif_size": period,
            #"consensus_motif_size": consensus_size,
            "num_repeats": copynumber,
            "repeat_sequence_length": indices[1] - indices[0],
            "motifs": motifs,
            "motifs_and_start_indices": motifs_and_start_indices,
            "score": score,
            #"motif_positions_with_interruptions": motif_positions_with_interruptions,
        })

    return results



"""
    # Selecting the best TRF result for a locus
    result_with_max_score = max(results, key=lambda x: x["trf_score"])
    max_score = result_with_max_score["trf_score"]

    result_with_min_consensus_motif_size = min([t for t in results if t["trf_score"] > 0.8 * max_score], key=lambda x: x["trf_consensus_motif_size"])
"""



if __name__ == "__main__":
    seq = "CAG"*2 + "CTG" + "CAG"*10 + "CTG"*2 + "CAG"*100

    r = TRFRunner("~/bin/trf409.macosx")
    records = list(r.run_TRF_and_get_motif_composition(
        nucleotide_sequence=seq,
        chromosome_name="chrN",
    ))

    print(f"Found {len(records)} TRF results in {seq}")
    for record in records:
        print(record)

