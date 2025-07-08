"""
Utils for using the TandemRepeatsFinder [Benson 1999] to find repeats in nucleotide sequences.
"""

import collections
import os
import re
import subprocess
import tempfile
import traceback
import unittest

from str_analysis.utils.dat_utils import parse_dat_file


class TRFRunner:
    """Class for running TandemRepeatFinder on a nucleotide sequence."""

    def __init__(self,
                 trf_executable_path,
                 match_score = 2,
                 mismatch_penalty = 7,
                 indel_penalty = 7,
                 pm = 80,
                 pi = 10,
                 minscore = 24, 
                 output_filename_prefix=None):
        """
        Args:
            trf_executable_path (str): path to the TRF executable
            match_score (int): see Tandem Repeats Finder docs
            mismatch_penalty (int): see Tandem Repeats Finder docs
            indel_penalty (int): see Tandem Repeats Finder docs
            pm (int): see Tandem Repeats Finder docs
            pi (int): see Tandem Repeats Finder docs
            minscore (int): see Tandem Repeats Finder docs
            output_filename_prefix (str): if specified, the TRF output files will have this prefix. If not specified, intermediate files will be deleted.
        """

        self.trf_executable_path = os.path.expanduser(trf_executable_path)
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.indel_penalty = indel_penalty
        self.pm = pm
        self.pi = pi
        self.minscore = minscore
        self.output_filename_prefix = output_filename_prefix

    def run_TRF(self, nucleotide_sequence, chromosome_name="chrN", min_motif_size=None, max_motif_size=None):
        """Run TRF on the given nucleotide sequence and return a generator of DatRecords"""
        temp_fasta_path = _write_sequences_to_temp_fasta(
            {chromosome_name : nucleotide_sequence},
            output_filename_prefix=self.output_filename_prefix)

        for dat_record in _run_tandem_repeats_finder(
                self.trf_executable_path,
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
                keep_intermediate_files=self.output_filename_prefix is not None,
        ):
            yield dat_record

    def run_TRF_and_parse_motif_composition(self, nucleotide_sequence, chromosome_name="chrN", min_motif_size=None, max_motif_size=None):
        """Run TRF on the given nucleotide sequence and return a generator of records with motif composition"""

        # Write the sequence to a temporary fasta file
        temp_fasta_path = _write_sequences_to_temp_fasta(
            {chromosome_name : nucleotide_sequence}, 
            output_filename_prefix=self.output_filename_prefix)

        # Run TRF and parse the motif composition
        for record in _run_tandem_repeats_finder(
            self.trf_executable_path,
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
            output_filename_prefix=self.output_filename_prefix,
            keep_intermediate_files=self.output_filename_prefix is not None,
        ):
            yield record


def _write_sequences_to_temp_fasta(sequences, output_filename_prefix=None):
    """Create a temporary .fasta file with the given nucleotide sequences.

    Args:
        sequences (dict): dictionary that maps chromosome names to nucleotide sequences

    Return:
        str: temp file path
    """
    if output_filename_prefix is None:
        temp_fasta_file = tempfile.NamedTemporaryFile("w+", suffix=".fasta", delete=False)
    else:
        temp_fasta_file = open(f"{output_filename_prefix}.fasta", "wt")

    for chromosome_name in sequences:
        temp_fasta_file.write(f">{chromosome_name}\n")
        temp_fasta_file.write(f"{sequences[chromosome_name]}\n")
    
    temp_fasta_file.close()

    print(f"Wrote {temp_fasta_file.name}")
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
        keep_intermediate_files=False,
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

        if not keep_intermediate_files:
            os.remove(output_dat_path)

    elif output_format == "html":
        subprocess.run(command, shell=True, stdout=subprocess.DEVNULL, check=False)

        # read html output
        trf_output_filename_prefix = f"{os.path.basename(input_fasta_path)}.{match_score}.{mismatch_penalty}.{indel_penalty}.{pm}.{pi}.{minscore}.{max_period}"
        i = 1
        html_file_path =  f"{trf_output_filename_prefix}.{i}.txt.html"
        while os.path.isfile(html_file_path):

            for record in _parse_trf_html_output(html_file_path, min_motif_size=min_motif_size, max_motif_size=max_motif_size):
                yield record

            if not keep_intermediate_files:
                # remove the TRF output html file
                os.remove(html_file_path)
                if os.path.isfile(f"{trf_output_filename_prefix}.{i}.html"):
                    # also remove the TRF output summary table html file
                    os.remove(f"{trf_output_filename_prefix}.{i}.html")

            i += 1
            html_file_path =  f"{trf_output_filename_prefix}.{i}.txt.html"



    else:
        raise ValueError(f"Invalid output_format: '{output_format}'. It must be either 'dat' or 'html'.")


def _parse_trf_html_output(html_file_path, min_motif_size=None, max_motif_size=None, debug=True):
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

        seq_start_index_0based = int(indices_score_match.group(1)) - 1
        seq_end_index = int(indices_score_match.group(2))
        score = int(indices_score_match.group(3))

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
        #motifs = []
        #motifs_and_start_indices = []
        motif_positions_with_interruptions = {}  # histogram of motif positions with interruptions
        total_mismatches = 0
        total_indels = 0

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
                consensus_motif_line = lines_in_block[-1]
                consensus_motif_match = re.search(r"(\d+)\s+([-ACGTN\s]+)", consensus_motif_line)
                if not consensus_motif_match:
                    raise ValueError(f"Could not parse consensus motif from line: {consensus_motif_line}")

                consensus_motif_sequence_string = consensus_motif_match.group(2)
                current_consensus_motif = consensus_motif_sequence_string.split(" ")[0]
                current_consensus_motif_without_deletions = current_consensus_motif.replace("-", "")
                if consensus_motif is None:
                    consensus_motif = current_consensus_motif_without_deletions
                elif consensus_motif != current_consensus_motif_without_deletions:
                    if i + 2 < len(lines):
                        print(f"WARNING: {html_file_path} consensus motif mismatch: {current_consensus_motif_without_deletions} != {consensus_motif} at line: {consensus_motif_line}")
                    continue

                repeat_sequence_line = lines_in_block[-2]
                repeat_sequence_match = re.search(r"(\d+)\s+([-ACGTN\s]+)", repeat_sequence_line)
                if not repeat_sequence_match:
                    raise ValueError(f"Could not parse motif from line: {repeat_sequence_line}")

                repeat_sequence_string = repeat_sequence_match.group(2)
                total_indels += repeat_sequence_string.count("-") + consensus_motif_sequence_string.count("-")

                #repeat_sequence_start_index_0based = int(repeat_sequence_match.group(1)) - 1

                interruptions_string = None
                if len(lines_in_block) > 2:
                    interruptions_line = lines_in_block[-3]

                    interruptions_line_offset, _ = repeat_sequence_match.span(2)  # offset of the repeat sequence in the interruptions line
                    interruptions_string = interruptions_line[interruptions_line_offset:]

                    total_mismatches += interruptions_string.count("*")
                    if debug:
                        print(interruptions_string)
                        print(repeat_sequence_string)
                        print(consensus_motif_sequence_string)

                interruptions_string_offset = 0
                for motif in repeat_sequence_string.split(" "):
                    #motifs.append(motif)
                    #motifs_and_start_indices.append((repeat_sequence_start_index_0based, motif))
                    #repeat_sequence_start_index_0based += len(motif)
    
                    if interruptions_string is not None:
                        motif_interruptions = interruptions_string[interruptions_string_offset:interruptions_string_offset + len(motif)]
                        if debug:
                            print(f"motif_interruptions: {motif_interruptions}")
                            print(f"motif              : {motif}")
                        if "*" in motif_interruptions:
                            for i, c in enumerate(motif_interruptions):
                                if c == "*":
                                    motif_positions_with_interruptions[i + 1] = motif_positions_with_interruptions.get(i + 1, 0) + 1

                    interruptions_string_offset += len(motif) + 1  # add +1 to account for the space between motifs

            except Exception as e:
                print(f"Error parsing {html_file_path} at lines {start_block_i} to {i}: {block_text}\n{e}")
                traceback.print_exc()
                continue

            block_counter += 1

        results.append({
            "motif": consensus_motif,
            "motif_size": len(consensus_motif),
            #"motif_size": period,
            #"consensus_motif_size": consensus_size,
            "num_repeats": copynumber,
            "repeat_sequence_length": seq_end_index - seq_start_index_0based,
            #"motifs": motifs,
            #"motifs_and_start_indices": motifs_and_start_indices,
            "total_mismatches": total_mismatches,
            "total_indels": total_indels,
            "score": score,
            "motif_positions_with_interruptions": motif_positions_with_interruptions,
        })

    return results



"""
    # Selecting the best TRF result for a locus
    result_with_max_score = max(results, key=lambda x: x["trf_score"])
    max_score = result_with_max_score["trf_score"]

    result_with_min_consensus_motif_size = min([t for t in results if t["trf_score"] > 0.8 * max_score], key=lambda x: x["trf_consensus_motif_size"])
"""



if __name__ == "__main__":
    #seq = "CAG"*2 + "CTG" + "CAG"*10 + "CTG"*2 + "CAG"*100
    seq = "CAG"*5 + "TAG" + "CTG" + "CTG" + "CTA" + "CAG"*5

    r = TRFRunner("~/bin/trf409.macosx")
    records = list(r.run_TRF_and_parse_motif_composition(
        nucleotide_sequence=seq,
        chromosome_name="chrN",
    ))

    print(f"Found {len(records)} TRF results in {seq}")
    for record in records:
        print(record)

