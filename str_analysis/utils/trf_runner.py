"""Utils for running TandemRepeatsFinder [Benson 1999] and parsing its output"""

import os
import re
import subprocess
import tempfile
import traceback

from str_analysis.utils.dat_utils import parse_dat_file


TRF_ALIGNMENT_LINE_REGEX = re.compile(r"(\d+)\s+([-ACGTN\s]+)")

class TRFRunner:
    """Class for running TandemRepeatFinder on a nucleotide sequence"""

    def __init__(self,
                 trf_executable_path,
                 match_score = 2,
                 mismatch_penalty = 7,
                 indel_penalty = 7,
                 pm = 80,
                 pi = 10,
                 minscore = 24, 
                 output_filename_prefix=None,
                 parse_motif_composition=False,
                 generate_motif_plot=False,
                 debug=False):
        """
        Args:
            trf_executable_path (str): path to the TRF executable
            match_score (int): see Tandem Repeats Finder docs
            mismatch_penalty (int): see Tandem Repeats Finder docs
            indel_penalty (int): see Tandem Repeats Finder docs
            pm (int): see Tandem Repeats Finder docs
            pi (int): see Tandem Repeats Finder docs
            minscore (int): see Tandem Repeats Finder docs
            output_filename_prefix (str): if specified, the TRF output files will have this prefix and will not be deleted.
                If not specified, the intermediate files will be deleted after they are parsed into in-memory data structures.
            parse_motif_composition (bool): if True, in addition to parsing the repeat coordinates and motif, also
                parse and return the motif composition from TRF's HTML output.
            generate_motif_plot (bool): if True, generate a logo (aka. motif) plot of the variability among different motifs observed in the repeat sequence.
            debug (bool): if True, print debug information
        """

        self.trf_executable_path = os.path.expanduser(trf_executable_path)
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.indel_penalty = indel_penalty
        self.pm = pm
        self.pi = pi
        self.minscore = minscore
        self.output_filename_prefix = output_filename_prefix
        self.parse_motif_composition = parse_motif_composition
        self.generate_motif_plot = generate_motif_plot
        self.debug = debug

        if self.generate_motif_plot and not self.parse_motif_composition:
            raise ValueError("generate_motif_plot requires parse_motif_composition to be True")


    def run_TRF(self, nucleotide_sequence, min_motif_size=None, max_motif_size=None):
        """Run TRF on the given nucleotide sequence and return a generator of records representing the TRF output for the given nucleotide sequence.
        This code assumes the goal is to find a single repeat that covers all or most of the input sequence.
        It is not designed to find multiple distict repeat loci within the sequence.

        Args:
            nucleotide_sequence (str): the nucleotide sequence to run TRF on
            min_motif_size (int): if specified, only return loci with this motif size or larger.
            max_motif_size (int): if specified, only return loci with this motif size or smaller.

        Yields:
            records representing the TRF output for the given nucleotide sequence.
        """

        # write the sequence to a temp FASTA file
        if self.output_filename_prefix is None:
            temp_fasta_file = tempfile.NamedTemporaryFile("w+", suffix=".fasta", delete=False)
        else:
            temp_fasta_file = open(f"{self.output_filename_prefix}.fasta", "wt")

        temp_fasta_file.write(f">seq_{len(nucleotide_sequence)}bp\n")
        temp_fasta_file.write(f"{nucleotide_sequence}\n")

        temp_fasta_file.close()

        if self.debug:
            print(f"Wrote {len(nucleotide_sequence):,d}bp sequence to {temp_fasta_file.name}")

        temp_fasta_path = temp_fasta_file.name

        keep_intermediate_files = self.output_filename_prefix is not None or self.debug
        
        # compose the TRF command
        # -l 6  = maximum TR length expected (in millions) (eg, -l 3 or -l=3 for 3 million). Human genome HG38 would need -l 6
        # -h = suppress html output
        # -ngs =  more compact .dat output on multisequence files, returns 0 on success. Output is printed to the screen.
        
        max_period = 2000  # maximum period size to report. Must be between 1 and 2000, inclusive
        command = f"{self.trf_executable_path} "
        command += f"{temp_fasta_path} "
        command += f"{self.match_score} "
        command += f"{self.mismatch_penalty} "
        command += f"{self.indel_penalty} "
        command += f"{self.pm} "
        command += f"{self.pi} "
        command += f"{self.minscore} "
        command += f"{max_period} -l 6 "

        redirect = subprocess.DEVNULL if not self.debug else None

        if self.parse_motif_composition:
            # run TRF
            subprocess.run(command, shell=True, stdout=redirect, stderr=redirect, check=False)

            # read html output
            trf_output_filename_prefix = f"{os.path.basename(temp_fasta_path)}"
            trf_output_filename_prefix += f".{self.match_score}"
            trf_output_filename_prefix += f".{self.mismatch_penalty}"
            trf_output_filename_prefix += f".{self.indel_penalty}"
            trf_output_filename_prefix += f".{self.pm}"
            trf_output_filename_prefix += f".{self.pi}"
            trf_output_filename_prefix += f".{self.minscore}"
            trf_output_filename_prefix += f".{max_period}"

            i = 1
            html_file_path = f"{trf_output_filename_prefix}.{i}.txt.html"
            while os.path.isfile(html_file_path):

                for record in self._parse_trf_html_output(nucleotide_sequence, html_file_path, min_motif_size=min_motif_size, max_motif_size=max_motif_size):
                    if self.generate_motif_plot:
                        self.create_motif_plot(record["repeats"], record["start_0based"], record["end_1based"], record["motif"])

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
            # run TRF
            output_dat_path = f"{temp_fasta_path}.dat"
            command += f"-h -ngs > {output_dat_path}"

            subprocess.check_output(command, shell=True, stderr=redirect)

            # parse output in .dat format
            for dat_record in parse_dat_file(output_dat_path):
                if min_motif_size is not None and dat_record.repeat_unit_length < min_motif_size:
                    continue
                if max_motif_size is not None and dat_record.repeat_unit_length > max_motif_size:
                    continue

                yield dat_record

            if not keep_intermediate_files:
                os.remove(output_dat_path)


    def _parse_trf_html_output(self, nucleotide_sequence, html_file_path, min_motif_size=None, max_motif_size=None):
        """Parse the HTML output of TRF and return a list of motifs"""

        with open(html_file_path, "rt") as f:
            html_content = f.read()

        results = []
        for locus_i, html_section in enumerate(html_content.split("<A NAME=")[1:]):
            #print(f"Parsing html section #{locus_i+1}")

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

            #period = int(period_copynumber_consensus_size_match.group(1))
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
            final_consensus_motif = None
            repeat_sequence = None
            repeats = []
            #repeats_and_start_indices = []
            motif_positions_with_interruptions = {}  # histogram of motif positions with interruptions
            total_mismatches = 0
            total_indels = 0

            line_i = 0
            while line_i < len(lines):
                line_i_block_start = line_i
                while line_i < len(lines) and lines[line_i].strip():
                    # find the end of this the block
                    line_i += 1

                lines_in_block = lines[line_i_block_start:line_i]
                line_i += 1

                if len(lines_in_block) < 2 or not any(l.strip() for l in lines_in_block):
                    continue

                consensus_motif_line = lines_in_block[-1]
                consensus_motif_match = TRF_ALIGNMENT_LINE_REGEX.search(consensus_motif_line)
                if not consensus_motif_match:
                    raise ValueError(f"Could not parse consensus motif from line: {consensus_motif_line}")

                # this is a continuation of the motif and alignment from the previous block
                is_continuation = consensus_motif_match.group(1) != "1"                    

                current_consensus_motif = consensus_motif_match.group(2).split(" ")[0].replace("-", "")
                if final_consensus_motif is None:
                    if consensus_motif is None:
                        assert not is_continuation, f"First block should not be a continuation in block: \n" + "\n".join(lines_in_block)
                        consensus_motif = current_consensus_motif
                    elif is_continuation:
                        consensus_motif += current_consensus_motif
                    else: 
                        final_consensus_motif = consensus_motif
                        consensus_motif = current_consensus_motif
                else:
                    if consensus_motif is None:
                        assert not is_continuation, f"First block should not be a continuation in block: \n" + "\n".join(lines_in_block)
                        consensus_motif = current_consensus_motif
                    elif is_continuation:
                        consensus_motif += current_consensus_motif
                    else: 
                        if consensus_motif != final_consensus_motif:
                            print(f"WARNING: {html_file_path} consensus motif mismatch: {consensus_motif} != {final_consensus_motif} in block: \n" + "\n".join(lines_in_block))
                        consensus_motif = current_consensus_motif


                repeat_sequence_line = lines_in_block[-2]
                repeat_sequence_match = TRF_ALIGNMENT_LINE_REGEX.search(repeat_sequence_line)
                if not repeat_sequence_match:
                    raise ValueError(f"Could not parse motif from line: {repeat_sequence_line}")


                total_indels += repeat_sequence_match.group(2).count("-") + consensus_motif_match.group(2).count("-")
                if not is_continuation:
                    repeat_sequence = repeat_sequence_match.group(2)
                else:
                    repeat_sequence += repeat_sequence_match.group(2)


                interruptions_string = None
                if len(lines_in_block) > 2:
                    interruptions_line = lines_in_block[-3]

                    interruptions_line_offset, _ = repeat_sequence_match.span(2)  # offset of the repeat sequence in the interruptions line
                    interruptions_string = interruptions_line[interruptions_line_offset:]

                    total_mismatches += interruptions_string.count("*")
                    if self.debug:
                        print()
                        print(interruptions_string)
                        print(repeat_sequence_match.group(2))
                        print(consensus_motif_match.group(2))

                interruptions_string_offset = 0
                for repeat_motif in repeat_sequence_match.group(2).split(" "):
                    repeats.append(repeat_motif)
                    #repeats_and_start_indices.append((repeat_sequence_start_index_0based, motif))
                    #repeat_sequence_start_index_0based += len(motif)

                    if interruptions_string is not None:
                        motif_interruptions = interruptions_string[interruptions_string_offset : interruptions_string_offset + len(repeat_motif)]
                        if self.debug:
                            print(f"motif_interruptions: {motif_interruptions}")
                            print(f"motif              : {repeat_motif}")

                        if "*" in motif_interruptions:
                            pos_offset = 0 if not is_continuation else int(consensus_motif_match.group(1)) - 1
                            for i, c in enumerate(motif_interruptions):
                                if c == "*":
                                    motif_positions_with_interruptions[i + 1 + pos_offset] = motif_positions_with_interruptions.get(i + 1 + pos_offset, 0) + 1

                    interruptions_string_offset += len(repeat_motif) + 1  # add +1 to account for the space between motifs

            if final_consensus_motif is None:
                final_consensus_motif = consensus_motif

            results.append({
                "start_0based": seq_start_index_0based,
                "end_1based": seq_end_index,
                "motif": final_consensus_motif,
                "motif_size": len(final_consensus_motif),
                #"motif_size": period,
                #"consensus_motif_size": consensus_size,
                "num_repeats": copynumber,
                "repeats": repeats,
                #"repeats_and_start_indices": repeats_and_start_indices,
                "total_mismatches": total_mismatches,
                "total_indels": total_indels,
                "alignment_length": seq_end_index - seq_start_index_0based,
                "alignment_score": score,
                "motif_positions_with_interruptions": motif_positions_with_interruptions,
            })

        
        # Filter out the short TRF results at a locus
        #metric = "alignment_length"
        result_with_max_alignment_score = max(results, key=lambda x: x["alignment_score"])
        max_alignment_score = result_with_max_alignment_score["alignment_score"]

        results = [t for t in results if t["alignment_score"] > 0.1 * max_alignment_score]
        #results = [t for t in results if t["start_0based"] < t["motif_size"] and t["end_1based"] > len(nucleotide_sequence) - t["motif_size"]]

        results.sort(key=lambda x: x["motif_size"])
        return results

    def create_motif_plot(self, motif_list, start_0based, end_1based, consensus_motif):
        if not motif_list:
            raise ValueError("No motifs to plot")

        filtered_motif_list = []
        for motif in motif_list:
            if len(motif) == len(consensus_motif):
                filtered_motif_list.append(motif)
            elif len(motif.replace("-", "")) == len(consensus_motif):
                # if the motif has deletions, we can still plot it
                filtered_motif_list.append(motif.replace("-", ""))
            else:
                print(f"WARNING: {motif} is not the same length as the consensus motif {consensus_motif}")

        # requires plotnine and plotnineseqsuite to be installed
        from plotnine import ggplot
        from plotnineseqsuite import geom_logo, theme_seq

        p = ggplot() + geom_logo(filtered_motif_list, method = 'probability' ) + theme_seq()
        p.save(f'{self.output_filename_prefix}_motif_plot_{start_0based}_{end_1based}_{consensus_motif}.png', "png", width=len(consensus_motif)*0.5, height=3, limitsize=False)




DEBUG=False

if __name__ == "__main__":
    #seq = "CAG"*2 + "CTG" + "CAG"*10 + "CTG"*2 + "CAG"*100
    import sys
    import pandas as pd
    df1 = pd.read_table("~/code/str-analysis/str_analysis/data/tests/CACNA1C_VNTR_sequences.tsv")
    df2 = pd.read_table("~/code/str-analysis/str_analysis/data/tests/ABCA7_VNTR_sequences.tsv")
    for parse_motif_composition in [True,]:

        #for df, column_name in [(df1, "CACNA1C_seq"),  (df2, "ABCA7_seq")]:
        #for df, column_name in [(df2, "ABCA7_seq")]:
        for df, column_name in [(df1, "CACNA1C_seq")]:
            for _, row in df.iterrows():
                #if row['sample_id'] != "NA18939":
                #    continue

                r = TRFRunner("~/bin/trf409.macosx", 
                              parse_motif_composition=parse_motif_composition, 
                              debug=DEBUG, 
                              output_filename_prefix=f"{row['sample_id']}_{column_name}",
                              generate_motif_plot=False,
                )

                # run TRF on the sequence
                seq = row[column_name]
                #print(f"Running TRF on {len(seq):,d}bp sequence from {column_name}")
                records = list(r.run_TRF(nucleotide_sequence=seq))
                records = [records[0]]
                sys.stdout.write(f"{row['sample_id']:<20} {column_name}: found {len(records):,d} TRF results in {len(seq):10,d}bp sequence")
                for trf_record_i, record in enumerate(records):
                    print(f" "*30 + f"TRF record #{trf_record_i + 1}: score = {record['alignment_score']:6,d}, start_diff = {record['start_0based']}, end_diff = {len(seq) - record['end_1based']}, cov = {(record['end_1based'] - record['start_0based'])/len(seq):0.1%}, motif = {record['motif']}, motif_size = {record['motif_size']}")
                    if DEBUG:
                        print(record)

