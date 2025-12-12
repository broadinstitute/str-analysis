"""Simple/naive method for taking a dictionary of motifs, and then parsing 1 or more nucleotide sequences to
extract a dictionary of counts showing how often the motifs in the dictionary were observed in the given
nucleotide sequences.
"""

import argparse
import collections
import gzip
import json
import os
import pandas as pd
import pysam
import re
from scipy.stats import binom

from str_analysis.utils.misc_utils import reverse_complement, parse_interval
from str_analysis.utils.cram_bam_utils import get_total_mapped_reads

MOTIF_INDEX_REGEXP = re.compile(r"[|][(\d+)][|]")
NUCLEOTIDE_SEQUENCE_REGEXP = re.compile(r"([ACGTRYMKSWHBVDN]+)", re.IGNORECASE)

SEPARATOR = "|"

MIN_MAPQ = 3


def _is_nucleotide_sequence(s):
    return "A" in s or "C" in s or "G" in s or "T" in s or "N" in s

def _parse_motif_ids_from_processed_sequence(sequence_with_motif_ids):
    return [int(s) for s in sequence_with_motif_ids.split(SEPARATOR) if s and not _is_nucleotide_sequence(s)]

class LocusParser:

    def __init__(self, reference_motif_frequency_dict):
        """Constructor.

        Args:
            reference_motif_frequency_dict (dict): dictionary that maps motifs to counts or relative frequencies.
                It represents the initial or expected set of frequencies that approximate what might be observed in the
                nucleotide sequences that will be parsed by this LocusParser.
        """
        if not reference_motif_frequency_dict:
            raise ValueError("reference_motif_frequency_dict not specified")

        self._reference_motif_frequency_dict = reference_motif_frequency_dict
        self._sorted_reference_motif_list = [
            (i, m) for i, m in
            enumerate(sorted(reference_motif_frequency_dict.keys(), key=lambda k: (reference_motif_frequency_dict[k], k), reverse=True))
        ]
        self._reference_motif_id_to_motif = dict(self._sorted_reference_motif_list)
        self._sorted_reference_motif_list = [
            (f"{SEPARATOR}{i}{SEPARATOR}", m) for i, m in self._sorted_reference_motif_list
        ]
        self._sorted_reference_motif_sizes = list(sorted({len(m) for m in reference_motif_frequency_dict}))
        self._min_reference_motif_size = self._sorted_reference_motif_sizes[0]
        self._reference_motif_sizes_set = list(sorted({len(m) for m in reference_motif_frequency_dict}))

        self._observed_motif_id_frequency_dict = collections.defaultdict(int)
        self._observed_motif_id_pair_frequency_dict = collections.defaultdict(int)
        self._observed_motif_id_triplet_frequency_dict = collections.defaultdict(int)
        self._novel_motif_frequency_dict = collections.defaultdict(int)

        # TODO use motif ids that are stable across different dictionaries?
        # For now, just use its position in the sorted motif list as the id.


    def _probability_of_result_by_chance(self, input_sequence, parsed_motif_ids):
        """Compute the approximate probability of seeing the parsed motif ids by chance in
        a sequence like the input sequence.

        Args:
             input_sequence (str): nucleotide sequence
             parsed_motif_ids (list): list of motif ids detected within the input_sequence

        Returns:
            float: the probability of chance occurance of these motifs in the given input sequence
        """

        # compute the relative frequency of each base in the input sequence
        total_bases = len(input_sequence)
        base_frequencies = {
            base: count/total_bases for base, count in collections.Counter(input_sequence).items()
        }

        total_bases_where_motif_could_occur = total_bases - self._min_reference_motif_size

        # compute the overall probability of observing any of the motifs
        probability_of_any_one_of_the_motifs_occurring = 0
        for motif_id in set(parsed_motif_ids):
            motif_sequence = self._reference_motif_id_to_motif[motif_id]
            motif_probability = 1
            for motif_base in motif_sequence:
                motif_probability *= base_frequencies.get(motif_base, 0.001)  # probability of base is 0.001 if not found
            probability_of_any_one_of_the_motifs_occurring += motif_probability

        total_motif_matches = len(parsed_motif_ids)
        pval = binom.sf(total_motif_matches - 1, total_bases_where_motif_could_occur, probability_of_any_one_of_the_motifs_occurring)

        return pval


    def convert_nucleotide_seq_to_motif_seq(
        self,
        sequence,
        check_reverse_complement=True,
        chance_occurrence_threshold=10**-3,
        record_reference_motif_counts=False,
        record_novel_motif_counts=False,
        record_motif_pair_counts=False,
        record_motif_triplet_counts=False,
    ):
        """Takes a nucleotide sequence consisting of A, C, G, T's and checks for the presence of each motif in the
        reference_motif_frequency_dict. If such motif(s) are found in the sequence, they are replaced by the string
        "|{motif id}|" while other nucleotides in the sequence are left unchanged.

        Args:
            sequence (str): nucleotide sequence consisting of A, C, G, T's
            check_reverse_complement (bool): if True, both the input sequence and its reverse complement will be checked
                for the presence of motifs from the reference_motif_frequency_dict. Whichever orientation contains
                the most motif occurrences will be considered the correct orientation. Checking the reverse complement is
                useful if the input sequence is from a short or long read since these are typically sequenced using an
                unstranded protocol where the orientation is arbitrary.
            chance_occurrence_threshold (float): if not None, this method will compute the probability that the identified
                motif occurrences were found by chance (taking into account the relative base pair frequencies within
                the input sequence and the sequence of the motifs in the reference_motif_frequency_dict). If
                this probability is higher than the threshold, the result will be considered spurious, and so the method
                will return None.
            record_reference_motif_counts (bool): if True, and motifs are found within the input sequence in a way
                that has a probability lower than the chance_occurrence_threshold, then the motifs detected in the
                input sequence will be counted in the internal observed_motif_frequency dictionary, allowing the
                later retrieval of observed motif frequencies across all parsed input sequences.
            record_novel_motif_counts (bool): if True, any nucleotide sequence fragments in the input sequence
                that remain at the ends, or between occurrences of motifs from the reference_motif_frequency_dict,
                and that have the same length as motifs in the reference_motif_frequency_dict, will be counted
                as novel motifs.
            record_motif_pair_counts (bool): if True, and consecutive repeats of motifs are found within the
                input sequence (such as 'AC|1||1|TG' or 'AC|1||0|TG') the pairs of consecutive motifs will
                be counted.
            record_motif_triplet_counts (bool): if True, and consecutive repeats of motifs are found within the
                input sequence (such as 'AC|1||1||0|TG'), these triplets of consecutive motifs will be counted.

        Return:
            str: nucleotide sequence with all known motifs replaced by "|{numerical motif id}|". For example:
                 "AC|0||0|AT|1||0|GT", meaning the input sequence contained 2 occurrences of the most frequent motif,
                 followed by 'AT', then 1 occurrence of the 2nd most frequent motif, then 1 more occurrence of the
                 most frequent motif. If no motifs were found in the input sequence or the probability of the motifs
                 occurring by chance is not lower than the chance_occurrence_threshold, None is returned.
        """
        sequence = sequence.upper()
        sequences_to_check = [sequence]
        if not set(sequence).issubset({"A", "C", "G", "T", "N"}):
            raise ValueError(f"The input sequence is not a valid nucleotide sequence: {sequence}")

        if check_reverse_complement:
            sequences_to_check.append(reverse_complement(sequence))

        best_remaining_nucleotide_count = None
        best_remaining_nucleotide_seqs = None
        best_parsed_seq = None
        best_original_seq = None
        for original_seq in sequences_to_check:
            # go through all motifs and convert the sequence
            seq = original_seq
            remaining_nucleotide_count = len(seq)
            remaining_nucleotide_seqs = [seq]
            maximum_parsing_achieved = False
            for motif_id, reference_motif in self._sorted_reference_motif_list:
                seq_before = seq
                seq = seq.replace(reference_motif, motif_id)
                if seq == seq_before:
                    # no motifs were replaced
                    continue

                # check how big the remaining nucleotide stretches are
                remaining_nucleotide_seqs = [s for s in seq.split(SEPARATOR) if len(s) > 0 and _is_nucleotide_sequence(s)]
                if not remaining_nucleotide_seqs:
                    remaining_nucleotide_count = 0
                    maximum_parsing_achieved = True
                    break  # done parsing

                remaining_nucleotide_seq_lengths = [len(s) for s in remaining_nucleotide_seqs]
                remaining_nucleotide_count = sum(remaining_nucleotide_seq_lengths)
                if max(remaining_nucleotide_seq_lengths) < self._min_reference_motif_size and len(remaining_nucleotide_seq_lengths) <= 2:
                    maximum_parsing_achieved = True
                    break  # done parsing


            if best_remaining_nucleotide_count is None or remaining_nucleotide_count < best_remaining_nucleotide_count:
                # chose whether to keep the results from the sequence or its reverse complement
                best_parsed_seq = seq
                best_original_seq = original_seq
                best_remaining_nucleotide_count = remaining_nucleotide_count
                best_remaining_nucleotide_seqs = remaining_nucleotide_seqs

            if maximum_parsing_achieved:
                chance_occurrence_threshold = None  # dont bother checking probability
                break

        if best_remaining_nucleotide_count == len(sequence):
            # no motifs found in sequence
            return None

        parsed_motif_ids = None
        if chance_occurrence_threshold or record_reference_motif_counts:
            parsed_motif_ids = _parse_motif_ids_from_processed_sequence(best_parsed_seq)

        if chance_occurrence_threshold and self._probability_of_result_by_chance(best_original_seq, parsed_motif_ids) >= chance_occurrence_threshold:
            # check probability of chance result
            return None

        if record_reference_motif_counts:
            for motif_id in parsed_motif_ids:
                self._observed_motif_id_frequency_dict[motif_id] += 1

        if record_novel_motif_counts and best_remaining_nucleotide_count < len(sequence):
            for remaining_seq in best_remaining_nucleotide_seqs:
                # if the remaining sequence is the same size as one of the motifs in the reference dict, add it as a new motif
                if len(remaining_seq) in self._reference_motif_sizes_set:
                    self._novel_motif_frequency_dict[remaining_seq] += 1

        if record_motif_pair_counts or record_motif_triplet_counts:
            previous_motif_ids = []
            for s in best_parsed_seq.split(SEPARATOR):
                if not s:
                    # skip empty strings between consecutive separators (eg. in "0||1") )
                    continue

                if _is_nucleotide_sequence(s):
                    # skip unparsed sequences and reset pair tracking
                    previous_motif_ids = []
                    continue

                # found parsed motif
                current_motif_id = int(s)
                if record_motif_pair_counts and len(previous_motif_ids) >= 1:
                    motif_id_pair = (previous_motif_ids[-1], current_motif_id)
                    self._observed_motif_id_pair_frequency_dict[motif_id_pair] += 1

                if record_motif_triplet_counts and len(previous_motif_ids) >= 2:
                    motif_id_triplet = (previous_motif_ids[-2], previous_motif_ids[-1], current_motif_id)
                    self._observed_motif_id_triplet_frequency_dict[motif_id_triplet] += 1
                
                previous_motif_ids.append(current_motif_id)

        return best_parsed_seq


    def process_sequences_in_alignment_file(self, alignment_file_path, genomic_intervals, include_low_quality_alignments=False, alignment_index_file_path=None):
        """Process sequences in an alignment file and count the number of occurrences of each motif.

        Args:
            alignment_file_path (str): path to a BAM/CRAM file
            genomic_intervals (list or str): genomic interval or list of intervals to extract from the alignment file
            include_low_quality_alignments (bool): whether to count low quality alignments (those with MAPQ < 3)

        Returns:
            dict: interval mapped to the average read depth of processed reads in that interval
        """

        if not alignment_file_path.startswith("gs://") and not alignment_file_path.startswith("s3://") and \
                not os.path.isfile(alignment_file_path):
            raise ValueError(f"Alignment file not found: {alignment_file_path}")
        if not alignment_file_path.endswith(".bam") and not alignment_file_path.endswith(".cram"):
            raise ValueError(f"Alignment file must end with .bam or .cram: {alignment_file_path}")
        if not genomic_intervals:
            raise ValueError("counted_region_list not specified")
        if isinstance(genomic_intervals, str):
            genomic_intervals = [genomic_intervals]

        interval_to_read_depth = collections.defaultdict(int)
        with pysam.AlignmentFile(alignment_file_path, index_filename=alignment_index_file_path) as f:
            for interval in genomic_intervals:
                chrom, start_0based, end_1based = parse_interval(interval)
                chrom = chrom.replace("chr", "")
                if f.references[0].startswith("chr"):
                    chrom = f"chr{chrom}"

                for read in f.fetch(chrom, start_0based, end_1based):
                    if not read.is_mapped or not read.query_alignment_sequence or (
                            not include_low_quality_alignments and read.mapq < MIN_MAPQ):
                        continue

                    # total number of aligned bases in this read that overlap the interval
                    bases_within_locus = min(end_1based, read.reference_end) - max(start_0based, read.reference_start)
                    interval_to_read_depth[interval] += bases_within_locus

                    self.convert_nucleotide_seq_to_motif_seq(
                        read.query_alignment_sequence,
                        check_reverse_complement=True,
                        record_reference_motif_counts=True,
                        record_novel_motif_counts=True,
                        record_motif_pair_counts=True,
                        record_motif_triplet_counts=True,
                    )

                interval_to_read_depth[interval] /= end_1based - start_0based   # divide by interval width

        return interval_to_read_depth

    def get_reference_motif_id_to_motif_dict(self):
        return self._reference_motif_id_to_motif

    def get_reference_motif_to_motif_id_dict(self):
        return {motif: motif_id for motif_id, motif in self._reference_motif_id_to_motif.items()}


    def get_observed_motif_frequency_dict(self, as_string=False):
        motif_frequency_dict = {
            self._reference_motif_id_to_motif[motif_id]: count for motif_id, count in self._observed_motif_id_frequency_dict.items()
        }
        if not as_string:
            return motif_frequency_dict

        items = []
        for motif, count in sorted(motif_frequency_dict.items(), key=lambda x: x[1], reverse=True):
            items.append(f"[{motif}]:{count}")
        return ",".join(items)


    def get_observed_motif_pair_frequency_dict(self, as_string=False):
        motif_pair_frequency_dict = {
            (self._reference_motif_id_to_motif[motif_id_pair[0]], self._reference_motif_id_to_motif[motif_id_pair[1]]): count
            for motif_id_pair, count in self._observed_motif_id_pair_frequency_dict.items()
        }
        if not as_string:
            return motif_pair_frequency_dict

        items = []
        for motif_pair, count in sorted(motif_pair_frequency_dict.items(), key=lambda x: x[1], reverse=True):
            items.append(f"[{motif_pair[0]}][{motif_pair[1]}]:{count}")
        return ",".join(items)


    def get_observed_motif_triplet_frequency_dict(self, as_string=False):
        motif_triplet_frequency_dict = {
            (
                self._reference_motif_id_to_motif[motif_id_triplet[0]],
                self._reference_motif_id_to_motif[motif_id_triplet[1]],
                self._reference_motif_id_to_motif[motif_id_triplet[2]],
            ): count for motif_id_triplet, count in self._observed_motif_id_triplet_frequency_dict.items()
        }
        if not as_string:
            return motif_triplet_frequency_dict

        items = []
        for motif_triplet, count in sorted(motif_triplet_frequency_dict.items(), key=lambda x: x[1], reverse=True):
            items.append(f"[{motif_triplet[0]}][{motif_triplet[1]}][{motif_triplet[2]}]:{count}")
        return ",".join(items)


    def get_novel_motif_frequency_dict(self, as_string=False):
        for motif in self._novel_motif_frequency_dict:
            if motif in self._reference_motif_frequency_dict:
                raise Exception(f"Motif {motif} is not novel: {self._reference_motif_frequency_dict}")

        if not as_string:
            return self._novel_motif_frequency_dict

        items = []
        for motif, count in sorted(self._novel_motif_frequency_dict.items(), key=lambda x: x[1], reverse=True):
            items.append(f"[{motif}]:{count}")
        return ",".join(items)


def parse_motif_composition_from_alignment_file(
    input_sequence_or_path,
    motif_frequency_dict,
    counted_region_list=None,
    other_region_list=None,
    include_low_quality_alignments=False,
    check_reverse_complement=True,
    output_json=None,
    verbose=False,
):
    """Returns a dictionary with results. See below for more details.
        
    Args:
        input_sequence_or_path (str): path to a BAM/CRAM file or a nucleotide sequence
        motif_frequency_dict (dict): Dictionary of population motif frequencies at this locus (such as from T2T assemblies)
        counted_region_list (list): list of genomic intervals to extract from the input BAM/CRAM file
        other_region_list (list): list of genomic intervals to extract from the input BAM/CRAM file for control region normalization
        include_low_quality_alignments (bool): whether to count low quality alignments (those with MAPQ < 3)
        check_reverse_complement (bool): whether to check the reverse complement of the input sequence
        output_json (str): path to a JSON file to write the results to
        verbose (bool): whether to print verbose output
    
    Returns:
        dict: a dictionary with the following fields:

        "input": BAM/CRAM file path or nucleotide input sequence
        "input_file_size": size of the input file in bytes or None if the input is a nucleotide sequence
        "total_mapped_reads": total number of mapped reads in the input BAM/CRAM file or None if the input is a nucleotide sequence
        "motif_frequency": a string representation of the motif frequencies observed in input reads or nucleotide sequence
        "novel_motif_frequency": a string representation of novel motif frequencies (ie. motifs not in the motif_frequency_table) 
            observed in input reads or nucleotide sequence
        "motif_pair_frequency": a string representation of motif pair frequencies (ie. consecutive motifs) observed in input reads or nucleotide sequence
        "motif_triplet_frequency": a string representation of motif triplet frequencies (ie. consecutive triplets of motifs) observed in input reads or nucleotide sequence
        "read_depth_counted_region": a dictionary with the read depth for each primary region specified by the counted_region_list
        "read_depth_other_region": a dictionary with the read depth for each control region specified by the other_region_list                
    """

    locus_parser = LocusParser(motif_frequency_dict)

    interval_read_depth_dict = collections.defaultdict(int)
    if os.path.isfile(input_sequence_or_path):
        input_is_file = True
        input_file = pysam.AlignmentFile(input_sequence_or_path)

        for interval in counted_region_list:
            interval_key = f"read_depth_counted_region_{interval}"
            chrom, start_0based, end_1based = parse_interval(interval)
            locus_width = end_1based - start_0based
            read_iterator = input_file.fetch(chrom, start_0based, end_1based)
            for read in read_iterator:
                if not read.is_mapped or not read.query_alignment_sequence or (
                        not include_low_quality_alignments and read.mapq < MIN_MAPQ):
                    continue

                # compute the total number of aligned bases in this read that overlap the interval
                read_bases_aligned_within_interval = min(end_1based, read.reference_end) - max(start_0based, read.reference_start)
                interval_read_depth_dict[interval_key] += read_bases_aligned_within_interval

                if read_bases_aligned_within_interval < read.reference_end - read.reference_start:
                    # don't count motifs in a read unless all of its aligned bases are within the locus interval
                    continue

                locus_parser.convert_nucleotide_seq_to_motif_seq(
                    read.query_sequence,   # process the entire sequence include soft-clips (eg. RFC1)
                    check_reverse_complement=True,
                    record_reference_motif_counts=True,
                    record_novel_motif_counts=True,
                    record_motif_pair_counts=True,
                    record_motif_triplet_counts=True,
                    chance_occurrence_threshold=10**-3,
                )

            interval_read_depth_dict[interval_key] /= locus_width

        if other_region_list:
            for interval in other_region_list:
                interval_key = f"read_depth_other_region_{interval}"

                chrom, start_0based, end_1based = parse_interval(interval)
                chrom = chrom.replace("chr", "")
                if input_file.references[0].startswith("chr"):
                    chrom = f"chr{chrom}"
                locus_width = end_1based - start_0based
                read_iterator = input_file.fetch(chrom, start_0based, end_1based)
                for read in read_iterator:
                    if not read.is_mapped or read.mapq < MIN_MAPQ:
                        continue

                    bases_within_locus = min(end_1based, read.reference_end) - max(start_0based, read.reference_start)  # total number of aligned bases in this read that overlap the interval
                    interval_read_depth_dict[interval_key] += bases_within_locus

                interval_read_depth_dict[interval_key] /= locus_width
    else:
        input_is_file = False
        unexpected_characters = set(input_sequence_or_path.upper()) - set("ATCGN")
        if unexpected_characters:
            raise ValueError("The input_sequence arg must be a nucleotide sequence consisting of A, C, G, T, and N. Invalid characters detected: " + ", ".join(unexpected_characters))

        parsed_sequence = locus_parser.convert_nucleotide_seq_to_motif_seq(
            input_sequence_or_path,
            check_reverse_complement=check_reverse_complement,
            record_reference_motif_counts=True,
            record_novel_motif_counts=True,
            record_motif_pair_counts=True,
            record_motif_triplet_counts=True,
            chance_occurrence_threshold=10**-3,
        )


    output_record = {
        "input":                    input_sequence_or_path,
        "input_file_size":          os.path.getsize(input_sequence_or_path) if input_is_file else None,
        "total_mapped_reads":       get_total_mapped_reads(input_sequence_or_path) if input_is_file else None,
        "motif_frequency":          locus_parser.get_observed_motif_frequency_dict(as_string=True),
        "motif_pair_frequency":     locus_parser.get_observed_motif_pair_frequency_dict(as_string=True),
        "motif_triplet_frequency":  locus_parser.get_observed_motif_triplet_frequency_dict(as_string=True),
        "novel_motif_frequency":    locus_parser.get_novel_motif_frequency_dict(as_string=True),
    }
    output_record.update(interval_read_depth_dict)
    if include_low_quality_alignments:
        output_record["include_low_quality_alignments"] = True
    else:
        output_record["min_mapq"] = MIN_MAPQ
        output_record["include_low_quality_alignments"] = False

    if output_json:
        if not output_json.endswith("gz"):
            output_json += ".gz"
        with gzip.open(output_json, "wt") as f:
            json.dump(output_record, f, indent=2)
        print(f"Wrote results to {output_json}")
    else:
        print("Done")

    if verbose:
        if not input_is_file:
            print("\nParsed motif sequence:")
            print(parsed_sequence)

        print("\nObserved motif frequencies:")
        for i, (motif, count) in enumerate(sorted(locus_parser.get_observed_motif_frequency_dict().items(), key=lambda x: x[1], reverse=True)):
            print(motif, count)

        print("\nObserved motif pair frequencies:")
        for i, (motif_pair, count) in enumerate(sorted(locus_parser.get_observed_motif_pair_frequency_dict().items(), key=lambda t: t[1], reverse=True)):
            print(f"{str(motif_pair):>10s}", count)

        print("\nNovel motif frequencies:")
        for i, (motif, count) in enumerate(sorted(locus_parser.get_novel_motif_frequency_dict().items(), key=lambda x: x[1], reverse=True)):
            print(motif, count)

    return output_record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--motif-frequency-table",
        help="A TSV table with two columns: 'motif' and 'number_of_times_seen_across_all_haplotypes'. "
             "The table should contain one or more expected motifs and their relative occurance counts in some "
             "reference set of sequences (such as haplotypes from the HPRC assemblies)",
        default="https://storage.googleapis.com/tandem-repeat-explorer/other/CACNA1C_from_dipcall.motifs.txt")

    parser.add_argument("--check-reverse-complement", action="store_true", help="Whether to also look for motif "
        "occurrences in the reverse complement of the input sequence. This is true by default (and therefore redundant) "
        "if the input is a BAM/CRAM path, but false by default if the input is a nucleotide sequence.")
    parser.add_argument("-L", "--counted-region", action="append", help="Genomic interval(s) to extract from the input BAM/CRAM. "
        "It should be specified in the format 'chrN:start-end' (0-based coordinates). This arg must be specified if the "
        "input_sequence is a BAM/CRAM file path. Example: chr1:12345-54321")
    parser.add_argument("--include-low-quality-alignments", action="store_true", help=f"Whether to count low quality "
                        "alignments (those with MAPQ < {MIN_MAPQ}).")
    parser.add_argument("--verbose", action="store_true", help="Print additional logging messages.")
    parser.add_argument("--output-json", help="If specified, output the results to this JSON file path.")
    parser.add_argument("--other-region", action="append", help="Optional genomic interval(s) from which to extract read depth for normalizing counts later. "
        "It should be specified in the format 'chrN:start-end' (0-based coordinates). This arg can only be used if the "
        "input_sequence is a BAM/CRAM file path. Example: chr1:12345-54321")
    parser.add_argument("input_sequence", help="Either a nucleotide sequence or the path of a BAM/CRAM file to parse")
    args = parser.parse_args()

    if os.path.isfile(args.input_sequence):
        if not args.counted_region:
            parser.error("Must specify --counted-region when the input is a BAM/CRAM file")
    else:
        unexpected_characters = set(args.input_sequence.upper()) - set("ATCGN")
        if unexpected_characters:
            if len(unexpected_characters) > 5:
                parser.error(f"File not found: {args.input_sequence}")
            else:
                parser.error("The input_sequence arg must be a nucleotide sequence consisting of A, C, G, T, and N. "
                             "Invalid characters detected: " + ", ".join(unexpected_characters))

        if args.counted_region:
            parser.error("The --counted-region arg is only supported if the input_sequence arg is a BAM/CRAM file")

    df = pd.read_table(args.motif_frequency_table)
    df["motif"] = df["motif"].str.strip().str.upper()
    motif_frequency_dict = dict(zip(df["motif"], df["number_of_times_seen_across_all_haplotypes"]))

    parse_motif_composition_from_alignment_file(
        input_sequence_or_path=args.input_sequence,
        motif_frequency_dict=motif_frequency_dict,
        counted_region_list=args.counted_region,
        other_region_list=args.other_region,
        include_low_quality_alignments=args.include_low_quality_alignments,
        check_reverse_complement=args.check_reverse_complement,
        output_json=args.output_json,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
