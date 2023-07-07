#!/usr/bin/env python3

import argparse
import logging
import math
import os

import pyfaidx

from str_analysis.utils.bam_utils import (compute_bam_stats, merge_bams,
                                          simulate_reads)
from str_analysis.utils.fasta_utils import get_reference_sequence
from str_analysis.utils.misc_utils import parse_interval, run

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')


def main():
    """Simulate STR expansions."""

    p = argparse.ArgumentParser(description="This script simulates STR expansions (either heterozygous or bi-allelic). "
        "It first creates a small synthetic reference sequence that contains the expanded STR + flanking sequence. "
        "Then it runs wgsim to simulate paired-end short reads from this reference, and then runs bwa to create a "
        "synthetic .bam file. This .bam file can then be passed to downstream STR calling tools to evaulate their "
        "accuracy. The main inputs to this script are: the human reference genome fasta, the STR locus coordinates, "
        "the repeat unit motif and number of copies to insert at this locus. "
        "This script requires the following programs to be installed and on PATH: bedtools, samtools, bwa, wgsim")
    p.add_argument("-R", "--ref-fasta", required=True, help="reference fasta file path. It should also have a "
                                                            "bwa index for running bwa mem.")


    p.add_argument("-p", "--padding", type=int, default=750, help="How many bases on either side of the reference "
        "repeat region to include in the synthetic reference from which reads are simulated")

    p.add_argument("-u1", "--new-repeat-unit1", help="The repeat unit to use for the short allele (aka. allele #1). "
        "If not specified, --new-repeat-unit2 will be used")
    p.add_argument("-n1", "--num-copies1", type=int, help="How many copies of --new-repeat-unit1 to insert for "
        "allele 1. Defaults to the number of copies in the reference genome")

    p.add_argument("--het", action="store_true", help="Simulate an individual that is heterozygous for "
                                                      "--num-copies of --new-repeat-unit")
    p.add_argument("--hom", action="store_true", help="Simulate an individual that is homozygous for "
                                                      "--num-copies of --new-repeat-unit")

    p.add_argument("-n2", "--num-copies2", type=int, help="How many copies of --new-repeat-unit2 to insert for "
        "allele 2")

    p.add_argument("--new-repeat-unit2", required=True, help="The repeat unit to use for the long allele "
        "(aka. allele #2)")

    p.add_argument("-t", "--num-copies2-max", type=int, help="If specified, multiple bams will be "
        "generated - one for each copy number between -n2 and -t, with step size -i.")
    p.add_argument("-i", "--num-copies2-increment", type=int, default=5, help="If -t is specified, multiple bams will "
        "be generated - one for each copy number between -n and -t, with step size -i.")
    p.add_argument("-l", "--num-copies2-list", type=str, help="As an alternative to -n2, -t and -i, this sets a "
                                                             "comma-seperated list of copy numbers to simulate.")
    p.add_argument("-f", "--force", action="store_true", help="Generated simulated .bam even if it already exists.")
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose output.")

    p.add_argument("--output-off-target-regions", action="store_true", help="Computes off-target regions and outputs "
        "them to a repeat_regions.bed output file.")
    p.add_argument("--min-off-target-reads", type=int, default=5, help="Minimum number of reads that has to map to "
        "a particular off-target region for it to be included in the output. Only relevant when "
        "--output-off-target-regions is used.")

    grp = p.add_argument_group("wgsim params")
    grp.add_argument("--coverage", type=float, help="Read death of simulated bams", default=30)
    grp.add_argument("--read-length", type=int, help="Read length of simulated reads", default=150)
    grp.add_argument("--fragment-length", type=float, help="Fragment length of simulated reads", default=345)
    grp.add_argument("--fragment-length-stddev", type=float, help="Fragment length standard deviation of simulated reads",
                   default=100)

    grp.add_argument("-e", "--wgsim-base-error-rate", type=float, help="wgsim -e arg (base error rate [0.020])")
    grp.add_argument("--wgsim-mutation-rate", type=float, help="wgsim -r arg (rate of mutations [0.0010])")
    grp.add_argument("--wgsim-fraction-indels", type=float, help="wgsim -R arg (fraction of indels [0.15])")
    grp.add_argument("--wgsim-p-indel-extension", type=float, help="wgsim -X arg (probability an indel is "
                                                                 "extended [0.30])")

    p.add_argument("--output-label", help="If specified, this label will be added to output filenames")
    p.add_argument("ref_repeat_coords", help="1-based coordinates in the reference genome (eg. chr1:12345-54321). "
        "The reference bases in this interval will be replaced with -n copies of the --new-repeat-unit sequence.")

    args = p.parse_args()

    # parse args
    het_or_hom_list = []
    if args.het or (not args.het and not args.hom):
        het_or_hom_list.append("het")
    if args.hom or (not args.het and not args.hom):
        het_or_hom_list.append("hom")

    # generate synthetic reference sequence
    chrom, start_1based, end_1based = parse_interval(args.ref_repeat_coords)
    padding = args.padding

    if args.new_repeat_unit1 is None:
        args.new_repeat_unit1 = args.new_repeat_unit2

    args.coverage = int(args.coverage)

    # validate paths
    if args.ref_fasta and not os.path.exists(args.ref_fasta):
        p.error(f"Path not found: {args.ref_fasta}")

    fasta_obj = pyfaidx.Fasta(args.ref_fasta, one_based_attributes=False, as_raw=True)

    output_filename_prefix = (
        (f"{args.output_label}_" if args.output_label else "") +
        f"{chrom}-{start_1based}-{end_1based}__rl{args.read_length}__{padding}bp_pad"
    ).replace(" ", "_")

    if args.wgsim_base_error_rate:
        output_filename_prefix += f"__wgsim_e-{args.wgsim_base_error_rate}"

    use_reference_allele = args.num_copies1 is None and args.new_repeat_unit1 is None

    if args.num_copies1 is None:
        args.num_copies1 = int(math.ceil((end_1based - start_1based + 1)/len(args.new_repeat_unit1 or args.new_repeat_unit2)))

    # if args.het, generate ALT allele1. This allele is the ref allele, or if num_copies1, then new_repeat_unit1 allele.
    if "het" in het_or_hom_list:
        if use_reference_allele:
            # calculate number of copies in reference
            logging.info(f"Getting allele1 from reference ({args.num_copies1}x{args.new_repeat_unit1})")
            synthetic_alt_allele1_reference_sequence = get_reference_sequence(
                fasta_obj,
                chrom,
                start_1based - padding,
                end_1based + padding)
        else:
            logging.info(f"Simulating allele1 reference with {args.num_copies1}x{args.new_repeat_unit1} repeats")
            synthetic_alt_allele1_reference_sequence = generate_synthetic_reference_sequence(
                fasta_obj,
                chrom,
                start_1based,
                end_1based,
                padding,
                args.new_repeat_unit1,
                args.num_copies1)

        # in HET case, divide coverage by 2 to generate ref and alt .bams with 1/2 original coverage, and then merge them.
        coverage = args.coverage / 2
        logging.info(f"Generating allele1 bam with {int(coverage)}x coverage")
        synthetic_alt_allele1_bam_path = simulate_reads(
            args.ref_fasta,
            synthetic_alt_allele1_reference_sequence,
            args.read_length,
            coverage,
            args.fragment_length,
            args.fragment_length_stddev,
            f"{output_filename_prefix}__{int(coverage)}x_cov__{args.num_copies1:3d}x{args.new_repeat_unit1}".replace(" ", "_"),
            wgsim_base_error_rate=args.wgsim_base_error_rate,
            wgsim_mutation_rate=args.wgsim_mutation_rate,
            wgsim_fraction_indels=args.wgsim_fraction_indels,
            wgsim_p_indel_extension=args.wgsim_p_indel_extension,
            force=args.force)

    # simulate ALT allele(s) 2, with het and/or hom genotype
    if args.num_copies2_list:
        num_copies2_list = [int(num_copies2) for num_copies2 in args.num_copies2_list.split(",")]
    elif args.num_copies2 and args.num_copies2_max:
        num_copies2_list = range(args.num_copies2, args.num_copies2_max+1, args.num_copies2_increment)
    elif args.num_copies2:
        num_copies2_list = [args.num_copies2]
    else:
        p.error("Must specify --num-copies or --num-copies-list")

    metadata_tsv_rows = []
    for het_or_hom in het_or_hom_list:
        simulated_bam_paths = []
        for num_copies2 in num_copies2_list:
            logging.info("-"*100)
            logging.info(f"Simulating allele2 reference with {num_copies2}x{args.new_repeat_unit2} repeats")
            synthetic_alt_allele2_reference_sequence = generate_synthetic_reference_sequence(
                fasta_obj,
                chrom,
                start_1based,
                end_1based,
                padding,
                args.new_repeat_unit2,
                num_copies2)

            # in HET case, divide coverage by 2 to generate ref and alt .bams with 1/2 original coverage, and then merge them.
            coverage = args.coverage / (2 if het_or_hom == "het" else 1)
            output_prefix = f"{output_filename_prefix}__{int(coverage)}x_cov__{num_copies2:3d}x{args.new_repeat_unit2}"
            output_prefix = output_prefix.replace(" ", "_")
            if het_or_hom == "hom":
                output_prefix += "__hom"
            logging.info(f"Generating allele2 bam with {int(coverage)}x coverage")
            synthetic_alt_allele2_bam_path = simulate_reads(
                args.ref_fasta,
                synthetic_alt_allele2_reference_sequence,
                args.read_length,
                coverage,
                args.fragment_length,
                args.fragment_length_stddev,
                output_prefix,
                force=args.force)

            if het_or_hom == "het":
                logging.info(f"Merging bams from allele1 ({args.num_copies1:3d}x{args.new_repeat_unit1}) "
                             f"and allele2 ({num_copies2:3d}x{args.new_repeat_unit2})")
                merged_bam_path = output_filename_prefix
                merged_bam_path += f"__{int(args.coverage)}x_cov"
                merged_bam_path += f"__allele1_{args.num_copies1:3d}x{args.new_repeat_unit1}"
                merged_bam_path += f"__allele2_{num_copies2:3d}x{args.new_repeat_unit2}"
                merged_bam_path += f"__het.bam"
                merged_bam_path = merged_bam_path.replace(" ", "_")
                merge_bams(
                    merged_bam_path,
                    synthetic_alt_allele1_bam_path,
                    synthetic_alt_allele2_bam_path,
                    force=args.force)

                #run(f"rm {synthetic_alt_bam_path} {synthetic_alt_bam_path}.bai")
            else:
                merged_bam_path = synthetic_alt_allele2_bam_path

            simulated_bam_paths.append(merged_bam_path)

            bam_stats = compute_bam_stats(
                merged_bam_path,
                args.ref_fasta,
                chrom,
                start_1based - padding,
                end_1based + padding)

            metadata_tsv_row = {
                "chrom": chrom,
                "start_1based": start_1based,
                "end_1based": end_1based,
                "padding": padding,
                "allele1_repeat_unit": args.new_repeat_unit1,
                "allele1_repeats": args.num_copies1,
                "allele2_repeat_unit": args.new_repeat_unit2,
                "allele2_repeats": num_copies2,
                "het_or_hom": het_or_hom,
                "simulated_bam_path": merged_bam_path,
            }
            metadata_tsv_row.update(bam_stats)

            metadata_tsv_rows.append(metadata_tsv_row)


            # for merged_bam_path, print bam stats and regions in other parts of the genome to which the simulated reads
            # mis-align (eg. off-target regions) - this is just for logging
            if args.verbose:
                compute_off_target_regions(
                    merged_bam_path,
                    args.ref_fasta,
                    interval=f"{chrom}:{start_1based - padding}-{end_1based + padding}",
                    verbose=True)

        # compute off-target regions based on all simulated bams together
        if args.output_off_target_regions:
            all_simulated_bams_merged_bam_path = f"{output_filename_prefix}__all_{het_or_hom}.bam"
            merge_bams(
                all_simulated_bams_merged_bam_path,
                *simulated_bam_paths,
                force=args.force)

            off_target_regions, _ = compute_off_target_regions(
                all_simulated_bams_merged_bam_path,
                args.ref_fasta,
                f"{chrom}:{start_1based - padding}-{end_1based + padding}",
                min_reads_threshold=args.min_off_target_reads, verbose=True)
            off_target_regions = [r for r, read_count in off_target_regions]

            repeat_regions_output_filename = f"{output_filename_prefix}__{int(args.coverage)}x_cov__repeat_regions.bed"
            logging.info(f"Writing out {repeat_regions_output_filename}")
            with open(repeat_regions_output_filename, "w") as repeat_region_bed_file:
                repeat_region_bed_file.write("\t".join(map(str, [
                    chrom,
                    start_1based - 1,
                    end_1based,
                    args.new_repeat_unit2,
                    num_copies2,  # max number of copies that were simulated above
                    ",".join(off_target_regions) if off_target_regions is not None else "",
                    #",".join(sorted(off_target_regions, key=genomic_order)),
                ])) + "\n")
            #run(f"rm {all_simulated_bams_merged_bam_path} {all_simulated_bams_merged_bam_path}.bai")

    metadata_output_filename = f"{output_filename_prefix}__{int(args.coverage)}x_cov__simulated_bam_metadata.tsv"
    logging.info(f"Writing out {metadata_output_filename}")
    header = list(metadata_tsv_rows[0].keys())
    with open(metadata_output_filename, "wt") as f:
        f.write("\t".join(header) + "\n")
        for row in metadata_tsv_rows:
            f.write("\t".join([str(row.get(c, "")) for c in header]) + "\n")

    logging.info(f"Done generating {het_or_hom_list} x {num_copies2_list} x {args.new_repeat_unit2}")


def generate_synthetic_reference_sequence(fasta_obj, chrom, start_1based, end_1based, padding_length, repeat_unit, num_copies):
    """Generates a nucleotide sequence that consists of {num_copies} of the {repeat_unit} surrounded by {padding_length}
    bases from the reference on either side of the interval given by {chrom}:{start_1based}-{end_1based}.
    """
    iupac_nucleotide_codes = {
        "N": "ACGT",
        "R": "AG",
        "Y": "CT",
        "S": "GC",
        "W": "AT",
        "K": "GT",
        "M": "AC",
    }

    if len(set(iupac_nucleotide_codes) & set(repeat_unit)) > 0:
        i = 0
        repeat_sequence = ""
        while len(repeat_sequence) < len(repeat_unit) * num_copies:
            next_repeat_unit = repeat_unit
            for iupac_code in iupac_nucleotide_codes:
                if iupac_code in next_repeat_unit:
                    possible_nucleotides = iupac_nucleotide_codes[iupac_code]
                    next_repeat_unit = next_repeat_unit.replace(iupac_code, possible_nucleotides[i % len(possible_nucleotides)])
            repeat_sequence += next_repeat_unit
            i += 1
    else:
        repeat_sequence = repeat_unit * num_copies

    if padding_length == 0:
        return repeat_sequence

    if chrom not in fasta_obj:
        raise ValueError(f"Invalid chromosome name: {chrom}")

    left_padding_start_1based = max(1, start_1based - padding_length)
    left_padding = get_reference_sequence(fasta_obj, chrom, left_padding_start_1based, start_1based - 1)

    chromosome_length = len(fasta_obj[chrom])
    right_padding_end_1based = min(end_1based + padding_length, chromosome_length)
    right_padding = get_reference_sequence(fasta_obj, chrom, end_1based + 1, right_padding_end_1based)

    return left_padding + repeat_sequence + right_padding


def compute_off_target_regions(merged_bam_path, ref_fasta, interval, min_reads_threshold=2, verbose=False):
    """Returns a list of genomic intervals that are outside the locus given by {chrom}:{start_1based}-{end_1based}
    but contain at least 'min_reads_threshold' aligned reads. This means that reads from the given locus can mis-align to
    these other regions.
    """

    #on_target_count = int(run(f"samtools view -c {merged_bam_path} {chrom}:{start_1based}-{end_1based}"))

    # compute off-target regions command
    off_target_reads_command = ""
    if interval is not None:
        chrom, start_1based, end_1based = parse_interval(interval)
        logging.info(f"Computing off-target regions for {chrom}:{start_1based}-{end_1based}")
        off_target_reads_command = f"echo '{chrom}\t{start_1based - 1}\t{end_1based}' "
        off_target_reads_command += f"| bedtools intersect -nonamecheck -v -wa -a {merged_bam_path} -b - "
        off_target_reads_command += f"| samtools view -F 4 -b - "  # exclude unmapped reads because this causes errors
    else:
        off_target_reads_command += f"samtools view -F 4 -b {merged_bam_path} "  # exclude unmapped reads because this causes errors

    logging.info(f"Minimum reads threshold for outputing an off-target region: {min_reads_threshold}")

    off_target_reads_command += f"| bedtools merge -d 100 "  # merge reads' intervals
    off_target_reads_command += f"| bedtools slop -b 150 -g {ref_fasta}.fai "
    off_target_reads_command += f"| bedtools sort "
    off_target_reads_command += f"| bedtools intersect -nonamecheck -wa -a - -b {merged_bam_path} -c "   # add counts of how many reads mapped to each off-target interval
    off_target_reads_command += f"| sort -n -k4 -r "   # sort by these counts, so that the 1st interval is the one with the most reads, etc.
    off_target_regions = [off_target_region for off_target_region in run(off_target_reads_command).strip().split("\n") if off_target_region]

    # print some stats
    total_read_count = int(run(f"samtools view -c {merged_bam_path}"))

    try:
        off_target_read_count = sum([int(off_target_region.split("\t")[-1]) for off_target_region in off_target_regions])
    except ValueError:
        off_target_read_count = 0

    if total_read_count > 0:
        fraction_of_reads_that_map_to_off_target_regions = off_target_read_count/total_read_count
    else:
        fraction_of_reads_that_map_to_off_target_regions = 0

    if verbose:
        logging.info(f"Found {len(off_target_regions)} off-target regions "
                 f"with {off_target_read_count} out of {total_read_count} reads "
                 f"({0 if not total_read_count else (100*off_target_read_count/total_read_count):0.1f}%) mismapping "
                 f"to them.")


        logging.info(f"Main interval: {interval}   {total_read_count - off_target_read_count}\n"
            + "\n".join(off_target_regions))

    off_target_region_list = []
    for off_target_region in off_target_regions:
        fields = off_target_region.strip("\n").split("\t")
        if len(fields) <= 3:
            logging.warning(f"Counts column missing for off_target_region: {off_target_region}. Skipping..." )
            continue

        # Skip loci with less than {coverage_threshold} reads in order to reduce noise
        num_reads_mapped_to_region = int(fields[3])
        if num_reads_mapped_to_region < min_reads_threshold:
            continue
        off_target_region_list.append((f"{fields[0]}:{max(1, int(fields[1]))}-{int(fields[2])}", num_reads_mapped_to_region))

    return off_target_region_list, fraction_of_reads_that_map_to_off_target_regions


if __name__ == "__main__":
    main()


# wgsim command line options
"""
wgsim [options] <in.ref.fa> <out.read1.fq> <out.read2.fq>

Options: -e FLOAT      base error rate [0.020]
         -d INT        outer distance between the two ends [500]
         -s INT        standard deviation [50]
         -N INT        number of read pairs [1000000]
         -1 INT        length of the first read [70]
         -2 INT        length of the second read [70]
         -r FLOAT      rate of mutations [0.0010]
         -R FLOAT      fraction of indels [0.15]
         -X FLOAT      probability an indel is extended [0.30]
         -S INT        seed for random generator [0, use the current time]
         -A FLOAT      discard if the fraction of ambiguous bases higher than FLOAT [0.05]
         -h            haplotype mode
         
"""

