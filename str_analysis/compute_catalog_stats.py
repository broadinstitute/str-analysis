import argparse
import collections
import ijson
import numpy as np
import os
import pandas as pd
import pysam
import re
import statistics
import tqdm
from ncls import NCLS

from str_analysis.utils.find_motif_utils import compute_repeat_purity
from str_analysis.utils.fasta_utils import create_normalize_chrom_function
from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.eh_catalog_utils import parse_motifs_from_locus_structure, get_variant_catalog_iterator
from str_analysis.utils.file_utils import file_exists

ACGT_REGEX = re.compile("^[ACGT]+$", re.IGNORECASE)
ACGTN_REGEX = re.compile("^[ACGTN]+$", re.IGNORECASE)

GENE_REGIONS = ["5utr", "3utr", "cds", "exon", "intergenic", "intron", "promoter"]


def build_ncls_trees(intervals_by_chrom):
    """Build NCLS interval trees from collected intervals.

    Args:
        intervals_by_chrom (dict): Dictionary mapping chromosome names to interval data
            with keys 'starts', 'ends', 'indices'

    Returns:
        dict: Dictionary mapping chromosome names to NCLS objects
    """
    ncls_trees = {}
    for chrom, intervals_data in intervals_by_chrom.items():
        if len(intervals_data['starts']) > 0:
            starts = np.array(intervals_data['starts'], dtype=np.int64)
            ends = np.array(intervals_data['ends'], dtype=np.int64)
            indices = np.array(intervals_data['indices'], dtype=np.int64)
            ncls_trees[chrom] = NCLS(starts, ends, indices)
    return ncls_trees


def find_overlapping_intervals(all_intervals, ncls_trees):
    """Find intervals that overlap each other by at least 2x the larger motif size.

    Args:
        all_intervals (list): List of interval dictionaries with keys 'chrom', 'start', 'end', 'motif_size'
        ncls_trees (dict): Dictionary mapping chromosome names to NCLS objects

    Returns:
        set: Set of tuples (chrom, start, end) representing overlapping intervals
    """
    overlapping_intervals = set()

    for interval_data in all_intervals:
        chrom = interval_data['chrom']
        start = interval_data['start']
        end = interval_data['end']
        motif_size = interval_data['motif_size']

        if chrom not in ncls_trees:
            continue

        overlaps = list(ncls_trees[chrom].find_overlap(start, max(end, start + 1)))

        for _, _, idx in overlaps:
            if idx < len(all_intervals):
                other = all_intervals[idx]
                # Skip self-overlap
                if other['start'] == start and other['end'] == end:
                    continue

                # Calculate overlap size
                overlap_start = max(start, other['start'])
                overlap_end = min(end, other['end'])
                overlap_size = max(0, overlap_end - overlap_start)

                larger_motif_size = max(motif_size, other['motif_size'])
                if overlap_size >= 2 * larger_motif_size:
                    overlapping_intervals.add((chrom, start, end))
                    overlapping_intervals.add((chrom, other['start'], other['end']))
                    break

    return overlapping_intervals


def count_merged_base_pairs(intervals_by_chrom):
    """Count total base pairs covered by merging overlapping intervals.

    Args:
        intervals_by_chrom (dict): Dictionary mapping chromosome names to interval data
            with keys 'starts', 'ends'

    Returns:
        int: Total number of base pairs covered by merged intervals
    """
    total_base_pairs = 0

    for chrom, intervals_data in intervals_by_chrom.items():
        if len(intervals_data['starts']) == 0:
            continue

        # Sort intervals by start position
        sorted_indices = np.argsort(intervals_data['starts'])
        starts = np.array(intervals_data['starts'])[sorted_indices]
        ends = np.array(intervals_data['ends'])[sorted_indices]

        # Merge overlapping intervals
        if len(starts) > 0:
            merged_start = starts[0]
            merged_end = ends[0]

            for i in range(1, len(starts)):
                if starts[i] <= merged_end:
                    # Overlapping, extend the current merged interval
                    merged_end = max(merged_end, ends[i])
                else:
                    # Non-overlapping, add the previous merged interval
                    total_base_pairs += merged_end - merged_start
                    merged_start = starts[i]
                    merged_end = ends[i]

            # Add the last merged interval
            total_base_pairs += merged_end - merged_start

    return total_base_pairs


def parse_args():
    parser = argparse.ArgumentParser(description="Compute and print stats for annotated repeat catalogs")
    parser.add_argument("-R", "--reference-fasta", help="Optional reference fasta file for computing reference repeat purity")
    parser.add_argument("--verbose", action="store_true", help="Print more information about what the script is doing")
    parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    parser.add_argument("-n", "--n-loci", type=int, help="Only process the first N loci in the catalog")
    parser.add_argument("variant_catalog_json_or_bed", nargs="+",
                        help="Repeat catalog in JSON or BED format. Repeat catalog(s) processed by "
                             "the annotate_and_filter_str_catalog.py script will have extra stats computed")

    args = parser.parse_args()

    for file_path in args.variant_catalog_json_or_bed:
        if not file_exists(os.path.expanduser(file_path)):
            parser.error(f"File not found: {file_path}")

    return args, parser



def compute_catalog_stats(catalog_name, records, reference_fasta_path=None, verbose=False, show_progress_bar=False, n_loci=None):
    """This script takes a TR catalog in BED format or in ExpansionHunter JSON format and outputs statistics about the
    distribution of the TRs in the catalog.

    Args:
        catalog_name (str): name of the catalog
        records (list of dicts): list of TRs in the catalog in ExpansionHunter JSON format
        reference_fasta_path (str): optional path to the reference fasta file

    Return:
        dict: dictionary of summary stats for this catalog
    """

    min_motif_size = min_locus_size = min_num_repeats_in_locus = 10**9
    max_motif_size = max_locus_size = max_num_repeats_in_locus = 0
    max_locus_size_reference_region = None
    max_locus_size_motif = None
    min_fraction_pure_bases = 1
    min_fraction_pure_bases_reference_region = None
    min_fraction_pure_bases_motif = None
    min_overall_mappability = 1
    min_overall_mappability_reference_region = None
    min_overall_mappability_motif = None

    if show_progress_bar:
        records = tqdm.tqdm(records, unit=" records", unit_scale=True, total=n_loci)

    # Collect intervals for each chromosome to build NCLS trees later
    intervals_by_chrom = collections.defaultdict(lambda: {'starts': [], 'ends': [], 'indices': [], 'motif_sizes': []})
    all_intervals = []  # Store all intervals with metadata for overlap checking
    overlapping_intervals = set()
    counters = collections.defaultdict(int)

    locus_sizes_by_motif_size = collections.defaultdict(list)  # used to compute the median locus sizes for each motif size
    reference_repeat_purity_by_motif_size = collections.defaultdict(list)  # used to compute the mean base purity for each motif size
    mappability_by_motif_size = collections.defaultdict(list)  # used to compute the mean mappability for each motif size

    ref_fasta = normalize_chrom = None
    if reference_fasta_path:
        ref_fasta = pysam.FastaFile(reference_fasta_path)
        does_chrom_start_with_chr = ref_fasta.references[0].startswith("chr")
        normalize_chrom = create_normalize_chrom_function(does_chrom_start_with_chr)


    has_reference_repeat_purity = False
    has_gene_annotations = False
    for record in records:
        counters["total"] += 1
        if n_loci is not None and counters["total"] > n_loci:
            break

        motifs = parse_motifs_from_locus_structure(record["LocusStructure"])
        if isinstance(record["ReferenceRegion"], list):
            reference_regions = record["ReferenceRegion"]
            variant_types = record["VariantType"]
            fraction_pure_bases_list = record.get("ReferenceRepeatPurity", [None]*len(reference_regions))

            counters["loci_with_adjacent_repeats"] += 1
        else:
            reference_regions = [record["ReferenceRegion"]]
            variant_types = [record["VariantType"]]
            fraction_pure_bases_list = [record.get("ReferenceRepeatPurity")]

        if "ReferenceRepeatPurity" not in record and ref_fasta is not None:
            fraction_pure_bases_list = []
            for reference_region, motif in zip(reference_regions, motifs):
                chrom, start_0based, end = parse_interval(reference_region)
                reference_fasta_sequence= ref_fasta.fetch(normalize_chrom(chrom), start_0based, end)
                fraction_pure_bases, _ = compute_repeat_purity(reference_fasta_sequence, motif, include_partial_repeats=True)
                fraction_pure_bases_list.append(fraction_pure_bases)
            has_reference_repeat_purity = True
        elif "ReferenceRepeatPurity" in record:
            has_reference_repeat_purity = True

        for motif, reference_region, variant_type, fraction_pure_bases in zip(
                motifs, reference_regions, variant_types, fraction_pure_bases_list):
            counters["total_repeat_intervals"] += 1

            chrom, start_0based, end = parse_interval(reference_region)
            if chrom.endswith("X"):
                counters["chrX"] += 1
            elif chrom.endswith("Y"):
                counters["chrY"] += 1
            elif chrom.upper().endswith("M") or chrom.upper().endswith("MT"):
                counters["chrM"] += 1
            elif len(chrom) > 5:
                counters["alt_contigs"] += 1

            motif_size = len(motif)
            locus_size = end - start_0based

            if motif_size <= 50:
                locus_sizes_by_motif_size[motif_size].append(locus_size)
                if fraction_pure_bases is not None and pd.notna(fraction_pure_bases):
                    reference_repeat_purity_by_motif_size[motif_size].append(fraction_pure_bases)
                    reference_repeat_purity_by_motif_size['all'].append(fraction_pure_bases)
                if "FlanksAndLocusMappability" in record:
                    mappability_by_motif_size[motif_size].append(record["FlanksAndLocusMappability"])

            min_motif_size = min(min_motif_size, motif_size)
            max_motif_size = max(max_motif_size, motif_size)
            min_locus_size = min(min_locus_size, locus_size)
            max_locus_size = max(max_locus_size, locus_size)
            num_repeats_per_locus = int(locus_size / motif_size)
            min_num_repeats_in_locus = min(min_num_repeats_in_locus, num_repeats_per_locus)
            max_num_repeats_in_locus = max(max_num_repeats_in_locus, num_repeats_per_locus)

            counters["sum_of_all_locus_interval_sizes"] += locus_size
            if locus_size == max_locus_size:
                max_locus_size_reference_region = reference_region
                max_locus_size_motif = motif

            if fraction_pure_bases is not None and pd.notna(fraction_pure_bases):
                min_fraction_pure_bases = min(min_fraction_pure_bases, fraction_pure_bases)
                if fraction_pure_bases == min_fraction_pure_bases:
                    min_fraction_pure_bases_reference_region = reference_region
                    min_fraction_pure_bases_motif = motif

            if motif_size == 1:
                counters["homopolymers"] += 1
            if locus_size % motif_size == 0:
                counters["trimmed"] += 1

            if not ACGT_REGEX.match(motif):
                counters["non_ACGT_motifs"] += 1
            if not ACGTN_REGEX.match(motif):
                counters["non_ACGTN_motifs"] += 1

            motif_size_bin = f"{motif_size}bp" if motif_size <= 6 else "7-24bp" if motif_size <= 24 else "25+bp"
            counters[f"motif_size:{motif_size_bin}"] += 1

            num_repeats_per_locus_bin = f"{num_repeats_per_locus}x" if num_repeats_per_locus <= 9 else "10-15x" if num_repeats_per_locus <= 15 else "16-25x" if num_repeats_per_locus <= 25 else "26-35x" if num_repeats_per_locus <= 35 else "36-50x" if num_repeats_per_locus <= 50 else "51+x"
            counters[f"num_repeats_per_locus:{num_repeats_per_locus_bin}"] += 1

            num_repeats_per_locus_detailed_bin = f"{num_repeats_per_locus}x" if num_repeats_per_locus <= 24 else "25-50x" if num_repeats_per_locus <= 50 else "51+x"
            #if num_repeats_per_locus <= 3:
            #    print(num_repeats_per_locus_bin, record["LocusId"], reference_region, motif)
            counters[f"num_repeats_per_locus_detailed:{num_repeats_per_locus_detailed_bin}"] += 1


            if record.get("GencodeGeneRegion"):
                has_gene_annotations = True
                gene_region = record["GencodeGeneRegion"].lower().replace("' ", "")
                if gene_region not in GENE_REGIONS:
                    raise ValueError(f"Unexpected gene region: {gene_region}")

                counters[f"gene_region:{gene_region}"] += 1
                if motif_size == 3:
                    counters[f"3bp_motif_gene_region:{gene_region}"] += 1

            if fraction_pure_bases is not None and pd.notna(fraction_pure_bases):
                fraction_pure_bases_bin = round(int(fraction_pure_bases*10)/10, 1)
                counters[f"fraction_pure_bases:{fraction_pure_bases_bin}"] += 1

            # Collect interval for later NCLS tree building and overlap checking
            interval_idx = len(all_intervals)
            intervals_by_chrom[chrom]['starts'].append(start_0based)
            intervals_by_chrom[chrom]['ends'].append(max(end, start_0based + 1))
            intervals_by_chrom[chrom]['indices'].append(interval_idx)
            intervals_by_chrom[chrom]['motif_sizes'].append(motif_size)
            all_intervals.append({
                'chrom': chrom,
                'start': start_0based,
                'end': end,
                'motif_size': motif_size,
            })

        if "FlanksAndLocusMappability" in record:
            min_overall_mappability = min(min_overall_mappability, record["FlanksAndLocusMappability"])
            if record["FlanksAndLocusMappability"] == min_overall_mappability:
                min_overall_mappability_reference_region = reference_region
                min_overall_mappability_motif = motif
            mappability_bin = round(int(record["FlanksAndLocusMappability"]*10)/10, 1)
            counters[f"mappability:{mappability_bin}"] += 1

    # Build NCLS trees for overlap detection and base pair counting
    if verbose:
        print("Building NCLS trees for overlap detection...")
    ncls_trees = build_ncls_trees(intervals_by_chrom)

    # Check for overlaps between intervals
    if verbose:
        print("Checking for overlapping intervals...")
    overlapping_intervals = find_overlapping_intervals(all_intervals, ncls_trees)

    # Count total base pairs covered by merging overlapping intervals
    if verbose:
        print("Counting total base pairs covered...")
    counters['total_base_pairs_covered_by_loci'] = count_merged_base_pairs(intervals_by_chrom)

    print("")
    print(f"Stats for {catalog_name}:")
    print(f"   {counters['total']:10,d} total loci")
    if counters["total_repeat_intervals"] == 0:
        return {
            "catalog": catalog_name,
            "total": f"{counters['total_repeat_intervals']:,d}",
        }

    total_genome_size = 3_088_286_401  # in hg38, including chrX, chrY, chrM
    print(f"   {counters['total_base_pairs_covered_by_loci']:10,d} base pairs covered by all loci ({counters['total_base_pairs_covered_by_loci']/total_genome_size:0.3%} of the genome)")
    print(f"   {counters['sum_of_all_locus_interval_sizes']:10,d} sum of all interval widths ({counters['sum_of_all_locus_interval_sizes']/counters['total_base_pairs_covered_by_loci']:0.3f}x overlap on average)")
    print(f"   {counters['loci_with_adjacent_repeats']:10,d} out of {counters['total']:10,d} ({counters['loci_with_adjacent_repeats']/counters['total']:6.1%}) loci define adjacent repeats")
    print(f"   {counters['total_repeat_intervals']:10,d} total repeat intervals")
    print(f"   {counters['trimmed']:10,d} out of {counters['total_repeat_intervals']:10,d} ({counters['trimmed']/counters['total_repeat_intervals']:6.1%}) repeat interval size is an integer multiple of the motif size (aka. trimmed)")
    print(f"   {counters['homopolymers']:10,d} out of {counters['total_repeat_intervals']:10,d} ({counters['homopolymers']/counters['total_repeat_intervals']:6.1%}) repeat intervals are homopolymers")
    print(f"   {len(overlapping_intervals):10,d} out of {counters['total_repeat_intervals']:10,d} ({len(overlapping_intervals)/counters['total_repeat_intervals']:6.1%}) repeat intervals overlap each other by at least two motif lengths")
    if counters["non_ACGT_motifs"]:
        print(f"   {counters['non_ACGT_motifs']:10,d} out of {counters['total_repeat_intervals']:10,d} ({counters['non_ACGT_motifs']/counters['total_repeat_intervals']:6.1%}) repeat intervals have non-ACGT motifs")
    if counters["non_ACGTN_motifs"]:
        print(f"   {counters['non_ACGTN_motifs']:10,d} out of {counters['total_repeat_intervals']:10,d} ({counters['non_ACGTN_motifs']/counters['total_repeat_intervals']:6.1%}) repeat intervals have non-ACGTN motifs")

    if len(overlapping_intervals) > 0:
        examples = list(overlapping_intervals)
        examples = ", ".join(f"{chrom}:{start}-{end}" for chrom, start, end in examples[:10])
        print(f"Examples of overlapping repeats: {examples}")

    print("")
    print("Ranges:")
    print(f"   Motif size range: {min_motif_size}-{max_motif_size}bp")
    print(f"   Locus size range: {min_locus_size}-{max_locus_size}bp")
    print(f"   Num repeats range: {min_num_repeats_in_locus}-{max_num_repeats_in_locus}x repeats")
    print("")
    print(f"   Max locus size = {max_locus_size:7,d}bp           @ {max_locus_size_reference_region} ({max_locus_size_motif})")
    if min_fraction_pure_bases_motif is not None:
        print(f"   Min reference repeat purity   = {min_fraction_pure_bases:5.2f}    @ {min_fraction_pure_bases_reference_region} ({min_fraction_pure_bases_motif})")
    if min_overall_mappability_motif is not None:
        print(f"   Min overall mappability   = {min_overall_mappability:5.2f}    @ {min_overall_mappability_reference_region} ({min_overall_mappability_motif})")
    if any([_ is not None for _ in reference_repeat_purity_by_motif_size['all']]) > 0:
        print(f"   Base-level   purity   median: {statistics.median(reference_repeat_purity_by_motif_size['all']):0.3f},  mean: {statistics.mean(reference_repeat_purity_by_motif_size['all']):0.3f}")
    print("")
    print(f"          chrX: {counters['chrX']:10,d} out of {counters['total_repeat_intervals']:10,d} ({counters['chrX']/counters['total_repeat_intervals']:6.1%}) repeat intervals")
    print(f"          chrY: {counters['chrY']:10,d} out of {counters['total_repeat_intervals']:10,d} ({counters['chrY']/counters['total_repeat_intervals']:6.1%}) repeat intervals")
    print(f"          chrM: {counters['chrM']:10,d} out of {counters['total_repeat_intervals']:10,d} ({counters['chrM']/counters['total_repeat_intervals']:6.1%}) repeat intervals")
    print(f"   alt contigs: {counters['alt_contigs']:10,d} out of {counters['total_repeat_intervals']:10,d} ({counters['alt_contigs']/counters['total_repeat_intervals']:6.1%}) repeat intervals")
    print("")
    print("Motif size distribution:")
    for motif_size_bin in "1bp", "2bp", "3bp", "4bp", "5bp", "6bp", "7-24bp", "25+bp":
        print(f"   {motif_size_bin:>10s}: {counters[f'motif_size:{motif_size_bin}']:10,d} out of {counters['total_repeat_intervals']:10,d} ({counters[f'motif_size:{motif_size_bin}']/counters['total_repeat_intervals']:6.1%}) repeat intervals")
    print("")
    print("Num repeats in reference:")
    for num_repeats_per_locus_bin in "1x", "2x", "3x", "4x", "5x", "6x", "7x", "8x", "9x", "10-15x", "16-25x", "26-35x", "36-50x", "51+x":
        print(f"   {num_repeats_per_locus_bin:>10s}: {counters[f'num_repeats_per_locus:{num_repeats_per_locus_bin}']:10,d} out of {counters['total_repeat_intervals']:10,d} ({counters[f'num_repeats_per_locus:{num_repeats_per_locus_bin}']/counters['total_repeat_intervals']:6.1%}) repeat intervals")

    if min_fraction_pure_bases_motif is not None:
        print("")
        print("Reference repeat purity distribution:")
        for fraction_pure_bases_bin in 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1:
            print(f"   {fraction_pure_bases_bin:10.1f}: {counters[f'fraction_pure_bases:{fraction_pure_bases_bin:.1f}']:10,d} out of {counters['total_repeat_intervals']:10,d} ({counters[f'fraction_pure_bases:{fraction_pure_bases_bin:.1f}']/counters['total_repeat_intervals']:6.1%}) repeat intervals")

    if min_overall_mappability_motif:
        print("")
        print("Mappability distribution:")
        for mappability_bin in 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1:
            print(f"   {mappability_bin:10.1f}: {counters[f'mappability:{mappability_bin:.1f}']:10,d} out of {counters['total']:10,d} ({counters[f'mappability:{mappability_bin:.1f}']/counters['total']:6.1%}) loci")

    result = {
        "catalog": catalog_name,
        "total": counters['total_repeat_intervals'],
        "count_chrX": counters['chrX'],
        "count_chrY": counters['chrY'],
        "count_chrM": counters["chrM"],
        "percent_of_genome_spanned_by_loci": f"{counters['total_base_pairs_covered_by_loci']/total_genome_size:0.3%}",
        "motif_size_range": f"{min_motif_size}-{max_motif_size}bp",
        "locus_size_range": f"{min_locus_size}-{max_locus_size}bp",
        "num_repeats_range": f"{min_num_repeats_in_locus}-{max_num_repeats_in_locus}x repeats",
        "percent_homopolymers": "%0.1f%%" % (100 * counters["motif_size:1bp"] / counters["total_repeat_intervals"]),
        "percent_2bp_motifs": "%0.1f%%" % (100 * counters["motif_size:1bp"] / counters["total_repeat_intervals"]),
        "percent_2bp_motifs": "%0.1f%%" % (100 * counters["motif_size:2bp"] / counters["total_repeat_intervals"]),
        "percent_3bp_motifs": "%0.1f%%" % (100 * counters["motif_size:3bp"] / counters["total_repeat_intervals"]),
        "percent_4bp_motifs": "%0.1f%%" % (100 * counters["motif_size:4bp"] / counters["total_repeat_intervals"]),
        "percent_5bp_motifs": "%0.1f%%" % (100 * counters["motif_size:5bp"] / counters["total_repeat_intervals"]),
        "percent_6bp_motifs": "%0.1f%%" % (100 * counters["motif_size:6bp"] / counters["total_repeat_intervals"]),
        "percent_7+bp_motifs": "%0.1f%%" % (100 * (counters["motif_size:7-24bp"] + counters["motif_size:25+bp"])/ counters["total_repeat_intervals"]),
        "percent_pure_bases": "%0.1f%%" % (100 * counters[f"fraction_pure_bases:1.0"] / counters["total_repeat_intervals"]),
        "percent_trimmed": "%0.1f%%" % (100 * counters["trimmed"] / counters["total_repeat_intervals"]),
        "percent_overlapping": "%0.1f%%" % (100 * len(overlapping_intervals) / counters['total_repeat_intervals']),
        "count_homopolymers": counters["motif_size:1bp"],
        "count_1bp_motifs": counters["motif_size:1bp"],
        "count_2bp_motifs": counters["motif_size:2bp"],
        "count_3bp_motifs":  counters["motif_size:3bp"],
        "count_4bp_motifs":  counters["motif_size:4bp"],
        "count_5bp_motifs":  counters["motif_size:5bp"],
        "count_6bp_motifs":  counters["motif_size:6bp"],
        "count_7+bp_motifs": counters["motif_size:7-24bp"] + counters["motif_size:25+bp"],
        "count_7-24bp_motifs": counters["motif_size:7-24bp"],
        "count_25+bp_motifs": counters["motif_size:25+bp"],
        "count_pure_bases": counters[f"fraction_pure_bases:1.0"],
        "count_trimmed": counters["trimmed"],
        "count_overlapping": len(overlapping_intervals),
        "min_motif_size": min_motif_size,
        "max_motif_size": max_motif_size,
        "min_locus_size": min_locus_size,
        "max_locus_size": max_locus_size,
    }

    if has_gene_annotations:
        for gene_region in GENE_REGIONS:
            result[f"count_gene_region_{gene_region}"] = counters[f"gene_region:{gene_region}"]
        for gene_region in GENE_REGIONS:
            result[f"count_3bp_motif_gene_region_{gene_region}"] = counters[f"3bp_motif_gene_region:{gene_region}"]

    print("")
    print("Locus sizes at each motif size:")
    for motif_size in list(sorted(locus_sizes_by_motif_size.keys()))[:10]:
        if motif_size not in locus_sizes_by_motif_size:
            continue
        min_size = int(min(locus_sizes_by_motif_size[motif_size]))
        median_size = int(statistics.median(locus_sizes_by_motif_size[motif_size]))
        max_size = int(max(locus_sizes_by_motif_size[motif_size]))
        if any(_ is not None for _ in reference_repeat_purity_by_motif_size[motif_size]):
            mean_reference_repeat_purity = statistics.mean(reference_repeat_purity_by_motif_size[motif_size])
        else:
            mean_reference_repeat_purity = float("nan")

        mean_mappability = statistics.mean(mappability_by_motif_size[motif_size]) if mappability_by_motif_size[motif_size] else None
        result[f"{motif_size}bp motifs: min locus size"] = min_size
        result[f"{motif_size}bp motifs: median locus size"] = median_size
        result[f"{motif_size}bp motifs: max locus size"] = max_size

        print(f"   {motif_size:3,d}bp motifs: locus size range:   {min_size:4,d} bp to {max_size:7,d} bp  (median: {int(median_size):4,d} bp) based on {len(locus_sizes_by_motif_size[motif_size]):10,d} loci. Mean base purity: {mean_reference_repeat_purity:0.2f}. ", f"Mean mappability: {mean_mappability:0.2f}" if mean_mappability is not None else "")
    if not has_reference_repeat_purity and not reference_fasta_path:
        print("WARNING: No reference repeat purity data was found in the catalog. Specify --reference-fasta to compute reference repeat purity.")

    for num_repeats_per_locus_detailed_bin in (
        "0x", "1x", "2x", "3x", "4x", "5x", "6x", "7x", "8x", "9x", "10x",
        "11x", "12x", "13x", "14x", "15x", "16x", "17x", "18x", "19x", "20x",
        "21x", "22x", "23x", "24x", "25-50x", "51+x"
    ):
        result[f"num_repeats_per_locus:{num_repeats_per_locus_detailed_bin}"] = counters[f"num_repeats_per_locus_detailed:{num_repeats_per_locus_detailed_bin}"]


    return result



def main():
    args, parser = parse_args()

    stat_table_rows = []
    for path in args.variant_catalog_json_or_bed:
        print("-"*50)
        print(f"Parsing {path}")
        catalog_name = os.path.basename(path)

        file_iterator = get_variant_catalog_iterator(path)
        stats = compute_catalog_stats(
            catalog_name, 
            file_iterator, 
            reference_fasta_path=args.reference_fasta, 
            verbose=args.verbose, 
            show_progress_bar=args.show_progress_bar,
            n_loci=args.n_loci)
        stat_table_rows.append(stats)

    if len(args.variant_catalog_json_or_bed) > 1:
        output_path = f"catalog_stats.{len(args.variant_catalog_json_or_bed)}_catalogs.tsv"
    else:
        output_path = re.sub("(.json|.bed)(.gz)?$", "", os.path.basename(args.variant_catalog_json_or_bed[0])) + ".catalog_stats.tsv"

    df = pd.DataFrame(stat_table_rows)
    columns_to_drop = ["motif_size_distribution", "fraction_pure_bases_distribution", "mappability_distribution"]
    for c in columns_to_drop:
        if c in df.columns:
            df.drop(columns=c, inplace=True)
    df.to_csv(output_path, sep="\t", index=False)
    print(f"Wrote {len(stat_table_rows):,d} rows to {output_path}")


if __name__ == "__main__":
    main()