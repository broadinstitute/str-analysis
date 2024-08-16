import argparse
import collections
import ijson
import os
import pandas as pd
import re
import tqdm
from intervaltree import IntervalTree, Interval
from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.eh_catalog_utils import parse_motifs_from_locus_structure, get_variant_catalog_iterator
from str_analysis.utils.file_utils import open_file, file_exists

ACGT_REGEX = re.compile("^[ACGT]+$", re.IGNORECASE)
ACGTN_REGEX = re.compile("^[ACGTN]+$", re.IGNORECASE)

def parse_args():
    parser = argparse.ArgumentParser(description="Compute and print stats for annotated repeat catalogs")
    parser.add_argument("--verbose", action="store_true", help="Print more information about what the script is doing")
    parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    parser.add_argument("variant_catalog_json_or_bed", nargs="+",
                        help="Repeat catalog in JSON or BED format. Repeat catalog(s) processed by "
                             "the annotate_and_filter_str_catalog.py script will have extra stats computed")

    args = parser.parse_args()

    for file_path in args.variant_catalog_json_or_bed:
        if not file_exists(os.path.expanduser(file_path)):
            parser.error(f"File not found: {file_path}")

    return args, parser



def compute_catalog_stats(catalog_name, records, verbose=False, show_progress_bar=False):
    """This script takes a TR catalog in BED format or in ExpansionHunter JSON format and outputs statistics about the
    distribution of the TRs in the catalog.

    Args:
        catalog_name (str): name of the catalog
        records (list of dicts): list of TRs in the catalog in ExpansionHunter JSON format

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
    min_fraction_pure_repeats = 1
    min_fraction_pure_repeats_reference_region = None
    min_fraction_pure_repeats_motif = None
    min_overall_mappability = 1
    min_overall_mappability_reference_region = None
    min_overall_mappability_motif = None

    if show_progress_bar:
        records = tqdm.tqdm(records, unit=" records", unit_scale=True)

    interval_trees = collections.defaultdict(IntervalTree)  # used to check for overlap between records in the catalog
    overlapping_intervals = set()
    counters = collections.defaultdict(int)
    for record in records:
        counters["total"] += 1
        motifs = parse_motifs_from_locus_structure(record["LocusStructure"])
        if isinstance(record["ReferenceRegion"], list):
            reference_regions = record["ReferenceRegion"]
            variant_types = record["VariantType"]
            fraction_pure_bases = record.get("FractionPureBases", [None]*len(reference_regions))
            fraction_pure_repeats = record.get("FractionPureRepeats", [None]*len(reference_regions))

            counters["loci_with_adjacent_repeats"] += 1
        else:
            reference_regions = [record["ReferenceRegion"]]
            variant_types = [record["VariantType"]]
            fraction_pure_bases = [record.get("FractionPureBases")]
            fraction_pure_repeats = [record.get("FractionPureRepeats")]

        for motif, reference_region, variant_type, fraction_pure_bases, fraction_pure_repeats in zip(
                motifs, reference_regions, variant_types, fraction_pure_bases, fraction_pure_repeats):
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
            min_motif_size = min(min_motif_size, motif_size)
            max_motif_size = max(max_motif_size, motif_size)
            min_locus_size = min(min_locus_size, locus_size)
            max_locus_size = max(max_locus_size, locus_size)
            motif_sizes_per_locus_size = int(locus_size / motif_size)
            min_num_repeats_in_locus = min(min_num_repeats_in_locus, motif_sizes_per_locus_size)
            max_num_repeats_in_locus = max(max_num_repeats_in_locus, motif_sizes_per_locus_size)

            counters["total_base_pairs_spanned_by_all_loci"] += locus_size
            if locus_size == max_locus_size:
                max_locus_size_reference_region = reference_region
                max_locus_size_motif = motif

            if fraction_pure_bases is not None:
                min_fraction_pure_bases = min(min_fraction_pure_bases, fraction_pure_bases)
                if fraction_pure_bases == min_fraction_pure_bases:
                    min_fraction_pure_bases_reference_region = reference_region
                    min_fraction_pure_bases_motif = motif
            if fraction_pure_repeats is not None:
                min_fraction_pure_repeats = min(min_fraction_pure_repeats, fraction_pure_repeats)
                if fraction_pure_repeats == min_fraction_pure_repeats:
                    min_fraction_pure_repeats_reference_region = reference_region
                    min_fraction_pure_repeats_motif = motif

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

            motif_sizes_per_locus_size_bin = f"{motif_sizes_per_locus_size}x" if motif_sizes_per_locus_size <= 9 else "10-15x" if motif_sizes_per_locus_size <= 15 else "16-25x" if motif_sizes_per_locus_size <= 25 else "26-35x" if motif_sizes_per_locus_size <= 35 else "36-50x" if motif_sizes_per_locus_size <= 50 else "51+x"
            counters[f"motif_sizes_per_locus_size:{motif_sizes_per_locus_size_bin}"] += 1

            motif_sizes_per_locus_size_detailed_bin = f"{motif_sizes_per_locus_size}x" if motif_sizes_per_locus_size <= 24 else "25-50x" if motif_sizes_per_locus_size <= 50 else "51+x"
            #if motif_sizes_per_locus_size <= 3:
            #    print(motif_sizes_per_locus_size_bin, record["LocusId"], reference_region, motif)
            counters[f"motif_sizes_per_locus_size_detailed:{motif_sizes_per_locus_size_detailed_bin}"] += 1

            if fraction_pure_bases is not None:
                fraction_pure_bases_bin = round(int(fraction_pure_repeats*10)/10, 1)
                counters[f"fraction_pure_bases:{fraction_pure_bases_bin}"] += 1

            # check for overlap
            for overlapping_interval in interval_trees[chrom].overlap(start_0based, end):
                overlapping_interval_motif_size = overlapping_interval.data
                larger_motif_size = max(motif_size, overlapping_interval_motif_size)
                if overlapping_interval.overlap_size(start_0based, end) >= 2*larger_motif_size:
                    overlapping_intervals.add((chrom, start_0based, end))
                    overlapping_intervals.add((chrom, overlapping_interval.begin, overlapping_interval.end))
                    break

            interval_trees[chrom].add(Interval(start_0based, end, data=len(motif)))

        if "EntireLocusMappability" in record:
            min_overall_mappability = min(min_overall_mappability, record["EntireLocusMappability"])
            if record["EntireLocusMappability"] == min_overall_mappability:
                min_overall_mappability_reference_region = reference_region
                min_overall_mappability_motif = motif
            mappability_bin = round(int(record["EntireLocusMappability"]*10)/10, 1)
            counters[f"mappability:{mappability_bin}"] += 1

    print("")
    print(f"Stats for {catalog_name}:")
    print(f"   {counters['total']:10,d} total loci")
    if counters["total_repeat_intervals"] == 0:
        return {
            "catalog": catalog_name,
            "total": f"{counters['total_repeat_intervals']:,d}",
        }

    total_genome_size = 3_088_286_401  # in hg38, including chrX, chrY, chrM
    print(f"   {counters['total_base_pairs_spanned_by_all_loci']:10,d} base pairs spanned by all loci ({counters['total_base_pairs_spanned_by_all_loci']/total_genome_size:0.3%} of the genome)")
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
    print(f"   Maximum locus size = {max_locus_size:5d}bp             @ {max_locus_size_reference_region} ({max_locus_size_motif})")
    if min_fraction_pure_bases_motif is not None:
        print(f"   Minimum fraction pure bases = {min_fraction_pure_bases:5.2f}      @ {min_fraction_pure_bases_reference_region} ({min_fraction_pure_bases_motif})")
    if min_fraction_pure_repeats_motif is not None:
        print(f"   Minimum fraction pure repeats = {min_fraction_pure_repeats:5.2f}    @ {min_fraction_pure_repeats_reference_region} ({min_fraction_pure_repeats_motif})")
    if min_overall_mappability_motif is not None:
        print(f"   Minimum overall mappability = {min_overall_mappability:5.2f}       @ {min_overall_mappability_reference_region} ({min_overall_mappability_motif})")
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
    for motif_sizes_per_locus_size_bin in "1x", "2x", "3x", "4x", "5x", "6x", "7x", "8x", "9x", "10-15x", "16-25x", "26-35x", "36-50x", "51+x":
        print(f"   {motif_sizes_per_locus_size_bin:>10s}: {counters[f'motif_sizes_per_locus_size:{motif_sizes_per_locus_size_bin}']:10,d} out of {counters['total_repeat_intervals']:10,d} ({counters[f'motif_sizes_per_locus_size:{motif_sizes_per_locus_size_bin}']/counters['total_repeat_intervals']:6.1%}) repeat intervals")

    if min_fraction_pure_bases_motif is not None:
        print("")
        print("Fraction pure bases distribution:")
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
        "percent_of_genome_spanned_by_loci": f"{counters['total_base_pairs_spanned_by_all_loci']/total_genome_size:0.3%}",
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
        "percent_pure_repeats": "%0.1f%%" % (100 * counters[f"fraction_pure_bases:1.0"] / counters["total_repeat_intervals"]),
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
        "count_pure_repeats": counters[f"fraction_pure_bases:1.0"],
        "count_trimmed": counters["trimmed"],
        "count_overlapping": len(overlapping_intervals),
        "min_motif_size": min_motif_size,
        "max_motif_size": max_motif_size,
        "min_locus_size": min_locus_size,
        "max_locus_size": max_locus_size,
    }

    for motif_sizes_per_locus_size_detailed_bin in (
        "0x", "1x", "2x", "3x", "4x", "5x", "6x", "7x", "8x", "9x", "10x",
        "11x", "12x", "13x", "14x", "15x", "16x", "17x", "18x", "19x", "20x",
        "21x", "22x", "23x", "24x", "25-50x", "51+x"
    ):
        result[f"motif_sizes_per_locus_size:{motif_sizes_per_locus_size_detailed_bin}"] = counters[f"motif_sizes_per_locus_size_detailed:{motif_sizes_per_locus_size_detailed_bin}"]


    return result



def main():
    args, parser = parse_args()

    stat_table_rows = []
    for path in args.variant_catalog_json_or_bed:
        print("-"*50)
        print(f"Parsing {path}")
        catalog_name = os.path.basename(path)

        file_iterator = get_variant_catalog_iterator(path)
        stats = compute_catalog_stats(catalog_name, file_iterator, verbose=args.verbose, show_progress_bar=args.show_progress_bar)
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