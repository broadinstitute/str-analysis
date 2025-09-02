import argparse
import collections
import ijson
import os
import pandas as pd
import re
import statistics
import tqdm
from intervaltree import IntervalTree, Interval
from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.eh_catalog_utils import parse_motifs_from_locus_structure, get_variant_catalog_iterator
from str_analysis.utils.file_utils import file_exists

ACGT_REGEX = re.compile("^[ACGT]+$", re.IGNORECASE)
ACGTN_REGEX = re.compile("^[ACGTN]+$", re.IGNORECASE)

GENE_REGIONS = ["5utr", "3utr", "cds", "exon", "intergenic", "intron", "promoter"]

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
    #min_fraction_pure_repeats = 1
    #min_fraction_pure_repeats_reference_region = None
    #min_fraction_pure_repeats_motif = None
    min_overall_mappability = 1
    min_overall_mappability_reference_region = None
    min_overall_mappability_motif = None

    if show_progress_bar:
        records = tqdm.tqdm(records, unit=" records", unit_scale=True)

    interval_trees = collections.defaultdict(IntervalTree)  # used to check for overlap between records in the catalog
    overlapping_intervals = set()
    counters = collections.defaultdict(int)

    locus_sizes_by_motif_size = collections.defaultdict(list)  # used to compute the median locus sizes for each motif size
    reference_repeat_purity_by_motif_size = collections.defaultdict(list)  # used to compute the mean base purity for each motif size
    mappability_by_motif_size = collections.defaultdict(list)  # used to compute the mean mappability for each motif size

    has_gene_annotations = False
    for record in records:
        counters["total"] += 1
        motifs = parse_motifs_from_locus_structure(record["LocusStructure"])
        if isinstance(record["ReferenceRegion"], list):
            reference_regions = record["ReferenceRegion"]
            variant_types = record["VariantType"]
            fraction_pure_bases = record.get("ReferenceRepeatPurity", [None]*len(reference_regions))
            #fraction_pure_repeats = record.get("FractionPureRepeats", [None]*len(reference_regions))

            counters["loci_with_adjacent_repeats"] += 1
        else:
            reference_regions = [record["ReferenceRegion"]]
            variant_types = [record["VariantType"]]
            fraction_pure_bases = [record.get("ReferenceRepeatPurity")]
            #fraction_pure_repeats = [record.get("FractionPureRepeats")]

        for motif, reference_region, variant_type, fraction_pure_bases in zip(
                motifs, reference_regions, variant_types, fraction_pure_bases):
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

            counters["total_base_pairs_spanned_by_all_loci"] += locus_size
            if locus_size == max_locus_size:
                max_locus_size_reference_region = reference_region
                max_locus_size_motif = motif

            if fraction_pure_bases is not None:
                min_fraction_pure_bases = min(min_fraction_pure_bases, fraction_pure_bases)
                if fraction_pure_bases == min_fraction_pure_bases:
                    min_fraction_pure_bases_reference_region = reference_region
                    min_fraction_pure_bases_motif = motif
            #if fraction_pure_repeats is not None:
            #    min_fraction_pure_repeats = min(min_fraction_pure_repeats, fraction_pure_repeats)
            #    if fraction_pure_repeats == min_fraction_pure_repeats:
            #        min_fraction_pure_repeats_reference_region = reference_region
            #        min_fraction_pure_repeats_motif = motif

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

            if fraction_pure_bases is not None:
                fraction_pure_bases_bin = round(int(fraction_pure_bases*10)/10, 1)
                counters[f"fraction_pure_bases:{fraction_pure_bases_bin}"] += 1

            # check for overlap with other loci in the catalog
            for overlapping_interval in interval_trees[chrom].overlap(start_0based, max(end, start_0based + 1)):
                overlapping_interval_motif_size = overlapping_interval.data
                larger_motif_size = max(motif_size, overlapping_interval_motif_size)
                if overlapping_interval.overlap_size(start_0based, end) >= 2*larger_motif_size:
                    overlapping_intervals.add((chrom, start_0based, end))
                    overlapping_intervals.add((chrom, overlapping_interval.begin, overlapping_interval.end))
                    break

            interval_trees[chrom].add(Interval(start_0based, max(end, start_0based + 1), data=len(motif)))

        if "FlanksAndLocusMappability" in record:
            min_overall_mappability = min(min_overall_mappability, record["FlanksAndLocusMappability"])
            if record["FlanksAndLocusMappability"] == min_overall_mappability:
                min_overall_mappability_reference_region = reference_region
                min_overall_mappability_motif = motif
            mappability_bin = round(int(record["FlanksAndLocusMappability"]*10)/10, 1)
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
    print(f"   Max locus size = {max_locus_size:7,d}bp           @ {max_locus_size_reference_region} ({max_locus_size_motif})")
    if min_fraction_pure_bases_motif is not None:
        print(f"   Min reference repeat purity   = {min_fraction_pure_bases:5.2f}    @ {min_fraction_pure_bases_reference_region} ({min_fraction_pure_bases_motif})")
    #if min_fraction_pure_repeats_motif is not None:
    #    print(f"   Min fraction pure repeats = {min_fraction_pure_repeats:5.2f}    @ {min_fraction_pure_repeats_reference_region} ({min_fraction_pure_repeats_motif})")
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
        print("Adding ")
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