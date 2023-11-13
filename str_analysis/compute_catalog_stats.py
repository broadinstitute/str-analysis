import argparse
import collections
import ijson
import os
import pandas as pd
import re
import tqdm
from intervaltree import IntervalTree, Interval
from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.file_utils import open_file, file_exists

ACGT_REGEX = re.compile("^[ACGT]+$", re.IGNORECASE)


def parse_args():
    parser = argparse.ArgumentParser(description="Compute and print stats for annotated repeat catalogs")
    parser.add_argument("--verbose", action="store_true", help="Print more information about what the script is doing")
    parser.add_argument("annotated_variant_catalog_json", nargs="+",
                        help="Repeat catalog(s) output by the annotate_and_filter_str_catalog.py script")

    args = parser.parse_args()

    for file_path in args.annotated_variant_catalog_json:
        if not file_exists(os.path.expanduser(file_path)):
            parser.error(f"File not found: {file_path}")

    return args, parser


def compute_catalog_stats(catalog_name, records):
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
    min_fraction_pure_bases = 1
    min_fraction_pure_bases_reference_region = None
    min_fraction_pure_bases_motif = None
    min_fraction_pure_repeats = 1
    min_fraction_pure_repeats_reference_region = None
    min_fraction_pure_repeats_motif = None
    min_overall_mappability = 1
    min_overall_mappability_reference_region = None
    min_overall_mappability_motif = None

    interval_trees = collections.defaultdict(IntervalTree)  # used to check for overlap between records in the catalog
    overlapping_intervals = set()
    counters = collections.defaultdict(int)
    for record in records:
        counters["total"] += 1
        if isinstance(record["ReferenceRegion"], list):
            motifs = []
            for m in record["LocusStructure"].strip("()*+").split(")"):
                if "(" in m:
                    m = m.split("(")[1]
                motifs.append(m.strip("()*+"))
            reference_regions = record["ReferenceRegion"]
            variant_types = record["VariantType"]
            fraction_pure_bases = record["FractionPureBases"]
            fraction_pure_repeats = record["FractionPureRepeats"]

            counters["loci_with_adjacent_repeats"] += 1
        else:
            motifs = [record["LocusStructure"].strip("()*+")]
            reference_regions = [record["ReferenceRegion"]]
            variant_types = [record["VariantType"]]
            fraction_pure_bases = [record["FractionPureBases"]]
            fraction_pure_repeats = [record["FractionPureRepeats"]]

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
            min_num_repeats_in_locus = min(min_num_repeats_in_locus, locus_size / motif_size)
            max_num_repeats_in_locus = max(max_num_repeats_in_locus, locus_size / motif_size)

            min_fraction_pure_bases = min(min_fraction_pure_bases, fraction_pure_bases)
            if fraction_pure_bases == min_fraction_pure_bases:
                min_fraction_pure_bases_reference_region = reference_region
                min_fraction_pure_bases_motif = motif
            min_fraction_pure_repeats = min(min_fraction_pure_repeats, fraction_pure_repeats)
            if fraction_pure_repeats == min_fraction_pure_repeats:
                min_fraction_pure_repeats_reference_region = reference_region
                min_fraction_pure_repeats_motif = motif

            if motif_size == 1:
                counters["homopolymers"] += 1
            if locus_size % motif_size == 0:
                counters["trimmed"] += 1

            if not ACGT_REGEX.match(motif):
                counters["non_acgt_motifs"] += 1

            motif_size_bin = f"{motif_size}bp" if motif_size <= 6 else "7-24bp" if motif_size <= 24 else "25+bp"
            fraction_pure_bases_bin = round(int(fraction_pure_repeats*10)/10, 1)

            counters[f"motif_size:{motif_size_bin}"] += 1
            counters[f"fraction_pure_bases:{fraction_pure_bases_bin}"] += 1

            # check for overlap
            for overlapping_interval in interval_trees[chrom].overlap(start_0based, end):
                larger_motif_size = max(motif_size, overlapping_interval.data["motif_size"])
                if overlapping_interval.overlap_size(start_0based, end) >= 2*larger_motif_size:
                    overlapping_intervals.add((chrom, start_0based, end))
                    overlapping_intervals.add((chrom, overlapping_interval.begin, overlapping_interval.end))
                    break

            interval_trees[chrom].add(Interval(start_0based, end, data={"motif_size": len(motif)}))

        min_overall_mappability = min(min_overall_mappability, record["OverallMappability"])
        if record["OverallMappability"] == min_overall_mappability:
            min_overall_mappability_reference_region = reference_region
            min_overall_mappability_motif = motif
        mappability_bin = round(int(record["OverallMappability"]*10)/10, 1)
        counters[f"mappability:{mappability_bin}"] += 1

    print("")
    print(f"Stats for {catalog_name}:")
    print(f"   {counters['total']:10,d} total loci")
    print(f"   {counters['loci_with_adjacent_repeats']:10,d} out of {counters['total']:10,d} ({counters['loci_with_adjacent_repeats']/counters['total']:6.1%}) loci have adjacent repeats")
    print(f"   {counters['total_repeat_intervals']:10,d} total repeat intervals")
    print(f"   {counters['trimmed']:10,d} out of {counters['total_repeat_intervals']:10,d} ({counters['trimmed']/counters['total_repeat_intervals']:6.1%}) repeat interval size is an integer multiple of the motif size (aka. trimmed)")
    print(f"   {counters['homopolymers']:10,d} out of {counters['total_repeat_intervals']:10,d} ({counters['homopolymers']/counters['total_repeat_intervals']:6.1%}) repeat intervals are homopolymers")
    print(f"   {len(overlapping_intervals):10,d} out of {counters['total_repeat_intervals']:10,d} ({len(overlapping_intervals)/counters['total_repeat_intervals']:6.1%}) repeat intervals overlap each other by at least two motif lengths")
    if counters["non_acgt_motifs"]:
        print(f"   {counters['non_acgt_motifs']:10,d} out of {counters['total_repeat_intervals']:10,d} ({counters['non_acgt_motifs']/counters['total_repeat_intervals']:6.1%}) repeat intervals have non-ACGT motifs")

    if len(overlapping_intervals) > 0:
        examples = list(overlapping_intervals)
        examples = ", ".join(f"{chrom}:{start}-{end}" for chrom, start, end in examples[:10])
        print(f"Examples of overlapping repeats: {examples}")

    print("")
    print("Ranges:")
    print(f"   Motif size range: {min_motif_size}-{max_motif_size}bp")
    print(f"   Locus size range: {min_locus_size}-{max_locus_size}bp")
    print(f"   Num repeats range: {min_num_repeats_in_locus:0.1f}-{max_num_repeats_in_locus:0.1f} repeats")
    print("")
    print(f"   Minimum fraction pure bases = {min_fraction_pure_bases:.2f}      @ {min_fraction_pure_bases_reference_region} ({min_fraction_pure_bases_motif})")
    print(f"   Minimum fraction pure repeats = {min_fraction_pure_repeats:.2f}    @ {min_fraction_pure_repeats_reference_region} ({min_fraction_pure_repeats_motif})")
    print(f"   Minimum overall mappability = {min_overall_mappability:.2f}       @ {min_overall_mappability_reference_region} ({min_overall_mappability_motif})")
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
    print("Fraction pure bases distribution:")
    for fraction_pure_bases_bin in 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1:
        print(f"   {fraction_pure_bases_bin:10.1f}: {counters[f'fraction_pure_bases:{fraction_pure_bases_bin:.1f}']:10,d} out of {counters['total_repeat_intervals']:10,d} ({counters[f'fraction_pure_bases:{fraction_pure_bases_bin:.1f}']/counters['total_repeat_intervals']:6.1%}) repeat intervals")

    print("")
    print("Mappability distribution:")
    for mappability_bin in 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1:
        print(f"   {mappability_bin:10.1f}: {counters[f'mappability:{mappability_bin:.1f}']:10,d} out of {counters['total']:10,d} ({counters[f'mappability:{mappability_bin:.1f}']/counters['total']:6.1%}) loci")

    result = {
        "Total": f"{counters['total_repeat_intervals']:,d}",
        "chrX": f"{counters['chrX']:,d} ({counters['chrX']/counters['total_repeat_intervals']:0.1%})",
        "chrY": f"{counters['chrY']:,d} ({counters['chrY']/counters['total_repeat_intervals']:0.1%})",
        "Motif size range": f"{min_motif_size}-{max_motif_size}bp",
        "Locus size range": f"{min_locus_size}-{max_locus_size}bp",
        "Num repeats range": f"{min_num_repeats_in_locus:0.1f}-{max_num_repeats_in_locus:0.1f} repeats",
        "Homopolymers": "%0.1f" % (100 * counters["motif_size:1bp"] / counters["total_repeat_intervals"]),
        "2bp motifs": "%0.1f" % (100 * counters["motif_size:2bp"] / counters["total_repeat_intervals"]),
        "3bp motifs": "%0.1f" % (100 * counters["motif_size:3bp"] / counters["total_repeat_intervals"]),
        "4bp motifs": "%0.1f" % (100 * counters["motif_size:4bp"] / counters["total_repeat_intervals"]),
        "5bp motifs": "%0.1f" % (100 * counters["motif_size:5bp"] / counters["total_repeat_intervals"]),
        "6bp motifs": "%0.1f" % (100 * counters["motif_size:6bp"] / counters["total_repeat_intervals"]),
        "7+bp motifs": "%0.1f" % (100 * (counters["motif_size:7-24bp"] + counters["motif_size:25+bp"])/ counters["total_repeat_intervals"]),
        "Pure repeats": "%0.1f" % (100 * counters[f"fraction_pure_bases:1.0"] / counters["total_repeat_intervals"]),
        "Trimmed": "%0.1f" % (100 * counters["trimmed"] / counters["total_repeat_intervals"]),
        "Overlapping": f"{len(overlapping_intervals):,d} ({(100 * len(overlapping_intervals) / counters['total_repeat_intervals']):0.1%})",

        "": "",
        "total_loci": counters["total"],
        "loci_with_adjacent_repeats": counters["loci_with_adjacent_repeats"],
        "total_repeats": counters["total_repeat_intervals"],
        "homopolymer_repeats": counters["homopolymers"],
        "trimmed_repeats": counters["trimmed"],
        "overlapping_repeat_intervals": len(overlapping_intervals),
        "chrX_repeats": counters["chrX"],
        "chrY_repeats": counters["chrY"],
        "chrM_repeats": counters["chrM"],
        "min_motif_size": min_motif_size,
        "max_motif_size": max_motif_size,
        "min_locus_size": min_locus_size,
        "max_locus_size": max_locus_size,
        "min_fraction_pure_bases": min_fraction_pure_bases,
        "min_fraction_pure_bases_reference_region": min_fraction_pure_bases_reference_region,
        "min_fraction_pure_bases_motif": min_fraction_pure_bases_motif,
        "min_fraction_pure_repeats": min_fraction_pure_repeats,
        "min_fraction_pure_repeats_reference_region": min_fraction_pure_repeats_reference_region,
        "min_fraction_pure_repeats_motif": min_fraction_pure_repeats_motif,
        "motif_size_distribution": {
            "1bp": counters["motif_size:1bp"],
            "2bp": counters["motif_size:2bp"],
            "3bp": counters["motif_size:3bp"],
            "4bp": counters["motif_size:4bp"],
            "5bp": counters["motif_size:5bp"],
            "6bp": counters["motif_size:6bp"],
            "7-24bp": counters["motif_size:7-24bp"],
            "25+bp": counters["motif_size:25+bp"],
        },
        "fraction_pure_bases_distribution": {},
        "mappability_distribution": {},
    }

    for fraction in 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1:
        result["fraction_pure_bases_distribution"][f"{fraction:.1f}"] = counters[f"fraction_pure_bases:{fraction:.1f}"]
        result["mappability_distribution"][f"{fraction:.1f}"] = counters[f"mappability:{fraction:.1f}"]

    return result


def main():
    args, parser = parse_args()

    stat_table_rows = []
    for path in args.annotated_variant_catalog_json:
        print("-"*50)
        print(f"Parsing {path}")
        catalog_name = os.path.basename(path).replace(".json", "")
        with open_file(path, "rt") as f:
            file_iterator = ijson.items(f, "item")
            if args.verbose:
                file_iterator = tqdm.tqdm(file_iterator, unit=" records", unit_scale=True)

            stats = compute_catalog_stats(catalog_name, file_iterator)
        stat_table_rows.append(stats)

    output_path = f"variant_catalog_stats.{len(args.annotated_variant_catalog_json)}_catalogs.tsv"
    df = pd.DataFrame(stat_table_rows)
    df.drop(columns=[
        "motif_size_distribution",
        "fraction_pure_bases_distribution",
        "mappability_distribution",
    ], inplace=True)
    df.to_csv(output_path, sep="\t", index=False)
    print(f"Wrote {len(stat_table_rows):,d} rows to {output_path}")


if __name__ == "__main__":
    main()