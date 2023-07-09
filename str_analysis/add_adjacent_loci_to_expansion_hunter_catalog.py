#!/usr/bin/env python3

"""This script takes an ExpansionHunter catalog and outputs a new catalog where adjacent loci have been combined
into a single spec with multiple adjacent repeats inorder to improve genotyping accuracy.
"""


import argparse
import collections
import intervaltree
import json
import os
import pybedtools
import pysam
import re
import tqdm

from str_analysis.utils.file_utils import open_file, file_exists
from str_analysis.utils.get_adjacent_repeats import get_adjacent_repeats, \
    MAX_DISTANCE_BETWEEN_REPEATS, MAX_TOTAL_ADJACENT_REGION_SIZE, MAX_OVERLAP_BETWEEN_ADJACENT_REPEATS
from str_analysis.utils.misc_utils import parse_interval


def process_input_record(
    record, reference_fasta, interval_tree,
    max_distance_between_adjacent_repeats=MAX_DISTANCE_BETWEEN_REPEATS,
    max_total_adjacent_region_size=MAX_TOTAL_ADJACENT_REGION_SIZE,
    max_overlap_between_adjacent_repeats=MAX_OVERLAP_BETWEEN_ADJACENT_REPEATS,
):
    """"Add adjacent loci to a record from the variant catalog

    Args:
        record (dict): A record from the input variant catalog. It must contain the following fields:
            - LocusId (string)
            - ReferenceRegion (a string representing 0-based reference coordinates of the repeat locus)
            - LocusStructure (a string containing only one repeat unit - ie. "(CAG)*" )
        reference_fasta (pysam.FastaFile): A pysam.FastaFile object for the reference genome
        interval_tree (dict): an IntervalTree containing 0-based intervals that represent the reference start and end
            coordinates of all TR loci found on the same chromosome as the given record. Adjacent repeats will be
            queried from this data structure.
        max_distance_between_adjacent_repeats (int)
        max_total_adjacent_region_size (int)
        max_overlap_between_adjacent_repeats (int)

    Return:
        dict: An copy of the input record, modified to include any adjacent repeats found in the interval_tree.
    """
    locus_id = record["LocusId"]
    locus_interval_0based = record["ReferenceRegion"]
    if not isinstance(locus_interval_0based, str):
        raise ValueError(f"ReferenceRegion is not a string: {record['LocusId']}")

    chrom, start_0based, end_1based = parse_interval(locus_interval_0based)
    repeat_unit = record["LocusStructure"].strip("()*")
    if record["LocusStructure"].count("(") > 1 or record["LocusStructure"].count(")") > 1 or record["LocusStructure"].count("*") > 1:
        raise ValueError(f"LocusStructure must contain only 1 repeat unit: {record['LocusId']}")

    output_record = record.copy()

    adjacent_repeats_left, adjacent_repeats_locus_structure, adjacent_repeats_right = get_adjacent_repeats(
        locus_interval_0based,
        repeat_unit=repeat_unit,
        pysam_fasta_file=reference_fasta,
        interval_tree_0based=interval_tree,
        max_distance_between_adjacent_repeats=max_distance_between_adjacent_repeats,
        max_total_adjacent_region_size=max_total_adjacent_region_size,
        max_overlap_between_adjacent_repeats=max_overlap_between_adjacent_repeats,
    )

    if adjacent_repeats_left or adjacent_repeats_right:
        reference_regions = []
        variant_ids = []
        variant_types = []
        for adjacent_repeat_label, adjacent_repeat in [
            (f"{locus_id}_ADJACENT_LEFT_{i+1}", adjacent_repeat_left) for i, adjacent_repeat_left in enumerate(adjacent_repeats_left)
        ] + [(locus_id, f"{chrom}:{start_0based}-{end_1based}-{repeat_unit}")] + [
            (f"{locus_id}_ADJACENT_RIGHT_{i+1}", adjacent_repeat_right) for i, adjacent_repeat_right in enumerate(adjacent_repeats_right)
        ]:
            adjacent_repeat_chrom, adjacent_repeat_start_0based, adjacent_repeat_end_1based, adjacent_repeat_unit = re.split("[:-]", adjacent_repeat)
            reference_regions.append(f"{adjacent_repeat_chrom}:{adjacent_repeat_start_0based}-{adjacent_repeat_end_1based}")
            variant_ids.append(adjacent_repeat_label)
            variant_types.append("Repeat")

        output_record["ReferenceRegion"] = reference_regions
        output_record["VariantId"] = variant_ids
        output_record["VariantType"] = variant_types
        output_record["LocusStructure"] = adjacent_repeats_locus_structure

    return output_record


def get_interval_tree_for_chrom(adjacent_loci_bed_file_path, chrom, start, end):
    """Load intervals from the given bed file into an IntervalTree for faster overlap queries. The bed file should
    contain all TR loci found in the reference genome.

    Args:
        adjacent_loci_bed_file_path (str): Path to a bed file containing adjacent loci
        chrom (str): Chromosome name
        start (int): Start coordinate of the chromosome region that spans all adjacent repeats that we want to load
        end (int): End coordinate of the chromosome region that spans all adjacent repeats that we want to load

    Return:
        intervaltree.IntervalTree: An IntervalTree containing 0-based intervals that represent the reference start and
            end coordinates of all TR loci found in the reference genome on the given chromosome.
    """
    interval_tree = intervaltree.IntervalTree()
    already_added_to_interval_tree = set() # avoid duplicates
    t = pybedtools.BedTool(os.path.expanduser(adjacent_loci_bed_file_path))
    print(f"Loading potential adjacent repeats from a {end-start+1:,d}bp region ({chrom}:{start:,d}-{end:,d}) in "
          f"{adjacent_loci_bed_file_path}...")
    for bed_record in tqdm.tqdm(t.all_hits(pybedtools.Interval(chrom, start, end)), unit=" bed records"):
        start_0based = bed_record.start
        end_1based = bed_record.end
        repeat_unit = bed_record.name

        key = (chrom, int(start_0based), int(end_1based), repeat_unit)
        if key not in already_added_to_interval_tree:
            interval_tree.append(intervaltree.Interval(int(start_0based), int(end_1based), data=repeat_unit))
            already_added_to_interval_tree.add(key)

    if not already_added_to_interval_tree:
        print(f"WARNING: no adjacent loci found in {adjacent_loci_bed_file_path} for chromosome {chrom}")

    return interval_tree

def get_min_and_max_coords_for_chrom(input_catalog, chrom):
    """Get the minimum and maximum coordinates of all TR loci found in the given chromosome in the input catalog

    Args:
        input_catalog (list): Input catalog as a list of dicts
        chrom (str): Chromosome name

    Return:
        2-tuple (min_coord, max_coord): the minimum and maximum coordinates of all TR loci found on the given chromosome
    """
    min_coord = 10**9
    max_coord = 0
    for record in input_catalog:
        reference_regions = record["ReferenceRegion"]
        if not isinstance(reference_regions, list):
            reference_regions = [reference_regions]
        for reference_region in reference_regions:
            if not reference_region.startswith(f"{chrom}:"):
                continue
            chrom, start, end = parse_interval(reference_region)
            min_coord = min(min_coord, start)
            max_coord = max(max_coord, end)

    return min_coord, max_coord


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref-fasta", help="Reference fasta path", required=True)
    parser.add_argument("-d", "--max-distance-between-adjacent-repeats",
                        default=MAX_DISTANCE_BETWEEN_REPEATS,
                        help="Loci no more than this many base pairs from each other will be combined into a single "
                             "locus, including any non-repetative bases between them")
    parser.add_argument("--max-overlap-between-adjacent-loci",
                        default=MAX_OVERLAP_BETWEEN_ADJACENT_REPEATS,
                        help="TandemRepeatFinder sometimes outputs loci that overlap by a few base pairs. For "
                             "ExpansionHunter catalogs that are based on TandemRepeatFinder output, this "
                             "controls how slightly overlapping loci are combined")
    parser.add_argument("--max-total-adjacent-region-size",
                        default=MAX_TOTAL_ADJACENT_REGION_SIZE,
                        help="The maximum size of the region spanning all adjacent repeats that we will load from "
                             "the adjacent loci bed file")
    parser.add_argument("--source-of-adjacent-loci",
                        default="gs://str-truth-set/hg38/ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_6bp.bed.gz",
                        help="The local or Google Storage (gs://) path of a .bed file that specifies repeat intervals "
                             "in the reference genome. It should have 0-based coordinates, and the name field (column 4) "
                             "should contain the repeat motif. This can be a large catalog of loci generated using "
                             "TandemRepeatFinder. This .bed file will be used as the source "
                             "for adjacent loci.", required=True)
    parser.add_argument("-o", "--output-catalog", help="Path where to write the output catalog. If not specified, "
                                                       "it will be based on the input catalog path")
    parser.add_argument("input_catalog", help="ExpansionHunter catalog json file path")
    args = parser.parse_args()

    for path in args.ref_fasta, args.input_catalog:
        if not file_exists(path):
            parser.error(f"{path} not found")

    ref_fasta = pysam.FastaFile(os.path.abspath(os.path.expanduser(args.ref_fasta)))
    ref_fasta_chromosomes = set(ref_fasta.references)

    with open_file(args.input_catalog) as f:
        input_catalog = json.load(f)

    # sort the input catalog by its ReferenceRegion chromosome
    input_catalog.sort(key=lambda record: (
        record["ReferenceRegion"][0].split(":")[0]
        if isinstance(record["ReferenceRegion"], list) else
        record["ReferenceRegion"].split(":")[0]))

    output_catalog = []
    already_added_to_output_catalog = set()  # avoid duplicates
    current_chrom = None
    current_interval_tree = None
    counters = collections.defaultdict(int)
    for record_i, record in enumerate(input_catalog):
        counters["total"] += 1
        missing_keys = {"ReferenceRegion", "LocusStructure", "LocusId"} - set(record.keys())
        if missing_keys:
            print(f"WARNING: record {record} is missing keys: {missing_keys}. Skipping...")
            continue

        if isinstance(record["ReferenceRegion"], list):
            counters["already had adjacent loci"] += 1
            print(f"WARNING: record {record} already has adjacent intervals specified in the input variant catalog. "
                  f"Will add it to the output variant catalog without modification...")
            output_catalog.append(record)
            continue

        chrom = record["ReferenceRegion"].split(":")[0]
        if chrom not in ref_fasta_chromosomes:
            print(f"ERROR: chrom. '{chrom}' from ReferenceRegion '{record['ReferenceRegion']}' in input variant catalog"
                  f" record #{record_i + 1} doesn't match any of the chromosomes in {args.ref_fasta}. Skipping...")
            continue

        if current_chrom != chrom:
            current_chrom = chrom
            min_coord, max_coord = get_min_and_max_coords_for_chrom(input_catalog, chrom)
            min_coord -= args.max_total_adjacent_region_size
            max_coord += args.max_total_adjacent_region_size
            current_interval_tree = get_interval_tree_for_chrom(args.source_of_adjacent_loci, chrom, min_coord, max_coord)

        output_record = process_input_record(record, ref_fasta, current_interval_tree)
        if isinstance(output_record["ReferenceRegion"], list):
            counters["added adjacent loci"] += 1
            output_record_key = output_record["LocusStructure"] + ";" + ";".join(output_record["ReferenceRegion"])
        else:
            output_record_key = output_record["LocusStructure"] + ";" + output_record["ReferenceRegion"]

        if output_record_key not in already_added_to_output_catalog:
            already_added_to_output_catalog.add(output_record_key)
            output_catalog.append(output_record)

    print(f"Added adjacent loci to {counters['added adjacent loci']:,d} out of {counters['total']:,d} records ("
          f"{counters['added adjacent loci'] / counters['total']:.1%})")
    if counters["already had adjacent loci"]:
        print(f"{counters['already had adjacent loci']:,d}  out of {counters['total']:,d} records already had adjacent "
              f"loci specified in the input catalog, and so were added to the output catalog without modification")

    if not args.output_catalog:
        args.output_catalog = re.sub(".json$", "", args.input_catalog) + ".combined_adjacent_loci.json"

    with open(args.output_catalog, "wt") as f:
        json.dump(output_catalog, f, indent=2)

    print(f"Wrote {len(output_catalog):,d} total records to {args.output_catalog}")


if __name__ == '__main__':
  main()