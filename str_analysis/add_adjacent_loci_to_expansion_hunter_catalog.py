#!/usr/bin/env python3

"""This script takes an ExpansionHunter catalog and outputs a new catalog where adjacent loci have been combined
into a single spec with multiple adjacent repeats inorder to improve genotyping accuracy.
"""


import argparse
import collections
import gzip
import intervaltree
import simplejson as json
import ijson
import os
import pybedtools
import pysam
import re

from str_analysis.utils.file_utils import open_file, file_exists
from str_analysis.utils.get_adjacent_repeats import get_adjacent_repeats, \
    MAX_DISTANCE_BETWEEN_REPEATS, MAX_TOTAL_ADJACENT_REGION_SIZE, MAX_OVERLAP_BETWEEN_ADJACENT_REPEATS
from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.file_utils import download_local_copy
from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif


def process_input_record(
    record, reference_fasta, interval_tree,
    max_distance_between_adjacent_repeats=MAX_DISTANCE_BETWEEN_REPEATS,
    max_total_adjacent_region_size=MAX_TOTAL_ADJACENT_REGION_SIZE,
    max_overlap_between_adjacent_repeats=MAX_OVERLAP_BETWEEN_ADJACENT_REPEATS,
    max_adjacent_repeats=None,
    add_extra_info_to_locus_id=False,
    add_extra_fields=False,
    set_variant_ids_to_locus_id_plus_motif=True,
):
    """Add adjacent loci to a record from the variant catalog

    Args:
        record (dict): A record from the input variant catalog. It must contain the following fields:
            - LocusId (string)
            - ReferenceRegion (a string representing 0-based reference coordinates of the repeat locus)
            - LocusStructure (a string containing only one repeat unit - ie. "(CAG)*" )
        reference_fasta (pysam.FastaFile): A pysam.FastaFile object for the reference genome
        interval_tree (intervaltree.IntervalTree): an IntervalTree containing 0-based intervals that represent the
            reference start and end coordinates of all TR loci found on the same chromosome as the given record.
            Adjacent repeats will be queried from this data structure.
        max_distance_between_adjacent_repeats (int)
        max_total_adjacent_region_size (int)
        max_overlap_between_adjacent_repeats (int)
        max_adjacent_repeats (int) Max adjacent repeats to add.
        add_extra_info_to_locus_id (bool) If True, add info about adjacent loci to the locus id in the output variant
        add_extra_fields (bool) If True, add extra fields to the output record that record the number of base pairs
            between the main repeat and any adjacent repeats.
        set_variant_ids_to_locus_id_plus_motif (bool) If True, set the variant id of adjacent loci to the locus id
            followed by the motif. This matches the naming convention in the official Illumina catalog. For
            example: FXN_A is the adjacent poly-A repeat next to the main FXN locus. If False, use a longer, more
            descriptive variant id like "FXN_ADJACENT_LEFT_1_A".

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
        max_adjacent_repeats=max_adjacent_repeats,
    )

    if adjacent_repeats_left or adjacent_repeats_right:
        reference_regions = []
        variant_ids = []
        variant_types = []
        if add_extra_fields:
            previous_reference_region_end = None
            distances_between_reference_regions = []
            reference_region_motifs = []

        for variant_id, adjacent_repeat in [
            (f"{locus_id}_ADJACENT_LEFT_{i+1}", adjacent_repeat_left) for i, adjacent_repeat_left in enumerate(adjacent_repeats_left)
        ] + [(locus_id, f"{chrom}:{start_0based}-{end_1based}-{repeat_unit}")] + [
            (f"{locus_id}_ADJACENT_RIGHT_{i+1}", adjacent_repeat_right) for i, adjacent_repeat_right in enumerate(adjacent_repeats_right)
        ]:
            adjacent_repeat_chrom, adjacent_repeat_start_0based, adjacent_repeat_end_1based, adjacent_repeat_unit = re.split("[:-]", adjacent_repeat)
            reference_regions.append(f"{adjacent_repeat_chrom}:{adjacent_repeat_start_0based}-{adjacent_repeat_end_1based}")
            if variant_id == locus_id:
                variant_ids.append(locus_id)
            elif set_variant_ids_to_locus_id_plus_motif:
                variant_ids.append(f"{locus_id}_{adjacent_repeat_unit}")
            else:
                variant_ids.append(f"{variant_id}_{adjacent_repeat_unit}")

            variant_types.append("Repeat")
            if add_extra_info_to_locus_id or add_extra_fields:
                reference_region_motifs.append(adjacent_repeat_unit)
                if previous_reference_region_end is not None:
                    distances_between_reference_regions.append(
                        int(adjacent_repeat_start_0based) - int(previous_reference_region_end))
                previous_reference_region_end = int(adjacent_repeat_end_1based)

        output_record["ReferenceRegion"] = reference_regions
        output_record["VariantId"] = variant_ids
        output_record["VariantType"] = variant_types
        output_record["LocusStructure"] = adjacent_repeats_locus_structure

        if add_extra_info_to_locus_id:
            output_record["LocusId"] += f";{len(reference_regions) - 1}_adjacent_repeats"
            if len(reference_regions) > 1:
                output_record["LocusId"] += f";max_dist_{max(distances_between_reference_regions)}bp"
        if add_extra_fields and distances_between_reference_regions:
            output_record["DistancesBetweenReferenceRegions"] = distances_between_reference_regions
            output_record["ReferenceRegionMotifs"] = reference_region_motifs
            output_record["MaxDistanceBetweenReferenceRegions"] = max(distances_between_reference_regions)

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
    if adjacent_loci_bed_file_path.startswith("gs://"):
        adjacent_loci_bed_file_path = download_local_copy(adjacent_loci_bed_file_path)

    interval_tree = intervaltree.IntervalTree()
    already_added_to_interval_tree = set()  # avoid duplicates
    t = pybedtools.BedTool(os.path.expanduser(adjacent_loci_bed_file_path))
    print(f"Parsing potential adjacent repeats from a {end-start+1:,d}bp region ({chrom}:{start:,d}-{end:,d}) in "
          f"{adjacent_loci_bed_file_path}...")
    for bed_record in t.all_hits(pybedtools.Interval(chrom, start, end)):
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


def toggle_chromosome_prefix(variant_catalog_record):
    """Toggle the ReferenceRegion field in the given variant catalog record between having and not having a 'chr' prefix

    Args:
        variant_catalog_record (dict): A variant catalog record
    """
    if isinstance(variant_catalog_record["ReferenceRegion"], list):
        for i in range(len(variant_catalog_record["ReferenceRegion"])):
            if variant_catalog_record["ReferenceRegion"][i].startswith("chr"):
                variant_catalog_record["ReferenceRegion"][i] = variant_catalog_record["ReferenceRegion"][i].replace("chr", "")
            else:
                variant_catalog_record["ReferenceRegion"][i] = "chr" + variant_catalog_record["ReferenceRegion"][i]

        chrom = variant_catalog_record["ReferenceRegion"][0].split(":")[0]
    else:
        if variant_catalog_record["ReferenceRegion"].startswith("chr"):
            variant_catalog_record["ReferenceRegion"] = variant_catalog_record["ReferenceRegion"].replace("chr", "")
        else:
            variant_catalog_record["ReferenceRegion"] = "chr" + variant_catalog_record["ReferenceRegion"]

        chrom = variant_catalog_record["ReferenceRegion"].split(":")[0]

    return chrom


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
    parser.add_argument("-n", type=int, help="Process at most this many loci. Useful for testing.")
    parser.add_argument("-d", "--max-distance-between-adjacent-repeats",
                        type=int,
                        default=MAX_DISTANCE_BETWEEN_REPEATS,
                        help="Loci no more than this many base pairs from each other will be combined into a single "
                             "locus, including any non-repetitive bases between them")
    parser.add_argument("--max-overlap-between-adjacent-repeats",
                        type=int,
                        default=MAX_OVERLAP_BETWEEN_ADJACENT_REPEATS,
                        help="TandemRepeatFinder sometimes outputs loci that overlap by a few base pairs. For "
                             "ExpansionHunter catalogs that are based on TandemRepeatFinder output, this "
                             "controls how slightly overlapping loci are combined")
    parser.add_argument("--max-total-adjacent-region-size",
                        type=int,
                        default=MAX_TOTAL_ADJACENT_REGION_SIZE,
                        help="Look for adjacent repeats that are no more than this many base pairs to the left or to "
                             "the right of the main repeat.")
    parser.add_argument("--max-adjacent-repeats",
                        type=int,
                        help="Add no more than this many of the closest adjacent repeats on either side. For example, "
                             "if set to 1, no more than 1 adjacent repeat will be added to the left and to the right "
                             "of the main locus.")
    parser.add_argument("--source-of-adjacent-loci",
                        required=True,
                        default="gs://str-truth-set/hg38/ref/other/"
                                "repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_6bp.bed.gz",
                        help="The local or Google Storage (gs://) path of a .bed file that specifies repeat intervals "
                             "in the reference genome. It should have 0-based coordinates and the name field (column 4)"
                             " should contain the repeat motif. This can be a large catalog of loci generated using "
                             "TandemRepeatFinder. This .bed file will be used as the source for adjacent loci.")
    parser.add_argument("--add-extra-info-to-locus-id",
                        action="store_true",
                        help="Add info about adjacent loci to the locus id in the output variant catalog. This will "
                             "include the number of adjacent loci and the max distance between them.")
    parser.add_argument("--add-extra-fields",
                        action="store_true",
                        help="Add extra fields to each locus in the output variant catalog that record the number of "
                             "base pairs between any adjacent loci and the main locus in an easy-to-parse format")
    parser.add_argument("-x", "--dont-add-adjacent-loci-to-this-locus-id", action="append", help="Don't compute "
                        "adjacent repeats for the given locus id, and instead just add it to the output catalog "
                        "without modification. This option can be specified more than once.")
    parser.add_argument("--output-dir",
                        help="Output directory. If not specified, the input catalog's directory will be used.")
    parser.add_argument("-o", "--output-catalog",
                        help="Path where to write the output catalog. If not specified, it will be based on the input "
                             "catalog path")
    parser.add_argument("input_catalog_paths",
                        nargs="+",
                        help="ExpansionHunter catalog json file path. If more than one catalog path is specified, the "
                             "catalogs will be processed sequentially, and the output catalog arg will be used as "
                             "an output filename suffix.")
    args = parser.parse_args()

    # validate command-line args
    for path in [args.ref_fasta] + args.input_catalog_paths:
        if not file_exists(path):
            parser.error(f"{path} not found")

    ref_fasta = pysam.FastaFile(os.path.abspath(os.path.expanduser(args.ref_fasta)))
    ref_fasta_chromosomes = set(ref_fasta.references)

    # parse the input catalog(s)
    input_catalogs = collections.defaultdict(dict)
    for input_catalog_path in args.input_catalog_paths:
        print(f"Loading input catalog {input_catalog_path}")
        with open_file(input_catalog_path) as f:
            if args.n:
                input_catalog = []
                for record in ijson.items(f, "item"):
                    if len(input_catalog) > args.n:
                        break
                    input_catalog.append(record)
            else:
                input_catalog = json.load(f)

        input_catalogs[input_catalog_path]["original_locus_id_sort_order"] = [
            record["LocusId"] for record in input_catalog]

        # sort the input catalog by its ReferenceRegion chromosome
        input_catalog.sort(key=lambda record: (
            record["ReferenceRegion"][0].split(":")[0]
            if isinstance(record["ReferenceRegion"], list) else
            record["ReferenceRegion"].split(":")[0]))

        input_catalogs[input_catalog_path]["input_catalog_records"] = input_catalog
        print(f"Loaded {len(input_catalog):,d} records from {input_catalog_path}")

    print(f"Processing {len(input_catalogs)} input variant catalog" + ("s" if len(input_catalogs) > 1 else f": {args.input_catalog_paths[0]}"))

    # validate input catalog records and get the min and max coordinates of the TR loci for each chromosome
    min_max_coords_for_chrom = collections.defaultdict(lambda: (10**10, 0))
    for input_catalog in input_catalogs.values():
        for record_i, record in enumerate(input_catalog["input_catalog_records"]):
            missing_keys = {"ReferenceRegion", "LocusStructure", "LocusId"} - set(record.keys())
            if missing_keys:
                print(f"WARNING: record {record} is missing keys: {missing_keys}. Skipping...")
                continue

            if isinstance(record["ReferenceRegion"], list):
                chrom = record["ReferenceRegion"][0].split(":")[0]
            else:
                chrom = record["ReferenceRegion"].split(":")[0]

            if chrom not in ref_fasta_chromosomes:
                chrom = toggle_chromosome_prefix(record)
                if chrom not in ref_fasta_chromosomes:
                    raise ValueError(f"chromosome '{chrom}' not found in reference fasta {args.ref_fasta} which has "
                                     f"chromosomes: " + ", ".join(sorted(ref_fasta_chromosomes)))

            if isinstance(record["ReferenceRegion"], list):
                chrom_start_end_tuples = [parse_interval(ref_region) for ref_region in record["ReferenceRegion"]]
            else:
                chrom_start_end_tuples = [parse_interval(record["ReferenceRegion"])]

            for current_chrom, current_start, current_end in chrom_start_end_tuples:
                current_start -= args.max_total_adjacent_region_size
                current_start = max(0, current_start)
                current_end += args.max_total_adjacent_region_size
                min_max_coords_for_chrom[current_chrom] = (
                    min(current_start, min_max_coords_for_chrom[current_chrom][0]),
                    max(current_end, min_max_coords_for_chrom[current_chrom][1]))

    # parse args.source_of_adjacent_loci within the genomic regions spanned by the input catalog(s)
    interval_trees_by_chrom = collections.defaultdict(intervaltree.IntervalTree)
    for chrom, (min_start, max_end) in min_max_coords_for_chrom.items():
        interval_trees_by_chrom[chrom] = get_interval_tree_for_chrom(
            args.source_of_adjacent_loci, chrom, min_start, max_end)

    # process the records from the input catalog(s) to add adjacent loci
    counters = collections.defaultdict(int)
    locus_ids_with_existing_adjacent_repeats = []
    locus_ids_with_new_adjacent_repeats = []
    for input_catalog_path, input_catalog in input_catalogs.items():
        output_catalog = []
        already_added_to_output_catalog = set()  # avoid duplicate records in the output catalog
        for record_i, record in enumerate(input_catalog["input_catalog_records"]):
            if args.n and record_i >= args.n:
                break
            counters["total records"] += 1

            if args.dont_add_adjacent_loci_to_this_locus_id and record["LocusId"] in args.dont_add_adjacent_loci_to_this_locus_id:
                print(f"Adding {record['LocusId']} without modification as requested by -x arg")
                output_catalog.append(record)
                counters["excluded locus id"] += 1
                continue

            # handle input catalog records that already have adjacent loci
            if isinstance(record["ReferenceRegion"], list):
                print(f"WARNING: record {record} already has adjacent intervals specified in the input variant catalog. "
                      f"Will add it to the output variant catalog without modification...")
                output_catalog.append(record)
                counters["already had adjacent loci"] += 1
                locus_ids_with_existing_adjacent_repeats.append(record["LocusId"])
                continue

            # add any adjacent loci to the current record
            chrom = record["ReferenceRegion"].split(":")[0]
            output_record = process_input_record(
                record,
                ref_fasta,
                interval_trees_by_chrom[chrom],
                max_distance_between_adjacent_repeats=args.max_distance_between_adjacent_repeats,
                max_total_adjacent_region_size=args.max_total_adjacent_region_size,
                max_overlap_between_adjacent_repeats=args.max_overlap_between_adjacent_repeats,
                max_adjacent_repeats=args.max_adjacent_repeats,
                add_extra_info_to_locus_id=args.add_extra_info_to_locus_id,
                add_extra_fields=args.add_extra_fields,
            )

            # check if adjacent loci were added
            if isinstance(output_record["ReferenceRegion"], list):
                locus_ids_with_new_adjacent_repeats.append(output_record["LocusId"])
                counters["added adjacent loci"] += 1
                counters["total reference regions"] += len(output_record["ReferenceRegion"])
                counters["total spacers"] += len(output_record["ReferenceRegion"]) - 1
                output_record_key = output_record["LocusStructure"] + ";" + ";".join(output_record["ReferenceRegion"])
            else:
                output_record_key = output_record["LocusStructure"] + ";" + output_record["ReferenceRegion"]
                counters["total reference regions"] += 1

            # if this is not a duplicate record, add it to the output catalog
            if output_record_key not in already_added_to_output_catalog:
                already_added_to_output_catalog.add(output_record_key)
                output_catalog.append(output_record)

            # compute stats about distances:
            if "DistancesBetweenReferenceRegions" in output_record:
                for distance in output_record["DistancesBetweenReferenceRegions"]:
                    counters[f"  {distance:3d}bp spacer"] += 1
            if "ReferenceRegionMotifs" in output_record:
                num_reference_regions = len(output_record["ReferenceRegion"])
                num_unique_motifs = len(set(compute_canonical_motif(m) for m in output_record["ReferenceRegionMotifs"]))
                same_motif = "same motif" if num_unique_motifs == 1 else "different motifs"
                counters[f"adjacent repeats count: {num_reference_regions:3d} reference regions with {same_motif}"] += 1

        # write the results to the output catalog
        if args.output_catalog:
            if len(input_catalogs) > 1:
                output_catalog_path = re.sub(".json$", "", input_catalog_path) + f".{args.output_catalog}.json"
            else:
                output_catalog_path = args.output_catalog
        else:
            output_catalog_path = re.sub(".json(.gz)?$", "", input_catalog_path) + ".with_adjacent_loci.json"
            if input_catalog_path.endswith("gz"):
                output_catalog_path += ".gz"

        if args.output_dir:
            output_catalog_path = os.path.join(args.output_dir, os.path.basename(output_catalog_path))

        # restore original sort order to match  loci in the input catalog
        output_catalog.sort(key=lambda record:
            input_catalog["original_locus_id_sort_order"].index(record["LocusId"]))

        fopen = gzip.open if output_catalog_path.endswith("gz") else open
        with fopen(output_catalog_path, "wt") as f:
            json.dump(output_catalog, f, indent=4)

        if counters["already had adjacent loci"]:
            print(f"Added {counters['already had adjacent loci']:,d} out of {counters['total records']:,d} records "
                  f"to the output catalog without modification since they already had adjacent repeats provided in the "
                  f"input catalog: ",
                  ", ".join(locus_ids_with_existing_adjacent_repeats[:10]),
                  "" if len(locus_ids_with_existing_adjacent_repeats) <= 10 else "...")
        else:
            print("None of the records in the input catalog(s) had adjacent loci specified")

        print(f"Added adjacent loci to {counters['added adjacent loci']:,d} out of {counters['total records']:,d} records ("
              f"{counters['added adjacent loci'] / counters['total records']:.1%}): ",
              ", ".join(locus_ids_with_new_adjacent_repeats[:10]),
              "" if len(locus_ids_with_new_adjacent_repeats) <= 10 else "...")
        if args.dont_add_adjacent_loci_to_this_locus_id:
            print(f"Added {counters['excluded locus id']:,d} out of {counters['total records']:,d} records ("
                  f"{counters['excluded locus id'] / counters['total records']:.1%}) without modification as "
                  f"requested by the -x arg")

        for key, value in sorted(counters.items(), key=lambda x: x[0]):
            if key.startswith("total") or key in {"already had adjacent loci", "added adjacent loci"}:
                continue
            if "adjacent repeats count" in key:
                label = "reference regions"
            elif "spacer" in key:
                label = "spacers"
            else:
                continue

            total = counters[f"total {label}"]
            print(f"{value:8,d} out of {total:8,d} {label} ({value / total:6.1%}) had {key}")

        print(f"Wrote {len(output_catalog):,d} total records to {output_catalog_path}")



if __name__ == '__main__':
    main()
