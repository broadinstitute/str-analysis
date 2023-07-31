import argparse
import collections
import ijson
import intervaltree
import json
import os
import re
import tqdm

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.file_utils import open_file, file_exists
from str_analysis.utils.misc_utils import parse_interval

REQUIRED_OUTPUT_FIELDS = {"LocusId", "ReferenceRegion", "LocusStructure", "VariantType"}

def parse_args():
    parser = argparse.ArgumentParser(description="Combines two or more repeat catalogs into a single catalog. Loci "
        "that overlap by more than ")
    parser.add_argument("-f", "--overlap-fraction", default=0.66, type=float, help="The minimum overlap for two loci "
        "to be considered as the same locus (assuming they are specified as having the same normalized motif). "
        "This is similar to the -f argument for 'bedtools intersect'.")
    parser.add_argument("--add-source-field", action="store_true", help="If specified, then the source field will be "
        "added to the output catalog. This requires the output format to be set to JSON.")
    parser.add_argument("--add-extra-fields-from-input-catalogs", action="store_true", help="If specified, then "
        "extra fields from the input catalogs will be added to the output catalog. This requires the output format "
        "to be set to JSON.")
    parser.add_argument("--output-format", choices=("JSON", "BED"), default="BED", help="Output file format.")
    parser.add_argument("--output-path", help="Output file path")
    parser.add_argument("--verbose", action="store_true", help="If specified, then print more stats")
    parser.add_argument("variant_catalog_json_or_bed", nargs="+", help="Paths of two or more repeat catalogs "
        "in JSON or BED format. For BED format, the chrom, start, and end should represent the repeat "
        "interval in 0-based coordinates, and the name field (column #4) should be the repeat unit. The order in which "
        "the catalogs are specified is important: All loci from the 1st catalog will be included in the output catalog."
        " This script will then parse subsequent catalogs sequentially, adding only loci that don't substantially "
        "overlap with loci in the previous catalog(s).")
    args = parser.parse_args()

    if len(args.variant_catalog_json_or_bed) < 2:
        parser.error("At least 2 input catalogs must be specified")

    if args.output_format == "BED":
        if args.add_source_field:
            parser.error("The --add-source-field option requires --output-format JSON to be specified")
        if args.add_extra_fields_from_input_catalogs:
            parser.error("The --add-extra-fields-from-input-catalogs option requires --output-format JSON to be specified")

    paths = parse_variant_catalog_paths_arg(args.variant_catalog_json_or_bed, parser)

    return args, paths


def parse_variant_catalog_paths_arg(variant_catalog_json_or_bed, argparser):
    paths = []
    for path in variant_catalog_json_or_bed:
        if not file_exists(path):
            argparser.error(f"{path} not found")

        if any(path.endswith(suffix) for suffix in [".bed", ".bed.gz", ".bed.bgz"]):
            paths.append((path, "BED"))
        elif any(path.endswith(suffix) for suffix in [".json", ".json.gz"]):
            paths.append((path, "JSON"))
        else:
            argparser.error(f"Unrecognized file extension: {path}")

    return paths


def parse_motifs_from_locus_structure(locus_structure):
    """Takes an ExpansionHunter LocusStructure like "(CAG)*AGAC(GCC)+" and returns a list of motifs ["CAG", "GCC"]"""
    return [
        locus_structure_part.split("(")[-1] for locus_structure_part in locus_structure.split(")")[:-1]
    ]

def get_variant_catalog_iterator(
        variant_catalog_json_or_bed,
        file_type,
        add_extra_fields_from_input_catalogs=False,
        verbose=False):
    """Takes the path of a JSON or BED file and returns an iterator over variant catalog records parsed from that file.

    Args:
        intervaltrees (dict): a dictionary that maps chromosome names to IntervalTrees
        variant_catalog_json_or_bed (str): path to a JSON or BED file containing variant catalog records
        file_type (str): either "JSON" or "BED"
        add_extra_fields_from_input_catalogs (bool): If False, then only the required fields will be kept from the
            input variant catalog and any extra fields will be discarded.
        verbose (bool): If True, add a progress bar
    """
    if verbose:
        print(f"Parsing {variant_catalog_json_or_bed}")

    if file_type == "JSON":
        with open_file(variant_catalog_json_or_bed, "rt") as f:
            file_iterator = ijson.items(f, "item")
            if verbose:
                file_iterator = tqdm.tqdm(file_iterator, unit=" records", unit_scale=True)

            for i, record in enumerate(file_iterator):
                missing_keys = REQUIRED_OUTPUT_FIELDS - record.keys()
                if missing_keys:
                    raise ValueError(f"Record #{i+1} in {variant_catalog_json_or_bed} is missing required field(s): "
                                     f"{missing_keys}")

                if not add_extra_fields_from_input_catalogs:
                    record = {k: record[k] for k in REQUIRED_OUTPUT_FIELDS}

                if isinstance(record["ReferenceRegion"], list):
                    reference_region_count = len(record["ReferenceRegion"])
                    start_0based, end_1based = None, None
                    for reference_region in record["ReferenceRegion"]:
                        unmodified_chrom, current_start_0based, current_end_1based = parse_interval(reference_region)
                        if start_0based is None or current_start_0based < start_0based:
                            start_0based = current_start_0based
                        if end_1based is None or current_end_1based > end_1based:
                            end_1based = current_end_1based
                else:
                    reference_region_count = 1
                    unmodified_chrom, start_0based, end_1based = parse_interval(record["ReferenceRegion"])

                motifs = parse_motifs_from_locus_structure(record["LocusStructure"])
                if len(motifs) != reference_region_count:
                    print(f"ERROR: {variant_catalog_json_or_bed} record {record_i+1:,d}: locus structure "
                          f"{record['LocusStructure']} contains a different number of motifs ({len(motifs)}) than the "
                          f"number of reference regions ({reference_region_count}): {record['ReferenceRegion']}. "
                          f"Skipping...")
                    continue

                yield unmodified_chrom, start_0based, end_1based, record

    elif file_type == "BED":
        with open_file(variant_catalog_json_or_bed, "rt") as file_iterator:
            if verbose:
                file_iterator = tqdm.tqdm(file_iterator, unit=" records", unit_scale=True)

            for line_i, line in enumerate(file_iterator):
                fields = line.strip().split("\t")
                unmodified_chrom = fields[0]
                chrom = unmodified_chrom.replace("chr", "")
                start_0based = int(fields[1])
                end_1based = int(fields[2])
                if "(" in fields[3] or ")" in fields[3]:
                    motifs = parse_motifs_from_locus_structure(fields[3])
                    if len(motifs) != 1:
                        filename = os.path.basename(variant_catalog_json_or_bed)
                        print(f"ERROR: {filename} line #{line_i+1:,d}: {chrom}:{start_0based }-{end_1based} "
                              f"locus structure {fields[3]} contains more than one motif. Skipping...")
                        continue

                    motif = motifs[0]
                else:
                    motif = fields[3]


                record = {
                    "LocusId": f"{chrom}-{start_0based + 1}-{end_1based}-{motif}",
                    "ReferenceRegion": f"{unmodified_chrom}:{start_0based + 1}-{end_1based}",
                    "LocusStructure": f"({motif})*",
                    "VariantType": "Repeat",
                }

                yield unmodified_chrom, start_0based, end_1based, record
    else:
        raise ValueError(f"Unrecognized file type: {file_type} for {variant_catalog_json_or_bed}")


def add_variant_catalog_to_interval_trees(
        interval_trees,
        variant_catalog_json_or_bed,
        file_type,
        min_overlap_fraction=0.01,
        add_source_field=False,
        add_extra_fields_from_input_catalogs=False,
        verbose=False):
    """Parses the the given input variant catalog and adds any new unique records to the IntervalTrees.

    Args:
        interval_trees (dict): a dictionary that maps chromosome names to IntervalTree objects for overlap detection
        variant_catalog_json_or_bed (str): path to a JSON or BED file containing input variant catalog records
        file_type (str): either "JSON" or "BED"
        add_source_field (bool): If True, then the source file path will be added to each record as a new "Source" field
        add_extra_fields_from_input_catalogs (bool): If False, then only the required fields will be kept in each input
            variant catalog record and any extra fields will be discarded.
        verbose (bool): If True, then print more stats about the number of records added to the output catalog
    """
    if verbose:
        print("- "*60)

    variant_catalog_filename = os.path.basename(variant_catalog_json_or_bed)
    counters = collections.defaultdict(int)
    for chrom, start_0based, end_1based, new_record in get_variant_catalog_iterator(
            variant_catalog_json_or_bed, file_type,
            add_extra_fields_from_input_catalogs=add_extra_fields_from_input_catalogs,
            verbose=verbose):

        chrom = chrom.replace("chr", "")
        # check for overlap with existing loci
        counters["total"] += 1
        new_record_canonical_motif = None
        found_matching_existing_record = False
        for overlapping_interval in interval_trees[chrom].overlap(start_0based, end_1based):
            overlap_size = overlapping_interval.overlap_size(start_0based, end_1based)

            # if the new record overlaps an existing record by less than the minimum overlap fraction, then it's not
            # considered a duplicate of the existing record
            if overlap_size < min_overlap_fraction * (end_1based - start_0based) \
                    and overlap_size < min_overlap_fraction * (overlapping_interval.end - overlapping_interval.begin):
                continue

            existing_record = overlapping_interval.data
            if new_record["LocusStructure"] == existing_record["LocusStructure"]:
                counters[f"overlapped an exising locus by at least {100*min_overlap_fraction}% " \
                         f"and had the exact same LocusStructure"] += 1
                found_matching_existing_record = True
                break
            elif isinstance(new_record["ReferenceRegion"], list):
                # since loci with multiple adjacent repeats may have multiple motifs, duplicate detection for these loci
                # is based on exact equality between their LocusStructures. The previous if statement ruled this out.
                continue
            else:
                # check if the existing record's canonical motif is the same as the canonical motif of the new record
                if new_record_canonical_motif is None:
                    motif = new_record["LocusStructure"].strip("()*+")
                    try:
                        new_record_canonical_motif = compute_canonical_motif(motif)
                    except Exception as e:
                        raise ValueError(f"Error computing canonical motif for {new_record}: {e}")
                existing_record_motif = new_record["LocusStructure"].strip("()*+")
                existing_record_canonical_motif = compute_canonical_motif(existing_record_motif)
                if new_record_canonical_motif == existing_record_canonical_motif:
                    counters[f"overlapped an existing locus by at least {100*min_overlap_fraction}% " \
                             f"and had the same canonical motif"] += 1
                    found_matching_existing_record = True
                    break

        if found_matching_existing_record:
            # don't add the current record to the output catalog
            continue

        if add_source_field:
            new_record["Source"] = variant_catalog_filename

        counters["added"] += 1
        interval_trees[chrom].add(intervaltree.Interval(start_0based, end_1based, data=new_record))

    print(f"Added {counters['added']:,d} out of {counters['total']:,d} ({counters['added']/counters['total']:6.1%}) "
          f"records from {variant_catalog_filename} to the output catalog")
    if verbose:
        for k, v in sorted(counters.items()):
            if k not in {"added", "total"}:
                print(" "*3, f"Discarded {counters[k]:7,d} out of {counters['total']:7,d} "
                      f"({counters[k]/counters['total']:6.1%}) duplicate records since they {k}")

def convert_interval_trees_to_output_records(interval_trees):
    """Converts the IntervalTrees to a generator for output catalog records."""
    counter = 0
    motif_size_stats = collections.defaultdict(int)

    for interval_tree in interval_trees.values():
        for interval in sorted(interval_tree, key=lambda i: (i.begin, i.end)):
            counter += 1
            yield interval.data


def write_output_catalog(interval_trees, output_path, output_format):
    """Writes the output catalog to a file.

    Args:
        interval_trees (dict): a dictionary that maps chromosome names to IntervalTree objects for overlap detection
        output_path (str): path of the output file
        output_format (str): either "JSON" or "BED"
    """

    output_catalog_record_generator = convert_interval_trees_to_output_records(interval_trees)

    # write the output catalog to the output file in the requested format
    if output_format == "JSON":
        open_file = gzip.open if output_path.endswith("gz") else open
        with open_file(output_path, "wt") as output_catalog:
            output_records_list = list(output_catalog_record_generator)
            json.dump(output_records_list, output_catalog, indent=1)
        print(f"Done writing {len(output_records_list):,d} output records to {output_path}")

    elif output_format == "BED":
        output_path = re.sub(".b?gz$", "", output_path)

        with open(output_path, "wt") as output_catalog:
            total = 0
            for record in output_catalog_record_generator:
                locus_structure = record["LocusStructure"]
                reference_regions = record["ReferenceRegion"]
                if not isinstance(reference_regions, list):
                    reference_regions = [reference_regions]

                motifs = parse_motifs_from_locus_structure(locus_structure)
                if len(motifs) != len(reference_regions):
                    raise ValueError(f"Number of motifs ({len(motifs)}) != number of reference regions "
                                     f"({len(reference_regions)}) in record: {record}")

                for reference_region, motif in zip(reference_regions, motifs):
                    chrom, start_0based, end_1based = parse_interval(reference_region)
                    total += 1
                    output_catalog.write("\t".join([
                        chrom,
                        str(start_0based),
                        str(end_1based),
                        motif,
                        ".",
                    ]) + "\n")

        os.system(f"bgzip -f {output_path}")
        os.system(f"tabix -f -p bed {output_path}.gz")
        print(f"Done writing {total:,d} output records to {output_path}.gz")



def print_catalog_stats(interval_trees, has_source_field=False):
    total = 0
    motif_counters = collections.defaultdict(int)
    source_counters = collections.defaultdict(int)
    for interval_tree in interval_trees.values():
        for inerval in interval_tree:
            total += 1
            motifs = parse_motifs_from_locus_structure(inerval.data["LocusStructure"])
            for motif in motifs:
                motif_label = f"{len(motif)}bp" if len(motif) <= 6 else f"7+bp"
                motif_counters[motif_label] += 1
            if has_source_field:
                source_counters[inerval.data["Source"]] += 1

    if has_source_field:
        print("Source of loci in output catalog:")
        for label, count in sorted(source_counters.items(), key=lambda x: x[1], reverse=True):
            print(f"   {count:9,d} out of {total:9,d} ({count/total:5.1%}) {label}")

    print("Motif sizes:")
    for label, count in sorted(motif_counters.items(), key=lambda x: x[1], reverse=True):
        print(f"   {count:9,d} out of {total:9,d} ({count/total:5.1%}) {label}")

def main():
    args, paths = parse_args()

    # parse each catalog and add each new unique record to this IntervalTrees dictionary
    interval_trees = collections.defaultdict(intervaltree.IntervalTree)
    for path, file_type in paths:
        add_variant_catalog_to_interval_trees(
            interval_trees, path, file_type,
            min_overlap_fraction=args.overlap_fraction,
            add_source_field=args.add_source_field,
            add_extra_fields_from_input_catalogs=args.add_extra_fields_from_input_catalogs,
            verbose=args.verbose)

    # write the output catalog to a file
    if not args.output_path:
        args.output_path = f"combined.{len(paths)}_catalogs.{args.output_format.lower()}"

    if args.verbose:
        print(f"Writing combined catalog to {args.output_path}")
    write_output_catalog(interval_trees, args.output_path, args.output_format)

    if args.verbose:
        print_catalog_stats(interval_trees, has_source_field=args.add_source_field)

if __name__ == "__main__":
    main()
