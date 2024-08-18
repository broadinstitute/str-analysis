import argparse
import collections
import ijson
import intervaltree
import gzip
import simplejson as json
import os
import re
import tqdm

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.eh_catalog_utils import parse_motifs_from_locus_structure, convert_json_records_to_bed_format_tuples
from str_analysis.utils.file_utils import open_file, file_exists
from str_analysis.utils.misc_utils import parse_interval

REQUIRED_OUTPUT_FIELDS = {"LocusId", "ReferenceRegion", "LocusStructure", "VariantType"}
SEPARATOR_FOR_MULTIPLE_SOURCES = " ||| "

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
    parser.add_argument("--overlapping-loci-action",
        choices=("keep-first", "keep-last", "keep-both", "keep-narrow", "keep-wider", "merge"),
        default="keep-first", help="When two loci overlap and have the same canonical motif, this option specifies "
        "which locus to keep in the output catalog. 'keep-first' keeps the first locus, 'keep-last' keeps the last "
        "locus, 'keep-both' keeps both loci, 'keep-narrow' keeps the locus with the narrower interval, 'keep-wider' "
        "keeps the locus with the wider interval, and 'merge' merges the two loci into a single locus.")
    parser.add_argument("-m", "--merge-adjacent-loci-with-same-motif", action="store_true",
        help="A repeat catalog generated by using a tool like TandemRepeatFinder to find all pure repeats in a "
        "reference genome provides a good starting point for a repeat catalog. However, the presence of a single base "
        "pair interruption in a reference repeat sequence will cause that repeat locus to be split into 2 adjacent "
        "loci with the same motif (under cyclic shift). This option enables detection of these adjacent loci (separated "
        "by 1 base pair and having the same motif) and merges them into a single locus in the output catalog.")
    parser.add_argument("--output-format", choices=("JSON", "BED"), help="Output file format. If not specified, both "
                                                                         "a JSON and a BED file will be generated.")
    parser.add_argument("--output-prefix", help="Output filename prefix")
    parser.add_argument("--verbose", action="store_true", help="If specified, then print more stats")
    parser.add_argument("--verbose-overlaps", action="store_true", help="If specified, print out overlapping definitions"
                        "that have similar motifs but different boundaries")
    parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    parser.add_argument("--write-bed-files-with-new-loci", action="store_true", help="If specified, then for every "
                        "input catalog except the first one, this script will output a BED file that contains the new "
                        "loci introduced by that catalog. This is useful for troubleshooting catalogs and "
                        "understanding the differences between them.")
    parser.add_argument("variant_catalog_json_or_bed", nargs="+", help="Paths of two or more repeat catalogs "
        "in JSON or BED format. For BED format, the chrom, start, and end should represent the repeat "
        "interval in 0-based coordinates, and the name field (column #4) should be the repeat unit. The order in which "
        "the catalogs are specified is important: All loci from the 1st catalog will be included in the output catalog."
        " This script will then parse subsequent catalogs sequentially, adding only loci that don't substantially "
        "overlap with loci in the previous catalog(s).")
    args = parser.parse_args()

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


def get_variant_catalog_iterator(
        variant_catalog_json_or_bed,
        file_type,
        add_extra_fields_from_input_catalogs=False,
        verbose=False,
        show_progress_bar=False):
    """Takes the path of a JSON or BED file and returns an iterator over variant catalog records parsed from that file.

    Args:
        intervaltrees (dict): a dictionary that maps chromosome names to IntervalTrees
        variant_catalog_json_or_bed (str): path to a JSON or BED file containing variant catalog records
        file_type (str): either "JSON" or "BED"
        add_extra_fields_from_input_catalogs (bool): If False, then only the required fields will be kept from the
            input variant catalog and any extra fields will be discarded.
        verbose (bool): If True, add a progress bar
        show_progress_bar (bool): If True, then show a progress bar
    """
    if verbose:
        print(f"Parsing {variant_catalog_json_or_bed}")

    if file_type == "JSON":
        with open_file(variant_catalog_json_or_bed, is_text_file=True) as f:
            file_iterator = ijson.items(f, "item")
            if show_progress_bar:
                file_iterator = tqdm.tqdm(file_iterator, unit=" records", unit_scale=True)

            for record_i, record in enumerate(file_iterator):
                missing_keys = REQUIRED_OUTPUT_FIELDS - record.keys()
                if missing_keys:
                    raise ValueError(f"Record #{record_i+1} in {variant_catalog_json_or_bed} is missing required "
                                     f"field(s): {missing_keys}")

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
        with open_file(variant_catalog_json_or_bed, is_text_file=True) as file_iterator:
            if show_progress_bar:
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
                        print(f"WARNING: {filename} line #{line_i+1:,d}: {chrom}:{start_0based }-{end_1based} "
                              f"locus structure {fields[3]} contains more than one motif. Skipping...")
                        continue

                    motif = motifs[0]
                else:
                    motif = fields[3]

                record = {
                    "LocusId": f"{chrom}-{start_0based}-{end_1based}-{motif}",
                    "ReferenceRegion": f"{unmodified_chrom}:{start_0based}-{end_1based}",
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
        overlapping_loci_action="keep-first",
        min_overlap_fraction=0.01,
        add_source_field=False,
        add_extra_fields_from_input_catalogs=False,
        verbose=False,
        verbose_overlaps=False,
        show_progress_bar=False,
        write_bed_files_with_new_loci=False
):
    """Parses the the given input variant catalog and adds any new unique records to the IntervalTrees.

    Args:
        interval_trees (dict): a dictionary that maps chromosome names to IntervalTree objects for overlap detection
        variant_catalog_json_or_bed (str): path to a JSON or BED file containing input variant catalog records
        file_type (str): either "JSON" or "BED"
        add_source_field (bool): If True, then the source file path will be added to each record as a new "Source" field
        add_extra_fields_from_input_catalogs (bool): If False, then only the required fields will be kept in each input
            variant catalog record and any extra fields will be discarded.
        verbose (bool): If True, then print more stats about the number of records added to the output catalog
        verbose_overlaps (bool): If True, print info about loci that overlap and have similar motifs, but different bondaries.
        show_progress_bar (bool): Whether to show a progress bar
        write_bed_files_with_new_loci (bool): If True, then write a BED file of loci not seen in previous catalogs.
    """
    if verbose:
        print("- "*60)

    variant_catalog_filename = os.path.basename(variant_catalog_json_or_bed)
    counters = collections.defaultdict(int)
    unique_loci = set()
    for unmodified_chrom, start_0based, end_1based, new_record in get_variant_catalog_iterator(
            variant_catalog_json_or_bed, file_type,
            add_extra_fields_from_input_catalogs=add_extra_fields_from_input_catalogs,
            verbose=verbose, show_progress_bar=show_progress_bar):

        chrom = unmodified_chrom.replace("chr", "")
        # check for overlap with existing loci
        counters["total"] += 1
        new_record_canonical_motif = None
        found_existing_record_with_similar_motif = False
        for overlapping_interval in interval_trees[chrom].overlap(start_0based, end_1based):
            overlap_size = overlapping_interval.overlap_size(start_0based, end_1based)

            # if the new record overlaps an existing record by less than the minimum overlap fraction, then it's not
            # considered a duplicate of the existing record
            if overlap_size < min_overlap_fraction * (end_1based - start_0based) \
                    and overlap_size < min_overlap_fraction * (overlapping_interval.end - overlapping_interval.begin):
                continue

            # handle overlapping interval
            existing_record = overlapping_interval.data
            if new_record["LocusStructure"] == existing_record["LocusStructure"]:
                counters[f"overlapped an exising locus by at least {100*min_overlap_fraction}% " \
                         f"and had the exact same LocusStructure"] += 1
                found_existing_record_with_similar_motif = True
                break
            elif not isinstance(new_record["ReferenceRegion"], list) and not isinstance(existing_record["ReferenceRegion"], list):
                # check if the existing record's canonical motif is the same as the canonical motif of the new record
                existing_record_motifs = parse_motifs_from_locus_structure(existing_record["LocusStructure"])
                if len(existing_record_motifs) != 1:
                    raise ValueError(f"Unexpected LocusStructure in {existing_record}.")

                try:
                    existing_record_canonical_motif = compute_canonical_motif(existing_record_motifs[0], include_reverse_complement=False)
                except Exception as e:
                    raise ValueError(f"Error computing canonical motif for {existing_record}: {e}")

                if new_record_canonical_motif is None:
                    # compute the new record's canonical motif if it hasn't been computed already
                    new_record_motifs = parse_motifs_from_locus_structure(new_record["LocusStructure"])
                    if len(new_record_motifs) != 1:
                        raise ValueError(f"Unexpected LocusStructure in {new_record}.")

                    try:
                        new_record_canonical_motif = compute_canonical_motif(new_record_motifs[0], include_reverse_complement=False)
                    except Exception as e:
                        raise ValueError(f"Error computing canonical motif for {new_record}: {e}")

                if existing_record_canonical_motif == new_record_canonical_motif:
                    counters[f"overlapped an existing locus by at least {100*min_overlap_fraction}% " \
                             f"and had the same canonical motif"] += 1
                    found_existing_record_with_similar_motif = True
                    break

                # check if one of the motifs contains the other
                if len(existing_record_canonical_motif) != len(new_record_canonical_motif):
                    if len(existing_record_canonical_motif) < len(new_record_canonical_motif):
                        short_motif, long_motif = existing_record_canonical_motif, new_record_canonical_motif
                    elif len(existing_record_canonical_motif) > len(new_record_canonical_motif):
                        short_motif, long_motif = new_record_canonical_motif, existing_record_canonical_motif

                    expanded_motif = short_motif * (1 + len(long_motif)//len(short_motif))
                    if expanded_motif[:len(long_motif)] == long_motif:
                        counters[f"overlapped an existing locus by at least {100*min_overlap_fraction}% " \
                                 f"and one motif was contained within the other"] += 1
                        found_existing_record_with_similar_motif = True
                        break

        if add_source_field:
            new_record["Source"] = variant_catalog_filename

        if found_existing_record_with_similar_motif:
            existing_record = overlapping_interval.data

            discard_new = False
            remove_existing = False
            if overlapping_loci_action == "keep-first":
                discard_new = True
            elif overlapping_loci_action == "keep-last":
                interval_trees[chrom].remove(overlapping_interval)
            elif overlapping_loci_action == "keep-both":
                pass
            elif overlapping_loci_action == "keep-narrow":
                if (end_1based - start_0based) >= (overlapping_interval.end - overlapping_interval.begin):
                    discard_new = True
                else:
                    remove_existing = True
            elif overlapping_loci_action == "keep-wider":
                if (end_1based - start_0based) <= (overlapping_interval.end - overlapping_interval.begin):
                    discard_new = True
                else:
                    remove_existing = True
            elif overlapping_loci_action == "merge":
                remove_existing = True
                min_start_0based = min(start_0based, overlapping_interval.begin)
                max_end_1based = max(end_1based, overlapping_interval.end)
                existing_record_locus_structure = existing_record["LocusStructure"]
                existing_record_motifs = parse_motifs_from_locus_structure(existing_record_locus_structure)
                if len(current_motifs) != 1:
                    raise ValueError(f"Unexpected LocusStructure in {existing_record}.")
                existing_record_motif  = existing_record_motifs[0]
                new_record = {
                    "LocusId": f"{chrom}-{min_start_0based}-{max_end_1based}-{existing_record_motif}",
                    "ReferenceRegion": f"{unmodified_chrom}:{min_start_0based}-{max_end_1based}",
                    "LocusStructure": existing_record_locus_structure,
                    "VariantType": existing_record["VariantType"],
                }
                if add_source_field:
                    new_record["Source"] = existing_record["Source"] + SEPARATOR_FOR_MULTIPLE_SOURCES + variant_catalog_filename
            else:
                raise ValueError(f"Unexpected overlapping_loci_action value: {overlapping_loci_action}")

            if verbose_overlaps:
                if new_record["ReferenceRegion"] != existing_record["ReferenceRegion"]:
                    print("="*100)
                    _, existing_start, existing_end = parse_interval(existing_record["ReferenceRegion"])
                    _, new_start, new_end = parse_interval(new_record["ReferenceRegion"])
                    existing_size = existing_end - existing_start
                    new_size = new_end - new_start
                    size_comparison = f"includes {existing_size - new_size} more repeats than" if existing_size > new_size \
                        else "is the same size as" if existing_size == new_size else \
                        f"includes {new_size - existing_size} fewer repeats than"
                    print(f"Existing record:", existing_record["ReferenceRegion"], existing_record["LocusStructure"], size_comparison, "the new record")
                    print(f"     New record:", new_record["ReferenceRegion"], new_record["LocusStructure"])
                    print(f"         Action:", "Replacing existing record with new record " if not discard_new and remove_existing else (
                        "Discarding new record" if discard_new else "Adding new record"
                    ))
            if remove_existing:
                interval_trees[chrom].remove(overlapping_interval)
            if discard_new:
               continue

        if write_bed_files_with_new_loci:
            motifs = parse_motifs_from_locus_structure(new_record["LocusStructure"])
            if len(motifs) > 1:
                print(f"Unexpected LocusStructure in {new_record}: {new_record['LocusStructure']}. Will not save to unique loci...")
            else:
                unique_loci.add((unmodified_chrom, start_0based, end_1based, motifs[0]))

        counters["added"] += 1
        interval_trees[chrom].add(intervaltree.Interval(start_0based, end_1based, data=new_record))

    print(f"Added {counters['added']:,d} out of {counters['total']:,d} ({counters['added']/(counters['total'] or 1):6.1%}) "
          f"records from {variant_catalog_filename} to the output catalog")
    if verbose:
        for k, v in sorted(counters.items()):
            if k not in {"added", "total"}:
                print(" "*3, f"Discarded {counters[k]:7,d} out of {counters['total']:7,d} "
                      f"({counters[k]/counters['total']:6.1%}) records since they {k}")

    if write_bed_files_with_new_loci:
        unqiue_loci_bed_prefix = re.sub("(.json|.bed)(.gz)?$", "", os.path.basename(variant_catalog_json_or_bed))
        unique_loci_bed_filename = f"{unqiue_loci_bed_prefix}.unique_loci.bed"
        with open(unique_loci_bed_filename, "wt") as unique_loci_bed:
            for chrom, start_0based, end_1based, motif in sorted(unique_loci):
                unique_loci_bed.write("\t".join(map(str, [
                    chrom,
                    start_0based,
                    end_1based,
                    motif,
                    ".",
                ])) + "\n")

        os.system(f"bgzip -f {unique_loci_bed_filename}")
        os.system(f"tabix -f {unique_loci_bed_filename}.gz")
        print(f"Wrote {counters['added']:,d} unique loci from {variant_catalog_filename} to {unique_loci_bed_filename}.gz")

def check_whether_to_merge_adjacent_loci(previous_interval, current_interval, add_source_field=False):
    """Checks whether the two loci should be merged into a single locus.

    Args:
        previous_interval (intervaltree.Interval): The previous locus
        current_interval (intervaltree.Interval): The current locus
        add_source_field (bool): see the --add-source-field option

    Return:
        bool: True if the two loci should be merged into a single locus
        intervaltree.Interval: The merged locus if the two loci should be merged, otherwise None
    """

    # check whether the repeats are within 1bp of each other
    if previous_interval.end + 1 != current_interval.begin:
        return False, None

    # check if their motif is the same.
    previous_motifs = parse_motifs_from_locus_structure(previous_interval.data["LocusStructure"])
    if len(previous_motifs) != 1:
        # can't merge loci with multiple motifs
        return False, None

    current_motifs = parse_motifs_from_locus_structure(current_interval.data["LocusStructure"])
    if len(current_motifs) != 1:
        # can't merge loci with multiple motifs
        return False, None

    previous_motif = previous_motifs[0]
    previous_motif_shifted = compute_canonical_motif(previous_motifs[0], include_reverse_complement=False)
    current_motif_shifted = compute_canonical_motif(current_motifs[0], include_reverse_complement=False)
    if previous_motif_shifted != current_motif_shifted:
        return False, None

    # the loci have the same motif and are 1bp appart, so create a merged Interval
    unmodified_chrom, _, _ = parse_interval(previous_interval.data["ReferenceRegion"])
    chrom = unmodified_chrom.replace("chr", "")
    start_0based = previous_interval.begin
    end_1based = current_interval.end
    locus_structure = previous_interval.data["LocusStructure"]
    merged_data = {
        "LocusId": f"{chrom}-{start_0based}-{end_1based}-{previous_motif}",
        "ReferenceRegion": f"{unmodified_chrom}:{start_0based}-{end_1based}",
        "LocusStructure": locus_structure,
        "VariantType": "Repeat",
    }

    if add_source_field:
        if not previous_interval.data.get("Source"):
            raise ValueError(f"'Source' field not found in record {previous_interval.data}")
        if not current_interval.data.get("Source"):
            raise ValueError(f"'Source' field not found in record {current_interval.data}")

        if previous_interval.data["Source"] == current_interval.data["Source"]:
            merged_data["Source"] = previous_interval.data["Source"]
        else:
            merged_data["Source"] = previous_interval.data["Source"] + SEPARATOR_FOR_MULTIPLE_SOURCES + current_interval.data["Source"]

    merged_interval = intervaltree.Interval(
        previous_interval.begin,
        current_interval.end,
        data=merged_data,
    )

    #print("Merged",
    #      previous_interval.data["ReferenceRegion"],  previous_interval.data["LocusStructure"], "and",
    #      current_interval.data["ReferenceRegion"], current_interval.data["LocusStructure"], "into",
    #      merged_interval.data["ReferenceRegion"], merged_interval.data["LocusStructure"])

    return True, merged_interval


def convert_interval_trees_to_output_records(
    interval_trees,
    merge_adjacent_loci_with_same_motif=False,
    add_source_field=False,
):
    """Converts the IntervalTrees to a generator for output catalog records.

    Args:
        interval_trees (dict): a dictionary that maps chromosome names to IntervalTree objects for overlap detection
        merge_adjacent_loci_with_same_motif (bool): see the --merge-adjacent-loci-with-same-motif option
        add_source_field (bool): see the --add-source-field option

    Yield:
        dict: A dictionary representing a record in the output catalog
    """
    counter = 0

    if not merge_adjacent_loci_with_same_motif:
        for interval_tree in interval_trees.values():
            for interval in sorted(interval_tree, key=lambda i: (i.begin, i.end)):
                counter += 1
                yield interval.data
    else:
        merged_counter = 0
        for interval_tree in interval_trees.values():
            previous_interval = None
            for interval in sorted(interval_tree, key=lambda i: (i.begin, i.end)):
                counter += 1
                if previous_interval is None:
                    previous_interval = interval
                    continue

                should_merge, merged_interval = check_whether_to_merge_adjacent_loci(
                    previous_interval, interval, add_source_field=add_source_field)
                if should_merge:
                    merged_counter += 1
                    previous_interval = merged_interval
                    continue

                yield previous_interval.data
                previous_interval = interval

            if previous_interval is not None:
                yield previous_interval.data

        print(f"Merged {merged_counter:,d} out of {counter:,d} ({merged_counter/counter:6.1%}) adjacent loci that "
              f"were within 1bp of each other and had the same motif after cyclic shift")


def write_output_catalog(output_catalog_record_iter, output_path, output_format):
    """Writes the output catalog to a file.

    Args:
        output_catalog_record_iter (iter): An iterator (in sorted genomic order) over records to include in the
            output catalog. Each record is a dictionary in the ExpansionHunter variant catalog format.
        output_path (str): path of the output file
        output_format (str): either "JSON" or "BED"
    """

    # write the output catalog to the output file in the requested format
    if output_format == "JSON":
        fopen = gzip.open if output_path.endswith("gz") else open
        with fopen(output_path, "wt") as output_catalog:
            output_records_list = list(output_catalog_record_iter)
            json.dump(output_records_list, output_catalog, indent=4, ignore_nan=True)
        print(f"Done writing {len(output_records_list):,d} output records to {output_path}")

    elif output_format == "BED":
        output_path = re.sub(".b?gz$", "", output_path)
        with open(output_path, "wt") as output_catalog:
            total = 0
            for bed_record in sorted(convert_json_records_to_bed_format_tuples(output_catalog_record_iter)):
                total += 1
                output_catalog.write("\t".join(map(str, bed_record)) + "\n")

        os.system(f"bgzip -f {output_path}")
        os.system(f"tabix -f -p bed {output_path}.gz")
        print(f"Done writing {total:,d} output records to {output_path}.gz")


def print_catalog_stats(interval_trees, has_source_field=False):
    total = 0
    motif_counters = collections.defaultdict(int)
    source_counters = collections.defaultdict(int)
    chrom_counters = collections.defaultdict(int)
    for chrom, interval_tree in interval_trees.items():
        for interval in interval_tree:
            chrom_type = chrom.upper() if any(c in chrom.upper() for c in "XYM") else "autosomes"
            chrom_counters[chrom_type] += 1
            
            total += 1
            motifs = parse_motifs_from_locus_structure(interval.data["LocusStructure"])
            for motif in motifs:
                motif_label = f"{len(motif)}bp" if len(motif) <= 6 else f"7+bp"
                motif_counters[motif_label] += 1
            if has_source_field:
                source_counters[interval.data["Source"]] += 1

    if has_source_field:
        print("Source of loci in output catalog:")
        for label, count in sorted(source_counters.items(), key=lambda x: x[1], reverse=True):
            print(f"   {count:9,d} out of {total:9,d} ({count/total:5.1%}) {label}")

    print("Motif sizes:")
    for label, count in sorted(motif_counters.items(), key=lambda x: x[1], reverse=True):
        print(f"   {count:9,d} out of {total:9,d} ({count/total:5.1%}) {label}")

    print("Chromsomes:")
    for label, count in sorted(chrom_counters.items(), key=lambda x: x[1], reverse=True):
        print(f"   {count:9,d} out of {total:9,d} ({count/total:5.1%}) are on {label}")

def fix_source_field_for_merged_adjacent_loci_with_multiple_sources(iterator):
    for record in iterator:
        if "Source" not in record:
            raise ValueError(f"Unexpected absence of 'Source' field in record: {record}")

        if SEPARATOR_FOR_MULTIPLE_SOURCES in record["Source"]:
            source = ", ".join(sorted(set(record["Source"].split(SEPARATOR_FOR_MULTIPLE_SOURCES))))
            record["Source"] = f"merged adjacent loci from: {source}"

        yield record

def main():
    args, paths = parse_args()

    # parse each catalog and add each new unique record to this IntervalTrees dictionary
    interval_trees = collections.defaultdict(intervaltree.IntervalTree)
    for i, (path, file_type) in enumerate(paths):
        add_variant_catalog_to_interval_trees(
            interval_trees, path, file_type,
            overlapping_loci_action="keep-first",
            min_overlap_fraction=args.overlap_fraction,
            add_source_field=args.add_source_field,
            add_extra_fields_from_input_catalogs=args.add_extra_fields_from_input_catalogs,
            verbose=args.verbose,
            verbose_overlaps=args.verbose_overlaps,
            show_progress_bar=args.show_progress_bar,
            write_bed_files_with_new_loci=args.write_bed_files_with_new_loci and i > 0,
        )

    # write the output catalog to a file
    output_formats = [args.output_format] if args.output_format else ["JSON", "BED"]
    if not args.output_prefix:
        args.output_prefix = f"combined.{len(paths)}_catalogs"

    for output_format in output_formats:
        output_catalog_record_generator = convert_interval_trees_to_output_records(
            interval_trees,
            merge_adjacent_loci_with_same_motif=args.merge_adjacent_loci_with_same_motif,
            add_source_field=args.add_source_field)

        if args.add_source_field and args.merge_adjacent_loci_with_same_motif:
            # convert the SEPARATOR_FOR_MULTIPLE_SOURCES to commas
            output_catalog_record_generator = fix_source_field_for_merged_adjacent_loci_with_multiple_sources(
                output_catalog_record_generator)

        output_path = f"{args.output_prefix}.{output_format.lower()}"
        if output_format == "JSON":
            output_path += ".gz"

        if args.verbose:
            print(f"Writing combined catalog to {output_path}")
        write_output_catalog(output_catalog_record_generator, output_path, output_format)

    if args.verbose:
        print_catalog_stats(interval_trees, has_source_field=args.add_source_field)


if __name__ == "__main__":
    main()
