#!/usr/bin/env python3

"""This script takes an ExpansionHunter catalog and outputs a new catalog where adjacent loci have been combined
into a single spec with multiple adjacent repeats inorder to improve genotyping accuracy.
"""


import argparse
import collections
import importlib.util
import intervaltree
import simplejson as json
import ijson
import os
import pysam
import re
import tempfile
import tqdm

from str_analysis.utils.eh_catalog_utils import JsonListFileWriter, parse_motifs_from_locus_structure
from str_analysis.utils.file_utils import open_file, file_exists
from str_analysis.utils.get_adjacent_repeats import get_adjacent_repeats, \
    MAX_DISTANCE_BETWEEN_REPEATS, MAX_TOTAL_ADJACENT_REGION_SIZE, MAX_OVERLAP_BETWEEN_ADJACENT_REPEATS
from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.file_utils import download_local_copy
from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif


EXTRA_FIELD_CHOICES = {
    "TRsInRegion",
    "DistancesBetweenReferenceRegions",
    "ReferenceRegionMotifs",
    "MaxDistanceBetweenReferenceRegions",
}

# the summary printed for each catalog lists at most this many locus ids of each kind
NUM_EXAMPLE_LOCUS_IDS_TO_PRINT = 10

def process_input_record(
    record, reference_fasta, interval_tree,
    max_distance_between_adjacent_repeats=MAX_DISTANCE_BETWEEN_REPEATS,
    max_total_adjacent_region_size=MAX_TOTAL_ADJACENT_REGION_SIZE,
    max_overlap_between_adjacent_repeats=MAX_OVERLAP_BETWEEN_ADJACENT_REPEATS,
    max_adjacent_repeats=None,
    add_extra_info_to_locus_id=False,
    add_extra_field=None,
    only_add_extra_fields=False,
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
        add_extra_field (set) Set of extra fields to add to the output record (see --add-extra-field).
        only_add_extra_fields (bool) If True, only add the extra field(s) to the output variant catalog. Don't add
            adjacent loci.
        set_variant_ids_to_locus_id_plus_motif (bool) If True, set the variant id of adjacent loci to the locus id
            followed by the motif. This matches the naming convention in the official Illumina catalog. For
            example: FXN_A is the adjacent poly-A repeat next to the main FXN locus. If False, use a longer, more
            descriptive variant id like "FXN_ADJACENT_LEFT_1_A".

    Return:
        2-tuple (output_record, num_adjacent_repeats_found):
            output_record (dict): A copy of the input record, modified to include any adjacent repeats found in the
                interval_tree (or, with only_add_extra_fields, only the requested extra fields).
            num_adjacent_repeats_found (int): The number of adjacent repeats found to the left and right of the locus.
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
        reference_region_motifs = []

        primary_repeat = f"{chrom}:{start_0based}-{end_1based}-{repeat_unit}"
        adjacent_left = [
            (f"{locus_id}_ADJACENT_LEFT_{i+1}", adjacent_repeat_left) for i, adjacent_repeat_left in enumerate(adjacent_repeats_left)
        ]
        adjacent_right = [
            (f"{locus_id}_ADJACENT_RIGHT_{i+1}", adjacent_repeat_right) for i, adjacent_repeat_right in enumerate(adjacent_repeats_right)
        ]

        for variant_id, adjacent_repeat in adjacent_left + [(locus_id, primary_repeat)] + adjacent_right:
            adjacent_repeat_chrom, adjacent_repeat_start_0based, adjacent_repeat_end_1based, adjacent_repeat_unit = re.split("[:-]", adjacent_repeat)
            reference_regions.append(f"{adjacent_repeat_chrom}:{adjacent_repeat_start_0based}-{adjacent_repeat_end_1based}")
            if variant_id == locus_id:
                variant_ids.append(locus_id)
            elif set_variant_ids_to_locus_id_plus_motif:
                variant_ids.append(f"{locus_id}_{adjacent_repeat_unit}")
            else:
                variant_ids.append(f"{variant_id}_{adjacent_repeat_unit}")

            # keep the main locus's own VariantType (for example, RareRepeat)
            variant_types.append(record.get("VariantType", "Repeat") if variant_id == locus_id else "Repeat")
            reference_region_motifs.append(adjacent_repeat_unit)

        if not only_add_extra_fields:
            output_record["ReferenceRegion"] = reference_regions
            output_record["LocusStructure"] = adjacent_repeats_locus_structure
            output_record["VariantId"] = variant_ids
            output_record["VariantType"] = variant_types

        if add_extra_info_to_locus_id:
            output_record["LocusId"] += (
                f";{len(reference_regions) - 1}_adjacent_repeats"
                f";max_dist_{max(compute_distances_between_reference_regions(reference_regions))}bp")
    else:
        reference_regions = [locus_interval_0based]
        reference_region_motifs = [repeat_unit]

    output_record.update(compute_extra_fields(reference_regions, reference_region_motifs, add_extra_field))

    return output_record, len(adjacent_repeats_left) + len(adjacent_repeats_right)


def compute_distances_between_reference_regions(reference_regions):
    """Returns the number of base pairs between each pair of consecutive reference regions.

    Args:
        reference_regions (list): reference regions in order, as strings like "chr1:100-150" (0-based start)

    Return:
        list: one distance per pair of consecutive regions (empty for a single region)
    """
    intervals = [parse_interval(reference_region) for reference_region in reference_regions]
    return [
        start_0based - previous_end_1based
        for (_, _, previous_end_1based), (_, start_0based, _) in zip(intervals, intervals[1:])
    ]


def compute_extra_fields(reference_regions, reference_region_motifs, add_extra_field):
    """Computes the requested extra fields (see --add-extra-field) for a locus made up of the given repeats.

    Args:
        reference_regions (list): the locus's reference regions in order, as strings like "chr1:100-150"
        reference_region_motifs (list): the motif of each of the reference regions
        add_extra_field (list): the names of the extra fields to compute, or None

    Return:
        dict: extra field name => value. TRsInRegion is included for every locus. The other fields are included only
            for a locus with more than one reference region (that is, with adjacent repeats).
    """
    extra_fields = {}
    if not add_extra_field:
        return extra_fields

    distances_between_reference_regions = compute_distances_between_reference_regions(reference_regions)
    if distances_between_reference_regions:
        if "DistancesBetweenReferenceRegions" in add_extra_field:
            extra_fields["DistancesBetweenReferenceRegions"] = distances_between_reference_regions
        if "ReferenceRegionMotifs" in add_extra_field:
            extra_fields["ReferenceRegionMotifs"] = reference_region_motifs
        if "MaxDistanceBetweenReferenceRegions" in add_extra_field:
            extra_fields["MaxDistanceBetweenReferenceRegions"] = max(distances_between_reference_regions)

    if "TRsInRegion" in add_extra_field:
        extra_fields["TRsInRegion"] = len(reference_regions)

    return extra_fields


def generate_bed_records_in_region(bed_file_path, tabix_index_path, chrom, start, end):
    """Yields the records of the given bed file that overlap the given region. When the bed file's tabix index exists,
    only the region is read. Otherwise, pybedtools.BedTool.all_hits loads the whole bed file into memory, which takes
    seconds and GBs of memory for a genome-wide bed file.

    Args:
        bed_file_path (str): local path of a bed file whose 4th column is the repeat motif
        tabix_index_path (str): local path of the bed file's tabix index (.tbi), which may not exist
        chrom (str): chromosome name
        start (int): 0-based start coordinate of the region
        end (int): end coordinate of the region

    Yields:
        a bed record with start, end and name (motif) attributes
    """
    if os.path.isfile(tabix_index_path):
        with pysam.TabixFile(bed_file_path, index=tabix_index_path) as tabix_file:
            if chrom in tabix_file.contigs:
                yield from tabix_file.fetch(chrom, start, end, parser=pysam.asBed())
    else:
        # pybedtools isn't one of str_analysis's required packages, so it's only needed for bed files without an index
        import pybedtools
        yield from pybedtools.BedTool(bed_file_path).all_hits(pybedtools.Interval(chrom, start, end))


def get_interval_tree_for_chrom(adjacent_loci_bed_file_path, chrom, start, end, verbose=False):
    """Load intervals from the given bed file into an IntervalTree for faster overlap queries. The bed file should
    contain all TR loci found in the reference genome.

    Args:
        adjacent_loci_bed_file_path (str): Path to a bed file containing adjacent loci
        chrom (str): Chromosome name
        start (int): Start coordinate of the chromosome region that spans all adjacent repeats that we want to load
        end (int): End coordinate of the chromosome region that spans all adjacent repeats that we want to load
        verbose (bool): If True, show progress bar

    Return:
        intervaltree.IntervalTree: An IntervalTree containing 0-based intervals that represent the reference start and
            end coordinates of all TR loci found in the reference genome on the given chromosome.
    """
    tabix_index_path = f"{adjacent_loci_bed_file_path}.tbi"
    if adjacent_loci_bed_file_path.startswith("gs://"):
        # download the bed file's index too, if it has one. Its local copy gets a name of its own rather than the local
        # bed file's name plus .tbi, so its path is passed on separately.
        if file_exists(tabix_index_path):
            tabix_index_path = download_local_copy(tabix_index_path)
        adjacent_loci_bed_file_path = download_local_copy(adjacent_loci_bed_file_path)
    adjacent_loci_bed_file_path = os.path.expanduser(adjacent_loci_bed_file_path)
    tabix_index_path = os.path.expanduser(tabix_index_path)

    interval_tree = intervaltree.IntervalTree()
    print(f"Parsing potential adjacent repeats from a {end-start+1:,d}bp region ({chrom}:{start:,d}-{end:,d}) in "
          f"{adjacent_loci_bed_file_path}...")

    bed_iterator = generate_bed_records_in_region(adjacent_loci_bed_file_path, tabix_index_path, chrom, start, end)
    if verbose:
        bed_iterator = tqdm.tqdm(bed_iterator, unit=" rows", unit_scale=True)

    num_rows_without_reference_bases = 0
    for bed_record in bed_iterator:
        if bed_record.end <= bed_record.start:
            # a row with no reference bases, like the insertions that filter_vcf_to_tandem_repeats.py outputs, can't be
            # an adjacent reference repeat, and IntervalTree doesn't allow empty intervals
            num_rows_without_reference_bases += 1
            continue

        # IntervalTree.add skips an interval (start, end, motif) that's already in the tree, so duplicate bed rows are
        # stored once
        interval_tree.add(intervaltree.Interval(bed_record.start, bed_record.end, data=bed_record.name))

    if num_rows_without_reference_bases:
        print(f"WARNING: skipped {num_rows_without_reference_bases:,d} rows whose end isn't after their start in "
              f"{adjacent_loci_bed_file_path} for chromosome {chrom}")

    if not interval_tree:
        print(f"WARNING: no adjacent loci found in {adjacent_loci_bed_file_path} for chromosome {chrom}")

    return interval_tree


def convert_region_to_reference_chromosome_naming(region, reference_chromosomes):
    """Returns the given region (e.g. "chr1:100-200") with its chromosome named the way the reference names it. If the
    reference has the region's chromosome, the region is returned unchanged. Otherwise, its 'chr' prefix is removed (or
    added, if it has none) when the reference has that name. If the reference has neither name, the region is returned
    unchanged.

    Args:
        region (str): a region like "chr1:100-200"
        reference_chromosomes (set): the chromosome names in the reference

    Return:
        str: the region with the reference's chromosome naming
    """
    if region.split(":")[0] in reference_chromosomes:
        return region
    toggled_region = region[len("chr"):] if region.startswith("chr") else f"chr{region}"
    return toggled_region if toggled_region.split(":")[0] in reference_chromosomes else region


def convert_record_to_reference_chromosome_naming(variant_catalog_record, reference_chromosomes):
    """Gives the ReferenceRegion, MainReferenceRegion and OfftargetRegions of the given variant catalog record the
    reference's chromosome naming (with or without a 'chr' prefix). Each region is converted on its own, since some
    catalogs name their off-target regions differently from their reference regions.

    Args:
        variant_catalog_record (dict): A variant catalog record. Its ReferenceRegion can be a string or a list. It's
            modified in place.
        reference_chromosomes (set): the chromosome names in the reference

    Return:
        str: the chromosome of the (first) ReferenceRegion after the conversion
    """
    reference_regions = variant_catalog_record["ReferenceRegion"]
    is_list = isinstance(reference_regions, list)
    converted_reference_regions = [
        convert_region_to_reference_chromosome_naming(reference_region, reference_chromosomes)
        for reference_region in (reference_regions if is_list else [reference_regions])
    ]
    variant_catalog_record["ReferenceRegion"] = converted_reference_regions if is_list else converted_reference_regions[0]

    if "MainReferenceRegion" in variant_catalog_record:
        variant_catalog_record["MainReferenceRegion"] = convert_region_to_reference_chromosome_naming(
            variant_catalog_record["MainReferenceRegion"], reference_chromosomes)

    if "OfftargetRegions" in variant_catalog_record:
        variant_catalog_record["OfftargetRegions"] = [
            convert_region_to_reference_chromosome_naming(region, reference_chromosomes)
            for region in variant_catalog_record["OfftargetRegions"]]

    return converted_reference_regions[0].split(":")[0]


def format_example_locus_ids(example_locus_ids, total_num_locus_ids):
    """Returns the given example locus ids joined by commas, followed by '...' if they're only some of the
    total_num_locus_ids locus ids."""
    return ", ".join(example_locus_ids) + (" ..." if total_num_locus_ids > len(example_locus_ids) else "")


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
                        help="The local or Google Storage (gs://) path of a .bed file that specifies repeat intervals "
                             "in the reference genome. It should have 0-based coordinates and the name field (column 4)"
                             " should contain the repeat motif. This can be a large catalog of loci generated using "
                             "TandemRepeatFinder. This .bed file will be used as the source for adjacent loci.")
    parser.add_argument("--add-extra-info-to-locus-id",
                        action="store_true",
                        help="Add info about adjacent loci to the locus id in the output variant catalog. This will "
                             "include the number of adjacent loci and the max distance between them.")
    parser.add_argument("-a", "--add-extra-field",
                        choices=EXTRA_FIELD_CHOICES,
                        action="append",
                        help="Add this extra field to the loci in the output variant catalog. TRsInRegion: the number "
                             "of TRs in the locus, counting the main repeat and any adjacent repeats (added to every "
                             "locus except those excluded with -x). DistancesBetweenReferenceRegions: the number of base pairs between each pair of "
                             "consecutive repeats. MaxDistanceBetweenReferenceRegions: the largest of those distances. "
                             "ReferenceRegionMotifs: the motif of each repeat. These last three are added only to loci "
                             "that have adjacent repeats. This option can be specified more than once.")
    parser.add_argument("--only-add-extra-fields",
                        action="store_true",
                        help="Only add the extra field(s) to the output variant catalog. Don't add adjacent loci.")
    parser.add_argument("-k", "--also-output-original-definitions",
                        action="store_true",
                        help="When adding adjacent loci to locus definions in the input catalog, output the original "
                             "definition without adjacent loci in addition to the new definition with adjacent loci."
                             "This can be used to compare whether adding adjacent loci to a locus definition has any "
                             "effect on the resulting genotypes.")
    parser.add_argument("-x", "--dont-add-adjacent-loci-to-this-locus-id", action="append", help="Don't compute "
                        "adjacent repeats for the given locus id, and instead just add it to the output catalog "
                        "without modification. This option can be specified more than once.")
    parser.add_argument("--output-dir",
                        help="Output directory. If not specified, the input catalog's directory will be used.")
    parser.add_argument("--temp-dir",
                        help="Directory for the temporary files that hold the input catalog's records, split by "
                             "chromosome, while they're processed. They take about as much space as the input catalog "
                             "uncompressed and without indentation: about 11 GB for the 5.66 million records of the "
                             "annotated TRExplorer v2.1 hg38 catalog. If not specified, the system's temporary "
                             "directory (TMPDIR) is used.")
    parser.add_argument("-o", "--output-catalog",
                        help="Path where to write the output catalog. If not specified, it will be based on the input "
                             "catalog path")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print more logging output")
    parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    parser.add_argument("input_catalog_paths",
                        nargs="+",
                        help="ExpansionHunter catalog json file path. If more than one catalog path is specified, the "
                             "catalogs will be processed sequentially, and the output catalog arg will be used as "
                             "an output filename suffix.")
    args = parser.parse_args()

    # validate command-line args
    for path in [args.ref_fasta, args.source_of_adjacent_loci] + args.input_catalog_paths:
        if not file_exists(path):
            parser.error(f"{path} not found")

    if args.only_add_extra_fields and not args.add_extra_field:
        parser.error("--only-add-extra-fields specified but no extra fields requested via --add-extra-field")

    if (not args.source_of_adjacent_loci.endswith(".bed") and
        not args.source_of_adjacent_loci.endswith(".bed.gz") and
        not args.source_of_adjacent_loci.endswith(".bed.bgz")):
        parser.error(f"source_of_adjacent_loci must be a .bed or .bed.gz file: {args.source_of_adjacent_loci}")

    # a bed file without a tabix index is read with pybedtools, which isn't one of str_analysis's required packages
    if not file_exists(f"{args.source_of_adjacent_loci}.tbi") and importlib.util.find_spec("pybedtools") is None:
        parser.error(f"{args.source_of_adjacent_loci} has no tabix index (.tbi), and reading it without one needs "
                     f"pybedtools, which isn't installed. Please bgzip and tabix-index it, or install pybedtools.")

    # name each input catalog's output catalog, and check that its directory exists now rather than after the slow first
    # pass below
    output_catalog_path_by_input_catalog_path = {}
    for input_catalog_path in args.input_catalog_paths:
        if args.output_catalog and len(args.input_catalog_paths) == 1:
            output_catalog_path = args.output_catalog
        else:
            # name the output after the input catalog, with --output-catalog as a suffix, and keep its gzip compression
            output_catalog_path = re.sub(r"\.json(\.b?gz)?$", "", input_catalog_path) + \
                f".{args.output_catalog or 'with_adjacent_loci'}.json"
            if input_catalog_path.endswith((".gz", ".bgz")):
                output_catalog_path += ".gz"

        if args.output_dir:
            output_catalog_path = os.path.join(args.output_dir, os.path.basename(output_catalog_path))

        output_catalog_dir = os.path.dirname(output_catalog_path)
        if output_catalog_dir and not os.path.isdir(output_catalog_dir):
            parser.error(f"Directory {output_catalog_dir} not found for output catalog {output_catalog_path}")
        output_catalog_path_by_input_catalog_path[input_catalog_path] = output_catalog_path

    ref_fasta =pysam.FastaFile(os.path.abspath(os.path.expanduser(args.ref_fasta)))
    ref_fasta_chromosomes = set(ref_fasta.references)

    # Read the input catalog(s) one record at a time, since a whole annotated catalog takes tens of GB of memory, and
    # write each catalog's records to one temporary file per chromosome, named by the chromosome of the record's
    # ReferenceRegion in the reference's naming. The records are then processed and written out one chromosome at a
    # time, in the order of these chromosome names, keeping their order within each chromosome.
    temp_dir = tempfile.TemporaryDirectory(prefix="add_adjacent_loci_", dir=args.temp_dir)
    temp_file_paths_by_chrom_by_input_catalog_path = {}
    num_records_by_input_catalog_path = {}
    min_max_coords_for_chrom = collections.defaultdict(lambda: (10**10, 0))
    for input_catalog_i, input_catalog_path in enumerate(args.input_catalog_paths):
        print(f"Loading input catalog {input_catalog_path}")
        temp_files_by_chrom = {}
        num_records = 0
        with open_file(input_catalog_path, is_text_file=True) as f:
            records = ijson.items(f, "item", use_float=True)
            if args.show_progress_bar:
                records = tqdm.tqdm(records, unit=" loci", unit_scale=True)

            for record in records:
                if args.n and num_records >= args.n:
                    break

                missing_keys = {"ReferenceRegion", "LocusStructure", "LocusId"} - set(record.keys())
                if missing_keys:
                    print(f"WARNING: record {record} is missing keys: {missing_keys}. Skipping...")
                    continue
                num_records += 1

                # give the record's regions the reference's chromosome naming, and get the min and max coordinates of
                # the TR loci for each chromosome
                chrom = convert_record_to_reference_chromosome_naming(record, ref_fasta_chromosomes)
                if chrom not in ref_fasta_chromosomes:
                    raise ValueError(f"chromosome '{chrom}' not found in reference fasta {args.ref_fasta} which "
                                     f"has chromosomes: " + ", ".join(sorted(ref_fasta_chromosomes)))

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

                if chrom not in temp_files_by_chrom:
                    temp_files_by_chrom[chrom] = open(os.path.join(
                        temp_dir.name, f"catalog{input_catalog_i}.chrom{len(temp_files_by_chrom)}.jsonl"), "wt")
                temp_files_by_chrom[chrom].write(json.dumps(record) + "\n")

        for temp_file in temp_files_by_chrom.values():
            temp_file.close()
        temp_file_paths_by_chrom_by_input_catalog_path[input_catalog_path] = {
            chrom: temp_file.name for chrom, temp_file in temp_files_by_chrom.items()}
        num_records_by_input_catalog_path[input_catalog_path] = num_records
        print(f"Loaded {num_records:,d} records from {input_catalog_path}")

    print(f"Processing {len(args.input_catalog_paths)} input variant catalog" + ("s" if len(args.input_catalog_paths) > 1 else f": {args.input_catalog_paths[0]}"))

    # chromosome (in the reference's naming) => IntervalTree of the potential adjacent loci in the region of that
    # chromosome that the input catalog(s) span. Parsing args.source_of_adjacent_loci takes seconds per chromosome. With
    # one input catalog, whose records are processed one chromosome at a time, only the current chromosome's tree is
    # kept, to save memory. With several input catalogs, every chromosome's tree is kept for the rest of the run, so that
    # each chromosome is parsed only once rather than once per catalog.
    interval_trees_by_chrom = {}

    # process the records from the input catalog(s) to add adjacent loci
    for input_catalog_path, temp_file_paths_by_chrom in temp_file_paths_by_chrom_by_input_catalog_path.items():
        print(f"Processing {input_catalog_path}")
        counters = collections.defaultdict(int)
        # the first NUM_EXAMPLE_LOCUS_IDS_TO_PRINT locus ids of each kind, for the summary
        example_locus_ids_with_existing_adjacent_repeats = []
        example_locus_ids_with_new_adjacent_repeats = []
        example_locus_ids_of_duplicate_records = []

        output_catalog_path = output_catalog_path_by_input_catalog_path[input_catalog_path]
        output_writer = JsonListFileWriter(output_catalog_path, use_gzip=output_catalog_path.endswith((".gz", ".bgz")))
        progress_bar = tqdm.tqdm(total=num_records_by_input_catalog_path[input_catalog_path], unit=" loci",
                                 unit_scale=True, disable=not args.show_progress_bar)
        for chrom in sorted(temp_file_paths_by_chrom):
            # avoid duplicate records in the output catalog. Duplicates have the same ReferenceRegion, so they're on the
            # same chromosome.
            already_added_to_output_catalog = set()

            with open(temp_file_paths_by_chrom[chrom]) as temp_file:
                for line in temp_file:
                    progress_bar.update(1)
                    record = json.loads(line)
                    counters["total records"] += 1

                    if args.dont_add_adjacent_loci_to_this_locus_id and record["LocusId"] in args.dont_add_adjacent_loci_to_this_locus_id:
                        print(f"Adding {record['LocusId']} without modification as requested by -x arg")
                        output_writer.write(record)
                        counters["excluded locus id"] += 1
                        continue

                    # handle input catalog records that already have adjacent loci
                    if isinstance(record["ReferenceRegion"], list):
                        print(f"WARNING: record {record} already has adjacent intervals specified in the input variant catalog. "
                              f"Will add it to the output variant catalog without adding adjacent loci to it...")
                        record.update(compute_extra_fields(
                            record["ReferenceRegion"], parse_motifs_from_locus_structure(record["LocusStructure"]),
                            args.add_extra_field))
                        output_writer.write(record)
                        counters["already had adjacent loci"] += 1
                        if len(example_locus_ids_with_existing_adjacent_repeats) < NUM_EXAMPLE_LOCUS_IDS_TO_PRINT:
                            example_locus_ids_with_existing_adjacent_repeats.append(record["LocusId"])
                        continue

                    # parse args.source_of_adjacent_loci within the region of this chromosome that the input catalog(s)
                    # span, unless it's already loaded
                    if chrom not in interval_trees_by_chrom:
                        if len(args.input_catalog_paths) == 1:
                            interval_trees_by_chrom.clear()
                        min_start, max_end = min_max_coords_for_chrom[chrom]
                        interval_trees_by_chrom[chrom] = get_interval_tree_for_chrom(
                            args.source_of_adjacent_loci, chrom, min_start, max_end, verbose=args.verbose)

                    # add any adjacent loci to the current record
                    output_record, num_adjacent_repeats_found = process_input_record(
                        record,
                        ref_fasta,
                        interval_trees_by_chrom[chrom],
                        max_distance_between_adjacent_repeats=args.max_distance_between_adjacent_repeats,
                        max_total_adjacent_region_size=args.max_total_adjacent_region_size,
                        max_overlap_between_adjacent_repeats=args.max_overlap_between_adjacent_repeats,
                        max_adjacent_repeats=args.max_adjacent_repeats,
                        add_extra_info_to_locus_id=args.add_extra_info_to_locus_id,
                        add_extra_field=args.add_extra_field,
                        only_add_extra_fields=args.only_add_extra_fields,
                        set_variant_ids_to_locus_id_plus_motif=False,
                    )

                    adjacent_loci_were_added = isinstance(output_record["ReferenceRegion"], list)
                    if adjacent_loci_were_added:
                        output_record_key = output_record["LocusStructure"] + ";" + ";".join(output_record["ReferenceRegion"])
                    else:
                        output_record_key = output_record["LocusStructure"] + ";" + output_record["ReferenceRegion"]

                    # with -k, also write the record's original definition (without adjacent loci), unless an
                    # identical original definition was already written. This is done before the duplicate check
                    # below, since a record whose definition with adjacent loci matches an earlier record's can still
                    # have its own original definition.
                    if adjacent_loci_were_added and args.also_output_original_definitions:
                        original_record_key = record["LocusStructure"] + ";" + record["ReferenceRegion"]
                        if original_record_key not in already_added_to_output_catalog:
                            already_added_to_output_catalog.add(original_record_key)
                            output_writer.write({
                                **record,
                                "LocusId": record["LocusId"] + "_ORIGINAL",
                                **compute_extra_fields([record["ReferenceRegion"]],
                                                       parse_motifs_from_locus_structure(record["LocusStructure"]),
                                                       args.add_extra_field),
                            })

                    # skip this record if an earlier record already added the same locus definition to the output
                    # catalog, so the summary and stats below count only the records that are written
                    if output_record_key in already_added_to_output_catalog:
                        counters["duplicate records"] += 1
                        if len(example_locus_ids_of_duplicate_records) < NUM_EXAMPLE_LOCUS_IDS_TO_PRINT:
                            example_locus_ids_of_duplicate_records.append(record["LocusId"])
                        continue
                    already_added_to_output_catalog.add(output_record_key)

                    if num_adjacent_repeats_found > 0:
                        counters["found adjacent loci"] += 1
                        if len(example_locus_ids_with_new_adjacent_repeats) < NUM_EXAMPLE_LOCUS_IDS_TO_PRINT:
                            example_locus_ids_with_new_adjacent_repeats.append(record["LocusId"])

                    if adjacent_loci_were_added and args.also_output_original_definitions:
                        output_record["LocusId"] = output_record["LocusId"] + "_MAIN"
                    output_writer.write(output_record)

                    # compute stats about distances. These use the extra fields rather than ReferenceRegion, since
                    # --only-add-extra-fields leaves ReferenceRegion as a single region.
                    if "DistancesBetweenReferenceRegions" in output_record:
                        counters["total spacers"] += len(output_record["DistancesBetweenReferenceRegions"])
                        for distance in output_record["DistancesBetweenReferenceRegions"]:
                            counters[f"  {distance:3d}bp spacer"] += 1
                    if "ReferenceRegionMotifs" in output_record:
                        num_reference_regions = len(output_record["ReferenceRegionMotifs"])
                        num_unique_motifs = len(set(compute_canonical_motif(m) for m in output_record["ReferenceRegionMotifs"]))
                        same_motif = "same motif" if num_unique_motifs == 1 else "different motifs"
                        counters[f"adjacent repeats count: {num_reference_regions:3d} reference regions with {same_motif}"] += 1

        output_writer.close()
        progress_bar.close()

        # an empty input catalog has no records, so divide by at least 1 when computing percentages
        total_records_for_percentages = max(counters["total records"], 1)

        if counters["already had adjacent loci"]:
            print(f"Added {counters['already had adjacent loci']:,d} out of {counters['total records']:,d} records "
                  f"to the output catalog without adding adjacent loci since they already had adjacent repeats provided "
                  f"in the input catalog: ",
                  format_example_locus_ids(
                      example_locus_ids_with_existing_adjacent_repeats, counters["already had adjacent loci"]))
        else:
            print("The input catalog(s) did not have any adjacent loci specified")

        # --only-add-extra-fields records the adjacent loci it finds in the extra fields without adding them
        found_or_added = "Found adjacent loci for" if args.only_add_extra_fields else "Added adjacent loci to"
        print(f"{found_or_added} {counters['found adjacent loci']:,d} out of {counters['total records']:,d} records ("
              f"{counters['found adjacent loci'] / total_records_for_percentages:.1%}): ",
              format_example_locus_ids(example_locus_ids_with_new_adjacent_repeats, counters["found adjacent loci"]))
        if counters["duplicate records"]:
            print(f"Skipped {counters['duplicate records']:,d} out of {counters['total records']:,d} records "
                  f"whose locus definition was the same as that of an earlier record: ",
                  format_example_locus_ids(example_locus_ids_of_duplicate_records, counters["duplicate records"]))
        if args.dont_add_adjacent_loci_to_this_locus_id:
            print(f"Added {counters['excluded locus id']:,d} out of {counters['total records']:,d} records ("
                  f"{counters['excluded locus id'] / total_records_for_percentages:.1%}) without modification as "
                  f"requested by the -x arg")

        for key, value in sorted(counters.items()):
            if key.startswith("total"):
                continue
            if "adjacent repeats count" in key:
                label = "records"
            elif "spacer" in key:
                label = "spacers"
            else:
                continue

            total = counters[f"total {label}"]
            print(f"{value:8,d} out of {total:8,d} {label} ({value / total:6.1%}) had {key}")

        print(f"Wrote {output_writer.num_records_written:,d} total records to {output_catalog_path}")

    temp_dir.cleanup()


if __name__ == '__main__':
    main()
