"""This script takes an ExpansionHunter variant catalog and adds off-target regions"""

import argparse
import collections
import simplejson as json
import re
import sqlite3

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.file_utils import open_file, download_local_copy

MIN_FRACTION_SIMULATED_READS_ABSORBED = 0.01
MIN_SIMULATED_READS_ABSORBED = 1
MAX_OFFTARGET_REGIONS = 10

OFFTARGET_DB_PATHS = {
    "37": "gs://str-analysis/hg37/ref/repeats_db/offtargets_GRCh37.db",
    "38": "gs://str-analysis/hg38/ref/repeats_db/offtargets_GRCh38.db",
}


def add_offtarget_regions(
        record,
        offtarget_db_connection,
        max_offtarget_regions=MAX_OFFTARGET_REGIONS,
        min_simulated_reads_absorbed=MIN_SIMULATED_READS_ABSORBED,
        min_fraction_of_simulated_reads_absorbed=MIN_FRACTION_SIMULATED_READS_ABSORBED,
        verbose=False):
    """Add a list of off-target regions to the given record.

    Args:
        record (dict): An ExpansionHunter variant catalog record. This record will be modified in-place.
        offtarget_db_connection: sqlite3 connection to the STR database that contains precomputed off-target regions
            for all motifs.
        max_offtarget_regions (int): The maximum number of off-target regions to add to the record.
        min_simulated_reads_absorbed (int): The minimum number of simulated reads that were absorbed by an off-target
            region in order for it to be included in the off-target region list.
        min_fraction_of_simulated_reads_absorbed (float): The minimum fraction of simulated reads that were absorbed
            by an off-target region in order for it to be included in the off-target region list.
    """

    if "LocusStructure" not in record:
        raise ValueError(f"LocusStructure not found in record: {record}")
    if record["LocusStructure"].count("(") != 1 or record["LocusStructure"].count(")") != 1:
        raise ValueError(f"Unexpected LocusStructure format. Expected exactly one repeat unit in parentheses: {record}")

    repeat_unit = record["LocusStructure"].strip("()*+?")
    if set(repeat_unit) - set("ACGTN"):
        raise ValueError(f"LocusStructure repeat unit contains non-ACGT characters: {record}")

    normalized_repeat_unit = compute_canonical_motif(repeat_unit, include_reverse_complement=True)
    # retrieve off-target regions from database
    chrom, start_0based, end_1based = re.split("[:-]", record["ReferenceRegion"])
    chrom = chrom.replace("chr", "")

    results = offtarget_db_connection.execute(
        """SELECT chrom, start_1based, end_1based, read_count, fraction_of_reads
           FROM offtargets 
           WHERE normalized_repeat_unit = ?
           AND ((chrom != ? AND chrom != ?) OR start_1based > ? OR end_1based < ?)
           AND read_count >= ?
           AND fraction_of_reads >= ?
           ORDER BY read_count DESC""",
        (
            normalized_repeat_unit, f"chr{chrom}", chrom, end_1based, start_0based,
            min_simulated_reads_absorbed, min_fraction_of_simulated_reads_absorbed,
        )).fetchall()

    offtarget_regions = []
    for i, offtarget_region in enumerate(results):
        if i < max_offtarget_regions:
            offtarget_regions.append(f"{offtarget_region[0]}:{offtarget_region[1] - 1}-{offtarget_region[2]}")

    if offtarget_regions:
        record = record.copy()
        record["LocusId"] = f"{record['LocusId']}_WITH_OFFTARGETS"
        record["OfftargetRegions"] = offtarget_regions
        record["VariantType"] = "RareRepeat"

    if verbose:
        print(f"Added {len(offtarget_regions)} off-target regions to {record['LocusId']} " +
              (f"where the last off-target region accounted for {results[-1][3]} simulated read(s) "
               f"(which was {results[-1][4]:0.1%} of all simulated reads)" if offtarget_regions else ""))


    return record if offtarget_regions else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome-version",
                        required=True,
                        help="Reference genome version (37 or 38). This is used to determine the path to the "
                             "This should downloaded from gs://str-analysis/hg38/ref/repeats_db/offtargets_GRCh38.db ",
                        choices=["37", "38"])
    parser.add_argument("--offtarget-db",
                        help="Local path of sqlite database that contains the catalog of off-target regions. "
                             "If not specified, it will be downloaded from " + " or ".join(OFFTARGET_DB_PATHS.values()))
    parser.add_argument("--max-offtarget-regions",
                        type=int,
                        default=MAX_OFFTARGET_REGIONS,
                        help="Max off-target list size. The newly added off-target lists will include at most this "
                             "many of the top regions (prioritized by the number of synthetic reads they absorb).")
    parser.add_argument("--min-simulated-reads-absorbed",
                        type=int,
                        default=MIN_SIMULATED_READS_ABSORBED,
                        help="Minimum number of simulated reads that were absorbed by an off-target region in order "
                             "for it to be included in the off-target region list.")
    parser.add_argument("--min-fraction-of-simulated-reads-absorbed",
                        type=float,
                        default=MIN_FRACTION_SIMULATED_READS_ABSORBED,
                        help="Minimum fraction of simulated reads that were absorbed by an off-target region in order "
                             "for it to be included in the off-target region list.")
    parser.add_argument("-k", "--keep-original-locus-definitions", action="store_true", help="If specified, each locus "
                              "definition in the input variant catalog will be included in the output catalog twice - "
                              "once with and and once without off-target regions.")
    parser.add_argument("-o", "--output-catalog",
                        help="Path where to write the output catalog. If not specified, it will be based on the input "
                             "catalog path")
    parser.add_argument("--verbose", action="store_true", help="Print details about off-target regions at each locus")
    parser.add_argument("input_catalog_path",
                        help="ExpansionHunter catalog json file path. A list of off-target regions will be added to "
                             "all entries that don't already have one. Entries where multiple adjacent repeats "
                             "with different motifs are specified will be skipped since it's not clear which of the "
                             "adjacent repeats should be used as the main repeat.")
    args = parser.parse_args()
        
    # parse the input catalog(s)
    try:
        with open_file(args.input_catalog_path) as f:
            input_catalog = json.load(f)
    except Exception as e:
        parser.error(f"Unable to parse {args.input_catalog_path}: {e}")

    print(f"Processing {len(input_catalog):,d} loci from {args.input_catalog_path}")

    output_catalog = []

    # connect to the sqlite database that contains the off-target regions
    offtarget_db_path = OFFTARGET_DB_PATHS[args.genome_version]
    if args.offtarget_db:
        local_db_path = args.offtarget_db
    else:
        local_db_path = download_local_copy(offtarget_db_path)
    offtarget_db_connection = sqlite3.connect(local_db_path)

    counters = collections.defaultdict(int)
    for record in input_catalog:
        if isinstance(record.get("ReferenceRegion"), list):
            counters["adjacent loci filter"] += 1
            continue

        if len(record.get("LocusStructure", "N").strip("(ACGT)*+?")) > 0:
            counters["non-ACGT bases filter"] += 1
            continue

        if args.keep_original_locus_definitions:
            output_catalog.append(record)

        record_with_offtarget_regions = add_offtarget_regions(
            record,
            offtarget_db_connection,
            args.max_offtarget_regions,
            args.min_simulated_reads_absorbed,
            args.min_fraction_of_simulated_reads_absorbed,
            verbose=args.verbose
        )

        if record_with_offtarget_regions is not None:
            output_catalog.append(record_with_offtarget_regions)

    total = len(input_catalog)
    if counters["adjacent loci filter"] > 0:
        print(f"Skipped {total - counters['adjacent loci filter']:,d} out of {total:,d} records "
              f"that have adjacent loci specified")
        total -= counters["adjacent loci filter"]

    if counters["non-ACGT bases filter"] > 0:
        print(f"Skipped {total - counters['non-ACGT bases filter']:,d} out of {total:,d} records "
              f"because they have an unexpected LocusStructure or a repeat unit that includes non-ACGT nucleotides")
        total -= counters["non-ACGT bases filter"]

    # write the results to the output catalog
    if args.output_catalog:
        output_catalog_path = args.output_catalog
    else:
        output_catalog_path = re.sub(".json$", "", args.input_catalog_path) + ".with_offtarget_regions.json"

    with open(output_catalog_path, "wt") as f:
        json.dump(output_catalog, f, indent=4)

    print(f"Wrote {len(output_catalog):,d} records to {output_catalog_path}")


if __name__ == '__main__':
    main()

