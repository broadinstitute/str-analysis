"""This script converts a bed file (with 0-based start coordinates) to an ExpansionHunter variant catalog. It requires
the 'name' field (column 4) to contain the repeat motif.
"""

import argparse
import collections
import gzip
import simplejson as json
import os
import re
import tqdm
from str_analysis.utils.misc_utils import parse_interval

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.eh_catalog_utils import get_variant_catalog_iterator, parse_motifs_from_locus_structure


def main():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--batch-size", type=int, help="Optionally, split the output into many variant catalogs with "
                                                  "at most this many loci per catalog")
    p.add_argument("-o", "--output-path", help="JSON variant catalog output path")
    p.add_argument("--trim", action="store_true", help="Trim loci to be a multiple of the repeat unit size")
    p.add_argument("-v", "--verbose", action="store_true")
    p.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    p.add_argument("bed_or_json_path", help="Input BED or JSON file path")
    args = p.parse_args()

    if not os.path.isfile(args.bed_or_json_path):
        p.error(f"{args.bed_or_json_path} file not found")

    # parse bed file and convert to variant catalog json format
    print(f"Parsing {args.bed_or_json_path}")
    json_records = parse_file(args.bed_or_json_path, trim=args.trim, verbose=args.verbose,
                              show_progress_bar=args.show_progress_bar)

    # sort records by normalized motif to maximize cache hit rate in the optimized version of ExpansionHunter
    print("Sorting records by normalized motif")
    def get_canonical_motif(record):
        motif = record["LocusStructure"].strip("()*")
        return compute_canonical_motif(motif, include_reverse_complement=True)

    json_records = sorted(json_records, key=get_canonical_motif)

    # write json records to output files
    if not args.output_path:
        output_path_prefix = re.sub(".bed(.gz)?$", "", args.bed_or_json_path)
    else:
        output_path_prefix = re.sub(".json(.gz)?$", "", args.output_path)

    if not args.batch_size:
        if not args.output_path:
            args.output_path = output_path_prefix + ".variant_catalog.json"

        fopen = gzip.open if args.output_path.endswith("gz") else open
        with fopen(args.output_path, "wt") as f:
            json.dump(json_records, f, indent=4)
        print(f"Wrote {len(json_records):,d} to {args.output_path}")
    else:
        output_path_prefix += f".{len(json_records)}_loci"
        for i in range(0, len(json_records), args.batch_size):
            output_path = output_path_prefix + f".batch{i // args.batch_size:03d}.json"
            with open(output_path, "wt") as f:
                json.dump(json_records[i:i + args.batch_size], f, indent=4)
        print(f"Wrote {len(json_records):,d} to {output_path_prefix}*.json")


def parse_file(bed_or_json_path, trim=True, verbose=False, show_progress_bar=False):
    json_records = []
    existing_locus_ids = set()
    counter = collections.defaultdict(int)
    fopen = gzip.open if bed_or_json_path.endswith("gz") else open
    for i, record in enumerate(get_variant_catalog_iterator(
        bed_or_json_path, show_progress_bar=show_progress_bar)):
        if isinstance(record["ReferenceRegion"], list):
            print(f"WARNING: Locus definition with adjacent repeats is not supported. {record}. Skipping...")
            continue

        chrom, start_0based, end_1based = parse_interval(record["ReferenceRegion"])
        repeat_units = parse_motifs_from_locus_structure(record["LocusStructure"])
        if len(repeat_units) != 1:
            print(f"WARNING: Locus with multiple motifs is not supported. {record}. Skipping...")
            continue
        repeat_unit = repeat_units[0]
        if not repeat_unit or (len(set(repeat_unit) - set("ACGTN")) > 0):
            raise ValueError(f"Unexpected characters in repeat unit in row #{i + 1}: {repeat_unit}. Record: {record}")

        counter["total input loci"] += 1
        if trim:
            trim_bp = (end_1based - start_0based) % len(repeat_unit)
            if trim_bp != 0:
                counter["trimmed locus"] += 1
                if verbose:
                    print(f"WARNING: {chrom}:{start_0based}-{end_1based} interval has size {end_1based - start_0based} "
                          f"which is not a multiple of the repeat unit {repeat_unit} (size {len(repeat_unit)}). "
                          f"Changing it to {chrom}:{start_0based}-{end_1based - trim_bp}")
                end_1based -= trim_bp
                assert (end_1based - start_0based) % len(repeat_unit) == 0

        locus_id = f"{chrom}-{start_0based}-{end_1based}-{repeat_unit}"
        if locus_id in existing_locus_ids:
            counter["skipped duplicate"] += 1
            if verbose:
                print(f"WARNING: skipping duplicate locus id: {locus_id}")
            continue

        existing_locus_ids.add(locus_id)

        json_records.append({
            "LocusId": locus_id,
            "ReferenceRegion": f"{chrom}:{start_0based}-{end_1based}",
            "LocusStructure": f"({repeat_unit})*",
            "VariantType": "Repeat",
        })

    return json_records

if __name__ == "__main__":
    main()
