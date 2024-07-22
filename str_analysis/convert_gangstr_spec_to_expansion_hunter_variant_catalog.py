"""This script converts a GangSTR repeat spec to an ExpansionHunter variant catalog. This simplifies the process of
switching from GangSTR to ExpansionHunter to genotype a set of loci previously genotyped using GangSTR.
"""

import argparse
import collections
import gzip
import simplejson as json
import os
from pprint import pformat
import re
import tqdm


def main():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-o", "--output-file", help="JSON file output path")
    p.add_argument("-v", "--verbose", action="store_true")
    p.add_argument("-t", "--trim", action="store_true", help="Trim locus intervals so that they contain exact multiples of the repeat unit")
    p.add_argument("gangstr_spec", help="path of the GangSTR repeat spec BED file")
    args = p.parse_args()

    if not args.output_file:
        args.output_file = re.sub(".bed(.gz)?$", "", args.gangstr_spec) + ".variant_catalog.json"

    if not os.path.isfile(args.gangstr_spec):
        p.error(f"{args.gangstr_spec} file not found")

    process_variant_catalog(args.gangstr_spec, args.output_file, trim_loci=args.trim, verbose=args.verbose)


def process_variant_catalog(gangstr_spec_path, output_file_path, trim_loci=False, verbose=False):
    print(f"Parsing {gangstr_spec_path}")
    json_records = []
    existing_locus_ids = set()
    counter = collections.defaultdict(int)
    fopen = gzip.open if gangstr_spec_path.endswith("gz") else open
    with fopen(gangstr_spec_path, "rt") as f:
        for i, row in tqdm.tqdm(enumerate(f), unit=" records", unit_scale=True):
            fields = row.strip("\n").split("\t")
            chrom = fields[0]
            start_0based = int(fields[1]) - 1
            end_1based = int(fields[2])
            repeat_unit = fields[4]

            if not repeat_unit or (len(set(repeat_unit) - set("ACGTN")) > 0):
                raise ValueError(f"Invalid repeat unit in row #{i + 1}: {repeat_unit}. Line: {row}")

            counter["total input loci"] += 1

            trim_bp = (end_1based - start_0based) % len(repeat_unit)
            if trim_bp != 0:
                counter["untrimmed locus"] += 1
            else:
                counter["trimmed locus"] += 1


            if trim_loci and trim_bp != 0:
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
            # TODO add support for off-target regions

    print(f"Writing records to {output_file_path}")
    fopen = gzip.open if output_file_path.endswith("gz") else open
    with fopen(output_file_path, "wt") as f:
        json.dump(json_records, f, indent=4, ignore_nan=True)

    print("Stats:", pformat(dict(counter)))
    print(f"Wrote {len(json_records):,d} records to {output_file_path}")


if __name__ == "__main__":
    main()
