"""This script converts an ExpansionHunter variant catalog to a TRGT catalog, including proper handling of
 variant catalog entries that include adjacent repeats (eg. with LocusStructure like "(A)*(ACG)*").
"""

import argparse
import gzip
import simplejson as json
import re
import tqdm

from str_analysis.utils.misc_utils import parse_interval


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-o", "--output-file", help="BED file output path")
    p.add_argument("expansion_hunter_catalog", help="ExpansionHunter variant catalog in JSON format")
    args = p.parse_args()

    if not args.output_file:
        args.output_file = re.sub(".json(.gz)?$", "", args.expansion_hunter_catalog) + ".trgt.bed"

    process_expansion_hunter_catalog(args.expansion_hunter_catalog, args.output_file)


def process_expansion_hunter_catalog(expansion_hunter_catalog_path, output_file_path):
    print(f"Parsing {expansion_hunter_catalog_path}")
    fopen = gzip.open if expansion_hunter_catalog_path.endswith("gz") else open
    with fopen(expansion_hunter_catalog_path, "rt") as f:
        expansion_hunter_catalog = json.load(f)

    print(f"Parsed {len(expansion_hunter_catalog):,d} records from {expansion_hunter_catalog_path}")

    fopen = gzip.open if output_file_path.endswith("gz") else open
    with fopen(output_file_path, "wt") as f:
        previous_chrom = None
        output_rows = []
        counter = 0
        for i, record in enumerate(tqdm.tqdm(expansion_hunter_catalog, unit=" variant catalog records")):
            locus_id = record["LocusId"]
            locus_structure = record["LocusStructure"]
            motifs = re.findall("[(]([A-Z]+)[)]", locus_structure)
            if not motifs:
                raise ValueError(f"Unable to parse LocusStructure '{locus_structure}' in variant catalog "
                                 f"record #{i+1}: {record}")

            reference_regions = record["ReferenceRegion"]
            if not isinstance(reference_regions, list):
                reference_regions = [reference_regions]

            if len(motifs) != len(reference_regions):
                raise ValueError(f"LocusStructure elements != # of entries in the list of ReferenceRegions in "
                                 f"variant catalog record #{i+1}: {record}")

            chrom = locus_start_0based = locus_end_1based = None
            for motif, reference_region in zip(motifs, reference_regions):
                chrom, start_0based, end_1based = parse_interval(reference_region)

                locus_start_0based = min(locus_start_0based, start_0based) if locus_start_0based is not None else start_0based
                locus_end_1based = max(locus_end_1based, end_1based) if locus_end_1based is not None else end_1based

            if previous_chrom is None:
                previous_chrom = chrom

            motif_string = ",".join(motifs)
            struc = "".join([f"({motif})n" for motif in motifs])

            output_rows.append([
                chrom,
                locus_start_0based,
                locus_end_1based,
                f"ID={locus_id};MOTIFS={motif_string};STRUC={struc}",
            ])

            if chrom != previous_chrom:
                for output_row in sorted(output_rows):
                    counter += 1
                    f.write("\t".join(map(str, output_row)) + "\n")

                output_rows = []
                previous_chrom = chrom

        for output_row in output_rows:
            counter += 1
            f.write("\t".join(map(str, output_row)) + "\n")

    print(f"Wrote {counter:,d} rows to {output_file_path}")


if __name__ == "__main__":
    main()
