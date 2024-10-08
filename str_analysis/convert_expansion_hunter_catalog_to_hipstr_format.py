"""This script converts an ExpansionHunter variant catalog to a GangSTR spec. Variant catalog entries
that include adjacent repeats (eg. with LocusStructure like "(A)*(ACG)*") are split into multiple GangSTR
specs - one per repeat.
"""

import argparse
import gzip
import ijson
import re
import tqdm
from str_analysis.utils.eh_catalog_utils import parse_motifs_from_locus_structure

from str_analysis.utils.misc_utils import parse_interval



def main():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--verbose", action="store_true", help="Print verbose output")
    p.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    p.add_argument("-o", "--output-file", help="bed file output path")
    p.add_argument("variant_catalog")
    args = p.parse_args()

    if not args.output_file:
        args.output_file = re.sub(".json(.gz)?$", "", args.variant_catalog) + ".hipstr_spec.bed"

    process_variant_catalog(args.variant_catalog, args.output_file, verbose=args.verbose, show_progress_bar=args.show_progress_bar)


def process_variant_catalog(variant_catalog_path, output_file_path, verbose=False, show_progress_bar=False):
    print(f"Parsing {variant_catalog_path}")
    fopen = gzip.open if variant_catalog_path.endswith("gz") else open
    with fopen(variant_catalog_path, "rt") as f:
        iterator = ijson.items(f, "item", use_float=True)
        if show_progress_bar:
            iterator = tqdm.tqdm(iterator, unit=" variant catalog records", unit_scale=True)
        with (gzip.open if output_file_path.endswith("gz") else open)(output_file_path, "wt") as f2:
            skipped_locus = 0
            for record in iterator:
                locus_id = record["LocusId"]
                motifs = parse_motifs_from_locus_structure(record["LocusStructure"])
                if not motifs:
                    raise ValueError(f"Unable to parse LocusStructure '{record['LocusStructure']}' in variant catalog record: {record}")

                reference_regions = record["ReferenceRegion"]
                if not isinstance(reference_regions, list):
                    reference_regions = [reference_regions]

                if len(motifs) != len(reference_regions):
                    raise ValueError(f"len(motifs) != len(reference_regions) in variant catalog record: {record}")

                for reference_region, motif in zip(reference_regions, motifs):
                    chrom, start_0based, end_1based = parse_interval(reference_region)
                    if start_0based >= end_1based:
                        skipped_locus += 1
                        print(f"WARNING: Skipping record #{skipped_locus} {locus_id} because the interval has width {end_1based - start_0based}bp")
                        continue

                    if len(motif) > 9 or (end_1based - start_0based) <= 1:
                        # HipSTR doesn't support motifs longer than 9bp or loci where start_1based == stop_1based
                        # (https://github.com/tfwillems/HipSTR/blob/master/src/region.cpp#L33-L35)
                        if verbose and len(motif) > 9:
                            skipped_locus += 1
                            print(f"Skipping record #{skipped_locus} {locus_id} because the motif > 9bp")
                        continue

                    motif_count = round((end_1based - start_0based)/len(motif), 3)
                    f2.write("\t".join(map(str, [
                        chrom,
                        start_0based + 1,
                        end_1based,
                        len(motif),
                        motif_count,
                        locus_id,
                        motif,
                    ])) + "\n")

    print(f"Wrote out {output_file_path}")


if __name__ == "__main__":
    main()
