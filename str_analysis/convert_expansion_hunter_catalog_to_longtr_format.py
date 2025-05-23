"""This script converts an ExpansionHunter variant catalog to a GangSTR spec. Variant catalog entries
that include adjacent repeats (eg. with LocusStructure like "(A)*(ACG)*") are split into multiple GangSTR
specs - one per repeat.
"""

import argparse
import collections
import gzip
import ijson
import re
import tqdm
from str_analysis.utils.eh_catalog_utils import parse_motifs_from_locus_structure

from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif


def main():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--verbose", action="store_true", help="Print verbose output")
    p.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    p.add_argument("-o", "--output-file", help="bed file output path")
    p.add_argument("variant_catalog")
    args = p.parse_args()

    if not args.output_file:
        args.output_file = re.sub(".json(.gz)?$", "", args.variant_catalog) + ".longtr_spec.bed"

    process_variant_catalog(args.variant_catalog, args.output_file, verbose=args.verbose, show_progress_bar=args.show_progress_bar)


def process_variant_catalog(variant_catalog_path, output_file_path, verbose=False, show_progress_bar=False):
    print(f"Parsing {variant_catalog_path}")
    fopen = gzip.open if variant_catalog_path.endswith("gz") else open
    with fopen(variant_catalog_path, "rt") as f:
        iterator = ijson.items(f, "item", use_float=True)
        if show_progress_bar:
            iterator = tqdm.tqdm(iterator, unit=" variant catalog records", unit_scale=True)
        with (gzip.open if output_file_path.endswith("gz") else open)(output_file_path, "wt") as f2:
            counters = collections.Counter()
            previously_seen_loci = set()  # dedup loci with identical chrom/start/end/canonical-motif
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

                counters["loci"] += 1
                for reference_region, motif in zip(reference_regions, motifs):
                    counters["regions"] += 1
                    chrom, start_0based, end_1based = parse_interval(reference_region)
                    if start_0based + 1 >= end_1based:
                        counters["skipped"] += 1
                        if verbose:
                            print(f"WARNING: Skipping record #{counters['skipped']} {locus_id} @ {chrom}:{start_0based}-{end_1based} because the interval has width {end_1based - start_0based}bp")
                        continue

                    counters["output"] += 1
                    f2.write("\t".join(map(str, [
                        chrom,
                        start_0based + 1,  # LongTR beds are 1-based
                        end_1based,
                        len(motif),
                        round((end_1based - start_0based)/len(motif), 3),
                        locus_id,
                        motif,
                    ])) + "\n")

                    canonical_motif = compute_canonical_motif(motif)
                    key = (chrom, start_0based, end_1based, canonical_motif)
                    if key in previously_seen_loci:
                        counters["skipped"] += 1
                        continue

    print(f"Parsed {counters['loci']:,d} locus definitions containing {counters['regions']:,d} ReferenceRegions from {variant_catalog_path}")
    print(f"Skipped {counters['skipped']:,d} records:")
    print(f"Wrote out {counters['output']:,d} records to {output_file_path}")


if __name__ == "__main__":
    main()
