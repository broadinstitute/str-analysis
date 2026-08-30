"""This script converts an ExpansionHunter variant catalog to a BED file.

The 5th column is a placeholder '.' by default. --motif-size-column puts the motif's length there
instead, which is the 5-column layout ATaRVa (https://github.com/SowpatiLab/ATaRVa) expects of its
--regions file: chromosome, start, end, motif, motif length. ATaRVa reads one motif per row, so
combine that option with --split-adjacent-repeats.
"""

import argparse
import gzip
import ijson
import os
import re
import tqdm

from str_analysis.utils.misc_utils import parse_interval


def main():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-s", "--split-adjacent-repeats", action="store_true", help="If a locus is defined "
                   "as having several adjacent repeats in the ExpansionHunter catalog, split it into separate "
                   "BED file rows. Otherwise, the entire locus will be represented as a single BED file row.")
    p.add_argument("--motif-size-column", action="store_true", help="Write the motif's length in the "
                   "5th column rather than a '.' placeholder. This is the layout ATaRVa expects of its "
                   "--regions file. Use together with --split-adjacent-repeats, since a row covering "
                   "several adjacent repeats has no single motif length to report.")
    p.add_argument("-o", "--output-file", help="BED file output path")
    p.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    p.add_argument("expansion_hunter_catalog", help="ExpansionHunter variant catalog in JSON format")
    args = p.parse_args()

    if not args.output_file:
        args.output_file = re.sub(".json(.gz)?$", "", args.expansion_hunter_catalog) + ".bed.gz"

    process_expansion_hunter_catalog(args.expansion_hunter_catalog, args.output_file, args.split_adjacent_repeats,
                                     show_progress_bar=args.show_progress_bar,
                                     motif_size_column=args.motif_size_column)


def process_expansion_hunter_catalog(expansion_hunter_catalog_path, output_file_path, split_adjacent_repeats=False,
                                     show_progress_bar=False, motif_size_column=False):
    print(f"Parsing {expansion_hunter_catalog_path}")
    fopen = gzip.open if expansion_hunter_catalog_path.endswith((".gz", ".bgz")) else open
    with fopen(expansion_hunter_catalog_path, "rt") as f:
        iterator = ijson.items(f, "item", use_float=True)
        if show_progress_bar:
            iterator = tqdm.tqdm(iterator, unit=" variant catalog records", unit_scale=True)
        with (gzip.open if output_file_path.endswith((".gz", ".bgz")) else open)(output_file_path, "wt") as f2:
            previous_chrom = None
            output_rows = []
            counter = 0
            for i, record in enumerate(iterator):
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

                if locus_start_0based >= locus_end_1based:
                    print(f"WARNING: Skipping locus {locus_id} because its ReferenceRegion "
                          f"{chrom}:{locus_start_0based}-{locus_end_1based} has a width = "
                          f"{locus_end_1based - locus_start_0based}bp")
                    continue

                if "|" in locus_structure:
                    print(f"WARNING: Skipping locus {locus_id} @ {chrom}:{locus_start_0based}-{locus_end_1based} because "
                          f"its LocusStructure {locus_structure} contains a sequence swap operation '|' which is not "
                          f"supported by TRGT.")
                    continue

                if split_adjacent_repeats:
                    for motif, reference_region in zip(motifs, reference_regions):
                        chrom, start_0based, end_1based = parse_interval(reference_region)

                        output_rows.append([
                            chrom,
                            start_0based,
                            end_1based,
                            motif,
                            len(motif) if motif_size_column else ".",
                        ])
                else:
                    motif_string = ",".join(motifs)

                    output_rows.append([
                        chrom,
                        locus_start_0based,
                        locus_end_1based,
                        motif_string,
                        # A row spanning adjacent repeats has no single motif length, so the
                        # placeholder stands even when a motif size column was asked for.
                        ".",
                    ])

                if chrom != previous_chrom:
                    for output_row in sorted(output_rows):
                        counter += 1
                        f2.write("\t".join(map(str, output_row)) + "\n")

                    output_rows = []
                    previous_chrom = chrom

            for output_row in sorted(output_rows):
                counter += 1
                f2.write("\t".join(map(str, output_row)) + "\n")

    bgzip_step = "| bgzip" if output_file_path.endswith((".gz", ".bgz")) else ""
    os.system(f"bedtools sort -i {output_file_path} {bgzip_step} > {output_file_path}.sorted")
    os.system(f"mv {output_file_path}.sorted {output_file_path}")
    print(f"Wrote {counter:,d} rows to {output_file_path}")
    if output_file_path.endswith((".gz", ".bgz")):
        os.system(f"tabix -f {output_file_path}")
        print(f"Added {output_file_path}.tbi index")


if __name__ == "__main__":
    main()
