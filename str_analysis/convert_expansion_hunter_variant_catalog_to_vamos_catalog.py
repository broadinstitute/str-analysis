"""This script converts an ExpansionHunter variant catalog to a BED file."""

import argparse
import gzip
import ijson
import os
import re
import tqdm

from str_analysis.utils.misc_utils import parse_interval

COLUMNS = [
    "chrom",
    "start_1based",
    "end",
    "motifs",
    "version",
    "STR_or_VNTR",
    "motif_size",
    "score1",
    "score2",
    "segdup",
    "coding",
]
"""
TODO: look up locus definitions in the official VAMOS catalog
df = pd.read_table(
    f"https://zenodo.org/records/11625069/files/vamos.motif.hg38.{version}.e0.1.tsv.gz?download=1",
    compression="gzip",
    names=COLUMNS,
)
"""

def main():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    p.add_argument("-o", "--output-file", help="TSV file output path")
    p.add_argument("expansion_hunter_catalog", help="ExpansionHunter variant catalog in JSON format")
    args = p.parse_args()

    if not args.output_file:
        args.output_file = re.sub(".json(.gz)?$", "", args.expansion_hunter_catalog) + ".tsv.gz"

    process_expansion_hunter_catalog(args.expansion_hunter_catalog, args.output_file, show_progress_bar=args.show_progress_bar)


def process_expansion_hunter_catalog(expansion_hunter_catalog_path, output_file_path, show_progress_bar=False):
    print(f"Parsing {expansion_hunter_catalog_path}")
    output_rows = []
    counter = 0

    fopen = gzip.open if expansion_hunter_catalog_path.endswith("gz") else open
    with fopen(expansion_hunter_catalog_path, "rt") as f:
        iterator = ijson.items(f, "item")
        if show_progress_bar:
            iterator = tqdm.tqdm(iterator, unit=" variant catalog records", unit_scale=True)

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

            if locus_start_0based + 1 >= locus_end_1based:
                print(f"WARNING: Skipping locus {locus_id} because its ReferenceRegion "
                      f"{chrom}:{locus_start_0based+1}-{locus_end_1based} has a width = "
                      f"{locus_end_1based - locus_start_0based - 1}bp")
                continue

            if "|" in locus_structure:
                print(f"WARNING: Skipping locus {locus_id} @ {chrom}:{locus_start_0based+1}-{locus_end_1based} because "
                      f"its LocusStructure {locus_structure} contains a sequence swap operation '|' which is not "
                      f"supported by TRGT.")
                continue

            for motif, reference_region in zip(motifs, reference_regions):
                chrom, start_0based, end_1based = parse_interval(reference_region)

                output_rows.append({
                    "chrom": chrom,
                    "start_1based":  start_0based + 1,
                    "end": end_1based,
                    "motifs": motif,
                    "version": ".",
                    "STR_or_VNTR": "STR" if len(motif) <= 6 else "VNTR",
                    "motif_size": len(motif),
                    "score1": 0,  # not sure what these scores mean
                    "score2": 0,
                    "segdup": "",  # should be set to "segDup" or empty string
                    "coding": "",  # should be set to "coding" or empty string
                })

    df = pd.DataFrame(output_rows)
    df = df[COLUMNS]
    df.to_csv(output_file_path, sep="\t", index=False, header=False)
    print(f"Wrote {counter:,d} rows to {output_file_path}")


if __name__ == "__main__":
    main()
