"""This script converts an ExpansionHunter variant catalog to a GangSTR spec. Variant catalog entries
that include adjacent repeats (eg. with LocusStructure like "(A)*(ACG)*") are split into multiple GangSTR
specs - one per repeat.
"""

import argparse
import gzip
import ijson
import re
import tqdm

from str_analysis.utils.misc_utils import parse_interval



def main():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--verbose", action="store_true", help="Print verbose output")
    p.add_argument("-o", "--output-file", help="bed file output path")
    p.add_argument("variant_catalog")
    args = p.parse_args()

    if not args.output_file:
        args.output_file = re.sub(".json(.gz)?$", "", args.variant_catalog) + ".gangstr_spec.bed"

    process_variant_catalog(args.variant_catalog, args.output_file, verbose=args.verbose)


def process_variant_catalog(variant_catalog_path, output_file_path, verbose=False):
    print(f"Parsing {variant_catalog_path}")
    with (gzip.open if variant_catalog_path.endswith("gz") else open)(variant_catalog_path, "rt") as f:
        with (gzip.open if output_file_path.endswith("gz") else open)(output_file_path, "wt") as f2:
            skipped_locus = 0
            for record in tqdm.tqdm(ijson.items(f, "item"), unit=" variant catalog records"):
                locus_structure = record["LocusStructure"]
                motifs = re.findall("[(]([A-Z]+)[)]", locus_structure)
                if not motifs:
                    raise ValueError(f"Unable to parse LocusStructure '{locus_structure}' in variant catalog record: {record}")

                reference_regions = record["ReferenceRegion"]
                if not isinstance(reference_regions, list):
                    reference_regions = [reference_regions]

                variant_types = record.get("VariantType", "Repeat")
                if not isinstance(variant_types, list):
                    variant_types = [variant_types]

                if len(motifs) != len(reference_regions):
                    raise ValueError(f"len(motifs) != len(reference_regions) in variant catalog record: {record}")

                if len(motifs) != len(variant_types):
                    raise ValueError(f"len(motifs) != len(variant_types) in variant catalog record: {record}")

                offtarget_regions = record.get("OfftargetRegions", [])
                for motif, variant_type, reference_region in zip(motifs, variant_types, reference_regions):
                    chrom, start_0based, end_1based = parse_interval(reference_region)

                    if start_0based + 1 >= end_1based:
                        print(f"WARNING: Skipping {motif} locus @ {chrom}:{start_0based+1}-{end_1based} because "
                              f"the interval has a width = {end_1based - start_0based - 1}bp")
                        continue

                    if len(motif) > 9 or (end_1based - start_0based) <= 1:
                        # HipSTR doesn't support motifs longer than 9bp or loci where start_1based == stop_1based
                        # (https://github.com/tfwillems/HipSTR/blob/master/src/region.cpp#L33-L35)
                        if verbose and len(motif) > 9:
                            skipped_locus += 1
                            print(f"Skipping locus #{skipped_locus} with motif size > 9bp:  {chrom}-{start_0based}-{end_1based}-{motif}")

                        continue

                    output_fields = [
                        chrom, start_0based + 1, end_1based, len(motif), int((end_1based - start_0based)/len(motif)),
                        f"{chrom}-{start_0based}-{end_1based}-{motif}"
                    ]
                    f2.write("\t".join(map(str, [
                        chrom,
                        start_0based + 1,  # GangSTR beds are 1-based
                        end_1based,
                        len(motif),
                        motif,
                        ",".join(offtarget_regions) if variant_type == "RareRepeat" else "",
                    ])) + "\n")

    print(f"Wrote out {output_file_path}")


if __name__ == "__main__":
    main()
