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
    p.add_argument("-o", "--output-file", help="bed file output path")
    p.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    p.add_argument("variant_catalog")
    args = p.parse_args()

    if not args.output_file:
        args.output_file = re.sub(".json(.gz)?$", "", args.variant_catalog) + ".gangstr_spec.bed"

    process_variant_catalog(args.variant_catalog, args.output_file, show_progress_bar=args.show_progress_bar)


def process_variant_catalog(variant_catalog_path, output_file_path, show_progress_bar=False):
    print(f"Parsing {variant_catalog_path}")

    fopen = gzip.open if variant_catalog_path.endswith("gz") else open
    with fopen(variant_catalog_path, "rt") as f:
        iterator = ijson.items(f, "item")
        if show_progress_bar:
            iterator = tqdm.tqdm(iterator, unit=" variant catalog records", unit_scale=True)

        with (gzip.open if output_file_path.endswith("gz") else open)(output_file_path, "wt") as f2:

            for record in iterator:
                locus_id = record["LocusId"]
                locus_structure = record["LocusStructure"]
                repeat_units = re.findall("[(]([A-Z]+)[)]", locus_structure)
                if not repeat_units:
                    raise ValueError(f"Unable to parse LocusStructure '{locus_structure}' in variant catalog record: {record}")

                reference_regions = record["ReferenceRegion"]
                if not isinstance(reference_regions, list):
                    reference_regions = [reference_regions]

                variant_types = record.get("VariantType", "Repeat")
                if not isinstance(variant_types, list):
                    variant_types = [variant_types]

                if len(repeat_units) != len(reference_regions):
                    raise ValueError(f"len(repeat_units) != len(reference_regions) in variant catalog record: {record}")

                if len(repeat_units) != len(variant_types):
                    raise ValueError(f"len(repeat_units) != len(variant_types) in variant catalog record: {record}")

                offtarget_regions = record.get("OfftargetRegions", [])
                for repeat_unit, variant_type, reference_region in zip(repeat_units, variant_types, reference_regions):
                    chrom, start_0based, end_1based = parse_interval(reference_region)
                    if start_0based + 1 >= end_1based:
                        print(f"WARNING: Skipping locus {locus_id} @ {chrom}:{start_0based+1}-{end_1based} because "
                              f"the interval has a width = {end_1based - start_0based - 1}bp")
                        continue

                    f2.write("\t".join(map(str, [
                        chrom,
                        start_0based + 1,  # GangSTR beds are 1-based
                        end_1based,
                        len(repeat_unit),
                        int((end_1based - start_0based)/len(repeat_unit)),
                        locus_id,
                    ])) + "\n")

    print(f"Wrote out {output_file_path}")


if __name__ == "__main__":
    main()
