import argparse
import collections
import gzip
import simplejson as json
import logging
import os
import pandas as pd
from pprint import pformat
import sys

from str_analysis.combine_str_json_to_tsv import convert_expansion_hunter_json_to_tsv_columns, parse_json_file

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

def parse_args(args_list=None):
    """Parse command line args and return the argparse args object"""
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--discard-hom-ref", action="store_true", help="Discard hom-ref calls")
    p.add_argument("--output-format", choices=("JSON", "TSV"), default="JSON", help="Output format")
    p.add_argument("-o", "--output-prefix", help="Combined table output filename prefix")
    p.add_argument("-v", "--verbose", action="store_true", help="Print additional logging messages")
    p.add_argument("json_path", help="EpxansionHunter output json file.")

    args = p.parse_args(args=args_list)

    if not os.path.isfile(args.json_path):
        p.error(f"File not found: {args.json_path}")

    return args, p


def main():
    args, parser = parse_args()

    if not args.output_prefix:
        args.output_prefix = os.path.basename(args.json_path).replace(".json", "").replace(".gz", "")

    if args.discard_hom_ref:
        args.output_prefix += ".without_hom_ref"
    else:
        args.output_prefix += ".postprocessed"

    if args.output_format == "JSON":
        output_filename = f"{args.output_prefix}.json.gz"
    elif args.output_format == "TSV":
        output_filename = f"{args.output_prefix}.tsv.gz"
    else:
        parser.error(f"Unknown output format: {args.output_format}")

    try:
        json_contents = parse_json_file(args.json_path)
    except Exception as e:
        parser.error(f"Unable to parse json: {e}")

    if "LocusResults" not in json_contents:
        parser.error(f"No 'LocusResults' key found in {args.json_path}")

    total_loci = len(json_contents["LocusResults"])

    output_records = []
    for variant_record in convert_expansion_hunter_json_to_tsv_columns(
        json_contents,
        json_file_path=args.json_path,
        yield_allele_records=False,
        include_extra_expansion_hunter_fields=False,
        discard_hom_ref=args.discard_hom_ref,
    ):
        output_records.append(variant_record)

    if args.output_format == "TSV":
        pd.DataFrame(output_records).to_csv(output_filename, sep="\t", index=False)
    elif args.output_format == "JSON":
        with gzip.open(output_filename, "wt") as output_file:
            json.dump(output_records, output_file, indent=4)

    logging.info(f"Wrote {len(output_records):,d} out of {total_loci:,d} records "
                 f"({len(output_records) / total_loci:.1%}) to {output_filename}")


if __name__ == "__main__":
    main()
