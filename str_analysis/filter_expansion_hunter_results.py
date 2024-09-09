#!/usr/bin/env python3

import argparse
import collections
import gzip
import simplejson as json
import logging
import os
import pathlib
import re
import sys

import numpy as np
import pandas as pd
from pprint import pformat

from str_analysis.utils.misc_utils import parse_interval
from str_analysis.combine_str_json_to_tsv import convert_expansion_hunter_json_to_tsv_columns, parse_json_file
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
ALREADY_WARNED_ABOUT = set()  # used for logging


def parse_args(args_list=None):
    """Parse command line args and return the argparse args object"""
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--discard-hom-ref", action="store_true", help="Discard hom-ref calls")
    p.add_argument("-o", "--output-prefix", help="Combined table output filename prefix")
    p.add_argument("-v", "--verbose", action="store_true", help="Print additional logging messages")
    p.add_argument("json_path", help="EpxansionHunter output json file.")

    args = p.parse_args(args=args_list)

    if not os.path.isfile(args.json_path):
        p.error(f"File not found: {args.json_path}")

    return args, p


def main():
    """Main"""

    args, parser = parse_args()

    if not args.discard_hom_ref:
        print("No filters to apply. Exiting..")
        sys.exit(0)

    if not args.output_prefix:
        args.output_prefix = os.path.basename(args.json_path).replace(".json", "").replace(".gz", "")

    if args.discard_hom_ref:
        output_filename = f"{args.output_prefix}.without_hom_ref.json.gz"
    else:
        output_filename = f"{args.output_prefix}.filtered.json.gz"

    try:
        json_contents = parse_json_file(args.json_path)
    except Exception as e:
        parser.error(f"Unable to parse json: {e}")

    if "LocusResults" not in json_contents:
        parser.error(f"No 'LocusResults' key found in {args.json_path}")

    total_loci = len(json_contents["LocusResults"])

    variant_records_counter = 0
    with gzip.open(output_filename, "wt") as variant_output_file:
        variant_output_file.write("[")
        for variant_record in convert_expansion_hunter_json_to_tsv_columns(
            json_contents,
            json_file_path=args.json_path,
            yield_allele_records=False,
            include_extra_expansion_hunter_fields=False,
            discard_hom_ref=args.discard_hom_ref,
        ):
            if variant_records_counter > 0:
                variant_output_file.write(", ")
            variant_records_counter += 1
            json.dump(variant_record, variant_output_file, indent=4)
        variant_output_file.write("]")
    logging.info(f"Wrote {variant_records_counter:,d} out of {total_loci:,d} records "
                 f"({variant_records_counter / total_loci:.1%}) to {output_filename}")


if __name__ == "__main__":
    main()
