#!/usr/bin/env python3

"""Combines json files into a single .tsv table for analysis."""


import argparse
import json
import logging
import os
import pandas as pd
import pathlib
import sys

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")


def parse_args(args_list=None):
    p = argparse.ArgumentParser()
    p.add_argument(
        "-m",
        "--sample-metadata",
        help="Table of sample annotations. If specified, all columns from this table will be added to the output table."
             " Both the json file(s) and this table should have a 'sample_id' field or column. Variations of this are"
             "also accepted - including 'SampleId', 'sample', and 'ParticipantId'",
    )
    p.add_argument(
        "-o",
        "--output-prefix",
        help="Combined table output filename prefix",
    )
    p.add_argument(
        "json_paths",
        help="EpxansionHunter output json path(s). If not specified, this script will retrieve all json files in the current directory and subdirectories",
        type=pathlib.Path,
        nargs="*"
    )
    args = p.parse_args(args=args_list)

    if not args.json_paths:
        args.json_paths = [p for p in pathlib.Path(".").glob("**/*.json")]
        print(f"Found {len(args.json_paths)} .json files under {os.getcwd()}")
        if len(args.json_paths) == 0:
            sys.exit(0)

    return args


SAMPLE_ID_COLUMN_ALAISES = {
    "sampleid",
    "sample_id",
    "sample",
    "participantid",
}


def get_sample_id_column_index(df):
    for i, column in enumerate(df.columns):
        if column.lower() in SAMPLE_ID_COLUMN_ALAISES:
            return i
    return -1


class ParseError(Exception):
    pass


def parse_json_files(json_paths):
    """Takes json file paths and yields the contents of each one as a dictionary or list"""
    for json_path in json_paths:
        if not os.path.isfile(json_path):
            raise ValueError(f"{json_path} not found")

        with open(json_path, "rt", encoding="UTF-8") as f:
            try:
                yield json.load(f)
            except Exception as e:
                raise ParseError(f"Unable to parse {json_path}: {e}")


def main():
    args = parse_args()

    output_prefix = args.output_prefix or f"combined.{len(args.json_paths)}_json_files"

    df = pd.DataFrame(parse_json_files(args.json_paths))

    if args.sample_metadata:
        sample_id_column_idx1 = get_sample_id_column_index(df)
        if sample_id_column_idx1 is -1:
            raise ValueError(f"'sample_id' field not found in json files. The fields found were: {df.columns}")

        metadata_df = pd.read_table(args.sample_metadata)

        sample_id_column_idx2 = get_sample_id_column_index(metadata_df)
        if sample_id_column_idx1 is -1:
            raise ValueError(f"'sample_id' column not found in {args.sample_metadata}. The columns found were: {metadata_df.columns}")

        df = pd.merge(df, metadata_df, how="left", left_on=df.columns[sample_id_column_idx1], right_on=metadata_df.columns[sample_id_column_idx2])

    output_filename = f"{output_prefix}.tsv"
    df.to_csv(output_filename, index=False, header=True, sep="\t")
    print(f"Wrote {len(df)} rows to {output_filename}")


if __name__ == "__main__":
    main()

