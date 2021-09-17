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
        "-d",
        "--add-dirname-column",
        action="store_true",
        help="Add Dirname column containing the relative path of directory containing the json file."
    )
    p.add_argument(
        "-f",
        "--add-filename-column",
        action="store_true",
        help="Add Filename column containing the json filename."
    )
    p.add_argument(
        "-m",
        "--sample-metadata",
        help="Table of sample annotations. If specified, all columns from this table will be added to the output table."
    )
    p.add_argument(
        "--json-sample-id-key",
        help="The json field that contains a sample id to use for joining with the --sample-metadata table. "
             "If not specified, 'sample_id' and other variations like 'SampleId', 'sample', and 'ParticipantId' "
             "will be tried."
    )
    p.add_argument(
        "--sample-metadata-key",
        help="The column name in the --sample-metdata table that contains a sample id. "
             "If not specified, 'sample_id' and other variations like 'SampleId', 'sample', and 'ParticipantId' "
             "will be tried."
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

    if args.json_paths:
        # check if files exist
        for json_path in args.json_paths:
            if not os.path.isfile(json_path):
                p.error(f"File not found: {json_path}")
    else:
        # find all .json files underneath the current directory
        args.json_paths = [p for p in pathlib.Path(".").glob("**/*.json")]
        print(f"Found {len(args.json_paths)} .json files under {os.getcwd()}")
        if len(args.json_paths) == 0:
            sys.exit(0)

    if args.sample_metadata_key and not args.sample_metadata:
        p.error("--sample-metadata-key should only be specified along with --sample-metadata")

    if args.json_sample_id_key and not args.sample_metadata:
        p.error("--json-sample-id-key should only be specified along with --sample-metadata")

    return args


SAMPLE_ID_COLUMN_ALAISES = {
    "sampleid",
    "sample_id",
    "sample",
    "participantid",
}


def get_sample_id_column_index(df, column_name=None):
    """Try to find a column in df that contains sample ids

    Args:
         df (DataFrame): pandas DataFrame
         column_name (str): if specified this function will get the index of this specific column name

    Return:
        Index of sample id column, or -1 if not found
    """
    for i, column in enumerate(df.columns):
        if column_name is not None:
            if column == column_name:
                return i
        elif column.lower() in SAMPLE_ID_COLUMN_ALAISES:
            return i

    return -1


class ParseError(Exception):
    pass


def parse_json_files(json_paths, add_dirname_column=False, add_filename_column=False):
    """Takes json file paths and yields the contents of each one as a dictionary or list"""
    for json_path in json_paths:
        if not os.path.isfile(json_path):
            raise ValueError(f"{json_path} not found")

        with open(json_path, "rt", encoding="UTF-8") as f:
            try:
                json_contents = json.load(f)
            except Exception as e:
                raise ParseError(f"Unable to parse {json_path}: {e}")

            if isinstance(json_contents, dict):
                if add_dirname_column:
                    json_contents["Dirname"] = os.path.dirname(json_path)
                if add_filename_column:
                    json_contents["Filename"] = os.path.basename(json_path)

                json_contents_excluding_complex_values = {
                    key: value for key, value in json_contents.items() if isinstance(value, (int, str, bool, float, tuple))
                }

                yield json_contents_excluding_complex_values


def main():
    args = parse_args()

    output_prefix = args.output_prefix or f"combined.{len(args.json_paths)}_json_files"

    df = pd.DataFrame(parse_json_files(
        args.json_paths,
        add_dirname_column=args.add_dirname_column,
        add_filename_column=args.add_filename_column))

    if args.sample_metadata:

        sample_id_column_idx1 = get_sample_id_column_index(df, column_name=args.json_sample_id_key)
        if sample_id_column_idx1 is -1:
            raise ValueError(f"'sample_id' field not found in json files. The fields found were: {df.columns}")

        metadata_df = pd.read_table(args.sample_metadata)

        sample_id_column_idx2 = get_sample_id_column_index(metadata_df, column_name=args.sample_metadata_key)
        if sample_id_column_idx2 is -1:
            raise ValueError(f"'sample_id' column not found in {args.sample_metadata}. The columns found were: {metadata_df.columns}")

        print(f"Doing a LEFT JOIN with {args.sample_metadata} using keys "
              f"{df.columns[sample_id_column_idx1]} and {metadata_df.columns[sample_id_column_idx2]}")
        df = pd.merge(df, metadata_df, how="left", left_on=df.columns[sample_id_column_idx1], right_on=metadata_df.columns[sample_id_column_idx2])

    output_filename = f"{output_prefix}.tsv"
    df.to_csv(output_filename, index=False, header=True, sep="\t")
    print(f"Wrote {len(df)} rows to {output_filename}")


if __name__ == "__main__":
    main()

