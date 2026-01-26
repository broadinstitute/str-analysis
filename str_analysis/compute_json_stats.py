#!/usr/bin/env python3
"""Compute statistics on fields in a JSON file containing a list of dictionaries.

Streams through the JSON file using ijson to handle large files efficiently,
counting how many times each field has a non-empty/non-null value.
"""

import argparse
import gzip
import ijson
import os
from collections import Counter
from tqdm import tqdm


class ProgressFileWrapper:
    """Wrapper around a file object that updates a tqdm progress bar as bytes are read."""

    def __init__(self, fh, progress_bar):
        self.fh = fh
        self.progress_bar = progress_bar

    def read(self, size=-1):
        data = self.fh.read(size)
        self.progress_bar.update(len(data))
        return data


def is_non_empty(value):
    """Check if a value is non-empty and non-null."""
    if value is None:
        return False
    if isinstance(value, str) and value.strip() == "":
        return False
    if isinstance(value, (list, dict)) and len(value) == 0:
        return False
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Compute statistics on fields in a JSON file containing a list of dictionaries."
    )
    parser.add_argument(
        "json_path",
        help="Path to JSON or JSON.gz file containing a list of dictionaries"
    )
    args = parser.parse_args()

    field_counts = Counter()
    total_rows = 0

    file_size = os.path.getsize(args.json_path)

    with tqdm(total=file_size, unit="B", unit_scale=True, desc="Reading") as progress_bar:
        # Open file, handling gzip if needed
        if args.json_path.endswith(".gz"):
            raw_fh = open(args.json_path, "rb")
            wrapped_fh = ProgressFileWrapper(raw_fh, progress_bar)
            fh = gzip.GzipFile(fileobj=wrapped_fh)
        else:
            raw_fh = open(args.json_path, "rb")
            fh = ProgressFileWrapper(raw_fh, progress_bar)

        # Stream through each item in the top-level array
        for record in ijson.items(fh, "item"):
            total_rows += 1
            for field, value in record.items():
                if is_non_empty(value):
                    field_counts[field] += 1

        raw_fh.close()

    # Print stats sorted by field name
    for field in sorted(field_counts.keys()):
        count = field_counts[field]
        percentage = 100.0 * count / total_rows if total_rows > 0 else 0
        print(f"{count:8,d} out of {total_rows:8,d} ({percentage:6.1f}%) entries have {field}")


if __name__ == "__main__":
    main()
