"""Compute average Q per sample"""

import argparse
import collections
import os
import pandas as pd

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("-o", "--output-prefix", help="Output TSV path")
	parser.add_argument("--sample-id-column", default="SampleId", help="Column name for sample ID")
	parser.add_argument("input_table", help="Input table with genotypes to filter")

	args = parser.parse_args()

	if not args.input_table.endswith(".tsv") and not args.input_table.endswith(".tsv.gz"):
		parser.error("Input table must have a .tsv or .tsv.gz extension")

	if not os.path.isfile(args.input_table):
		parser.error(f"Input table not found: {args.input_table}")

	if not args.output_prefix:
		args.output_prefix = args.input_table.replace(".tsv", "").replace(".gz", "")

	return args


def main():
	args = parse_args()
	print(f"Parsing {args.input_table}")

	output_path = f"{args.output_prefix}.average_Q_per_sample.tsv"
	df = pd.read_table(args.input_table)

	output_records = []
	for sample_id, sample_df in df.groupby(args.sample_id_column):
		output_records.append({
			"SampleId": sample_id,
			"AverageQ": (sample_df["Q: Allele 1"].mean() + sample_df["Q: Allele 2"].mean())/2,
		})

	output_df = pd.DataFrame(output_records)
	output_df.to_csv(output_path, sep="\t", index=False)
	print(f"Wrote {len(output_df):,d} rows to {output_path}")

if __name__ == "__main__":
	main()