"""This script takes a combined table of TR genotypes (generated by the combine_str_json_to_tsv script) and
filters it to the N most expanded genotypes per locus, based on either the long allele, or on the short allele
(for recessive searches). The output is a table with the same columns as the input, but with only the N most expanded
genotypes per locus.
The input table must have a "LocusId" column and a "Num Repeats: Allele 1" and "Num Repeats: Allele 2" column. If
--min-purity is specified, it should also have an "AllelePurity" column with comma-separated purity values for each
allele.
"""

import argparse
import collections
import os
import pandas as pd

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("-n", type=int, default=50, help="Number of genotypes to keep per locus")
	parser.add_argument("-o", "--output-prefix", help="Output TSV path")
	parser.add_argument("--min-purity", type=float, help="Optionally apply a purity threshold before filtering")
	group = parser.add_mutually_exclusive_group()
	group.add_argument("--by-long-allele", action="store_true", help="Filter by long allele")
	group.add_argument("--by-short-allele", action="store_true", help="Filter by short allele")
	parser.add_argument("input_table", help="Input table with genotypes to filter")

	args = parser.parse_args()

	if not args.input_table.endswith(".tsv") and not args.input_table.endswith(".tsv.gz"):
		parser.error("Input table must have a .tsv or .tsv.gz extension")

	if not os.path.isfile(args.input_table):
		parser.error(f"Input table not found: {args.input_table}")

	if not args.output_prefix:
		args.output_prefix = args.input_table.replace(".tsv", "").replace(".gz", "")

	return args

def filter_by_purity(min_purity_threshold, both_alleles=False):
	def filter_func(allele_purity_values):
		if both_alleles:
			return all(p != "." and float(p) >= min_purity_threshold for p in allele_purity_values.split(","))
		else:
			allele_purity_values = allele_purity_values.split(",")
			p = allele_purity_values[-1]
			return p != "." and float(p) >= min_purity_threshold

	return filter_func


def convert_allele_histogram_dict_to_string(allele_histogram_dict):
    data = sorted(allele_histogram_dict.items())
    return ",".join(f"{repeat_number}x:{allele_count}" for repeat_number, allele_count in data)


def get_stdev_of_allele_histogram_dict(allele_histogram_dict):
    total = sum(allele_histogram_dict.values())
    mean = sum(repeat_number * count for repeat_number, count in allele_histogram_dict.items()) / total
    return (sum((repeat_number - mean) ** 2 * count for repeat_number, count in allele_histogram_dict.items()) / total) ** 0.5


def add_allele_histogram(df, histogram_key, stdev_key):
	locus_id_to_histogram_dict = collections.defaultdict(collections.Counter)
	for locus_id, locus_df in df.groupby("LocusId"):
		for c in ["Num Repeats: Allele 1", "Num Repeats: Allele 2"]:
			for allele_size in locus_df[c]:
				if pd.isna(allele_size):
					continue
				locus_id_to_histogram_dict[locus_id][int(allele_size)] += 1

	for locus_id, histogram_dict in locus_id_to_histogram_dict.items():
		locus_id_selector = df["LocusId"] == locus_id
		df.loc[locus_id_selector, histogram_key] = convert_allele_histogram_dict_to_string(histogram_dict)
		df.loc[locus_id_selector, stdev_key] = get_stdev_of_allele_histogram_dict(histogram_dict)


def process_table(df, args, by_long_allele=True):
	total = len(df)

	if args.min_purity:
		add_allele_histogram(df, "AlleleHistBeforePurityFilter", "AlleleStdevBeforePurityFilter")
		if by_long_allele:
			print(f"Filtering to genotypes with long allele purity of at least {args.min_purity}")
			df = df[df["AllelePurity"].apply(filter_by_purity(args.min_purity, both_alleles=False))]
		else:
			print(f"Filtering to genotypes where both alleles have purity of at least {args.min_purity}")
			df = df[df["AllelePurity"].apply(filter_by_purity(args.min_purity, both_alleles=True))]
		print(f"Kept {len(df):,d} out of {total:,d} ({len(df) / total:.2%}) rows after filtering by purity")
		add_allele_histogram(df, "AlleleHist", "AlleleStdev")
	else:
		add_allele_histogram(df, "AlleleHist", "AlleleStdev")

	print(f"Sorting by {'long' if by_long_allele else 'short'} allele")
	if by_long_allele:
		df.sort_values(["Num Repeats: Allele 2", "Num Repeats: Allele 1"], ascending=False, inplace=True)
	else:
		df.sort_values(["Num Repeats: Allele 1", "Num Repeats: Allele 2"], ascending=False, inplace=True)

	print(f"Filtering to top {args.n} genotypes per locus")
	df = df.groupby("LocusId").head(args.n)
	print(f"Kept {len(df):,d} out of {total:,d} ({len(df) / total:.2%}) rows after filtering to top {args.n} genotypes per locus")

	output_path = args.output_prefix
	if args.min_purity:
		output_path += f".purity_{args.min_purity}"
	output_path += f".top_{args.n}_by_{'long' if by_long_allele else 'short'}_allele.tsv.gz"
	df.to_csv(output_path, sep="\t", index=False)
	print(f"Wrote {len(df):,d} rows to {output_path}")


def main():
	args = parse_args()
	print(f"Parsing {args.input_table}")
	df = pd.read_table(args.input_table)
	print(f"Parsed {len(df):,d} rows from {args.input_table}")

	if args.by_long_allele:
		by_long_allele_settings = [False]
	elif args.by_short_allele:
		by_long_allele_settings = [True]
	else:
		by_long_allele_settings = [True, False]

	for by_long_allele in by_long_allele_settings:
		if len(by_long_allele_settings) > 1:
			print(f"="*80)
			print(f"Processing {args.input_table} with sorting by {'long' if by_long_allele else 'short'} allele")
		process_table(df, args, by_long_allele)

if __name__ == "__main__":
	main()