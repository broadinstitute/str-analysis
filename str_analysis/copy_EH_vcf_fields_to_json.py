"""This script copies useful fields from the ExpansionHunter output VCF to the ExpansionHunter output json since
some info (like the 'LowDepth' filter) are only available in the VCF and not in the json.
"""

import argparse
import json
import os


def main():
	p = argparse.ArgumentParser(description="Copy fields from the ExpansionHunter output VCF to the ExpansionHunter json")
	p.add_argument("-o", "--output-json", help="Output json path. If not specified, the input json file will be updated in place")
	p.add_argument("expansion_hunter_output_vcf", help="ExpansionHunter output VCF path")
	p.add_argument("expansion_hunter_output_json", help="ExpansionHunter output json path")
	args = p.parse_args()

	for path in args.expansion_hunter_output_vcf, args.expansion_hunter_output_json:
		if not os.path.isfile(path):
			p.error(f"{path} file not found")

	if not args.output_json:
		args.output_json = args.expansion_hunter_output_json

	with open(args.expansion_hunter_output_json, "rt") as f:
		expansion_hunter_output_json = json.load(f)

	vcf_row_counter = 0
	vcf_row_lookup = {}
	with open(args.expansion_hunter_output_vcf, "rt") as f:
		for line in f:
			if line.startswith("#"):
				continue
			vcf_row_counter += 1
			fields = line.strip().split("\t")

			filter_field = fields[6]
			info_field = fields[7]
			info_fields = info_field.split(";")

			# example: END=16327723;REF=30;RL=90;RU=TGC;VARID=ATXN1;REPID=ATXN1
			info_fields_dict = dict([key_value.split("=") for key_value in info_fields if "=" in key_value])
			if "VARID" not in info_fields_dict:
				raise ValueError(f"VARID key not found in info field of line: {line.strip()}")

			var_id = info_fields_dict["VARID"]
			if var_id in vcf_row_lookup:
				raise ValueError(f"Duplicate variant id found in VCF: {line.strip()}")

			genotype_format = fields[8].split(":")
			genotype_values = fields[9].split(":")
			genotype_dict = dict(zip(genotype_format, genotype_values))

			vcf_row_lookup[var_id] = (filter_field, info_fields_dict, genotype_dict)

	print(f"Parsed {vcf_row_counter:,d} rows from {args.expansion_hunter_output_vcf}")

	variant_counter = 0
	all_locus_results = expansion_hunter_output_json.get("LocusResults", [])
	for locus_id, locus_results in all_locus_results.items():
		for variant_id, variant_results in locus_results["Variants"].items():
			variant_counter += 1
			if variant_id not in vcf_row_lookup:
				raise ValueError(f"Variant id {variant_id} not found in VCF")

			filter_field, info_fields_dict, genotype_dict = vcf_row_lookup[variant_id]

			del vcf_row_lookup[variant_id]

			variant_results["Filter"] = filter_field
			for key in "SO", "ADSP", "ADFL", "ADIR", "LC":
				variant_results[key] = genotype_dict[key]

	if len(vcf_row_lookup) > 0:
		raise ValueError(f"VCF variant id(s) not found in the ExpansionHunter json: {', '.join(sorted(vcf_row_lookup.keys()))}")

	with open(args.output_json, "wt") as f:
		json.dump(expansion_hunter_output_json, f, indent=4)

	print(f"Wrote {variant_counter:,d} variants at {len(all_locus_results):,d} loci to {args.output_json}")


if __name__ == "__main__":
	main()