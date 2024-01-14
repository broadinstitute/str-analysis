"""This script takes an ExpansionHunter variant catalog as well as the ExpansionHunter output json, vcf, and bam file
for a given sample, and then runs REViewer on loci whose genotype exceeded a certain locus-specific threshold.
"""

import argparse
import collections
import gzip
import simplejson as json
import os
from pprint import pformat
import re

from str_analysis.utils.file_utils import file_exists


def parse_args():
	parser = argparse.ArgumentParser(description="Annotate and filter an STR variant catalog.")
	parser.add_argument("--samtools-path", default="samtools", help="Path of samtools executable")
	parser.add_argument("--reviewer-path", default="REViewer", help="Path of REViewer executable")
	parser.add_argument("-r", "--reference-fasta", required=True, help="Reference fasta file")
	parser.add_argument("-c", "--catalog", required=True, help="Variant catalog in JSON format")
	parser.add_argument("--region-extension-length", type=int, default=1000, help="Length of region extracted "
					    "around each locus. Must match the setting used when running ExpansionHunter")
	parser.add_argument("-j", "--json", required=True, help="ExpansionHunter output json file")
	parser.add_argument("-b", "--bam", required=True, help="ExpansionHunter output bam file")
	parser.add_argument("-v", "--vcf", required=True, help="ExpansionHunter output vcf file")
	parser.add_argument("-o", "--output-prefix", help="Output prefix. If not specified, the input json filename will be used.")
	parser.add_argument("--verbose", action="store_true", help="Print verbose output")

	parser.add_argument("--is-above-threshold", choices=[
		"short-allele-CI-lower-bound", "short-allele", "short-allele-CI-upper-bound",
		"long-allele-CI-lower-bound", "long-allele", "long-allele-CI-lower-bound"], required=True,
		help="Which part of the ExpansionHunter genotype must be greater than or equal to the threshold in order to "
			 "run REViewer? The choices are the allele size point estimate or the lower or upper bounds of the "
			 "confidence interval")

	g = parser.add_mutually_exclusive_group(required=True)
	g.add_argument("--use-threshold", type=int, help="Specify a single threshold on the command line to be used for all loci")
	g.add_argument("--use-custom-threshold-field", help="Specify the name of the variant catalog field that contains the threshold for each locus")
	g.add_argument("--use-pathogenic-threshold", action="store_true", help="Parse thresholds from the variant catalog's "
				   "Diseases > PathogenicMin fields for each locus. See https://github.com/broadinstitute/str-analysis/blob/main/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json "
				   "for an example catalog where these fields are defined.")

	parser.add_argument("--run-reviewer-if-threshold-not-available", help="By default, loci will be skipped if a "
						"threshold is not provided in the variant catalog. If this flag is specified, the "
					    "script will instead run REViewer for such loci.")
	args = parser.parse_args()

	for file_path in [args.reference_fasta, args.catalog, args.json, args.bam, args.vcf]:
		if not file_exists(os.path.expanduser(file_path)):
			parser.error(f"File not found: {file_path}")

	if not args.output_prefix:
		args.output_prefix = re.sub(".json(.gz)?$", "", os.path.basename(args.json))

	return args, parser


def run(c):
	print(c)
	os.system(c)

def compute_loci_to_process(args):
	loci_to_process = []

	# populate the threshold_lookup and main_reference_region_lookup from fields in the provided variant catalog
	with open(args.catalog) as f:
		variant_catalog = json.load(f)

	# parse threshold for each locus from variant catalog. This code expects the extra fields
	# "MainReferenceRegion" and "Diseases" to be defined in the catalog for each locus as in
	# https://github.com/broadinstitute/str-analysis/blob/main/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json
	threshold_lookup = {}
	main_reference_region_lookup = {}
	for record in variant_catalog:
		if "LocusId" not in record:
			print(f"WARNING: 'LocusId' key not found in variant catalog record for {pformat(record)}")
			continue

		locus_id = record["LocusId"]
		if "MainReferenceRegion" in record:
			# see https://github.com/broadinstitute/str-analysis/blob/main/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json
			main_reference_region_lookup[locus_id] = record["MainReferenceRegion"]

		if args.use_threshold:
			# no need to parse the threshold from the variant catalog since it was explicitly set on the command-line
			continue
		elif args.use_custom_threshold_field:
			try:
				threshold_lookup[locus_id] = int(record[args.use_custom_threshold_field])
			except Exception as e:
				print(f"WARNING: Unable to parse '{args.use_custom_threshold_field}' field for {locus_id}: {e}")
				continue
		elif args.use_pathogenic_threshold:
			try:
				pathogenic_min_thresholds = [
					int(disease_record["PathogenicMin"]) for disease_record in record.get("Diseases", []) if "PathogenicMin" in disease_record
				]

				if not pathogenic_min_thresholds:
					print(f"WARNING: Did not find any 'Diseases' records with a 'PathogenicMin' key for {locus_id}")
					continue

				threshold_lookup[locus_id] = min(pathogenic_min_thresholds)
			except Exception as e:
				print(f"WARNING: Unable to parse pathogenic threshold for {locus_id}: {e}")
				continue
		else:
			raise ValueError("Either --use-threshold or --use-pathogenic-threshold or --use-custom-threshold-field must be specified")

	# parse ExpansionHunter output json and determine which loci have genotypes above the threshold
	with open(args.json) as f:
		expansion_hunter_output_json = json.load(f)

	if "LocusResults" not in expansion_hunter_output_json:
		raise ValueError(f"'LocusResults' key not found in {args.json}")

	for locus_id, locus_results in expansion_hunter_output_json["LocusResults"].items():
		if args.use_threshold:
			threshold = args.use_threshold
		else:
			# check if the variant catalog had a threshold for this locus
			if locus_id not in threshold_lookup:
				if args.run_reviewer_if_threshold_not_available:
					loci_to_process.append(locus_id)
					continue
				else:
					if args.verbose:
						print(f"Skipping {locus_id} because no threshold was provided in the variant catalog")
					continue

			threshold = threshold_lookup[locus_id]

		# get the main variant results dictionary for this locus
		if len(locus_results.get("Variants", [])) == 1:
			variant_results = list(locus_results["Variants"].values())[0]

		elif len(locus_results.get("Variants", [])) > 1:
			if main_reference_region_lookup.get(locus_id) is None or not any(v.get("ReferenceRegion") == main_reference_region_lookup.get(locus_id) for v in locus_results["Variants"].values()):
				if args.verbose:
					print(f"WARNING: Skipping {locus_id} because it consists of {len(locus_results['Variants'])} adjacent repeats, and it's not clear which one to apply the threshold to.")
				continue
			for variant_id, v in locus_results["Variants"].items():
				if v.get("ReferenceRegion") == main_reference_region_lookup.get(locus_id):
					variant_results = v
					break
		else:
			print(f"WARNING: 'Variants' key not found for {locus_id} in {args.json}")
			continue

		# get the genotype
		if "Genotype" not in variant_results:
			print(f"WARNING: 'Genotype' key not found for {variant_id} in {args.json}")
			continue
		genotype = variant_results["Genotype"]
		genotype = [int(x) for x in genotype.split("/")]

		if "GenotypeConfidenceInterval" not in variant_results:
			print(f"WARNING: 'GenotypeConfidenceInterval' key not found for {variant_id} in {args.json}")
			continue
		genotype_confidence_interval = variant_results["GenotypeConfidenceInterval"]
		genotype_confidence_interval = [[int(boundary) for boundary in ci.split("-")] for ci in genotype_confidence_interval.split("/")]

		# Compare the genotype to the threshold. If it's a haploid genotype (ie. male on chrX), then apply the threshold
		# to the 1 allele regardless of whether the short or long allele was specified via --is-above-threshold
		if args.is_above_threshold == "short-allele-CI-lower-bound":
			if genotype_confidence_interval[0][0] >= threshold:
				loci_to_process.append(locus_id)
		elif args.is_above_threshold == "short-allele":
			if genotype[0] >= threshold:
				loci_to_process.append(locus_id)
		elif args.is_above_threshold == "short-allele-CI-upper-bound":
			if genotype_confidence_interval[0][1] >= threshold:
				loci_to_process.append(locus_id)
		elif args.is_above_threshold == "long-allele-CI-lower-bound":
			if genotype_confidence_interval[-1][0] >= threshold:
				loci_to_process.append(locus_id)
		elif args.is_above_threshold == "long-allele":
			if genotype[-1] >= threshold:
				loci_to_process.append(locus_id)
		elif args.is_above_threshold == "long-allele-CI-upper-bound":
			if genotype_confidence_interval[-1][1] >= threshold:
				loci_to_process.append(locus_id)
		else:
			raise ValueError(f"Unexpected value for --is-above-threshold: {args.is_above_threshold}")

	return loci_to_process


def main():
	args, parser = parse_args()

	if args.verbose:
		print("Args:")
		for key, value in sorted(args.__dict__.items(), key=lambda x: str(x[1])):
			if value is None:
				continue
			key += " = "
			if isinstance(value, list):
				print(f"   {key:30s}")
				for v in value:
					print(f"       {v}")
			else:
				print(f"   {key:30s} {str(value)}")

	# parse ExpansionHunter output json
	loci_to_process = compute_loci_to_process(args)

	if not loci_to_process:
		print("No loci need REViewer to be run since all genotypes are below thresholds. Exiting...")
		return


	print(f"Running REViewer for {len(loci_to_process)} loci")
	run(f"{args.samtools_path} sort {args.bam} -o {args.output_prefix}.sorted.bam")
	run(f"{args.samtools_path} index {args.output_prefix}.sorted.bam")

	for locus in loci_to_process:
		print("-"*20)
		print(f"Running REViewer for {locus}")
		run(" ".join([
			args.reviewer_path, 
			f"--reads {args.output_prefix}.sorted.bam",
			f"--vcf {args.vcf}",
			f"--reference {args.reference_fasta}",
			f"--catalog {args.catalog}",
			f"--locus {locus}",
			f"--region-extension-length {args.region_extension_length}",
			f"--output-prefix {args.output_prefix}",
		]))
	print("Done")

if __name__ == "__main__":
	main()

