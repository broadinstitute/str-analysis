"""This script takes an ExpansionHunter variant catalog and splits any locus definitions that contain adjacent repeats
 (eg. with LocusStructure like "(A)*(ACG)*"). It outputs each locus as a separate record.
"""

import argparse
import gzip
import ijson
import json
import os
import re
import tqdm

from str_analysis.utils.eh_catalog_utils import get_variant_catalog_iterator
from str_analysis.utils.misc_utils import parse_interval


def main():
	p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	p.add_argument("-o", "--output-file", help="JSON file output path")
	p.add_argument("expansion_hunter_catalog", help="ExpansionHunter variant catalog in JSON format")
	args = p.parse_args()

	if not args.output_file:
		args.output_file = re.sub(".json(.gz)?$", "", os.path.basename(args.expansion_hunter_catalog))
		args.output_file += ".split.json"
		if args.expansion_hunter_catalog.endswith("gz"):
			args.output_file += ".gz"

	process_expansion_hunter_catalog(args.expansion_hunter_catalog, args.output_file)


def process_expansion_hunter_catalog(expansion_hunter_catalog_path, output_file_path):
	print(f"Parsing {expansion_hunter_catalog_path}")
	input_file_iterator = get_variant_catalog_iterator(expansion_hunter_catalog_path)

	input_records_split_counter = 0
	output_records_split_counter = 0
	output_records = []
	for i, record in enumerate(tqdm.tqdm(input_file_iterator, unit=" variant catalog records", unit_scale=True)):
		locus_id = record["LocusId"]
		locus_structure = record["LocusStructure"]
		if "|" in locus_structure:
			print(f"WARNING: Skipping locus {locus_id} @ {chrom}:{locus_start_0based+1}-{locus_end_1based} because "
				  f"its LocusStructure {locus_structure} contains a sequence swap operation '|'")
			continue

		motifs = re.findall("[(]([A-Z]+)[)]", locus_structure)
		if not motifs:
			raise ValueError(f"Unable to parse LocusStructure '{locus_structure}' in variant catalog "
							 f"record #{i+1}: {record}")

		reference_regions = record["ReferenceRegion"]
		if not isinstance(reference_regions, list):
			reference_regions = [reference_regions]

		if len(motifs) != len(reference_regions):
			raise ValueError(f"LocusStructure elements != # of entries in the list of ReferenceRegions in "
							 f"variant catalog record #{i+1}: {record}")

		if len(reference_regions) == 1:
			output_records.append(record)
			continue

		input_records_split_counter += 1
		variant_ids = record.get("VariantId")
		if not variant_ids:
			raise ValueError(f"VariantId field missing in record: {record}")

		if len(variant_ids) != len(reference_regions):
			raise ValueError(f"Number of VariantIds ({len(variant_ids)}) != number of reference regions "
							 f"({len(reference_regions)}) in record: {record}")

		variant_types = record.get("VariantType")
		if not variant_types:
			raise ValueError(f"VariantType field missing in record: {record}")
		if len(variant_types) != len(reference_regions):
			raise ValueError(f"Number of VariantTypes ({len(variant_types)}) != number of reference regions "
							 f"({len(reference_regions)}) in record: {record}")

		for variant_id, variant_type, motif, reference_region in zip(variant_ids, variant_types, motifs, reference_regions):
			chrom, start_0based, end_1based = parse_interval(reference_region)

			output_record = record.copy()
			output_record["LocusId"] = variant_id
			del output_record["VariantId"]
			output_record["VariantType"] = variant_type
			output_record["ReferenceRegion"] = reference_region
			output_record["LocusStructure"] = f"({motif})*"
			output_records.append(output_record)
			output_records_split_counter += 1

	print(f"Split {input_records_split_counter:,d} loci into {output_records_split_counter:,d} output records")
	fopen = gzip.open if output_file_path.endswith("gz") else open
	with fopen(output_file_path, "wt") as f:
		json.dump(output_records, f, indent=2)

	print(f"Wrote {len(output_records):,d} total records to {output_file_path}")


if __name__ == "__main__":
	main()
