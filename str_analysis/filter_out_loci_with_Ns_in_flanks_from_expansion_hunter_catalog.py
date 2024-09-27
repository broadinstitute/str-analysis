"""This script takes an ExpansionHunter variant catalog and filters out loci that would trigger this ExpansionHunter
error: "[Error loading locus 1-10001-10108-TAACCC: Flanks can contain at most 5 characters N but found 1000 Ns]"
"""

import argparse
import gzip
import json
import os
import pyfaidx
import re
import tqdm

from str_analysis.utils.misc_utils import parse_interval
from utils.eh_catalog_utils import get_variant_catalog_iterator

# based on https://github.com/Illumina/ExpansionHunter/blob/master/ehunter/io/LocusSpecDecoding.cpp#L70-L79
MAX_N_IN_FLANKS = 5

def main():
	p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	p.add_argument("-o", "--output-file", help="JSON catalog output file path")
	p.add_argument("-R", "--reference-fasta-path", help="Reference genome fasta path.", required=True)
	p.add_argument("--region-extension-length", type=int, default=1000, help="Length of flanking regions")
	p.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
	p.add_argument("-f", "--output-list-of-filtered-loci", help="Output a list of filtered loci to this TXT file path")
	p.add_argument("-v", "--verbose", action="store_true")
	p.add_argument("variant_catalog", help="Path of the ExpansionHunter variant catalog in JSON format")
	args = p.parse_args()

	if not os.path.isfile(args.reference_fasta_path):
		p.error(f"{args.reference_fasta_path} file not found")

	if not args.output_file:
		args.output_file = re.sub(".json(.gz)?$", "", args.variant_catalog) + ".without_loci_with_flanking_Ns.json.gz"

	if args.output_list_of_filtered_loci and not args.output_list_of_filtered_loci.endswith(".txt"):
		p.error(f"Output list of filtered loci must end with .txt")

	fasta_obj = pyfaidx.Fasta(
		args.reference_fasta_path, one_based_attributes=False, as_raw=True, sequence_always_upper=True)

	num_flanking_bases = args.region_extension_length

	print(f"Parsing {args.variant_catalog}")
	fopen = gzip.open if args.output_file.endswith("gz") else open
	output_file = fopen(args.output_file, "wt")

	output_list_of_filtered_loci_file = None
	if args.output_list_of_filtered_loci:
		output_list_of_filtered_loci_file = open(args.output_list_of_filtered_loci, "wt")

	output_file.write("[")
	output_record_counter = 0
	filtered_loci_counter = 0
	for record in get_variant_catalog_iterator(args.variant_catalog, show_progress_bar=args.show_progress_bar):
		if not "ReferenceRegion" in record:
			p.error(f"ReferenceRegion not found in variant catalog record #{i+1}: {record}")

		if isinstance(record["ReferenceRegion"], list):
			reference_regions = record["ReferenceRegion"]
			reference_region_chrom, reference_region_start_0based, _ = parse_interval(reference_regions[0])
			_, _, reference_region_end = parse_interval(reference_regions[-1])
		else:
			reference_region_chrom, reference_region_start_0based, reference_region_end = parse_interval(
				record["ReferenceRegion"])

		left_flank_start = max(reference_region_start_0based - num_flanking_bases, 0)
		right_flank_end = min(reference_region_end + num_flanking_bases, len(fasta_obj[reference_region_chrom]))

		left_flank = fasta_obj[reference_region_chrom][left_flank_start:reference_region_start_0based]
		right_flank = fasta_obj[reference_region_chrom][reference_region_end:right_flank_end]

		total_Ns_in_flanks = left_flank.count("N") + right_flank.count("N")
		if total_Ns_in_flanks > MAX_N_IN_FLANKS:
			locus_id = record["LocusId"]
			if args.verbose:
				print(f"Skipping locus {locus_id} because it has {total_Ns_in_flanks} Ns in the flanks")
			if output_list_of_filtered_loci_file:
				output_list_of_filtered_loci_file.write(f"{locus_id}\n")
			filtered_loci_counter += 1
			continue

		if output_record_counter > 0:
			output_file.write(", ")
		output_file.write(json.dumps(record, indent=4))
		output_record_counter += 1
	output_file.write("]")
	output_file.close()
	print(f"Wrote {output_record_counter:,d} records to {args.output_file}")
	
	if output_list_of_filtered_loci_file:
		print(f"Wrote {filtered_loci_counter:,d} loci to {args.output_list_of_filtered_loci}")
		output_list_of_filtered_loci_file.close()


if __name__ == "__main__":
	main()

