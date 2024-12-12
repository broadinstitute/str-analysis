"""This script converts a TRGT repeat catalog BED file to an ExpansionHunter variant catalog JSON file.
Currently, this script will skip any loci that contain interruptions in the reference repeat sequence.
"""

import argparse
import collections
import gzip
import simplejson as json
import os
from pprint import pformat
import pyfaidx
import re
import tqdm

from str_analysis.utils.eh_catalog_utils import parse_motifs_from_locus_structure
from str_analysis.utils.trgt_utils import convert_trgt_locus_to_expansion_hunter_format

def main():
	p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	p.add_argument("-r", "--reference-fasta", required=True, help="Reference FASTA file path")
	p.add_argument("-o", "--output-file", help="ExansionHunter variant catalog JSON output file path")
	p.add_argument("-v", "--verbose", action="store_true")
	p.add_argument("--set-locus-id", action="store_true", help="Ignore the ID field in the TRGT catalog and set the "
															   "locus ID to {chrom}-{start_0based}-{end}-{motif}")
	p.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
	p.add_argument("trgt_bed_path", help="path of the TRGT repeat catalog BED file")
	args = p.parse_args()

	if not args.output_file:
		args.output_file = re.sub(".bed(.gz)?$", "", args.trgt_bed_path) + ".variant_catalog.json"

	if not os.path.isfile(args.trgt_bed_path):
		p.error(f"{args.trgt_bed_path} file not found")

	process_variant_catalog(args.trgt_bed_path, args.reference_fasta, args.output_file,
							set_locus_id=args.set_locus_id,
							verbose=args.verbose,
							show_progress_bar=args.show_progress_bar)


def process_variant_catalog(trgt_bed_path_path, reference_fasta_path, output_file_path, set_locus_id=False,
							verbose=False, show_progress_bar=False):

	fasta_obj = pyfaidx.Fasta(reference_fasta_path, one_based_attributes=False, as_raw=True)

	print(f"Parsing {trgt_bed_path_path}")
	output_json_records = []
	counter = collections.defaultdict(int)
	fopen = gzip.open if trgt_bed_path_path.endswith("gz") else open
	with fopen(trgt_bed_path_path, "rt") as f:
		iterator = f
		if show_progress_bar:
			iterator = tqdm.tqdm(iterator, unit=" records", unit_scale=True)

		for i, row in enumerate(iterator):
			fields = row.strip("\n").split("\t")
			if len(fields) < 4:
				print(f"WARNING: skipping invalid row which has fewer than 4 columns: {fields}")
				continue
			chrom = fields[0]
			start_0based = int(fields[1])
			end_1based = int(fields[2])
			info_fields = fields[3]  # Example: ID=chr1_153777445_153777638;MOTIFS=GAAGAGGAAGTGAGA,GGGGACA;STRUC=GCAGGGAAGGACTGTGTCAGTCCAG(GAAGAGGAAGTGAGA)n(GGGGACA)nGCACGGCTTCCTGACGGTGTGACTC
			info_fields_dict = {}
			for key_value in info_fields.split(";"):
				key_value = key_value.split("=")
				if len(key_value) != 2:
					print(f"WARNING: skipping invalid key-value pair '{key_value}' in line {fields}")
					continue
				key, value = key_value
				info_fields_dict[key] = value

			locus_structure = info_fields_dict.get("STRUC", "")
			motifs_from_locus_structure = parse_motifs_from_locus_structure(locus_structure)
			motifs = info_fields_dict.get("MOTIFS")
			if not motifs:
				raise ValueError(f"MOTIFS field not found in row #{i + 1}: {fields}")
			motifs = motifs.split(",")
			if set(motifs) != set(motifs_from_locus_structure):
				print(f"WARNING: MOTIFS {tuple(motifs)} != motifs in locus structure {tuple(motifs_from_locus_structure)} in row #{i + 1}: {fields}")
				continue

			motifs = motifs_from_locus_structure
			for motif in motifs:
				if not motif or (len(set(motif) - set("ACGTN")) > 0):
					raise ValueError(f"Invalid repeat unit in row #{i + 1}: {motif}. Line: {row}")

			counter["total input loci"] += 1

			# get reference sequence
			if chrom not in fasta_obj:
				raise ValueError(f"Chromosome '{chrom}' not found in reference FASTA file which has chromosomes {list(fasta_obj.keys())}")

			if len(motifs) == 1:
				motif = motifs[0]
				if set_locus_id:
					locus_id = f"{chrom}-{start_0based}-{end_1based}-{motif}"
				else:
					locus_id = info_fields_dict.get("ID")
					if not locus_id:
						raise ValueError(f"ID field not found in row #{i + 1}: {fields}. One work-around is to rerun with --set-locus-id")

				#if verbose:
				#	print(f"Adding locus {locus_id} with motif {motif} and reference region {chrom}:{start_0based}-{end_1based}")

				output_json_records.append({
					"LocusId": locus_id,
					"ReferenceRegion": f"{chrom}:{start_0based}-{end_1based}",
					"LocusStructure": f"({motif})*",
					"VariantType": "Repeat",
				})

			else:
				if set_locus_id:
					locus_id = None
				else:
					locus_id = info_fields_dict.get("ID")
					if not locus_id:
						raise ValueError(
							f"ID field not found in row #{i + 1}: {fields}. One work-around is to rerun with --set-locus-id")

				output_records = list(convert_trgt_locus_to_expansion_hunter_format(
					fasta_obj, chrom, start_0based, end_1based, locus_structure, locus_id_prefix=locus_id, verbose=False))

				if len(output_records) > 0:
					output_json_records.extend(output_records)
				else:
					counter["skipped"] += 1
					#if verbose:
					#	print(
					#		f"WARNING: Skipped {info_fields_dict.get('ID')} due to issues with parsing the reference "
					#		f"sequence. This typically happens when reference repeat sequence contains interruptions. "
					#		f"TRGT catalog specified the motifs as: {', '.join(motifs)}")

	print(f"Parsed {counter['total input loci']:,d} loci from {trgt_bed_path_path}")
	print(f"Writing {len(output_json_records):,d} records to {output_file_path}")
	fopen = gzip.open if output_file_path.endswith("gz") else open
	with fopen(output_file_path, "wt") as f:
		json.dump(output_json_records, f, indent=4, ignore_nan=True)

	print("Stats:")
	for key, value in sorted(counter.items(), key=lambda x: x[1], reverse=True):
		print(f"{value:10,d} loci {key}")

	print(f"Wrote {len(output_json_records):,d} records to {output_file_path}")



if __name__ == "__main__":
	main()
