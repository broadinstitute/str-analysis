# This script checks concordance between the variant catalogs in this repo and disease-associated TR reference sites:
# https://stripy.org/database/ and  https://harrietdashnow.com/STRchive/

import argparse
import json
import os
import pandas as pd
from pprint import pformat
import pysam
import re
import requests

from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.combine_str_catalogs import parse_motifs_from_locus_structure

OTHER_CATALOG_NAME_STRCHIVE = "STRchive"
OTHER_CATALOG_NAME_STRIPY = "STRipy"


def get_json_from_url(url):
	"""Download and parse json from a url"""
	return requests.get(url, headers={'Cache-Control': 'no-cache'}).json()

def get_official_expansion_hunter_catalog_dict():
	"""Download the official ExpansionHunter catalog from the Illumina github repo

	Return:
		dict: locus_id -> locus_info
	"""
	print("Downloading official ExpansionHunter catalog from the Illumina GitHub repo")

	official_EH_catalog_data = get_json_from_url(
		"https://raw.githubusercontent.com/Illumina/ExpansionHunter/master/variant_catalog/grch38/variant_catalog.json"
	)

	official_EH_catalog_loci = {
		locus["LocusId"]: locus for locus in official_EH_catalog_data
	}

	return official_EH_catalog_loci


def get_stripy_dict():
	"""Download and parse relevant info from the STRipy website https://stripy.org/database

	Return:
		dict: locus_id -> locus_info
	"""

	print("Downloading loci from stripy.org")
	tables = pd.read_html("https://stripy.org/database/")

	assert len(tables) >= 1
	stripy_df = tables[0]
	assert "Locus" in stripy_df.columns

	stripy_lookup = {}
	for _, row in stripy_df.iterrows():
		info = row.to_dict()
		locus_id = info["Locus"]
		print(f"Downloading {locus_id} info from stripy.org")
		for key in "Locus", "Normal", "Intermediate", "Pathogenic", "Inheritance", "Disease":
			if key in info:
				del info[key]

		tables = pd.read_html(f"http://stripy.org/database/{locus_id}")

		locus_df = tables[0]
		locus_df.columns = ["Name", "Value"]
		for _, locus_df_row in locus_df.iterrows():
			value = locus_df_row["Value"]
			if locus_df_row["Name"] == "Reference region":
				for reference_region in value.split("  "):
					match = re.search("[(]h.*[)]", reference_region)
					if not match or match.group(0) not in ("(hg19)", "(hg38)", "(hs1)"):
						if locus_id != "XYLT1":
							raise ValueError(f"Unable to parse reference regions: {reference_region}")
						else:
							continue

					reference_genome = match.group(0).strip("() ")
					info[f"ReferenceRegion_{reference_genome}"] = reference_region.split(" ")[0]
			elif locus_df_row["Name"] == "Repeat unit":
				for motif in value.split(" "):
					info["Motif"] = motif
					break
			else:
				info[locus_df_row["Name"]] = value

		# replace weird chars
		for key in "Region", "Location":
			info[key] = info[key].replace("â\x80²", "'")

		info["CanonicalMotif"] = compute_canonical_motif(info["Motif"])

		disease_df = tables[1].astype(str)
		diseases = []
		for _, disease_df_row in disease_df.iterrows():
			normal_max = None
			if not disease_df_row["Normal"] or disease_df_row["Normal"].strip() == "-":
				pass
			elif "-" in disease_df_row["Normal"]:
				try:
					normal_max = int(disease_df_row["Normal"].split("-")[1])
				except ValueError:
					print(f"WARNING: Unable to parse normal max: {disease_df_row['Normal']}")
					input("?")
			else:
				try:
					normal_max = int(disease_df_row["Normal"].replace("&leq;", ""))
				except ValueError:
					print(f"WARNING: Unable to parse normal max: {disease_df_row['Normal']}")
					input("?")

			intermediate_range = disease_df_row["Intermediate"]
			if intermediate_range.strip() == "-":
				intermediate_range = None

			pathogenic_min = None
			if not disease_df_row["Pathogenic"] or disease_df_row["Pathogenic"].strip() == "-":
				pass
			elif "-" in disease_df_row["Pathogenic"]:
				try:
					pathogenic_min = int(disease_df_row["Pathogenic"].split("-")[0])
				except ValueError:
					print(f"WARNING: Unable to parse normal max: {disease_df_row['Pathogenic']}")
					input("?")

			else:
				try:
					pathogenic_min = int(disease_df_row["Pathogenic"].replace("&GreaterEqual;", ""))
				except ValueError:
					print(f"WARNING: Unable to parse normal max: {disease_df_row['Pathogenic']}")
					input("?")

			#print(disease_df_row["Normal"], normal_max)
			#print(disease_df_row["Intermediate"], intermediate_range)
			#print(disease_df_row["Pathogenic"], pathogenic_min)

			disease_name = disease_df_row["Disease"].replace("\n", " ")
			disease_symbol_match = re.search("[(](.*?)[)]", disease_df_row["Disease"])
			diseases.append({
				"Name": disease_name,
				"Symbol": disease_symbol_match.group(1) if disease_symbol_match else None,
				"Inheritance": disease_df_row["Inheritance"],
				"Onset": disease_df_row["Onset"],
				"NormalMax": normal_max,
				"IntermediateRange": intermediate_range,
				"PathogenicMin": pathogenic_min,
			})

		info["Diseases"] = diseases
		stripy_lookup[locus_id] = info

	#from pprint import pprint
	#pprint(info)
	#print("--")

	return stripy_lookup


def get_strchive_dict():
	"""Download the locus metadata json from the STRchive github repo which underlies https://harrietdashnow.com/STRchive

	Return:
		dict: locus_id -> locus_info
	"""
	print("Downloading loci from STRchive")
	strchive_data = get_json_from_url(
		"https://raw.githubusercontent.com/hdashnow/STRchive/main/data/STR-disease-loci.processed.json"
	)

	strchive_id_to_gnomad_map = {
		"OPDM_ABCD3": "ABCD3",
		"FRA2A_AFF3": "AFF3",
		"JBS_CBL": "CBL",
		"SCA27B_FGF14": "FGF14",
		#"CPEO_POLG": "POLG",
		"OPDM4_RILPL1": "RILPL1",
		"SCA_THAP11": "THAP11",
		"SCA4_ZFHX3": "ZFHX3",
		"FRA7A_ZNF713": "ZNF713",
	}

	previously_seen_gnomad_genes = set()
	for d in strchive_data:
		if not d["gnomAD_gene"] or not d["gnomAD_gene"].strip():
			if d["id"] in strchive_id_to_gnomad_map:
				d["gnomAD_gene"] = strchive_id_to_gnomad_map[d["id"]]
			else:
				print(f"WARNING: STRchive gnomAD_gene field not set for {d['id']}: {d['chrom']}:{d['start_hg38']}-{d['stop_hg38']}")
		if d["gnomAD_gene"] in previously_seen_gnomad_genes:
			print(f"WARNING: Duplicate gnomAD_gene field: {d['gnomAD_gene']}")
		previously_seen_gnomad_genes.add(d["gnomAD_gene"])

	strchive_lookup = {
		r["gnomAD_gene"]: r for r in strchive_data
	}

	for d in strchive_lookup.values():
		d["CanonicalMotif"] = compute_canonical_motif(d["reference_motif_reference_orientation"])


	return strchive_lookup


def get_gnomad_catalog():
	"""Download the gnomAD catalog from the str-analysis github repo"""

	print("Downloading the gnomAD catalog from the str-analysis github repo")
	gnomad_catalog = get_json_from_url(
		"https://raw.githubusercontent.com/broadinstitute/str-analysis/main/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json"
	)

	# check internal consistency
	for d in gnomad_catalog:
		if isinstance(d["ReferenceRegion"], list):
			motifs = parse_motifs_from_locus_structure(d["LocusStructure"])
			assert len(motifs) == len(d["ReferenceRegion"]), pformat(d)
			assert len(motifs) == len(d["VariantType"]), pformat(d)
			assert len(motifs) == len(d["VariantId"]), pformat(d)
		else:
			if d["LocusId"] != "RFC1":
				assert d["LocusStructure"] == "(" + d["RepeatUnit"] + ")*", pformat(d)
			assert d["ReferenceRegion"] == d["MainReferenceRegion"], pformat(d)

	return gnomad_catalog


def output(file, line):
	"""Print line to screen and write it to the file"""
	print(line)
	file.write(line + "\n")


def compare_catalogs(args, official_EH_catalog_loci, gnomad_catalog, stripy_lookup, strchive_lookup, compare_with=OTHER_CATALOG_NAME_STRIPY):
	"""Compare the gnomAD catalog with the other catalog"""

	fasta_file = pysam.FastaFile(os.path.expanduser(args.reference_fasta))

	output_path = f"gnomAD_{compare_with}_comparison.txt"
	output_file = open(output_path, "wt")
	print(f"Writing to {output_path}")

	output_stripy_variant_catalog = []
	output_strchive_variant_catalog = []

	gnomad_locus_ids = {d["LocusId"] for d in gnomad_catalog}
	stripy_locus_ids = set(stripy_lookup.keys())
	strchive_locus_ids = set(strchive_lookup.keys())
	other_catalog_name = compare_with
	other_catalog_locus_ids = strchive_locus_ids if compare_with == OTHER_CATALOG_NAME_STRCHIVE else stripy_locus_ids
	other_catalog_lookup = strchive_lookup if compare_with == OTHER_CATALOG_NAME_STRCHIVE else stripy_lookup

	# check which loci are only in one of the catalogs
	gnomad_records = [d for d in gnomad_catalog if d["LocusId"] in other_catalog_lookup]
	for locus_id in other_catalog_locus_ids - gnomad_locus_ids:
		output(output_file, f"{locus_id} is in {other_catalog_name} but not in gnomAD")
	output(output_file, "------")
	for locus_id in gnomad_locus_ids - other_catalog_locus_ids:
		output(output_file, f"{locus_id} is in gnomAD but not in {other_catalog_name}")

	output(output_file, "------")

	# compare motifs for loci that are in both catalogs
	for d in gnomad_records:
		locus_id = d["LocusId"]
		canonical_motif = compute_canonical_motif(d["RepeatUnit"])
		if canonical_motif != other_catalog_lookup[locus_id]["CanonicalMotif"] and len(canonical_motif) != 5:
			output(output_file, f"Canonical motifs for {locus_id:<6s} differ between gnomAD ({canonical_motif}) and {other_catalog_name} ({other_catalog_lookup[locus_id]['CanonicalMotif']})")

	output(output_file, "------")


	# compare other fields across loci that are in both catalogs
	counter = 0
	for d in gnomad_records:
		locus_id = d["LocusId"]
		motif_size = len(d["RepeatUnit"])
		gnomad_reference_region = d["MainReferenceRegion"]
		gnomad_chrom, gnomad_start, gnomad_end = parse_interval(gnomad_reference_region)
		gnomad_repeats = (gnomad_end - gnomad_start)/motif_size
		gnomad_motif = d["RepeatUnit"]
		gnomad_canonical_motif = compute_canonical_motif(gnomad_motif)
		gnomad_pathogenic_min = ", ".join(map(str, [disease_info.get("PathogenicMin") for disease_info in d.get("Diseases", [])]))
		gnomad_ref = fasta_file.fetch(gnomad_chrom, gnomad_start, gnomad_end)

		stripy_reference_region = None
		gnomad_differs_from_stripy = False
		if locus_id in stripy_lookup:
			stripy_info = stripy_lookup[locus_id]
			stripy_reference_region = stripy_info["ReferenceRegion_hg38"]
			stripy_chrom, stripy_start, stripy_end = parse_interval(stripy_info["ReferenceRegion_hg38"])
			stripy_repeats = (stripy_end - stripy_start)/motif_size
			stripy_motif = stripy_info["Motif"]
			stripy_canonical_motif = compute_canonical_motif(stripy_motif)
			stripy_pathogenic_min = ", ".join(map(str, [disease_info.get("PathogenicMin") for disease_info in stripy_info.get("Diseases", [])]))
			stripy_ref = fasta_file.fetch(stripy_chrom, stripy_start, stripy_end)

			gnomad_vs_stripy_distance = max(abs(gnomad_start - stripy_start), abs(gnomad_end - stripy_end))
			gnomad_differs_from_stripy = gnomad_reference_region != stripy_reference_region and (gnomad_vs_stripy_distance >= motif_size or gnomad_repeats != stripy_repeats)
			gnomad_differs_from_stripy = gnomad_differs_from_stripy or (len(gnomad_motif) != 5 and gnomad_canonical_motif != stripy_canonical_motif)

		strchive_reference_region = None
		gnomad_differs_from_strchive = False
		if locus_id in strchive_lookup:
			strchive_info = strchive_lookup[locus_id]
			strchive_reference_region = strchive_info["chrom"] + ":" + strchive_info["start_hg38"] + "-" + strchive_info["stop_hg38"]
			strchive_chrom, strchive_start, strchive_end = parse_interval(strchive_reference_region)
			strchive_repeats = (strchive_end - strchive_start)/motif_size
			strchive_pathogenic_min = strchive_info.get("pathogenic_min")
			strchive_motif = strchive_info["reference_motif_reference_orientation"]
			strchive_canonical_motif = compute_canonical_motif(strchive_motif)
			strchive_ref = fasta_file.fetch(strchive_chrom, strchive_start, strchive_end)

			gnomad_vs_strchive_distance = max(abs(gnomad_start - strchive_start), abs(gnomad_end - strchive_end))
			gnomad_differs_from_strchive = gnomad_reference_region != strchive_reference_region and (gnomad_vs_strchive_distance >= motif_size or gnomad_repeats != strchive_repeats)
			gnomad_differs_from_strchive = gnomad_differs_from_strchive or (len(gnomad_motif) != 5 and gnomad_canonical_motif != strchive_canonical_motif)

		if (gnomad_differs_from_stripy and compare_with == OTHER_CATALOG_NAME_STRIPY) or (gnomad_differs_from_strchive and compare_with == OTHER_CATALOG_NAME_STRCHIVE):
			counter += 1
			output(output_file, "------")
			whats_different = "ReferenceRegion" if (gnomad_reference_region != stripy_reference_region and compare_with == OTHER_CATALOG_NAME_STRIPY) or (gnomad_reference_region != strchive_reference_region and compare_with == OTHER_CATALOG_NAME_STRCHIVE) else "motif"
			output(output_file, f"{locus_id} hg38 {whats_different} isn't the same between gnomAD and {other_catalog_name}")

			if locus_id in official_EH_catalog_loci and not isinstance(official_EH_catalog_loci[locus_id]['ReferenceRegion'], list):
				official_EH_reference_region = official_EH_catalog_loci[locus_id]['ReferenceRegion']
				_, official_EH_start, official_EH_end = parse_interval(official_EH_reference_region)
				official_EH_motif = official_EH_catalog_loci[locus_id]['LocusStructure'].strip("()*+")
				official_EH_repeats = (official_EH_end - official_EH_start)/motif_size
				output(output_file, "%3s %s" % (" ", f"chr{official_EH_reference_region}  ({official_EH_repeats:4.1f} x {official_EH_motif}  {official_EH_end - official_EH_start:5d}bp) official catalog from EH repo"))

			min_left_pos = min(gnomad_start, stripy_start, strchive_start) if locus_id in strchive_lookup else min(gnomad_start, stripy_start)

			output(output_file, f"%3s %-100s%{gnomad_end-min_left_pos}s" % (
				" ", f"{gnomad_reference_region}  ({gnomad_repeats:4.1f} x {gnomad_motif}  {gnomad_end - gnomad_start:5d}bp) gnomAD    {gnomad_pathogenic_min} = pathogenic.min.  hg38 seq:", gnomad_ref))

			if locus_id in stripy_lookup:
				output(output_file, f"%3s %-100s%{stripy_end-min_left_pos}s" % (
					" ", f"{stripy_info['ReferenceRegion_hg38']}  ({stripy_repeats:4.1f} x {stripy_motif}  {stripy_end - stripy_start:5d}bp) STRipy    {stripy_pathogenic_min} = pathogenic.min.  hg38 seq:", stripy_ref))

				if gnomad_differs_from_stripy:
					output_stripy_variant_catalog.append({
						"LocusId": f"{locus_id}_STRIPY",
						"LocusStructure": "(" + stripy_info["Motif"] + ")*",
						"ReferenceRegion": stripy_info["ReferenceRegion_hg38"],
						"VariantType": "Repeat",
					})

			if locus_id in strchive_lookup:
				output(output_file, f"%3s %-100s%{strchive_end-min_left_pos}s" % (
					" ", f"{strchive_reference_region}  ({strchive_repeats:4.1f} x {strchive_motif}  {strchive_end - strchive_start:5d}bp) STRchive  {strchive_pathogenic_min} = pathogenic.min.  hg38 seq:", strchive_ref))

				if gnomad_differs_from_strchive:
					output_strchive_variant_catalog.append({
						"LocusId": f"{locus_id}_STRCHIVE",
						"LocusStructure": "(" + strchive_info["reference_motif_reference_orientation"] + ")*",
						"ReferenceRegion": strchive_info["chrom"] + ":" + strchive_info["start_hg38"] + "-" + strchive_info["stop_hg38"],
						"VariantType": "Repeat",
					})


	output(output_file, f"\n{counter:,d} out of {len(gnomad_records):,d} loci have a different hg38 ReferenceRegion or motif in gnomAD and {other_catalog_name}")

	output_file.close()

	if output_stripy_variant_catalog and compare_with == OTHER_CATALOG_NAME_STRIPY:
		with open("STRipy_definitions.json", "wt") as f:
			json.dump(output_stripy_variant_catalog, f, indent=4)
			print(f"Wrote {len(output_stripy_variant_catalog):,d} loci to STRipy_definitions.json")
	if output_strchive_variant_catalog and compare_with == OTHER_CATALOG_NAME_STRCHIVE:
		with open("STRchive_definitions.json", "wt") as f:
			json.dump(output_strchive_variant_catalog, f, indent=4)
			print(f"Wrote {len(output_strchive_variant_catalog):,d} loci to STRchive_definitions.json")


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-R", "--reference-fasta", help="hg38 reference fasta file path", default="~/hg38.fa")
	args = parser.parse_args()

	args.reference_fasta = os.path.expanduser(args.reference_fasta)

	if not os.path.isfile(args.reference_fasta):
		parser.error(f"Reference fasta file not found: {args.reference_fasta}")

	official_EH_catalog_loci = get_official_expansion_hunter_catalog_dict()
	gnomad_catalog = get_gnomad_catalog()
	strchive_lookup = get_strchive_dict()
	stripy_lookup = get_stripy_dict()

	for compare_with in OTHER_CATALOG_NAME_STRIPY, OTHER_CATALOG_NAME_STRCHIVE:
		print("="*100)
		print(f"Comparing gnomAD with {compare_with}")
		compare_catalogs(args, official_EH_catalog_loci, gnomad_catalog, stripy_lookup, strchive_lookup, compare_with=compare_with)


if __name__ == "__main__":
	main()
