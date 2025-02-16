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
import sys

from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.merge_loci import parse_motifs_from_locus_structure

OTHER_CATALOG_NAME_STRCHIVE = "STRchive"
OTHER_CATALOG_NAME_STRIPY = "STRipy"
OTHER_CATALOG_NAME_TRGT = "TRGT"


def get_json_from_url(url):
	"""Download and parse json from a url"""
	return requests.get(url, headers={'Cache-Control': 'no-cache'}).json()


def get_text_from_url(url):
	"""Download and parse json from a url"""
	return requests.get(url, headers={'Cache-Control': 'no-cache'}).text



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


def get_strchive_dict(gnomad_catalog):
	"""Download the locus metadata json from the STRchive github repo which underlies https://harrietdashnow.com/STRchive

	Return:
		dict: locus_id -> locus_info
	"""

	gnomad_loci = {d["LocusId"] for d in gnomad_catalog}

	print("Downloading loci from STRchive")
	strchive_data = get_json_from_url(
		"https://raw.githubusercontent.com/dashnowlab/STRchive/refs/heads/main/data/STRchive-loci.json"
	)

	strchive_id_to_gnomad_map = {
		"pre-MIR7-2_CHNG3": "PRE-MIR7-2",
		"HOXA13_1": "HOXA13",
		"HOXA13_2": "HOXA13",
		"HOXA13_3": "HOXA13",
		"C9orf72": "C9ORF72",
		"ARX_1": "ARX",
		"ARX_2": "ARX",
	}

	previously_seen_gnomad_genes = set()
	for d in strchive_data:
		if d["id"] in strchive_id_to_gnomad_map:
			d["gnomad"] = strchive_id_to_gnomad_map[d["id"]]
		elif d["gene"] in gnomad_loci:
			d["gnomad"] = d["gene"]
		elif len(d["gnomad"]) > 0:
			d["gnomad"] = d["gnomad"][0]
		else:
			print(f"WARNING: STRchive locus is absent from gnomAD: {d['id']}")
			continue

		if d["gnomad"]:
			if d["gnomad"] in previously_seen_gnomad_genes:
				print(f"WARNING: Duplicate gnomad field: {d['gnomad']}")

			previously_seen_gnomad_genes.add(d["gnomad"])

	strchive_lookup = {
		r["gnomad"]: r for r in strchive_data if r.get("gnomad")
	}

	for d in strchive_lookup.values():
		assert len(d["reference_motif_reference_orientation"]) == 1, d
		d["reference_motif_reference_orientation"] = d["reference_motif_reference_orientation"][0]
		d["CanonicalMotif"] = compute_canonical_motif(d["reference_motif_reference_orientation"])

		d["Diseases"] = [{
			"Name": d["disease"],
			"Symbol": d["disease_id"],
			"Inheritance": ",".join(d["inheritance"]),
			"Onset": d["age_onset"],
			"NormalMax": d["benign_max"],
			"IntermediateRange": f"{d['intermediate_min']}-{d['intermediate_max']}",
			"PathogenicMin": d["pathogenic_min"],
		}]

	return strchive_lookup


def get_trgt_catalog():

	print("Downloading the TRGT catalog from the trgt github repo")
	bed_contents = get_text_from_url(
		"https://raw.githubusercontent.com/PacificBiosciences/trgt/main/repeats/pathogenic_repeats.hg38.bed"
	)

	trgt_lookup = {}
	for row in bed_contents.strip().split("\n"):
		fields = row.split("\t")
		chrom = fields[0]
		start_0based = int(fields[1])
		end = int(fields[2])

		info = dict(key_value.split("=") for key_value in fields[3].split(";"))
		locus_id = info["ID"]
		motifs = info["MOTIFS"].split(",")
		trgt_lookup[locus_id] = {
			"chrom": chrom,
			"start_0based": start_0based,
			"end": end,
			"ReferenceRegion": f"{chrom}:{start_0based}-{end}",
			"CanonicalMotif": compute_canonical_motif(motifs[0]),
			"RepeatUnit": motifs[0],
			"LocusStructure": info["STRUC"].replace("n", "*"),
			"Diseases": [],
		}

	return trgt_lookup


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

def compare_threshold(output_file, gnomad_records, other_catalog_lookup, other_catalog_name,
					  threshold_field="PathogenicMin", aggregator_func=min):

	for d in gnomad_records:
		locus_id = d["LocusId"]
		gnomad_threshold = list(filter(None, [
			disease_info.get(threshold_field, None) for disease_info in d.get("Diseases", [])
		]))
		try:
			gnomad_threshold = aggregator_func(map(int, gnomad_threshold))
		except:
			pass

		other_catalog_threshold = list(filter(None, [
			disease_info.get(threshold_field, None) for disease_info in other_catalog_lookup[locus_id].get("Diseases", [])
		]))
		try:
			other_catalog_threshold = aggregator_func(map(int, other_catalog_threshold))
		except:
			pass

		if not gnomad_threshold and not other_catalog_threshold:
			continue
		elif not gnomad_threshold:
			output(output_file, f"{threshold_field} for {locus_id:<6s} is set to {other_catalog_threshold} in {other_catalog_name} but not in gnomAD")
		elif not other_catalog_threshold:
			output(output_file, f"{threshold_field} for {locus_id:<6s} is set to {gnomad_threshold} in gnomAD but not in {other_catalog_name}")
		elif gnomad_threshold != other_catalog_threshold:
			output(output_file, f"{threshold_field} for {locus_id:<6s} differs between gnomAD ({gnomad_threshold}) and {other_catalog_name} ({other_catalog_threshold})")



def compare_catalogs(args, official_EH_catalog_loci, gnomad_catalog, stripy_lookup, strchive_lookup, trgt_lookup,
					 compare_with=OTHER_CATALOG_NAME_STRIPY):
	"""Compare the gnomAD catalog with the other catalog"""

	fasta_file = pysam.FastaFile(os.path.expanduser(args.reference_fasta))

	output_path = f"gnomAD_{compare_with}_comparison.txt"
	output_file = open(output_path, "wt")
	print(f"Writing to {output_path}")

	output_stripy_variant_catalog = []
	output_strchive_variant_catalog = []
	output_trgt_variant_catalog = []

	gnomad_locus_ids = {d["LocusId"] for d in gnomad_catalog}
	other_catalog_name = compare_with
	if compare_with == OTHER_CATALOG_NAME_STRCHIVE:
		other_catalog_locus_ids = set(strchive_lookup.keys())
	elif compare_with == OTHER_CATALOG_NAME_STRIPY:
		other_catalog_locus_ids = set(stripy_lookup.keys())
	elif compare_with == OTHER_CATALOG_NAME_TRGT:
		other_catalog_locus_ids = set(trgt_lookup.keys())
	else:
		raise ValueError(f"Unknown catalog name: {compare_with}")

	if compare_with == OTHER_CATALOG_NAME_STRCHIVE:
		other_catalog_lookup = strchive_lookup
	elif compare_with == OTHER_CATALOG_NAME_STRIPY:
		other_catalog_lookup = stripy_lookup
	elif compare_with == OTHER_CATALOG_NAME_TRGT:
		other_catalog_lookup = trgt_lookup
	else:
		raise ValueError(f"Unknown catalog name: {compare_with}")

	# check which loci are only in one of the catalogs
	gnomad_records = [d for d in gnomad_catalog if d["LocusId"] in other_catalog_lookup]
	for locus_id in sorted(other_catalog_locus_ids - gnomad_locus_ids):
		output(output_file, f"{locus_id} is in {other_catalog_name} but not in gnomAD")
	output(output_file, "------")
	for locus_id in sorted(gnomad_locus_ids - other_catalog_locus_ids):
		output(output_file, f"{locus_id} is in gnomAD but not in {other_catalog_name}")
	output(output_file, "------")

	# compare motifs for loci that are in both catalogs
	for d in gnomad_records:
		locus_id = d["LocusId"]
		canonical_motif = compute_canonical_motif(d["RepeatUnit"])
		if canonical_motif != other_catalog_lookup[locus_id]["CanonicalMotif"] and len(canonical_motif) != 5:
			output(output_file, f"Canonical motifs for {locus_id:<6s} differ between gnomAD ({canonical_motif}) and {other_catalog_name} ({other_catalog_lookup[locus_id]['CanonicalMotif']})")

	output(output_file, "------")

	# compare inheritance mode
	if compare_with != OTHER_CATALOG_NAME_TRGT:
		for d in gnomad_records:
			locus_id = d["LocusId"]
			gnomad_inheritance = set([disease_info.get("Inheritance", "") for disease_info in d.get("Diseases", [])])
			other_catalog_inheritance = set([
				disease_info.get("Inheritance", "") for disease_info in other_catalog_lookup[locus_id].get("Diseases", [])
			])
			other_catalog_inheritance = {i for i_list in other_catalog_inheritance for i in i_list.split("/")}
			other_catalog_inheritance = {{
				"Autosomal dominant": "AD",
				"Autosomal recessive": "AR",
				"X-linked recessive": "XR",
				"X-linked dominant": "XD",
			}.get(i, i) for i in other_catalog_inheritance}

			if gnomad_inheritance != other_catalog_inheritance:
				output(output_file, f"Inheritance mode for {locus_id:<6s} differ between gnomAD ({','.join(sorted(gnomad_inheritance))}) and {other_catalog_name} ({','.join(sorted(other_catalog_inheritance))})")

		output(output_file, "------")

	# compare STRipy, STRchive age of onset
	if compare_with != OTHER_CATALOG_NAME_TRGT:
		output(output_file, "Age of Onset:")
		output(output_file, "%20s %35s       %-100s" % ("LocusId", "STRipy", "STRchive"))

		for d in gnomad_records:
			locus_id = d["LocusId"]
			if locus_id not in stripy_lookup or locus_id not in strchive_lookup:
				continue
			stripy_info = stripy_lookup[locus_id]
			strchive_info = strchive_lookup[locus_id]
			stripy_age_of_onset = ", ".join(set([disease_info.get("Onset", "") for disease_info in stripy_info.get("Diseases", [])]))
			strchive_age_of_onset = ", ".join(set([disease_info.get("Onset", "") for disease_info in strchive_info.get("Diseases", [])]))
			output(output_file, "%20s %35s       %-100s" % (locus_id, stripy_age_of_onset, strchive_age_of_onset))

		output(output_file, "------")

	# compare thresholds between catalogs (normal max, pathogenic min)
	if compare_with != OTHER_CATALOG_NAME_TRGT:
		compare_threshold(output_file, gnomad_records, other_catalog_lookup, other_catalog_name, threshold_field="PathogenicMin", aggregator_func=min)
		output(output_file, "------")
		compare_threshold(output_file, gnomad_records, other_catalog_lookup, other_catalog_name, threshold_field="NormalMax", aggregator_func=max)

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
			#print("Parsing STRchive locus", locus_id)
			strchive_info = strchive_lookup[locus_id]
			strchive_reference_region = f"{strchive_info['chrom']}:{strchive_info['start_hg38']}-{strchive_info['stop_hg38']}"
			strchive_chrom, strchive_start, strchive_end = parse_interval(strchive_reference_region)
			strchive_repeats = (strchive_end - strchive_start)/motif_size
			strchive_pathogenic_min = strchive_info.get("pathogenic_min")
			if isinstance(strchive_pathogenic_min, float):
				strchive_pathogenic_min = int(strchive_pathogenic_min)

			strchive_motif = strchive_info["reference_motif_reference_orientation"]
			strchive_canonical_motif = compute_canonical_motif(strchive_motif)
			strchive_ref = fasta_file.fetch(strchive_chrom, strchive_start, strchive_end)

			gnomad_vs_strchive_distance = max(abs(gnomad_start - strchive_start), abs(gnomad_end - strchive_end))
			gnomad_differs_from_strchive = gnomad_reference_region != strchive_reference_region and (gnomad_vs_strchive_distance >= motif_size or gnomad_repeats != strchive_repeats)
			gnomad_differs_from_strchive = gnomad_differs_from_strchive or (len(gnomad_motif) != 5 and gnomad_canonical_motif != strchive_canonical_motif)

		trgt_reference_region = None
		gnomad_differs_from_trgt = False
		if locus_id in trgt_lookup:
			#print("Parsing TRGT locus", locus_id)
			trgt_info = trgt_lookup[locus_id]
			trgt_chrom, trgt_start, trgt_end = trgt_info["chrom"], trgt_info["start_0based"], trgt_info["end"]
			trgt_reference_region = trgt_info["ReferenceRegion"]
			trgt_canonical_motif = trgt_info["CanonicalMotif"]
			trgt_repeats = (trgt_end - trgt_start)/motif_size
			trgt_ref = fasta_file.fetch(trgt_chrom, trgt_start, trgt_end)

			gnomad_vs_trgt_distance = max(abs(gnomad_start - trgt_start), abs(gnomad_end - trgt_end))
			gnomad_differs_from_trgt = gnomad_reference_region != trgt_reference_region and (gnomad_vs_trgt_distance >= motif_size or gnomad_repeats != trgt_repeats)
			gnomad_differs_from_trgt = gnomad_differs_from_trgt or (len(gnomad_motif) != 5 and gnomad_canonical_motif != trgt_canonical_motif)

		if (
			(gnomad_differs_from_stripy and compare_with == OTHER_CATALOG_NAME_STRIPY)
			or (gnomad_differs_from_strchive and compare_with == OTHER_CATALOG_NAME_STRCHIVE)
			or (gnomad_differs_from_trgt and compare_with == OTHER_CATALOG_NAME_TRGT)
		):
			counter += 1
			output(output_file, "------")
			whats_different = "ReferenceRegion" if (gnomad_reference_region != stripy_reference_region and compare_with == OTHER_CATALOG_NAME_STRIPY) or (gnomad_reference_region != strchive_reference_region and compare_with == OTHER_CATALOG_NAME_STRCHIVE) or (gnomad_reference_region != trgt_reference_region and compare_with == OTHER_CATALOG_NAME_TRGT) else "motif"
			output(output_file, f"{locus_id} hg38 {whats_different} isn't the same between gnomAD and {other_catalog_name}")

			if locus_id in official_EH_catalog_loci and not isinstance(official_EH_catalog_loci[locus_id]['ReferenceRegion'], list):
				official_EH_reference_region = official_EH_catalog_loci[locus_id]['ReferenceRegion']
				_, official_EH_start, official_EH_end = parse_interval(official_EH_reference_region)
				official_EH_motif = official_EH_catalog_loci[locus_id]['LocusStructure'].strip("()*+")
				official_EH_repeats = (official_EH_end - official_EH_start)/motif_size
				output(output_file, "%3s %s" % (" ", f"chr{official_EH_reference_region}  ({official_EH_repeats:4.1f} x {official_EH_motif}  {official_EH_end - official_EH_start:5d}bp) official catalog from EH repo"))

			if locus_id in strchive_lookup and locus_id in stripy_lookup:
				min_left_pos = min(gnomad_start, stripy_start, strchive_start)
			elif locus_id in strchive_lookup:
				min_left_pos = min(gnomad_start, strchive_start)
			elif locus_id in stripy_lookup:
				min_left_pos = min(gnomad_start, stripy_start)
			elif locus_id in trgt_lookup:
				min_left_pos = min(gnomad_start, trgt_start)
			else:
				raise ValueError(f"Unknown locus_id: {locus_id}")

			output(output_file, f"%3s %-100s%{min(400, gnomad_end-min_left_pos)}s" % (
				" ", f"{gnomad_reference_region}  ({gnomad_repeats:4.1f} x {gnomad_motif}  {gnomad_end - gnomad_start:5d}bp) gnomAD    {gnomad_pathogenic_min} = pathogenic.min.  hg38 seq:", gnomad_ref))

			if locus_id in strchive_lookup:
				output(output_file, f"%3s %-100s%{min(400, strchive_end-min_left_pos)}s" % (
					" ", f"{strchive_reference_region}  ({strchive_repeats:4.1f} x {strchive_motif}  {strchive_end - strchive_start:5d}bp) STRchive  {strchive_pathogenic_min} = pathogenic.min.  hg38 seq:", strchive_ref))

				if gnomad_differs_from_strchive:
					output_strchive_variant_catalog.append({
						"LocusId": f"{locus_id}_STRCHIVE",
						"LocusStructure": "(" + strchive_info["reference_motif_reference_orientation"] + ")*",
						"ReferenceRegion": f"{strchive_info['chrom']}:{strchive_info['start_hg38']}-{strchive_info['stop_hg38']}",
						"VariantType": "Repeat",
					})

			if locus_id in stripy_lookup:
				output(output_file, f"%3s %-100s%{min(400, stripy_end-min_left_pos)}s" % (
					" ", f"{stripy_info['ReferenceRegion_hg38']}  ({stripy_repeats:4.1f} x {stripy_motif}  {stripy_end - stripy_start:5d}bp) STRipy    {stripy_pathogenic_min} = pathogenic.min.  hg38 seq:", stripy_ref))

				if gnomad_differs_from_stripy:
					output_stripy_variant_catalog.append({
						"LocusId": f"{locus_id}_STRIPY",
						"LocusStructure": "(" + stripy_info["Motif"] + ")*",
						"ReferenceRegion": stripy_info["ReferenceRegion_hg38"],
						"VariantType": "Repeat",
					})

			#if locus_id in trgt_lookup:
			#	output(output_file, f"%3s %-100s%{trgt_end-min_left_pos}s" % (
			#		" ", f"{trgt_reference_region}  ({trgt_repeats:4.1f} x {trgt_canonical_motif}  {trgt_end - trgt_start:5d}bp) TRGT                        hg38 seq:", trgt_ref))

			#	if gnomad_differs_from_strchive:
			#		output_trgt_variant_catalog.append({
			#			"LocusId": f"{locus_id}_TRGT",
			#			"LocusStructure": trgt_info["LocusStructure"],
			#			"ReferenceRegion": trgt_info["ReferenceRegion"],
			#			"VariantType": "Repeat",
			#		})


	output(output_file, f"\n{counter:,d} out of {len(gnomad_records):,d} loci have a different hg38 ReferenceRegion or motif in gnomAD vs other catalogs")

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
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-R", "--reference-fasta", help="hg38 reference fasta file path", default="~/hg38.fa")
	args = parser.parse_args()

	args.reference_fasta = os.path.expanduser(args.reference_fasta)

	if not os.path.isfile(args.reference_fasta):
		parser.error(f"Reference fasta file not found: {args.reference_fasta}")

	official_EH_catalog_loci = get_official_expansion_hunter_catalog_dict()
	gnomad_catalog = get_gnomad_catalog()
	trgt_lookup = get_trgt_catalog()

	strchive_lookup = get_strchive_dict(gnomad_catalog)
	stripy_lookup = get_stripy_dict()

	for compare_with in OTHER_CATALOG_NAME_STRCHIVE,:  # OTHER_CATALOG_NAME_TRGT, OTHER_CATALOG_NAME_STRIPY,
		print("="*100)
		print(f"Comparing gnomAD with {compare_with}")
		compare_catalogs(args, official_EH_catalog_loci, gnomad_catalog, stripy_lookup, strchive_lookup, trgt_lookup,
						 compare_with=compare_with)


if __name__ == "__main__":
	main()
