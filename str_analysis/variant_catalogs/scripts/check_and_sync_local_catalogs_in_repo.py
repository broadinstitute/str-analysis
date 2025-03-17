"""This script checks internal consistency and fields for the STR catalogs in this repo.
Since many fields are duplicated between the GRCh37 & GRCh38 versions of the catalogs, this script
syncs the duplicated data between them and reports any inconsistencies.
"""

import json
import os

from str_analysis.utils.eh_catalog_utils import parse_motifs_from_locus_structure

#%%

json_paths = (
    "../variant_catalog_without_offtargets.GRCh37.json",
    "../variant_catalog_without_offtargets.GRCh38.json",
    "../variant_catalog_with_offtargets.GRCh37.json",
    "../variant_catalog_with_offtargets.GRCh38.json",
)

json_path_pairs = [
    (json_paths[0], json_paths[2]),
    (json_paths[1], json_paths[3]),
]

json_path_without_offtargets_GRCh38 = json_paths[1]

data = {}
for json_path in json_paths:
    with open(json_path, "rt") as f:
        data[json_path] = json.load(f)

    print(f"Parsed {len(data[json_path])} records from {json_path}")


with open("../../data/non_ref_motif.locus_info.json", "rt") as f:
    pathogenic_motif_info = json.load(f)

REQUIRED_KEYS = {
    "LocusId",
    "LocusStructure",
    "ReferenceRegion",
    "VariantType",
}

EXTRA_KEYS = {
    "RepeatUnit",
    "Gene",
    "GeneId",
    "GeneRegion",
    "Diseases",
    "MainReferenceRegion",
}

for json_path, catalog in data.items():
    changes_made = False
    variant_ids = set()
    for record in catalog:
        locus_id = record["LocusId"]
        missing_keys = (REQUIRED_KEYS | EXTRA_KEYS) - set(record.keys())
        if missing_keys:
            print(f"ERROR: Missing keys in {json_path} {locus_id}: {missing_keys}")
        #extra_keys = set(record.keys()) - (REQUIRED_KEYS | EXTRA_KEYS)
        #if extra_keys:
        #    print(f"WARNING: Extra keys in {json_path} {locus_id}: {extra_keys}")
        if isinstance(record["ReferenceRegion"], list):
            motifs = parse_motifs_from_locus_structure(record["LocusStructure"])
            if len(record["ReferenceRegion"]) != len(motifs):
                print(f"ERROR: ReferenceRegion length != LocusStructure length in {json_path} {locus_id}")
            if len(motifs) != len(record["VariantType"]):
                print(f"ERROR: VariantType length != LocusStructure length in {json_path} {locus_id}")
            if "VariantId" not in record or len(motifs) != len(record["VariantId"]):
                print(f"ERROR: VariantId length != LocusStructure length in {json_path} {locus_id}")
            if record["MainReferenceRegion"] not in set(record["ReferenceRegion"]):
                print(f"ERROR: MainReferenceRegion is not one of the ReferenceRegions in {json_path} {locus_id}")
        else:
            if record["MainReferenceRegion"] != record["ReferenceRegion"]:
                print(f"ERROR: MainReferenceRegion != ReferenceRegion in {json_path} {locus_id}")

        for variant_id in record.get("VariantId", []):
            if variant_id in variant_ids:
                print(f"ERROR: Duplicate variant id {variant_id} in {json_path} {locus_id}")
            variant_ids.add(variant_id)

        if locus_id in pathogenic_motif_info:
            motifs = pathogenic_motif_info[locus_id].get('Motifs', {})
            #print(f"Found pathogenic motif info for {json_path} {locus_id}: {motifs}")
            if "PathogenicMotif" in record:
                del record["PathogenicMotif"]
                changes_made = True
                
            if "BenignMotif" in record:
                del record["BenignMotif"]
                changes_made = True
                
            if set(record.get("PathogenicMotifs", [])) != set(motifs.get("PATHOGENIC", [])):
                record["PathogenicMotifs"] = list(sorted(motifs["PATHOGENIC"]))
                print(record.get("PathogenicMotifs", []), motifs["PATHOGENIC"])
                changes_made = True

            if set(record.get("BenignMotifs", [])) != set(motifs.get("BENIGN", [])):
                record["BenignMotifs"] = list(sorted(motifs["BENIGN"]))
                print(record.get("BenignMotifs", []), motifs["BENIGN"])
                changes_made = True

    if changes_made:
        with open(json_path, "wt") as f:
            json.dump(catalog, f, indent=4)
        print(f"Wrote {len(catalog)} records to {json_path}")

#%%

for json_path_without_offtargets, json_path_with_offtargets in json_path_pairs:
    catalog_without_offtargets = data[json_path_without_offtargets]
    catalog_with_offtargets = data[json_path_with_offtargets]
    print(f"Comparing {json_path_without_offtargets} to {json_path_with_offtargets}")
    for record_without_offtargets in catalog_without_offtargets:
        locus_id = record_without_offtargets["LocusId"]
        record_with_offtargets = next((r for r in catalog_with_offtargets if r["LocusId"] == locus_id), None)

        # sync "Inheritance" field
        #if "Inheritance" in record_without_offtargets:
        #    record_with_offtargets["Inheritance"]  = record_without_offtargets["Inheritance"]

        if not record_with_offtargets:
            print(f"ERROR: {locus_id} not found in {json_path_with_offtargets}")
            continue
        for key in record_without_offtargets.keys():
            if key not in record_with_offtargets:
                print(f"ERROR: {key} not found in {json_path_with_offtargets} {locus_id}")
                continue
            if key == "VariantType":
                if isinstance(record_without_offtargets.get(key, []), list):
                    if set(record_without_offtargets.get(key, [])) != {"Repeat"}:
                        print(f"ERROR: VariantType is not set to 'Repeat' in {json_path_without_offtargets} {locus_id}: {record_without_offtargets[key]}")
                        continue
                    if "RareRepeat" not in set(record_with_offtargets.get(key, [])):
                        print(f"ERROR: VariantType is not set to 'RareRepeat' in {json_path_with_offtargets} {locus_id}: {record_with_offtargets[key]}")
                        continue
                else:
                    if record_without_offtargets[key] != "Repeat":
                        print(f"ERROR: VariantType is not set to 'Repeat' in {json_path_without_offtargets} {locus_id}: {record_without_offtargets[key]}")
                        continue
                    if record_with_offtargets[key] != "RareRepeat":
                        print(f"ERROR: VariantType is not set to 'RareRepeat' in {json_path_with_offtargets} {locus_id}: {record_with_offtargets[key]}")
                        continue

            elif record_without_offtargets[key] != record_with_offtargets[key]:
                print(f"ERROR: {key} mismatch in {json_path_without_offtargets} {locus_id}: {record_without_offtargets[key]} != {record_with_offtargets[key]}")
                continue
    #with open(json_path_with_offtargets, "wt") as f:
    #    json.dump(catalog_with_offtargets, f, indent=4)
#%%

os.chdir("../gnomad_notes")

# generate gnomad notes files if they don't exist already
for record in data[json_path_without_offtargets_GRCh38]:
    locus_id = record["LocusId"]
    filename = f"{locus_id}.md"
    if os.path.isfile(filename):
        continue
    print(f"Creating {filename}")
    with open(filename, "wt") as f:
        f.write("\n")


#%%
