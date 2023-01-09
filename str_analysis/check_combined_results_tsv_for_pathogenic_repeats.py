"""This script makes it easier to review an ExpansionHunter (or other) STR callset at known pathogenic loci and
manually identify samples that have pathogenic expansions.

The script takes a .tsv file generated by the combine_str_json_to_tsv.py script (which
combines multiple ExpansionHunter .json output files into a single .tsv table). It's assumed that this input .tsv will
have at least the columns listed below. One way to add these Sample_* and VariantCatalog_* columns is to run
combine_str_json_to_tsv.py with the --sample-metadata and the --variant-catalog args.

The input .tsv must have at least these input columns:

    "LocusId", "SampleId", "Sample_affected", "Sample_sex",
    "Num Repeats: Allele 1", "Num Repeats: Allele 2", "CI end: Allele 1", "CI end: Allele 2",

This script then prints out a table for each locus with sample genotypes sorted by most expanded to least expanded.
It can print all sample genotypes, or set a cut-off based on the known pathogenic threshold, or alternatively, by
print all affected individuals that are more expanded than some user-defined number of unaffected individuals.
It also takes into account locus inheritance - for example, printing 2 separate tables (males, females) for X-linked
recessive loci.
"""

import argparse
import os
import subprocess
from ast import literal_eval

import pandas as pd
from tabulate import tabulate


REQUIRED_COLUMNS = [
    "LocusId",
    "SampleId",
    "Sample_affected",
    "Sample_sex",
    "Num Repeats: Allele 1",
    "Num Repeats: Allele 2",
    "CI end: Allele 1",
    "CI end: Allele 2",
]

OPTIONAL_COLUMNS = [
    "RepeatUnit",
    "VariantId",

    "Sample_analysis_status",
    "Sample_coded_phenotype",
    "Sample_phenotypes",
    "Sample_genome_version",

    "VariantCatalog_Inheritance",
    "VariantCatalog_Diseases",

    "Genotype",
    "GenotypeConfidenceInterval",
]

GNOMAD_STR_TABLE_PATH = "https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.3/tsv/gnomAD_STR_genotypes__2022_01_20.tsv.gz"


def run(command):
    print(command)
    subprocess.check_call(command, shell=True)


def parse_args():
    p = argparse.ArgumentParser()
    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument("--use-affected", action="store_true", help="Use affected status to determine which samples to "
        "include in the output for each locus. Only include those affected samples that have larger expansions than "
        "the top N unaffected samples (adjustsable by --max-n-unaffected).")
    grp.add_argument("--use-thresholds", action="store_true", help="Use known pathogenic thresholds to determine which "
        "samples to include in the output for each locus. All samples with expansions above the pathogenic threshold "
        "will be included.")
    p.add_argument("-n", "--max-rows", type=int, default=10000000, help="Limit the max number of samples to include in "
        "the output for each locus.")
    p.add_argument("--max-n-unaffected", type=int, default=10, help="After this many unaffected samples are "
        "encountered with genotypes above a particular expansion size at locus, all samples (affected or unaffected) "
        "that have smaller expansions at that locus will be ignored")
    p.add_argument("--use-gnomad", action="store_true", help="Include samples from the gnomAD v3.1 STR release"
        "@ https://gnomad.broadinstitute.org/downloads#v3-short-tandem-repeats")
    p.add_argument("-l", "--locus", action="append", help="If specified, only these locus ids will be processed")
    p.add_argument("--highlight-samples", nargs="*", help="If specified, this can be the path of a text file that "
        "contains sample ids (one per line) or just 1 or more sample ids listed on the commandline - one per line")
    p.add_argument("combined_tsv_path", nargs="+", help="Path of combined ExpansionHunter .tsv table generated by the "
        "combine_str_json_to_tsv.py script. It's assumed that combine_str_json_to_tsv.py "
        "was run with --sample-metadata and --variant-catalog args to add sample-level and locus-level metadata columns")

    # Parse and validate command-line args + read in the combined table(s) from the given command_tsv_path(s)
    args = p.parse_args()

    for combined_tsv_path in args.combined_tsv_path:
        if not os.path.isfile(combined_tsv_path):
            p.error(f"{combined_tsv_path} not found")

    return args


def load_gnomad_df():
    print("Loading gnomAD STR table...")
    gnomad_df = pd.read_table(GNOMAD_STR_TABLE_PATH, usecols=[
        "Id",
        "LocusId",
        "ReferenceRegion",
        "Motif",
        "IsAdjacentRepeat",
        "Sex",
        "Allele1",
        "Allele2",
        "Genotype",
        "GenotypeConfidenceInterval",
    ],
        dtype={
          "Id": str,
          "Sex": str,
          "Allele1": str,
          "Allele2": str,
          "GenotypeConfidenceInterval": str,
          "IsAdjacentRepeat": bool,
        },
    )
    gnomad_df = gnomad_df[~gnomad_df["IsAdjacentRepeat"]]
    print(f"Processing {len(gnomad_df)} gnomAD table rows...")
    for column_name, value in [
        ("Sample_genome_version", "38"),
        ("Sample_affected", "Not Affected"),
    ]:
        gnomad_df.loc[:, column_name] = value

    gnomad_df.loc[:, "Id"] = "gnomAD:" + gnomad_df["Id"]
    gnomad_df.loc[:, "VariantId"] = gnomad_df["LocusId"]

    gnomad_df.loc[:, ("CI1", "CI2")] = gnomad_df["GenotypeConfidenceInterval"].str.split("/", expand=True)
    gnomad_df.loc[:, "CI1"] = gnomad_df["CI1"].fillna("")
    gnomad_df.loc[:, "CI2"] = gnomad_df["CI2"].fillna("")
    gnomad_df.loc[:, ("CI start: Allele 1", "CI end: Allele 1")] = gnomad_df["CI1"].str.split("-", expand=True)
    gnomad_df.loc[:, ("CI start: Allele 2", "CI end: Allele 2")] = gnomad_df["CI2"].str.split("-", expand=True)
    gnomad_df.drop(columns=["CI1", "CI2", "CI start: Allele 1", "CI start: Allele 2"], inplace=True)
    gnomad_df.rename({
        "Id": "SampleId",
        "Sex": "Sample_sex",
        "Allele1": "Num Repeats: Allele 1",
        "Allele2": "Num Repeats: Allele 2",
        "Motif": "RepeatUnit",
    }, axis="columns", inplace=True)

    return gnomad_df


def load_results_tables(args):
    all_dfs = []
    for combined_tsv_path in args.combined_tsv_path:

        # read in table
        df = pd.read_table(combined_tsv_path, low_memory=False)
        df.loc[:, "Num Repeats: Allele 1"] = df["Num Repeats: Allele 1"].fillna(0).astype(int)
        df.loc[:, "Num Repeats: Allele 2"] = df["Num Repeats: Allele 2"].fillna(0).astype(int)

        df.loc[:, "Num Repeats: Min Allele 1, 2"] = df.apply(lambda row: (
            min(row["Num Repeats: Allele 1"], row["Num Repeats: Allele 2"])
            if
            row["Num Repeats: Allele 1"] > 0 and row["Num Repeats: Allele 2"] > 0
            else
            max(row["Num Repeats: Allele 1"], row["Num Repeats: Allele 2"])
        ), axis=1)

        df.loc[:, "Num Repeats: Max Allele 1, 2"] = df[["Num Repeats: Allele 1", "Num Repeats: Allele 2"]].max(axis=1)

        # check that all required columns are present
        missing_required_columns = set(REQUIRED_COLUMNS) - set(df.columns)
        if missing_required_columns:
            raise ValueError(f"{combined_tsv_path} is missing these required columns: {missing_required_columns}")

        # fill in values for missing optional columns
        missing_optional_columns = set(OPTIONAL_COLUMNS) - set(df.columns)
        if missing_optional_columns:
            print(f"WARNING: {combined_tsv_path} is missing these columns: {missing_optional_columns}. "
                  f"Filling them with None...")
            for c in missing_optional_columns:
                df.loc[:, c] = None

        for integer_column in ("CI end: Allele 1", "CI end: Allele 2"):
            df.loc[:, integer_column] = pd.to_numeric(df[integer_column], errors="coerce")

        all_dfs.append((df, combined_tsv_path))

    return all_dfs


def print_results_for_locus(args, locus_df, locus_id, highlight_sample_ids):
    """Prints the sorted table of results for a single locus"""

    for column_name in (
            "Sample_affected", "Sample_sex", "Sample_analysis_status", "Sample_coded_phenotype",
            "Sample_genome_version",
    ):
        # split values like a; b; b; and collapse to a; b
        locus_df.loc[:, column_name] = locus_df[column_name].apply(
            lambda s: "; ".join(set(s.split("; "))) if not pd.isna(s) else s)
        # truncate long text so it fits on screen
        locus_df.loc[:, column_name] = locus_df[column_name].str[0:50]

    unexpected_affected_status_values = set(locus_df[
        ~locus_df["Sample_affected"].isin({"Affected", "Not Affected", "Unknown"})
    ].Sample_affected)
    if unexpected_affected_status_values:
        raise ValueError(f"Some rows have unexpected affected value(s): {unexpected_affected_status_values}")

    # Sort
    locus_df = locus_df.sort_values(
        by=["Num Repeats: Allele 2", "Num Repeats: Allele 1", "Sample_affected"],
        ascending=False)

    # Get the 1st row and use it to look up metadata values which are the same across all rows for the locus
    # (ie. Inheritance Mode)
    first_row = locus_df.iloc[0].to_dict()
    inheritance_mode = first_row.get("VariantCatalog_Inheritance") or "Unknown"

    disease_info = first_row.get("VariantCatalog_Diseases")
    intermediate_threshold_min = None
    pathogenic_threshold_min = None
    if disease_info and not pd.isna(disease_info):
        try:
            disease_info = literal_eval(disease_info)
        except Exception as e:
            raise ValueError(f"Unable to parse {disease_info} as json: {e}")

        try:
            intermediate_threshold_min = min(int(float(d["NormalMax"]) + 1) for d in disease_info if d.get("NormalMax"))
        except ValueError as e:
            print(f"WARNING: {locus_id} couldn't parse NormalMax fields from disease_info {disease_info}: {e}")

        try:
            pathogenic_threshold_min = min(int(float(d["PathogenicMin"])) for d in disease_info if d.get("PathogenicMin"))
        except ValueError as e:
            print(f"WARNING: {locus_id} couldn't parse PathogenicMin fields from disease_info {disease_info}: {e}")

    reference_region = first_row["ReferenceRegion"]
    genome_version = f"GRCh{first_row['Sample_genome_version']}" if first_row.get('Sample_genome_version') else ""
    motif = first_row.get("RepeatUnit")
    locus_description = f"{locus_id} ({reference_region}: {genome_version})  https://stripy.org/database/{locus_id}"
    print("**Locus**: ", locus_description)
    print("**Disease**: ", str(disease_info))
    print("**Inheritance**: ", inheritance_mode)
    print("**Pathogenic Threshold**: >=", pathogenic_threshold_min, "x", motif)
    print("**Intermediate Threshold**: >=", intermediate_threshold_min, "x", motif)

    # replace NA with "Unknown" strings
    locus_df.loc[:, "Sample_affected"] = locus_df["Sample_affected"].fillna("Unknown")
    locus_df.loc[:, "Sample_sex"] = locus_df["Sample_sex"].fillna("Unknown")

    # create a list of dfs that are subsets of locus_df and where rows pass thresholds
    dfs_to_process = []
    if locus_id == "COMP":
        # COMP is a special case where contractions below 5 repeats are also pathogenic
        df_comp = locus_df[
            (locus_df["CI end: Allele 1"] < 5) |
            (locus_df["CI end: Allele 2"] < 5)
        ]
        dfs_to_process.append(df_comp)

    # add 2 rows to each table representing the pathogenic thresholds
    threshold_records = []
    threshold_records_for_male_samples = []
    for label, threshold in [
        ("Intermediate Threshold", intermediate_threshold_min),
        ("Pathogenic Threshold",   pathogenic_threshold_min),
    ]:
        if not threshold: continue

        threshold_record = {c: "---" for c in locus_df.columns}
        threshold_record.update({
            "SampleId": f"**{label}**", "LocusId": locus_id, "VariantId": locus_id,
            "Genotype": f">= {threshold}",
            "Num Repeats: Min Allele 1, 2": threshold,
            "Num Repeats: Max Allele 1, 2": threshold,
            "Num Repeats: Allele 1": threshold,
            "Num Repeats: Allele 2": 0,
        })

        threshold_records_for_male_samples.append(dict(threshold_record))

        threshold_record.update({
            "Num Repeats: Allele 2": threshold,
        })
        threshold_records.append(dict(threshold_record))

    threshold_records = pd.DataFrame(threshold_records)
    threshold_records_for_male_samples = pd.DataFrame(threshold_records_for_male_samples)
    args.max_rows += 2  # add 2 to allow for these threshold rows

    # create a list of dfs to print, filtered by the pathogenic thresholds and/or affected status
    threshold = (intermediate_threshold_min or pathogenic_threshold_min) if args.use_thresholds else 0

    if inheritance_mode == "XR":

        # split results by male/female
        male_df = locus_df[~locus_df["Genotype"].str.contains("/")]
        if args.use_thresholds:
            male_df = male_df[
                          (male_df["CI end: Allele 1"] >= threshold)
                      ].iloc[0:args.max_rows]
        male_df = pd.concat([threshold_records_for_male_samples, male_df], ignore_index=True)
        dfs_to_process.append(male_df)

        female_df = locus_df[locus_df["Genotype"].str.contains("/")]
        if args.use_thresholds:
            female_df = female_df[
                            (female_df["CI end: Allele 1"] >= threshold) &
                            (female_df["CI end: Allele 2"] >= threshold)
                            ].iloc[0:args.max_rows]

        female_df = pd.concat([threshold_records, female_df], ignore_index=True)
        dfs_to_process.append(female_df)

    else:
        if args.use_thresholds:
            if inheritance_mode == "AR":
                locus_df = locus_df[
                                        (locus_df["CI end: Allele 1"] >= threshold) &
                                        (locus_df["CI end: Allele 2"] >= threshold)
                                        ].iloc[0:args.max_rows]
            else:
                # for x-linked dominant, it's enough for one of the alleles to be above the threshold
                # (ie. male vs. female genotyes)
                locus_df = locus_df[
                                        (locus_df["CI end: Allele 1"] >= threshold) |
                                        (locus_df["CI end: Allele 2"] >= threshold)
                                        ].iloc[0:args.max_rows]

        locus_df = pd.concat([threshold_records, locus_df], ignore_index=True)
        dfs_to_process.append(locus_df)

    if args.use_affected:
        # filter by affected status
        filtered_dfs_list = []
        for df_to_process in dfs_to_process:
            unaffected_counter = 0
            idx = 0
            for affected_status in df_to_process["Sample_affected"]:
                idx += 1
                if affected_status and ("Not Affected" in affected_status or "Unknown" in affected_status):
                    unaffected_counter += 1

                if unaffected_counter >= args.max_n_unaffected:
                    break

            df_to_process = df_to_process.iloc[:idx]

            filtered_dfs_list.append(df_to_process)
        dfs_to_process = filtered_dfs_list

    for i, df_to_process in enumerate(dfs_to_process):
        print(f"Found {len(df_to_process)} samples passed filters in table {i+1} out of "
              f"{len(dfs_to_process)} for locus {locus_id}")

        if len(df_to_process) == 0:
            continue

        if inheritance_mode in {"AR", "XR"}:
            df_to_process = df_to_process.sort_values(
                by=["Num Repeats: Min Allele 1, 2", "Sample_affected"],
                ascending=False)
        elif inheritance_mode in {"AD", "XD", "Unknown"}:
            df_to_process = df_to_process.sort_values(
                by=["Num Repeats: Max Allele 1, 2", "Sample_affected"],
                ascending=False)
        else:
            raise ValueError(f"Unexpected inheritance mode: {inheritance_mode}")

        # Print the filtered results for this locus
        df_to_process = df_to_process[[
            "SampleId",
            "LocusId",
            #"VariantId",
            "Sample_affected",
            "Sample_sex",
            "Sample_genome_version",

            "Genotype",
            "GenotypeConfidenceInterval",
            #"Num Repeats: Min Allele 1, 2",

            "VariantCatalog_Inheritance",
            "RepeatUnit",

            "Sample_analysis_status",
            "Sample_coded_phenotype",
            "Sample_phenotypes",
        ]]

        if highlight_sample_ids:
            highlight_sample_ids = set(highlight_sample_ids)
            df_to_process.loc[:, "SampleId"] = df_to_process["SampleId"].apply(
                lambda s: (s if s not in highlight_sample_ids else f"==> {s}"))

        # Shorten some column names so more can fit on screen
        df_to_process.rename(columns={
            "Sample_affected": "affected",
            "Sample_sex": "sex",
            "Sample_analysis_status": "analysis_status",
            "Sample_coded_phenotype": "coded_phenotype",
            "Sample_phenotypes": "hpo",
            "Sample_genome_version": "genome",
            "GenotypeConfidenceInterval": "GenotypeCI",
            "VariantCatalog_Inheritance": "Mode",
            "RepeatUnit": "Motif",
        }, inplace=True)

        # Print the candidate pathogenic rows for this locus
        print(tabulate(df_to_process, headers="keys", tablefmt="psql", showindex=False))


def main():
    args = parse_args()

    highlight_sample_ids = []
    if args.highlight_samples:
        for h in args.highlight_samples:
            if os.path.isfile(h):
                with open(h, "rt") as f:
                    sample_ids_list = [s.strip() for s in f.readlines()]
                print(f"Read {len(sample_ids_list)} sample ids to highlight from: {h}")
                highlight_sample_ids.extend(sample_ids_list)
            else:
                highlight_sample_ids.append(h)

        print(f"Will highlight {len(highlight_sample_ids)} sample ids: {highlight_sample_ids}")
        
    all_dfs = load_results_tables(args)

    all_locus_ids = set()
    for df, _ in all_dfs:
        all_locus_ids |= set(df.LocusId)

    all_variant_ids = set()
    for df, _ in all_dfs:
        all_variant_ids |= set(df.VariantId)

    variant_ids_that_are_not_locus_ids = all_variant_ids - all_locus_ids
    if variant_ids_that_are_not_locus_ids:
        print("WARNING: discarding records with VariantIds:", ", ".join(variant_ids_that_are_not_locus_ids))
        for i, (df, df_source_path) in enumerate(all_dfs):
            df = df[~df.VariantId.isin(variant_ids_that_are_not_locus_ids)]
            all_dfs[i] = (df, df_source_path)

    if args.use_gnomad:
        gnomad_df = load_gnomad_df()

        for i, df in enumerate(all_dfs):
            print("gnomad_df columns: ", gnomad_df.columns)
            print("df columns: ", df.columns)
            df_with_gnomad = pd.concat([df, gnomad_df])
            df_with_gnomad.fillna("", inplace=True)
            all_dfs[i] = df_with_gnomad

    # Process each locus
    for locus_id in sorted(all_locus_ids):
        if args.locus and locus_id not in args.locus:
            continue

        # Filter to rows for the current locus
        for i, (full_df, df_source_path) in enumerate(all_dfs):
            print("="*100)  # print a divider
            print(f"** {locus_id} from {df_source_path} **")
            locus_df = full_df[full_df.LocusId == locus_id]
            if len(locus_df) == 0:
                return

            print_results_for_locus(args, locus_df, locus_id, highlight_sample_ids)


if __name__ == "__main__":
    main()
