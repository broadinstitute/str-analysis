"""This script takes an ExpansionHunter catalog and outputs a new catalog where adjacent loci have been combined.
"""
import simplejson as json
from pprint import pformat
import re

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif

# Example input table row
"""
$1                                   contig : chr15
$2                                    start : 71459006
$3                                      end : 71459475
$4                                    motif : AG
$5                          top_case_zscore : 1.0
$6                                   Source : test.outlier_locus.tsv
$7                                MotifSize : 2
$8                           CanonicalMotif : AG
$9       MatchedKnownDiseaseAssociatedLocus :
$10      MatchedKnownDiseaseAssociatedMotif : False
$11                       GencodeGeneRegion : intron
$12                         GencodeGeneName : THSD4
$13                           GencodeGeneId : ENSG00000187720
$14                     GencodeTranscriptId : ENST00000355327
$15                          ManeGeneRegion : intron
$16                            ManeGeneName : THSD4
$17                              ManeGeneId : ENSG00000187720
$18                        ManeTranscriptId : ENST00000355327
$19                      MatchedReferenceTR : 15-71459381-71459405-TC
$20  MatchedTRFrom:illumina_variant_catalog : 15-71459381-71459404-TC
$21      OverlapsWith:GRCh38GenomicSuperDup : False
$22                                SampleId : sample1
$23                         NormalizedCount : 1.0
$24                       SampleRankAtLocus : 1
$25                     TotalSamplesAtLocus : 1
"""


import argparse
import pandas as pd

VALID_GENE_REGIONS = {"CDS", "UTR", "5UTR", "3UTR", "promoter", "exon", "intron", "intergenic"}

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-motif-size", type=int, help="Minimum motif size to include in the output catalog")
    parser.add_argument("--max-motif-size", type=int, help="Maximum motif size to include in the output catalog")
    parser.add_argument("-ms", "--motif-size", type=int, action="append", help="Only include loci with these motif sizes")
    parser.add_argument("-m", "--motif", action="append", help="Only include loci whose canonical motif Matched this motif")

    parser.add_argument("--max-rank", type=int, help="Maximum sample rank to include in the output catalog")

    parser.add_argument("--gene-column-prefix", default="Gencode", help="Prefix for gene columns in the input table")
    parser.add_argument("-t", "--region-type", action="append", help="Gene region(s) to include in the output catalog",
                        choices=VALID_GENE_REGIONS)
    parser.add_argument("-xt", "--exclude-region-type", action="append", choices=VALID_GENE_REGIONS,
                        help="Gene region(s) to exclude from the output catalog",)

    parser.add_argument("-g", "--gene-name", action="append", help="Only include loci in this gene")
    parser.add_argument("-xg", "--exclude-gene-name", action="append", help="Exclude loci in this gene")

    parser.add_argument("--gene-id", action="append", help="Only include loci with this gene id")
    parser.add_argument("--exclude-gene-id", action="append", help="Exclude loci with this gene id")

    parser.add_argument("-s", "--sample-id", action="append", help="Sample IDs to include in the output catalog")

    parser.add_argument("-ok", "--only-known-disease-associated-loci", action="store_true",
                        help="Only include known disease-associated loci")
    parser.add_argument("-xk", "--exclude-known-disease-associated-loci", action="store_true",
                        help="Exclude known disease-associated loci")

    parser.add_argument("-om", "--only-known-disease-associated-motifs", action="store_true",
                        help="Only include loci that have known disease-associated motifs")

    parser.add_argument("-R", "--matched-reference-repeats", action="store_true",
                        help="Only include loci that matched the reference repeats")

    parser.add_argument("-l", "--locus-id", action="append", help="Only include the locus with this locus id")
    parser.add_argument("-xl", "--exclude-locus-id", action="append", help="Exclude the locus with this locus id")

    #parser.add_argument("-L", "--region", action="append", help="Only include loci in this genomic region")

    parser.add_argument("--matched-tr-from", action="append", help="Column name suffix after 'MatchedTRFrom:'", default=[])
    parser.add_argument("--exclude-matched-tr-from", action="append", help="Exclude loci that matched a TR in this column", default=[])

    parser.add_argument("--overlaps-with", action="append", help="Column name suffix after 'OverlapsWith:'", default=[])
    parser.add_argument("--exclude-overlaps-with", action="append", help="Exclude loci that overlapped with an interval in this column", default=[])

    parser.add_argument("--output-prefix", help="Optional output filename prefix. If not specified, the input_tsv filename "
                                                "(with out the .tsv) will be used as the prefix")

    parser.add_argument("input_tsv", nargs="+", help="Table(s) generated by annotate_EHdn_locus_outliers.py")
    args = parser.parse_args()

    df = None
    for input_tsv in args.input_tsv:
        print(f"Parsing {input_tsv}")
        current_df = pd.read_table(input_tsv, dtype={"contig": str, "start": int, "end": int})
        current_df = current_df.rename(columns={"motif": "Motif"})

        if f"{args.gene_column_prefix}GeneRegion" not in current_df.columns:
            parser.error(f"{args.gene_column_prefix}GeneRegion column not found in {input_tsv}")

        for suffix in args.matched_tr_from:
            if f"MatchedTRFrom:{suffix}" not in current_df.columns:
                parser.error(f"MatchedTRFrom:{suffix} column not found in {input_tsv}")

        for suffix in args.overlaps_with:
            if f"OverlapsWith:{suffix}" not in current_df.columns:
                parser.error(f"OverlapsWith:{suffix} column not found in {input_tsv}")

        expected_columns = [
            "contig",
            "start",
            "end",

            "Source",
            "Motif",
            "MotifSize",
            "CanonicalMotif",
            "MatchedReferenceTR",
            "MatchedKnownDiseaseAssociatedLocus",
            "MatchedKnownDiseaseAssociatedMotif",
            f"{args.gene_column_prefix}GeneRegion",
            f"{args.gene_column_prefix}GeneName",
            f"{args.gene_column_prefix}GeneId",
            f"{args.gene_column_prefix}TranscriptId",
            "SampleId",
            "SampleRankAtLocus",
            f"TotalSamplesAtLocus",
            "NormalizedCount",
        ]
        if args.overlaps_with:
            for matched_tr_from in args.matched_tr_from:
                expected_columns.append(f"MatchedTRFrom:{matched_tr_from}")

            for overlaps_with in args.overlaps_with:
                expected_columns.append(f"OverlapsWith:{overlaps_with}")


        missing_columns = set(expected_columns) - set(current_df.columns)
        if missing_columns:
            parser.error(f"{input_tsv} is missing expected column(s): {missing_columns}")

        current_df = current_df[expected_columns]
        if df is None:
            df = current_df
        else:
            df = pd.concat([df, current_df], axis=0)

        print(f"Loaded {len(current_df):,d} rows from {input_tsv}")

    return df, args


def group_by_locus(df, gene_column_prefix):
    df = df[[
        "MatchedReferenceTR",
        "Motif",
        "MotifSize",
        "CanonicalMotif",
        f"{gene_column_prefix}GeneRegion",
        f"{gene_column_prefix}GeneName",
        f"{gene_column_prefix}GeneId",
        f"{gene_column_prefix}TranscriptId",
        f"TotalSamplesAtLocus",
        "SampleId",
        "SampleRankAtLocus",
        "Source",
        "NormalizedCount",
    ]]

    df_loci = df.groupby([
        "MatchedReferenceTR",
        "Motif",
        "MotifSize",
        "CanonicalMotif",
        f"{gene_column_prefix}GeneRegion",
        f"{gene_column_prefix}GeneName",
        f"{gene_column_prefix}GeneId",
        f"{gene_column_prefix}TranscriptId",
    ]).agg({
        "SampleId": lambda x: ",".join(x),
        "SampleRankAtLocus": lambda x: ",".join(map(str, x)),
        "TotalSamplesAtLocus": lambda x: ",".join(map(str, x)),
        "Source": lambda x: ",".join(map(str, x)),
    }).reset_index()

    return df_loci


def main():
    df, args = parse_args()

    # Filter by motif size
    if args.min_motif_size is not None:
        before = len(df)
        df = df[df["MotifSize"] >= args.min_motif_size]
        print(f"Filter by min motif size: filtered out {before - len(df):,d} out of {before:,d} rows ({(before - len(df)) / (before or 1):.1%}%)")
    if args.max_motif_size is not None:
        before = len(df)
        df = df[df["MotifSize"] <= args.max_motif_size]
        print(f"Filter by max motif size: filtered out {before - len(df):,d} out of {before:,d} rows ({(before - len(df)) / (before or 1):.1%}%)")

    if args.motif_size is not None:
        before = len(df)
        df = df[df["MotifSize"].isin(set(args.motif_size))]
        print(f"Filter to specific motif sizes {list(sorted(args.motif_size))}: filtered out {before - len(df):,d} out of {before:,d} rows ({(before - len(df)) / (before or 1):.1%}%)")

    # Filter by motif
    if args.motif is not None:
        canonical_motifs = {compute_canonical_motif(m) for m in args.motif}
        before = len(df)
        df = df[df["CanonicalMotif"].isin(canonical_motifs)]
        print(f"Filter by motif: filtered out {before - len(df):,d} out of {before:,d} rows ({(before - len(df)) / (before or 1):.1%}%)")

    # Filter by region
    #if args.region is not None:
    #    before = len(df)
    #    # NOT YET IMPLEMENTED
    #    print(f"Filter by region: filtered out {before - len(df):,d} out of {before:,d} rows ({(before - len(df)) / (before or 1):.1%}%)")

    # Filter by region type
    if args.region_type is not None:
        region_types = set(args.region_type)
        for region_type in region_types:
            if region_type == "UTR" or region_type == "5UTR":
                region_types.add("5' UTR")
            if region_type == "UTR" or region_type == "3UTR":
                region_types.add("3' UTR")
        before = len(df)
        df = df[df[f"{args.gene_column_prefix}GeneRegion"].isin(region_types)]
        print(f"Filter to region type {region_types}: filtered out {before - len(df):,d} out of {before:,d} rows ({(before - len(df)) / (before or 1):.1%}%)")

    # Filter by gene name
    if args.gene_name is not None:
        before = len(df)
        df = df[df[f"{args.gene_column_prefix}GeneName"].isin(args.gene_name)]
        print(f"Filter by gene name: filtered out {before - len(df):,d} out of {before:,d} rows ({(before - len(df)) / (before or 1):.1%}%)")

    # Filter by gene id
    if args.gene_id is not None:
        before = len(df)
        df = df[df[f"{args.gene_column_prefix}GeneId"].isin(args.gene_id)]
        print(f"Filter by gene id: filtered out {before - len(df):,d} out of {before:,d} rows ({(before - len(df)) / (before or 1):.1%}%)")

    # Filter by sample id
    if args.sample_id is not None:
        before = len(df)
        df = df[df["SampleId"].isin(args.sample_id)]
        print(f"Filter by sample id: filtered out {before - len(df):,d} out of {before:,d} rows ({(before - len(df)) / (before or 1):.1%}%)")

    # Filter by sample rank
    if args.max_rank is not None:
        before = len(df)
        df = df[df["SampleRankAtLocus"] <= args.max_rank]
        print(f"Filter by sample rank: filtered out {before - len(df):,d} out of {before:,d} rows ({(before - len(df)) / (before or 1):.1%}%)")

    # Filter by --matched-reference-repeats
    if args.matched_reference_repeats:
        before = len(df)
        df = df[~df["MatchedReferenceTR"].isna()]
        print(f"Filter by MatchedReferenceTR: filtered out {before - len(df):,d} out of {before:,d} rows ({(before - len(df)) / (before or 1):.1%}%)")

    # Filter by known disease-associated loci
    if args.only_known_disease_associated_loci:
        before = len(df)
        df = df[~df["MatchedKnownDiseaseAssociatedLocus"].isna()]
        print(f"Filter to known disease-associated loci: filtered out {before - len(df):,d} out of {before:,d} rows ({(before - len(df)) / (before or 1):.1%}%)")

    # Filter by known disease-associated motifs
    if args.only_known_disease_associated_motifs:
        before = len(df)
        df = df[df["MatchedKnownDiseaseAssociatedMotif"]]
        print(f"Filter to known disease-associated motifs: filtered out {before - len(df):,d} out of {before:,d} rows ({(before - len(df)) / (before or 1):.1%}%)")

    # Filter by locus id
    if args.locus_id is not None:
        before = len(df)
        df = df[df["MatchedReferenceTR"].isin(args.locus_id)]
        print(f"Filter by locus id: filtered out {before - len(df):,d} out of {before:,d} rows ({(before - len(df)) / (before or 1):.1%}%)")

    # Filter --exclude-locus-id
    if args.exclude_locus_id is not None:
        before = len(df)
        df = df[~df["MatchedReferenceTR"].isin(args.exclude_locus_id)]
        print(f"Filter by --exclude-locus-id: filtered out {before - len(df):,d} out of {before:,d} rows ({(before - len(df)) / (before or 1):.1%}%)")

    # Filter by MatchedTrFrom:
    for matched_tr_from in args.matched_tr_from:
        before = len(df)
        df = df[~df[f"MatchedTrFrom:{matched_tr_from}"].isna()]
        print(f"Filter by MatchedTrFrom:{matched_tr_from}: filtered out {before - len(df):,d} out of {before:,d} rows ({(before - len(df)) / (before or 1):.1%}%)")

    # Filter by OverlapsWith:
    for overlaps_with in args.overlaps_with:
        before = len(df)
        df = df[~df[f"OverlapsWith:{overlaps_with}"]]
        print(f"Filter by OverlapsWith:{overlaps_with}: filtered out {before - len(df):,d} out of {before:,d} rows ({(before - len(df)) / (before or 1):.1%}%)")

    # Filter --exclude-region-type
    if args.exclude_region_type is not None:
        exclude_region_types = set(args.exclude_region_type)
        for region_type in exclude_region_types:
            if region_type == "UTR" or region_type == "5UTR":
                exclude_region_types.add("5' UTR")
            if region_type == "UTR" or region_type == "3UTR":
                exclude_region_types.add("3' UTR")
        before = len(df)
        df = df[~df[f"{args.gene_column_prefix}GeneRegion"].isin(exclude_region_types)]
        print(f"Filter by exclude region type: filtered out {before - len(df):,d} out of {before:,d} rows ({(before - len(df)) / (before or 1):.1%}%)")

    # Filter by --exclude-gene-name
    if args.exclude_gene_name is not None:
        before = len(df)
        df = df[df[f"{args.gene_column_prefix}GeneName"].isin(args.exclude_gene_name)]
        print(f"Filter by exclude gene name: filtered out {before - len(df):,d} out of {before:,d} rows ({(before - len(df)) / (before or 1):.1%}%)")

    # Filter by --exclude-gene-id
    if args.exclude_gene_id is not None:
        before = len(df)
        df = df[df[f"{args.gene_column_prefix}GeneId"].isin(args.exclude_gene_id)]
        print(f"Filter by exclude gene id: filtered out {before - len(df):,d} out of {before:,d} rows ({(before - len(df)) / (before or 1):.1%}%)")

    # Filter by --exclude-known-disease-associated-loci
    if args.exclude_known_disease_associated_loci:
        before = len(df)
        for _, row in group_by_locus(df[~df["MatchedKnownDiseaseAssociatedLocus"].isna()], args.gene_column_prefix).iterrows():
            d = row.to_dict()
            d["SampleId"] = ", ".join(d["SampleId"].split(",")[0:7])
            d = {k: v for k, v in d.items() if k in {"CanonicalMotif", f"{args.gene_column_prefix}GeneName", f"{args.gene_column_prefix}GencodeGeneRegion", "SampleId", "MatchedReferenceTR", "MotifSize"}}
            print("Filter out known disease-associated locus: ", d)
        df = df[df["MatchedKnownDiseaseAssociatedLocus"].isna()]
        print(f"Filter out known disease-associated loci: filtered out {before - len(df):,d} out of {before:,d} rows ({(before - len(df)) / (before or 1):.1%}%)")

    # Filter by --exclude-matched-tr-from
    if args.exclude_matched_tr_from:
        for exclude_matched_tr_from in args.exclude_matched_tr_from:
            before = len(df)
            df = df[df[f"MatchedTrFrom:{exclude_matched_tr_from}"].isna()]
            print(f"Filter out MatchedTrFrom:{exclude_matched_tr_from} filtered out {before - len(df):,d} out of {before:,d} rows ({(before - len(df)) / (before or 1):.1%}%)")

    # Filter by --exclude-overlaps-with
    if args.exclude_overlaps_with:
        for exclude_overlaps_with in args.exclude_overlaps_with:
            before = len(df)
            df = df[df[f"OverlapsWith:{exclude_overlaps_with}"]]
            print(f"Filter out OverlapsWith:{exclude_overlaps_with} filtered out {before - len(df):,d} out of {before:,d} rows ({(before - len(df)) / (before or 1):.1%}%)")

    if args.output_prefix is None:
        args.output_prefix = re.sub(".tsv(.gz)?$", "", args.input_tsv[0])

    # Generate per-sample per-locus table
    output_tsv_path = f"{args.output_prefix}.filtered.tsv"
    df.to_csv(output_tsv_path, sep="\t", index=False)
    print(f"Wrote {len(df):,d} rows to {output_tsv_path}")

    # Filter to only rows with a MatchedReferenceTR
    before = len(df)
    df = df[~df["MatchedReferenceTR"].isna()]
    print(f"Filter out rows without a MatchedReferenceTR: filtered out {before - len(df):,d} out of {before:,d} rows ({(before - len(df)) / (before or 1):.1%}%)")

    df_source_sample_ids = df[["Source", "SampleId"]].drop_duplicates()
    df_source_locus_id = df[["Source", "MatchedReferenceTR"]].drop_duplicates()
    df_canonical_motif_locus_id = df[["CanonicalMotif", "MatchedReferenceTR"]].drop_duplicates()
    print(f"These rows contain: ")
    print(f"   {len(df['MatchedReferenceTR'].unique()):,d} unique loci")
    print(f"   {len(df_source_sample_ids):,d} unique samples")
    print(f"   {len(df['Source'].unique()):,d} unique source files: ")

    for _, row in df_source_locus_id.groupby("Source").count()[["MatchedReferenceTR"]].sort_values(by="MatchedReferenceTR", ascending=False).iterrows():
        print(f"      {row['MatchedReferenceTR']:>8,d} loci are from {row.name}")

    print("   These canonical motifs have 5 or more loci:")
    for _, row in df_canonical_motif_locus_id.groupby("CanonicalMotif").count()[["MatchedReferenceTR"]].sort_values(by="MatchedReferenceTR", ascending=False).iterrows():
        if row['MatchedReferenceTR'] >= 5:
            print(f"      {row['MatchedReferenceTR']:>8,d} loci have canonical motif {row.name}")

    print("   This is the number of samples from each source that have at least one outlier locus:")
    for _, row in df_source_sample_ids.groupby("Source").count()[["SampleId"]].sort_values(by="SampleId", ascending=False).iterrows():
        print(f"      {row['SampleId']:>8,d} samples are from {row.name}")


    # Group by locus
    df_loci = group_by_locus(df, args.gene_column_prefix)

    # Generate per-locus table
    output_tsv_path = f"{args.output_prefix}.filtered_and_grouped_by_locus.tsv"
    df_loci.to_csv(output_tsv_path, sep="\t", index=False)
    print(f"Wrote {len(df_loci):,d} rows to {output_tsv_path}")

    # Generate ExpansionHunter variant catalog
    variant_catalog = []
    variant_catalog_output_path = f"{args.output_prefix}.variant_catalog.json"
    for _, row in df_loci.iterrows():
        locus_id = row["MatchedReferenceTR"]
        chrom, start_0based, end_1based, repeat_unit = locus_id.split("-")

        variant_catalog_record = {
            "LocusId": locus_id,
            "ReferenceRegion": f"{chrom}:{start_0based}-{end_1based}",
            "LocusStructure": f"({repeat_unit})*",
            "VariantType": "Repeat",

            # extra fields that are ignored by ExpansionHunter:
            "RepeatUnit": repeat_unit,
            "GeneName": row[f"{args.gene_column_prefix}GeneName"],
            "GeneId": row[f"{args.gene_column_prefix}GeneId"],
            "TranscriptId": row[f"{args.gene_column_prefix}TranscriptId"],
            "GeneRegion": row[f"{args.gene_column_prefix}GeneRegion"],
        }

        if len(row["SampleId"].split(",")) < 10:
            variant_catalog_record.update({
                "Samples": row["SampleId"].split(","),
                "SampleRanks": list(map(int, row["SampleRankAtLocus"].split(","))),
                "Sources": row["Source"].split(","),
                "TotalSamplesAtLocus": row["TotalSamplesAtLocus"],
            })
        variant_catalog.append(variant_catalog_record)


    with open(variant_catalog_output_path, "wt") as f:
        json.dump(variant_catalog, f, indent=4)

    print(f"Wrote {len(variant_catalog):,d} records to {variant_catalog_output_path}")


if __name__ == "__main__":
    main()
