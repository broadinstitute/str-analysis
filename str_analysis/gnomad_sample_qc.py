import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from pprint import pprint
import seaborn as sns
import sys

from str_analysis.generate_gnomad_json import load_data_df, process_sample_id

os.chdir(os.path.expanduser("~/code/str-analysis"))

#%%

args = argparse.Namespace(
    expansion_hunter_tsv="./local_files/gnomad_str_data/data/without_offtargets/combined_expansion_hunter.19243_json_files.variants.tsv",
    expansion_hunter_tsv_using_offtargets="./local_files/gnomad_str_data/data/with_offtargets/combined_expansion_hunter.19243_json_files.variants.tsv",
    non_ref_motif_tsv="./local_files/gnomad_str_data/data/without_offtargets/combined.173160_json_files.tsv",
    non_ref_motif_tsv_using_offtargets="./local_files/gnomad_str_data/data/with_offtargets/combined.173187_json_files.tsv",
    gnomad_metadata_tsv="~/code/sample_metadata/metadata/gnomad_v3.1_metadata_v3.1_with_read_lengths.tsv.gz",
    known_pathogenic_strs_tsv="./str_analysis/data/known_pathogenic_str_loci_metadata_hg38_hg19.tsv",
    existing_readviz_filename_list="./local_files/gnomad_str_data/all_renamed_svg_image_paths.txt",
    gnomad_cram_paths_tsv="~/code/sample_metadata/metadata/gnomAD_v3_and_v3_1__all_release_samples_cram_paths_and_metadata_with_150bp_read_lengths.tsv",
    hgdp_cram_paths="gs://ibd-external-datasets/HGDP/broad_reprocessed_crams/*.cram",
    include_all_age_and_pcr_info=True,
    no_readviz=True,
)

gnomad_str_calls_df = load_data_df(args)
pprint(list(gnomad_str_calls_df.columns))

"""
['SampleId',
 'LocusId',
 'VariantId',
 'ReferenceRegion',
 'Motif: Allele 1',
 'Motif: Allele 2',
 'Num Repeats: Allele 1',
 'Num Repeats: Allele 2',
 'Genotype',
 'GenotypeConfidenceInterval',
 'RepeatUnit',
 'ReadvizFilename',
 'Genotype__UsingOfftargetRegions',
 'GenotypeConfidenceInterval__UsingOfftargetRegions',
 'Num Repeats: Allele 1__UsingOfftargetRegions',
 'Num Repeats: Allele 2__UsingOfftargetRegions',
 'population_inference.pop',
 'sex_imputation.sex_karyotype',
 'age',
 'pcr_protocol',
 'read_length']
 """

#%%

gnomad_metadata_df = pd.read_table(args.gnomad_metadata_tsv)
gnomad_metadata_df["s"] = gnomad_metadata_df["project_meta.sample_id"].str.replace("-", "_").str.replace(" ", "_")
gnomad_metadata_df = gnomad_metadata_df[["s", "project_meta.neuro_case"]]

#%%
gnomad_str_calls_df = gnomad_str_calls_df.set_index("SampleId").join(
    gnomad_metadata_df.set_index("s"), how="left").reset_index().rename(columns={"index": "SampleId"})

#%%

# Downloaded from https://broadinstitute.slack.com/archives/C05DBQMJYSZ/p1687525486024799
sv_team_pcr_status_tsv = "./local_files/gnomAD-SV_v3.master_sample_metadata.pre_00c_qc.gatksv_sample_id.inferred_pcr_status"
sv_team_pcr_status_df = pd.read_table(sv_team_pcr_status_tsv)
sv_team_pcr_status_df = sv_team_pcr_status_df[["external_sample_id", "inferred_pcr_status"]]
sv_team_pcr_status_df["external_sample_id"] = sv_team_pcr_status_df.external_sample_id.apply(process_sample_id)
sv_team_pcr_status_df["inferred_pcr_status"] = sv_team_pcr_status_df["inferred_pcr_status"].replace({
    "PCRMINUS": "pcr_free",
    "PCRPLUS": "pcr_plus",
})

#%%

gnomad_str_calls_df = gnomad_str_calls_df.set_index("SampleId").join(
    sv_team_pcr_status_df.set_index("external_sample_id"),
    how="left").reset_index().rename(columns={"index": "SampleId"})

#%%
# mark HGDP samples as pcr_free, since they were sequenced with a PCR-free protocol
gnomad_str_calls_df.loc[gnomad_str_calls_df.SampleId.str.startswith("HGDP"), "pcr_protocol"] = "pcr_free"
# fill in missing PCR status with inferred PCR status
gnomad_str_calls_df.loc[gnomad_str_calls_df.pcr_protocol == "pcr_info_not_available", "pcr_protocol"] = np.nan
gnomad_str_calls_df.pcr_protocol.fillna(gnomad_str_calls_df.inferred_pcr_status, inplace=True)
# discard the 3 remaining samples without PCR status
assert len(set(gnomad_str_calls_df[gnomad_str_calls_df.pcr_protocol.isna()].SampleId)) < 5
gnomad_str_calls_df = gnomad_str_calls_df[~gnomad_str_calls_df.pcr_protocol.isna()]

#%%

# write out a gnomAD metadata table for PCR-free samples

gnomad_cram_paths_df = pd.read_table(args.gnomad_cram_paths_tsv)
gnomad_cram_paths_df = gnomad_cram_paths_df[["sample_id", "cram_path", "crai_path"]]
gnomad_cram_paths_df["sample_id"] = gnomad_cram_paths_df.sample_id.apply(process_sample_id)

gnomad_metadata_df = gnomad_str_calls_df[[
    'SampleId',
    'population_inference.pop', 'sex_imputation.sex_karyotype', 'age',
    'pcr_protocol', 'read_length', 'project_meta.neuro_case',
]].drop_duplicates()

gnomad_combined_metadata_df = gnomad_cram_paths_df.set_index("sample_id").join(
    gnomad_metadata_df.set_index("SampleId"), how="left").reset_index()

gnomad_combined_metadata_df = gnomad_combined_metadata_df.set_index("sample_id").join(
    sv_team_pcr_status_df.set_index("external_sample_id"), how="left").reset_index().rename(columns={"index": "sample_id"})

gnomad_combined_metadata_df.loc[gnomad_combined_metadata_df["sample_id"].str.startswith("HGDP"), "pcr_protocol"] = "pcr_free"

# fill in missing PCR status with inferred PCR status
gnomad_combined_metadata_df.loc[gnomad_combined_metadata_df["pcr_protocol"] == "pcr_info_not_available", "pcr_protocol"] = np.nan
gnomad_combined_metadata_df.pcr_protocol.fillna(gnomad_combined_metadata_df["inferred_pcr_status"], inplace=True)

# discard the 3 remaining samples without PCR status
assert len(set(gnomad_combined_metadata_df[gnomad_combined_metadata_df["pcr_protocol"].isna()]["sample_id"])) < 5
gnomad_combined_metadata_df = gnomad_combined_metadata_df[~gnomad_combined_metadata_df["pcr_protocol"].isna()]
gnomad_combined_metadata_df = gnomad_combined_metadata_df.drop(columns=["inferred_pcr_status"])
gnomad_combined_metadata_df["sex"] = gnomad_combined_metadata_df['sex_imputation.sex_karyotype'].apply(lambda x:
    "M" if x == "XY" else ("F" if x == "XX" else None))

gnomad_combined_metadata_df.to_csv(
    "~/code/sample_metadata/metadata/gnomad_v3.1_metadata_v3.1_with_inferred_pcr_status.tsv.gz",
    sep="\t", header=True, index=False)
gnomad_combined_metadata_df[gnomad_combined_metadata_df.pcr_protocol == "pcr_free"].to_csv(
    "~/code/sample_metadata/metadata/gnomad_v3.1_metadata_v3.1_with_inferred_pcr_status.only_pcr_free.tsv.gz",
    sep="\t", header=True, index=False)

#%%

import hailtop.fs as hfs
hgdp_crams = {x["path"] for x in hfs.ls(args.hgdp_cram_paths)}


#%%

loci_to_include_in_sample_qc = {
    #'AFF2',        # chrX
    #'AR',          # chrX
    #'ARX_1',       # contains N in motif
    #'ARX_2',       # contains N in motif
    'ATN1',
    'ATXN1',
    'ATXN10',
    'ATXN2',
    'ATXN3',
    'ATXN7',
    'ATXN8OS',
    ##'BEAN1',       # pentanucleotide locus with non-ref motifs
    'C9ORF72',
    'CACNA1A',
    'CNBP',
    'COMP',
    ###'CSTB',       # 12 nucleotide motif
    ##'DAB1',        # pentanucleotide locus with non-ref motifs
    'DIP2B',
    #'DMD',          # chrX
    'DMPK',
    ###'EIF4A3',     # 20 nucleotide motif
    #'FMR1',         # chrX
    #'FOXL2',         # contains N in motif
    'FXN',
    'GIPC1',
    'GLS',
    #'HOXA13_1',     # contains N in motif
    #'HOXA13_2',     # contains N in motif
    #'HOXA13_3',     # contains N in motif
    'HOXD13',
    'HTT',
    'JPH3',
    'LRP12',
    ##'MARCHF6',     # pentanucleotide locus with non-ref motifs
    'NIPA1',
    'NOP56',
    'NOTCH2NLC',
    'NUTM2B-AS1',
    'PABPN1',
    #'PHOX2B',      # contains N in motif
    'PPP2R2B',
    'PRDM12',
    #'PRNP',        # 24 nucleotide motif
    ##'RAPGEF2',   # pentanucleotide locus with non-ref motifs
    ##'RFC1',      # pentanucleotide locus with non-ref motifs
    'RILPL1',
    #'RUNX2',      # contains N in motif
    ##'SAMD12',    # pentanucleotide locus with non-ref motifs
    #'SOX3',       # contains N in motif
    ##'STARD7',    # pentanucleotide locus with non-ref motifs
    'TBP',
    #'TBX1',       # contains N in motif
    'TCF4',
    ##'TNRC6A',    # pentanucleotide locus with non-ref motifs
    ###'VWA1',     # 10 nucleotide motif
    ####'XYLT1',   # this locus has unique features that make it not work well
    ##'YEATS2',    # pentanucleotide locus with non-ref motifs
    #'ZIC2',       # contains N in motif
    #'ZIC3',       # contains N in motif
}

loci_to_include_in_sample_qc = {
    'ATN1',
    'ATXN1',
    'ATXN10',
    'ATXN2',
    'ATXN3',
    'ATXN7',
    'ATXN8OS',
    'C9ORF72',
    'CACNA1A',
    'CNBP',
    'COMP',
    'DIP2B',
    'DMPK',
    'FXN',
    'GIPC1',
    'GLS',
    'HOXD13',
    'HTT',
    'JPH3',
    'LRP12',
    'NIPA1',
    'NOP56',
    'NOTCH2NLC',
    'NUTM2B-AS1',
    'PABPN1',
    'PPP2R2B',
    'PRDM12',
    'RILPL1',
    'TBP',
    'TCF4',
}

print(len(loci_to_include_in_sample_qc), "loci to include in sample QC:", loci_to_include_in_sample_qc)


#%%

gnomad_str_calls_df.loc[:, "LociToIncludeInSampleQc"] = gnomad_str_calls_df.VariantId.isin(loci_to_include_in_sample_qc)

assert len(set(gnomad_str_calls_df[gnomad_str_calls_df.LociToIncludeInSampleQc].VariantId)) == len(loci_to_include_in_sample_qc)

#%%

# calculate genotype quality using approximately the formula in https://github.com/gymrek-lab/EnsembleTR/blob/main/ensembletr/utils.py#L36

gnomad_str_calls_df[['ConfidenceInterval: Allele 1', 'ConfidenceInterval: Allele 2']] = gnomad_str_calls_df['GenotypeConfidenceInterval'].str.split("/", expand=True)

gnomad_str_calls_df[['ConfidenceIntervalLowerBound: Allele 1', 'ConfidenceIntervalUpperBound: Allele 1']] = gnomad_str_calls_df['ConfidenceInterval: Allele 1'].str.split("-", expand=True).astype(float)
gnomad_str_calls_df[['ConfidenceIntervalLowerBound: Allele 2', 'ConfidenceIntervalUpperBound: Allele 2']] = gnomad_str_calls_df['ConfidenceInterval: Allele 2'].str.split("-", expand=True).astype(float)

gnomad_str_calls_df['ConfidenceIntervalSize: Allele 1'] = gnomad_str_calls_df['ConfidenceIntervalUpperBound: Allele 1'] - gnomad_str_calls_df['ConfidenceIntervalLowerBound: Allele 1']
gnomad_str_calls_df['ConfidenceIntervalSize: Allele 2'] = gnomad_str_calls_df['ConfidenceIntervalUpperBound: Allele 2'] - gnomad_str_calls_df['ConfidenceIntervalLowerBound: Allele 2']

gnomad_str_calls_df['GenotypeQuality: Allele 1'] = np.where(
    gnomad_str_calls_df['Num Repeats: Allele 1'].astype(float) > 0,
    1/np.exp(gnomad_str_calls_df['ConfidenceIntervalSize: Allele 1']/gnomad_str_calls_df['Num Repeats: Allele 1'].astype(float)),
    float('nan'))
gnomad_str_calls_df['GenotypeQuality: Allele 2'] = np.where(
    gnomad_str_calls_df['Num Repeats: Allele 2'].astype(float) > 0,
    1/np.exp(gnomad_str_calls_df['ConfidenceIntervalSize: Allele 2']/gnomad_str_calls_df['Num Repeats: Allele 2'].astype(float)),
    float('nan'))

gnomad_str_calls_df['GenotypeQuality'] = 0.8 * gnomad_str_calls_df[['GenotypeQuality: Allele 1', 'GenotypeQuality: Allele 2']].min(axis=1) + 0.2 * gnomad_str_calls_df[['GenotypeQuality: Allele 1', 'GenotypeQuality: Allele 2']].max(axis=1)

#%%
gnomad_str_calls_df["Num Repeats: Allele 1"] = gnomad_str_calls_df["Num Repeats: Allele 1"].astype(float)
gnomad_str_calls_df["Num Repeats: Allele 2"] = gnomad_str_calls_df["Num Repeats: Allele 2"].astype(float)
gnomad_str_calls_df['Num Repeats: Allele 1 + 2'] = gnomad_str_calls_df['Num Repeats: Allele 1'].fillna(0) + gnomad_str_calls_df['Num Repeats: Allele 2'].fillna(0)

#%%
# calculate quartile rank for each sample across loci that work

# compute how far each sample's 'Num Repeats: Allele 2' is from the mean of the population for that LocusId
df_for_qc = gnomad_str_calls_df
n_rows_in_df_for_qc = len(df_for_qc)
#df_for_qc = df_for_qc[df_for_qc.VariantId.isin(loci_to_include_in_sample_qc)]
print(f"Filtered out {n_rows_in_df_for_qc - len(df_for_qc):,d} out of {n_rows_in_df_for_qc:,d} rows "
      f"({(n_rows_in_df_for_qc - len(df_for_qc))/n_rows_in_df_for_qc*100:.2f}%): from loci that were not in loci_to_include_in_sample_qc")
n_rows_in_df_for_qc = len(df_for_qc)
df_for_qc = df_for_qc[~df_for_qc["Num Repeats: Allele 1"].isna() & ~df_for_qc["Num Repeats: Allele 2"].isna()]
print(f"Filtered out {n_rows_in_df_for_qc - len(df_for_qc):,d} out of {n_rows_in_df_for_qc:,d} rows "
      f"({(n_rows_in_df_for_qc - len(df_for_qc))/n_rows_in_df_for_qc*100:.2f}%): that were missing Num Repeats: "
      f"Allele 1 or Num Repeats: Allele 2")


#df_by_locus_id = df_without_na.reset_index().set_index("LocusId")
#df_without_na["Num Repeats: Allele 2: locus mean"] = df_grpby_locus_id["Num Repeats: Allele 2"].mean()
#df_without_na["Num Repeats: Allele 2: locus std"] = df_grpby_locus_id["Num Repeats: Allele 2"].std()
#df_by_locus_id["Num Repeats: Allele 2: locus median"] = df_grpby_locus_id["Num Repeats: Allele 2"].median()
#df_without_na["Num Repeats: Allele 2: norm"] = (df_without_na["Num Repeats: Allele 2"] - df_without_na["Num Repeats: Allele 2: locus mean"])/df_without_na["Num Repeats: Allele 2: locus std"]
#df_without_na["Num Repeats: Allele 2: outlier_by_stdev"] = df_without_na["Num Repeats: Allele 2: norm"] >= 7

# add pathogenic thresholds
df_pathogenic = pd.read_table(args.known_pathogenic_strs_tsv)
df_pathogenic = df_pathogenic[["LocusId", "PathogenicMin", "Gene Region", "RU", "InheritanceMode", "Disease"]].set_index("LocusId")
df_for_qc = df_for_qc.set_index("LocusId").join(df_pathogenic, how="left").reset_index()

# compute rank of each sample/locus pair for that locus
df_for_qc = df_for_qc.assign(allele2_percentile=df_for_qc.groupby("LocusId")['Num Repeats: Allele 2'].rank(pct=True).mul(100))
df_for_qc = df_for_qc.assign(allele1plus2_percentile=df_for_qc.groupby("LocusId")['Num Repeats: Allele 1 + 2'].rank(pct=True).mul(100))

df_for_qc["Num Repeats: Allele 2: outlier_by_percentile"] = df_for_qc.allele2_percentile > 99
df_for_qc["Num Repeats: Allele 2: above_pathogenic_threshold"] = df_for_qc["Num Repeats: Allele 2"] >= df_for_qc["PathogenicMin"]
df_for_qc["GenotypeQuality: low_genotype_quality"] = df_for_qc["GenotypeQuality"] < 0.2

# sum the percentile ranks across loci
df_for_qc_groupby_sample_id = df_for_qc.groupby("SampleId")
df_for_qc = df_for_qc.set_index("SampleId")
df_for_qc = df_for_qc.assign(allele2_percentile_sum_across_loci=df_for_qc_groupby_sample_id["allele2_percentile"].sum())
df_for_qc = df_for_qc.assign(allele1plus2_percentile_sum_across_loci=df_for_qc_groupby_sample_id["allele1plus2_percentile"].sum())
df_for_qc = df_for_qc.assign(genotype_quality_sum=df_for_qc_groupby_sample_id["GenotypeQuality"].sum())
df_for_qc = df_for_qc.assign(genotype_quality_mean=df_for_qc_groupby_sample_id["GenotypeQuality"].mean())
df_for_qc = df_for_qc.assign(genotype_quality_median=df_for_qc_groupby_sample_id["GenotypeQuality"].median())
df_for_qc = df_for_qc.assign(sample_is_outlier_at_this_many_loci=df_for_qc_groupby_sample_id["Num Repeats: Allele 2: outlier_by_percentile"].sum())
df_for_qc = df_for_qc.assign(sample_is_pathogenic_count=df_for_qc_groupby_sample_id["Num Repeats: Allele 2: above_pathogenic_threshold"].sum())
df_for_qc = df_for_qc.assign(low_genotype_quality_count=df_for_qc_groupby_sample_id["Num Repeats: Allele 2: above_pathogenic_threshold"].sum())
df_for_qc = df_for_qc.reset_index()


#%%
output_path = "gnomAD_sample_qc__2023_06_22.tsv"
df_for_qc.to_csv(output_path, sep="\t", index=False, header=True)
print(f"Wrote {len(df_for_qc):,d} rows for {len(set(df_for_qc.VariantId)):,d} loci and "
      f"{len(set(df_for_qc.SampleId)):,d} samples to {os.path.abspath(output_path)}")
print("Columns:", list(df_for_qc.columns))

#%%
sys.exit(0)

#df_for_qc[df_for_qc.pcr_protocol == "pcr_info_not_available"][["SampleId"]].sort_values(by="SampleId").drop_duplicates().to_csv("gnomad_samples_without_pcr_info.tsv", sep="\t", index=False, header=True)

#%%
#%%


#%%
plt.close()
ax = sns.histplot(data=gnomad_str_calls_df[gnomad_str_calls_df.VariantId == "PABPN1"], x='Num Repeats: Allele 2', hue="pcr_protocol", multiple="dodge");
ax.set_yscale("symlog", linthresh=10)
plt.show()

#%%

plt.close()
ax = sns.histplot(data=gnomad_str_calls_df[gnomad_str_calls_df.VariantId == "PABPN1"], x='Num Repeats: Allele 2', hue="pcr_protocol", multiple="dodge");
ax.set_yscale("symlog", linthresh=10)
plt.show()

#%%

plt.close()
ax = sns.histplot(data=gnomad_str_calls_df[gnomad_str_calls_df.VariantId == "PABPN1"], x='Num Repeats: Allele 2', hue="project_meta.project_id", multiple="dodge");
ax.set_yscale("symlog", linthresh=10)
plt.show()

"""
'project_meta.project_id',
 'project_meta.v2_project_name',
 'project_meta.product',
 'project_meta.pdo',
 'project_meta.title',
 'bam_metrics.median_coverage'
"""
#%%

"project_meta.project_id"
#%%

for column in [
    "pcr_protocol", # "read_length", "bam_metrics.median_coverage",
    "project_meta.project_id",
    "project_meta.v2_project_name",
    "project_meta.product",
    "project_meta.pdo",
    "project_meta.title",
]:
    column_values = sorted(list(set(gnomad_str_calls_df[column].fillna("NA").astype(str))))
    print(f"{len(column_values):3d} {column} values: ({100 * sum(gnomad_str_calls_df[column].isna()) / len(gnomad_str_calls_df):0.1f}% are N/A)", ", ".join(column_values[:5]))
    # print percent N/A

#%%

print(gnomad_str_calls_df.groupby("project_meta.title").count().sort_values("SampleId", ascending=False)["SampleId"] / len(set(gnomad_str_calls_df.VariantId)))

#%%

#
#%%

df_plot = df_without_na.reset_index()
plt.close()
#ax = sns.histplot(data=df_plot[df_plot.VariantId == "HTT"], x='Num Repeats: Allele 2', hue="Num Repeats: Allele 2: outlier_by_percentile", multiple="stack");
#ax = sns.histplot(data=df_plot[df_plot.VariantId == "PABPN1"], x='Num Repeats: Allele 2', hue="Num Repeats: Allele 2: outlier_by_percentile", multiple="stack");
ax = sns.histplot(data=df_plot[df_plot.VariantId == "ATXN8OS"], x='Num Repeats: Allele 2', hue="Num Repeats: Allele 2: outlier_by_percentile", multiple="stack");

ax.set_yscale("symlog", linthresh=10)
plt.show()

#%%


#%%