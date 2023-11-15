#!/usr/bin/env python3

"""
This script takes ExpansionHunter calls and metadata files and reprocesses them into tables and .json files
that can be loaded into the gnomAD browser, and also shared publicly through the gnomAD downloads page.

The script is intended to be run twice:
1. The 1st time, it is run without the --existing-readviz-filename-list arg so it can generate the
readviz_rename_list__...tsv file which maps regular readviz image filenames like NA12878.ATXN1.svg to encrypted ones
like  ffa0880117e0791d51b0ef85b56f3a54216.svg
2. The readviz_rename_list__...tsv output file can then be used to rename readviz images to their encrypted filenames
so that these image files can be made public without revealing sample ids in the filenames.
3. To handle any missing readviz images (such as those where REViewer crashed) all renamed image files can then be
listed in a text file - lets call it "all_renamed_readviz_files.txt"
4. This script can then be run a 2nd time, with  --existing-readviz-filenames-list "all_renamed_readviz_files.txt"
in order to exclude missing readviz image paths from output files, thereby letting the gnomAD browser know which
images are missing.

Output documentation:

======================================================
==     gnomAD_STR_genotypes__[datestamp].tsv.gz     ==
======================================================

This is a flat table that contains the ExpansionHunter genotypes from
all samples at each of the 59 disease-associated loci. Also, it contains Population, Sex, Age, and PcrProtocol
metadata columns, so the data in this table can be used to generate any of the plots displayed in the gnomAD
browser STR pages. It also contains results not currently available through the browser.
These are in the GenotypeUsingOfftargetRegions and GenotypeConfidenceIntervalUsingOfftargetRegions columns
which store ExpansionHunter calls generated using off-target regions. Finally, the ReadvizFilename column
links each row to a REViewer read visualization image that is available through the browser. This should allow users
to construct the full public url of the image and programmatically download specific images of interest.

Below is an example of column names and values from a typical row in this table:

                                              Id : PABPN1
                                         LocusId : PABPN1
                                 ReferenceRegion : chr14:23321472-23321490
                                           Chrom : chr14
                                    Start_0based : 23321472
                                             End : 23321490
                                           Motif : GCG
                                IsAdjacentRepeat : False
                                      Population : nfe
                                             Sex : XY
                                             Age : 20-25
                                     PcrProtocol : pcr_free
                                        Genotype : 6/13
                                         Allele1 : 6
                                         Allele2 : 13
                      GenotypeConfidenceInterval : 6-6/11-30
                   GenotypeUsingOfftargetRegions : 6/13
                    Allele1UsingOfftargetRegions : 6
                    Allele2UsingOfftargetRegions : 13
 GenotypeConfidenceIntervalUsingOfftargetRegions : 6-6/11-30
                                 ReadvizFilename : c82034cf2aad813a07a8523898d64c81148.svg

Id : This id is unique to each STR locus and repeat, meaning that it differs between repeats and any adjacent repeats
    at a locus. For example the main GAA repeat at the FXN locus has id "FXN" while the adjacent poly-A repeat has id
    "FXN_A". This id corresponds to the "VariantId" field in the ExpansionHunter variant catalogs @
    https://github.com/broadinstitute/str-analysis/tree/main/str_analysis/variant_catalogs
LocusId: This id is unique to each STR locus. It corresponds to the "LocusId" field in the ExpansionHunter variant
    catalogs @ https://github.com/broadinstitute/str-analysis/tree/main/str_analysis/variant_catalogs
    and can be used to look up reference information about each locus there - including the
    gene name, disease associations, mode of inheritance, and pathogenic threshold. For most but not all loci, the
    LocusId is identical to the name of the gene that contains the locus.
ReferenceRegion: Genomic interval delineating the exact boundaries of the STR repeat in the reference genome.
    The start coordinate is 0-based.
Chrom: The chromosome of the ReferenceRegion. This is provided as a separate column for convenience.
Start_0based: The 0-based start coordinate of the ReferenceRegion. This is provided as a separate column for convenience.
End: The end coordinate of the ReferenceRegion. This is provided as a separate column for convenience.
Motif: The repeat unit of the STR locus. For example this would be "GAA" at the FXN locus.
IsAdjacentRepeat: True or False depending on whether this row represents the main repeat at a locus or an adjacent repeat.
    Adjacent repeats are included for some loci either for technical reasons to improve ExpansionHunter accuracy, or
    due to research interest in the size of these adjacent repeats.
Population: The gnomAD ancestry group of the individual. Possible values are: "afr", "ami", "amr", "asj", "eas", "fin",
    "mid", "nfe", "oth", "sas"
Sex: The sex karyotype of the genotyped individual. Possible values are "XX" and "XY".
Age: The age of the individual at the time when they enrolled in one of the research studies underlying gnomAD.
    The values represent 5 year bins such as "20-25", as well as ">80" for individuals over 80 years old and "<20" for
    individuals younger than 20. For individuals with unknown age, the value is "age_not_available"
PcrProtocol: Possible values are "pcr_free", "pcr_plus" and "pcr_info_not_available"
Genotype: The ExpansionHunter genotype for this individual at this locus, generated using the variant catalog without
    off-target regions (see https://github.com/broadinstitute/str-analysis/tree/main/str_analysis/variant_catalogs).
    These are the genotypes used to generate all plots in the gnomAD browser STR pages.
Allele1: The shorter repeat size from the genotype. This is provided as a separate column for convenience.
Allele2: The longer repeat size from the genotype; empty in the special case of hemizygous genotypes (e.g., in male samples
    at loci on chrX). This is provided as a separate column for convenience.
GenotypeConfidenceInterval: The ExpansionHunter confidence intervals associated with the genotype in the "Genotype" column.

GenotypeUsingOfftargetRegions: Same meaning as the "Genotype" column, but generated using the variant catalog with off-target regions.
Allele1UsingOfftargetRegions: Same meaning as the "Allele1" column, but generated using the variant catalog with off-target regions.
Allele2UsingOfftargetRegions: Same meaning as the "Allele2" column, but generated using the variant catalog with off-target regions.
GenotypeConfidenceIntervalUsingOfftargetRegions: Same meaning as the "GenotypeConfidenceInterval" column, but generated using the
    variant catalog with off-target regions.

ReadvizFilename: The filename of the SVG image generated by REViewer based on the ExpansionHunter call reported in the "Genotype"
    column. This can be used to compute the public url and programatically retrieve this image from the gnomAD browser server.

"""


import argparse
import collections
from datetime import datetime
import gzip
import hashlib
import simplejson as json
import math
import os
import pandas as pd
import pkgutil
import pwd
import requests
import tqdm

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif

# Map STR locus ids to readable names for STR loci that are adjacent to the main known pathogenic loci
from str_analysis.utils.export_json import export_json
from str_analysis.utils.known_pathogenic_strs_tsv import parse_known_pathogenic_strs_tsv
from str_analysis.utils.misc_utils import parse_interval

ADJACENT_REPEAT_LABELS = {
    "ATXN7_GCC": "Adjacent Right STR",
    "ATXN8OS_CTA": "Adjacent Left STR",
    "HTT_CCG": "Adjacent Right STR",
    "FXN_A": "Adjacent Left Homopolymer",
    "CNBP_CAGA": "Adjacent Right STR #1",
    "CNBP_CA": "Adjacent Right STR #2",
    "NOP56_CGCCTG": "Adjacent Right STR",
    "PRNP_CCTCAGGGCGGTGGTGGCTGGGGGCAG": "Adjacent Left STR",
}

# Map gene name to Ensembl gene id for genes that contain known pathogenic STRs
GENE_NAME_TO_GENE_ID = {
    'ATXN8': 'ENSG00000230223',
    'AFF2': 'ENSG00000155966',
    'AR': 'ENSG00000169083',
    'ARX': 'ENSG00000004848',
    'ATN1': 'ENSG00000111676',
    'ATXN1': 'ENSG00000124788',
    'ATXN10': 'ENSG00000130638',
    'ATXN2': 'ENSG00000204842',
    'ATXN3': 'ENSG00000066427',
    'ATXN7': 'ENSG00000163635',
    'BEAN1': 'ENSG00000166546',
    'C9orf72': 'ENSG00000147894',
    'CACNA1A': 'ENSG00000141837',
    'CBL2': 'ENSG0000011039',
    'CNBP': 'ENSG00000169714',
    'COMP': 'ENSG00000105664',
    'CSTB': 'ENSG00000160213',
    'DAB1': 'ENSG00000173406',
    'DIP2B': 'ENSG00000066084',
    'DMD': 'ENSG00000198947',
    'DMPK': 'ENSG00000104936',
    'EIF4A3': 'ENSG00000141543',
    'FMR1': 'ENSG00000102081',
    'FOXL2': 'ENSG00000183770',
    'FXN': 'ENSG00000165060',
    'GIPC1': 'ENSG00000123159',
    'GLS': 'ENSG00000115419',
    'HOXA13': 'ENSG00000106031',
    'HOXD13': 'ENSG00000128714',
    'HTT': 'ENSG00000197386',
    'JPH3': 'ENSG00000154118',
    'LOC642361': 'ENSG00000272447',
    'LRP12': 'ENSG00000147650',
    'MARCHF6': 'ENSG00000145495',
    'NIPA1': 'ENSG00000170113',
    'NOP56': 'ENSG00000101361',
    'NOTCH2NLC': 'ENSG00000286219',
    'PABPN1': 'ENSG00000100836',
    'PHOX2B': 'ENSG00000109132',
    'PPP2R2B': 'ENSG00000156475',
    'PRDM12': 'ENSG00000130711',
    'PRNP': 'ENSG00000171867',
    'RAPGEF2': 'ENSG00000109756',
    'RILPL1': 'ENSG00000188026',
    'RFC1': 'ENSG00000035928',
    'RUNX2': 'ENSG00000124813',
    'SAMD12': 'ENSG00000177570',
    'SOX3': 'ENSG00000134595',
    'STARD7': 'ENSG00000084090',
    'TBP': 'ENSG00000112592',
    'TBX1': 'ENSG00000184058',
    'TCF4': 'ENSG00000196628',
    'TNRC6A': 'ENSG00000090905',
    'VWA1': 'ENSG00000179403',
    'XYLT1': 'ENSG00000103489',
    'YEATS2': 'ENSG00000163872',
    'ZIC2': 'ENSG00000043355',
    'ZIC3': 'ENSG00000156925',
}

# Discard samples with read length below the MIN_READ_LENGTH threshold (in base pairs)
MIN_READ_LENGTH = 150

# Round ages to the nearest N years so that they can be shared publicly without increasing identifiability
AGE_RANGE_SIZE = 5

# Truncate the age distribution at this lower and upper bound.
LOWER_AGE_CUTOFF = 20
UPPER_AGE_CUTOFF = 80

# Show age for not more than this many of the most expanded samples per locus for each sex/population bucket
MAX_AGES_PER_BUCKET_TO_DISPLAY_IN_THE_READVIZ_SECTION = 100

# Use this value instead of the age range for samples where age is not available or not shown.
AGE_NOT_AVAILABLE = "age_not_available"

PCR_INFO_NOT_AVAILABLE = "pcr_info_not_available"

# Show age only for the these larger sub-populations to avoid increasing identifiability in smaller populations
POPULATIONS_WITH_AGE_DISPLAYED_IN_READVIZ_SECTION = {"sas", "oth", "asj", "amr", "fin", "eas", "afr", "nfe", "mid"}

# Fraction of genotypes that can be missing for a locus before generating an error
MISSING_GENOTYPES_ERROR_THRESHOLD = 0.01

# Fraction of readviz images that can be missing for a locus before generating an error
MISSING_READVIZ_ERROR_THRESHOLD = 0.04

# Fraction of samples that can be missing age information before generating an error
MISSING_AGE_THRESHOLD = 0.5

# Fraction of samples that can be pcr-free/pcr-plus information before generating an error
MISSING_PCR_PROTOCOL_THRESHOLD = 0.25

# Expected number of known pathogenic repeats
EXPECTED_N_KNOWN_PATHOGENIC_REPEATS = 60

# gnomAD projects for which it's ok to release individual-level genotype data
PUBLIC_PROJECT_NAMES = {"1000 Genomes Project", "Human Genome Diversity Project"}

# Add this "salt" value to the sha512 hash to prevent dictionary attacks on the encrypted sample ids
salt = pwd.getpwuid(os.getuid()).pw_name


def parse_args():
    """Parse command-line args, perform basic validation, and then return the args object."""

    p = argparse.ArgumentParser()
    p.add_argument(
        "--expansion-hunter-tsv",
        help="Table generated by running python3 -m str_analysis.combine_expansionhunter_json_to_tsv on all samples "
             "called by ExpansionHunter."
    )
    p.add_argument(
        "--non-ref-motif-tsv",
        help="Table generated by running python3 -m str_analysis.combine_json_to_tsv on all loci called by "
             "str_analysis.call_non_ref_motifs.",
    )
    p.add_argument(
        "--expansion-hunter-tsv-using-offtargets",
        help="Table generated by running python3 -m str_analysis.combine_expansionhunter_json_to_tsv on all samples "
             "called by ExpansionHunter while using a variant catalog that includes off-target regions."
    )
    p.add_argument(
        "--non-ref-motif-tsv-using-offtargets",
        help="Table generated by running python3 -m str_analysis.combine_json_to_tsv on all loci called by "
             "str_analysis.call_non_ref_motifs with options set to use a variant catalog that "
             "includes off-target regions.",
    )
    p.add_argument(
        "--gnomad-metadata-tsv",
        default="~/code/sample_metadata/metadata/gnomad_v3.1_metadata_v3.1_with_read_lengths.tsv.gz",
        help="gnomAD metadata table path.",
    )
    p.add_argument(
        "--known-pathogenic-strs-tsv",
        default="./data/known_pathogenic_str_loci_metadata_hg38_hg19.tsv",
        help="Table of known pathogenic STRs.",
    )
    p.add_argument(
        "--inferred-pcr-status-table",
        default="~/code/str-analysis/local_files/gnomAD-SV_v3.master_sample_metadata.pre_00c_qc.gatksv_sample_id.inferred_pcr_status",
        help="Table of inferred PCR+/PCR- labels for all samples")

    grp = p.add_mutually_exclusive_group()
    grp.add_argument(
        "--existing-readviz-filename-list",
        help="A text file that lists all readviz .svg filenames that exist (one per line). These are the encrypted "
             "public filenames that don't contain sample ids - for example: ffa0880117e0791d51b0ef85b56f3a54216.svg",
    )
    grp.add_argument(
        "--no-readviz",
        action="store_true",
        help="Use this flag to indicate that readviz was not generated for this dataset.",
    )
    p.add_argument(
        "--include-all-age-and-pcr-info",
        action="store_true",
        help="If this is specified, don't drop age and/or pcr info values from the user-friendly genotypes table in "
             "the add_histograms_and_compute_readviz_paths(..) method."
    )
    p.add_argument(
        "--output-dir",
        default="gs://gnomad-browser/STRs",
        help="Where to write output files. Supports local and Google storage (gs://) paths.",
    )
    p.add_argument(
        "--output-filename-suffix",
        default="",
        help="An optional label to append to all output filenames",
    )
    p.add_argument(
        "--local-output-dir",
        help="Optionally specify local directory to write output files to",
    )
    args = p.parse_args()

    print("Args:")
    for key, value in sorted(args.__dict__.items()):
        print(f"     {key} = {value}")

    for path in args.expansion_hunter_tsv, args.non_ref_motif_tsv, args.gnomad_metadata_tsv, \
                args.known_pathogenic_strs_tsv:
        if path and not os.path.isfile(os.path.expanduser(path)):
            p.error(f"{path} file not found")

    return args

def process_sample_id(sample_id):
    """Utility method for normalizing a gnomAD sample id

    Args:
        sample_id (str): gnomAD sample id
    Return:
        str: normalized gnomAD sample id
    """
    sample_id = sample_id.replace("RP-1400::", "").replace("v3.1::", "")
    return sample_id.strip().replace(" ", "_").replace("-", "_").split(".")[0].split("_SM_")[0]


def parse_expansion_hunter_tsv_and_non_ref_motif_tsv(expansion_hunter_tsv, non_ref_motif_tsv, no_readviz_images=False):
    """Parses the expansion_hunter_tsv and non_ref_motif_tsv files and returns a combined table with records
    from both.

    Args:
        expansion_hunter_tsv (str): path of either --expansion-hunter-tsv or --expansion-hunter-tsv-using-offtargets
        non_ref_motif_tsv (str): path of either --non-ref-motif-tsv or --non-ref-motif-tsv-using-offtargets
        no_readviz_images (bool): If True, this method will assume the REViewer images were not generated for this
            dataset.
    """

    def split_by_forward_slash(expansion_hunter_call_repeat_unit):
        repeat_units = expansion_hunter_call_repeat_unit.split("/")
        return repeat_units[0].strip(), repeat_units[-1].strip()

    combined_table_columns = [
        "SampleId", "LocusId", "VariantId", "ReferenceRegion",
        "Motif: Allele 1", "Motif: Allele 2",
        "Num Repeats: Allele 1", "Num Repeats: Allele 2",
        "Genotype", "GenotypeConfidenceInterval",
        "RepeatUnit", "ReadvizFilename",
    ]

    df = pd.read_table(expansion_hunter_tsv, low_memory=False)
    df.loc[:, "SampleId"] = df.SampleId.apply(process_sample_id)
    df.loc[:, "Motif: Allele 1"] = df["RepeatUnit"]
    df.loc[:, "Motif: Allele 2"] = df["RepeatUnit"]
    df.loc[:, "ReadvizFilename"] = None if no_readviz_images else df["SampleId"] + "." + df["LocusId"] + ".svg"
    df = df[combined_table_columns]

    # Parse the args.non_ref_motif_tsv generated by call_non_ref_motifs
    if non_ref_motif_tsv:
        print(f"Loading {non_ref_motif_tsv}")
        non_ref_motifs_df = pd.read_table(non_ref_motif_tsv, low_memory=False)
        non_ref_motifs_df = non_ref_motifs_df[~non_ref_motifs_df["expansion_hunter_call_genotype"].isna()]

        non_ref_motifs_df["Motif: Allele 1"], non_ref_motifs_df["Motif: Allele 2"] = zip(
            *non_ref_motifs_df["expansion_hunter_call_repeat_unit"].apply(split_by_forward_slash))

        non_ref_motifs_df.loc[:, "Num Repeats: Allele 1"], non_ref_motifs_df.loc[:, "Num Repeats: Allele 2"] = zip(
            *non_ref_motifs_df["expansion_hunter_call_genotype"].apply(split_by_forward_slash))

        non_ref_motifs_df.loc[:, "SampleId"] = non_ref_motifs_df.sample_id.apply(process_sample_id)
        non_ref_motifs_df.loc[:, "LocusId"] = non_ref_motifs_df["locus_id"]
        non_ref_motifs_df.loc[:, "VariantId"] = non_ref_motifs_df["locus_id"]
        non_ref_motifs_df.loc[:, "ReferenceRegion"] = non_ref_motifs_df["locus_coords"]
        non_ref_motifs_df.loc[:, "Genotype"] = non_ref_motifs_df["expansion_hunter_call_genotype"]
        non_ref_motifs_df.loc[:, "GenotypeConfidenceInterval"] = non_ref_motifs_df["expansion_hunter_call_CI"]
        non_ref_motifs_df.loc[:, "RepeatUnit"] = None   # will be set later
        non_ref_motifs_df.loc[:, "ReadvizFilename"] = None if no_readviz_images else non_ref_motifs_df["expansion_hunter_call_reviewer_svg"]
        non_ref_motifs_df = non_ref_motifs_df[combined_table_columns]

        df = df[~df["LocusId"].isin(set(non_ref_motifs_df["LocusId"]))]
        df = pd.concat([df, non_ref_motifs_df])

    return df


def load_data_df(args):
    """Load the tables specified by args.expansion_hunter_tsv, args.non_ref_motif_tsv, and args.gnomad_metadata_tsv.
    Rename and select relevant columns, combine the tables, then return a single combined table.

    Args:
        args (argparse.Namespace): The argparse parsed arguments object.

    Return:
        pandas.DataFrame: The result of combining the 3 tables.
    """

    print(f"Loading {args.expansion_hunter_tsv}")

    # Parse ExpansionHunter tsv
    df = parse_expansion_hunter_tsv_and_non_ref_motif_tsv(
        args.expansion_hunter_tsv, args.non_ref_motif_tsv, no_readviz_images=args.no_readviz)

    if args.expansion_hunter_tsv and args.non_ref_motif_tsv_using_offtargets:
        df_using_offtargets = parse_expansion_hunter_tsv_and_non_ref_motif_tsv(
            args.expansion_hunter_tsv_using_offtargets, args.non_ref_motif_tsv_using_offtargets, no_readviz_images=True)

        join_on_columns = ["SampleId", "LocusId", "VariantId", "ReferenceRegion"]
        df_using_offtargets = df_using_offtargets.set_index(join_on_columns)
        df_using_offtargets = df_using_offtargets[
            ["Genotype", "GenotypeConfidenceInterval", "Num Repeats: Allele 1", "Num Repeats: Allele 2"]
        ]
        df = df.set_index(join_on_columns).join(df_using_offtargets, how="left", rsuffix="__UsingOfftargetRegions")
        df = df.reset_index()
        num_missing_genotypes_using_offtargets = sum(pd.isna(df["Genotype__UsingOfftargetRegions"]))
        if num_missing_genotypes_using_offtargets > 0:
            raise ValueError(f"Could not find {num_missing_genotypes_using_offtargets} out of {len(df)} "
                             f"({100*num_missing_genotypes_using_offtargets/len(df):0.1f}%) "
                             f"genotypes in the off-targets table: {df[pd.isna(df['Genotype__UsingOfftargetRegions'])]}")

    # Parse gnomAD metadata tsv
    print(f"Loading {args.gnomad_metadata_tsv}")
    gnomad_df = pd.read_table(args.gnomad_metadata_tsv, low_memory=False)
    gnomad_df = gnomad_df[gnomad_df.release]
    gnomad_df.loc[:, "age"] = gnomad_df["project_meta.age"].fillna(gnomad_df["project_meta.age_alt"])
    gnomad_df["age"].fillna(AGE_NOT_AVAILABLE, inplace=True)

    gnomad_df.loc[:, "pcr_free"] = gnomad_df["project_meta.product"].apply(
        lambda s: pd.NA if not s or pd.isna(s) else (True if "pcr-free" in s.lower() else False), convert_dtype="boolean")
    gnomad_df["pcr_free"].fillna(gnomad_df["project_meta.v2_pcr_free"].astype("boolean"), inplace=True)
    gnomad_df["pcr_free"].fillna(PCR_INFO_NOT_AVAILABLE, inplace=True)
    gnomad_df.loc[:, "pcr_protocol"] = gnomad_df["pcr_free"].replace({True: "pcr_free", False: "pcr_plus"})

    gnomad_df = gnomad_df[[
        "s", "project_meta.sample_id", "population_inference.pop", "sex_imputation.sex_karyotype",
        "age", "pcr_protocol", "read_length",
        "project_meta.title", "project_meta.neuro_case",
    ]]
    gnomad_df.loc[:, "s"] = gnomad_df["project_meta.sample_id"].apply(process_sample_id)

    is_public_project = gnomad_df["project_meta.title"].isin(PUBLIC_PROJECT_NAMES)
    gnomad_df.loc[~is_public_project, "project_meta.title"] = None
    gnomad_df.loc[is_public_project, "public_sample_id"] = gnomad_df[is_public_project].s
    gnomad_df.rename(columns={"project_meta.title": "public_project_id", "project_meta.neuro_case": "is_neuro_case"},
                     inplace=True)

    unknown_sample_ids = set(df.SampleId) - set(gnomad_df.s)
    if len(unknown_sample_ids) > 0:
        print(f"WARNING: Dropping {len(unknown_sample_ids)} sample ids in {args.expansion_hunter_tsv} that "
              f"were not found in the gnomAD metadata table, or were found but are not 'release': ", unknown_sample_ids)
        assert len(unknown_sample_ids) < 5, unknown_sample_ids

        df = df[~df.SampleId.isin(unknown_sample_ids)]

    # Merge the data frames
    print(f"Combining STR data tables with gnomAD metadata")
    df = pd.merge(left=df, right=gnomad_df, how="inner", left_on="SampleId", right_on="s").drop(columns="s")

    print(f"Found {len(set(df.SampleId))} gnomAD 'release' samples")
    df = df[df.read_length >= MIN_READ_LENGTH]
    print(f"Kept {len(set(df.SampleId))} gnomAD samples after filtering to ReadLength >= {MIN_READ_LENGTH}")


    # Add imputed PCR labels
    inferred_pcr_status_df = pd.read_table(args.inferred_pcr_status_table)
    inferred_pcr_status_df = inferred_pcr_status_df[["external_sample_id", "inferred_pcr_status"]]
    inferred_pcr_status_df["external_sample_id"] = inferred_pcr_status_df["external_sample_id"].apply(process_sample_id)
    inferred_pcr_status_df["inferred_pcr_status"] = inferred_pcr_status_df["inferred_pcr_status"].replace({
        "PCRMINUS": "pcr_free",
        "PCRPLUS": "pcr_plus",
    })

    df = df.set_index("SampleId").join(
        inferred_pcr_status_df.set_index("external_sample_id"),
        how="left").reset_index().rename(columns={"index": "SampleId"})

    # mark HGDP samples as pcr_free, since they are known to have been sequenced with a PCR-free protocol
    df.loc[df.SampleId.str.startswith("HGDP"), "pcr_protocol"] = "pcr_free"

    # fill in missing PCR status with inferred PCR status
    df.loc[df["pcr_protocol"] == PCR_INFO_NOT_AVAILABLE, "pcr_protocol"] = float('nan')
    df["pcr_protocol"] = df["pcr_protocol"].fillna(df["inferred_pcr_status"])
    df.drop(columns="inferred_pcr_status", inplace=True)

    # discard any remaining samples with unknown PCR status
    n_samples_before_filter = len(set(df.SampleId))
    df = df[~df.pcr_protocol.isna()]
    n_samples_after_filter = len(set(df.SampleId))
    if n_samples_after_filter < n_samples_before_filter:
        print(f"Dropping {n_samples_before_filter - n_samples_after_filter} samples with unknown PCR status")
    assert n_samples_before_filter - n_samples_after_filter < 5

    locus_id_with_max_samples = None
    max_num_samples = 0
    for locus_id in sorted(set(df.LocusId)):
        num_samples = len(set(df[df.LocusId == locus_id].SampleId))
        if num_samples > max_num_samples:
            locus_id_with_max_samples = locus_id
            max_num_samples = num_samples

        print(f"Found {num_samples} {locus_id} 'release' samples")

    if locus_id_with_max_samples is not None:
        missing_age_counter = sum(df[df.LocusId == locus_id_with_max_samples]["age"] == AGE_NOT_AVAILABLE)
        message = (f"missing age for {missing_age_counter} out of {max_num_samples} "
                   f"({100*missing_age_counter/max_num_samples:0.1f}%) individuals")
        if missing_age_counter/max_num_samples > MISSING_AGE_THRESHOLD:
            raise ValueError(f"ERROR: {message}")
        print(f"WARNING: {message}")

        missing_pcr_protocol_counter = sum(df[df.LocusId == locus_id_with_max_samples]["pcr_protocol"] == "pcr_plus")
        message = (f"missing pcr_protocol for {missing_pcr_protocol_counter} out of {max_num_samples} " 
                  f"({100*missing_pcr_protocol_counter/max_num_samples:0.1f}%) individuals")
        if missing_pcr_protocol_counter/max_num_samples > MISSING_PCR_PROTOCOL_THRESHOLD:
            raise ValueError(f"ERROR: {message}")
        print(f"WARNING: {message}")

    return df


def init_gnomad_json(df):
    """Compute an initial .json structure with a key for each STR locus. Initialize sub-dictionaries that will hold
    the allele count histogram and scatter plot counts for each locus.

    Args:
        df (pandas.DataFrame): Combined DataFrame generated by load_data_df(..)

    Return:
        dict: An initial version of the main .json structure being generated by this script.
    """

    # Compute the STR loci
    df = df[["LocusId", "VariantId", "ReferenceRegion", "RepeatUnit"]]
    df = df.drop_duplicates()

    # Init sub-dictionaries for each locus
    gnomad_json = {}
    for _, row in tqdm.tqdm(df.iterrows(), unit=" rows", total=len(df)):
        locus_id = row["LocusId"]
        variant_id = row["VariantId"]
        adjacent_repeat_label = ADJACENT_REPEAT_LABELS[variant_id] if variant_id in ADJACENT_REPEAT_LABELS else None

        if locus_id not in gnomad_json:
            gnomad_json[locus_id] = {
                "LocusId": locus_id,
            }

        repeat_specific_fields = {
            "ReferenceRegion": row["ReferenceRegion"],
            "ReferenceRepeatUnit": row.get("RepeatUnit"),
            "AlleleCountHistogram":  {},
            "AlleleCountScatterPlot": {},
            "AgeDistribution": {},
        }

        if adjacent_repeat_label is not None:
            if "AdjacentRepeats" not in gnomad_json[locus_id]:
                gnomad_json[locus_id]["AdjacentRepeats"] = {}

            if variant_id not in gnomad_json[locus_id]["AdjacentRepeats"]:
                gnomad_json[locus_id]["AdjacentRepeats"][adjacent_repeat_label] = repeat_specific_fields
        else:
            gnomad_json[locus_id].update(repeat_specific_fields)

    return gnomad_json


def add_gene_ids(gnomad_json):
    """Add the GeneId field to gnomad_json.

    Args:
        gnomad_json (dict): The main .json structure being generated by this script.
    """
    for locus_id in gnomad_json:
        if not gnomad_json[locus_id].get("GeneName"):
            continue

        gene_name = gnomad_json[locus_id]["GeneName"]
        if gene_name in GENE_NAME_TO_GENE_ID:
            gnomad_json[locus_id]["GeneId"] = GENE_NAME_TO_GENE_ID[gene_name]
            continue

        # Get gene id via the Ensembl API.
        response = None
        while response is None or not response.ok or not response.json():
            print(f"Getting gene id for {gene_name}")
            request_url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_name}"
            request_url += "?content-type=application/json;expand=1"
            response = requests.get(request_url)

        response_json = response.json()
        if not response_json.get('id'):
            print("Unable to get ensembl details for", gene_name)
            continue

        gene_id = response_json['id']
        gnomad_json[locus_id]["GeneId"] = gene_id


def add_known_pathogenic_STR_annotations(args, gnomad_json):
    """Load the args.known_pathogenic_strs_tsv table and add metadata from it to gnomad_json.

    Args:
        args (argparse.Namespace): The argparse parsed arguments object.
        gnomad_json (dict): The main .json structure being generated by this script.
    """

    known_pathogenic_strs_info = parse_known_pathogenic_strs_tsv(args.known_pathogenic_strs_tsv)
    if len(known_pathogenic_strs_info) != EXPECTED_N_KNOWN_PATHOGENIC_REPEATS:
        raise ValueError(f"{args.known_pathogenic_strs_tsv} contains {len(known_pathogenic_strs_info)} pathogenic loci."
                         f" Expected {EXPECTED_N_KNOWN_PATHOGENIC_REPEATS} loci.")
    locus_ids_without_annotations = set(gnomad_json.keys()) - set(known_pathogenic_strs_info)
    if locus_ids_without_annotations:
        print(f"WARNING: These LocusIds were not found in known pathogenic STRs spreadsheet: {locus_ids_without_annotations}")

    # Compute STRipy urls
    for locus_id in gnomad_json:
        stripy_name = locus_id
        stripy_url = f"https://stripy.org/database/{stripy_name}"
        r = requests.get(stripy_url)
        if r.ok and "invalid locus" not in r.content.decode("UTF-8").lower():
            known_pathogenic_strs_info[locus_id]["STRipyName"] = stripy_name
        else:
            print(f"WARNING: STRipy page not found for {locus_id}")

    # Add the metadata to gnomad_json
    for locus_id in gnomad_json:
        if locus_id in known_pathogenic_strs_info:
            gnomad_json[locus_id].update(known_pathogenic_strs_info[locus_id])


def compute_most_common_motif_lookup_dict(df):
    """Create a lookup dictionary that maps (LocusId, canonical motif) pairs to the most common non-canonical motif
    among observed motifs that share this same canonical motif. This allows converting motif rearrangements such as
    AAAAG, AAAGA, AAGAA, etc. at the RFC1 locus into "AAAAG" which is the rearrangement that is seen most frequently
    in the general population. Similarly, for the HTT locus, "AGC", "CAG", and "GCA" would get converted to "CAG" since
    that's the only rearrangement that's seen in practice.

    Args:
         df (pandas.DataFrame): Combined DataFrame generated by load_data_df(..)

    Return:
         dict: A dictionary of the form {("RFC1", "AAAAG"): "AAAAG", ...}
    """

    # First, create a dictionary that maps each (LocusId, Motif) pair to the number of times it occurs in df.
    # Example entries:  ('RFC1', 'GAAAG'): 805,  ('RFC1', 'AAAGG'): 774, etc.
    motif_counts = pd.concat([
        df[["LocusId", "Motif: Allele 1"]].rename(columns={"Motif: Allele 1": "Motif"}),
        df[["LocusId", "Motif: Allele 2"]].rename(columns={"Motif: Allele 2": "Motif"}),
    ]).value_counts().to_dict()

    # Create a new dictionary that maps (LocusId, canonical motif) pairs to the most common non-canonical motif
    # observed among motifs that share the same canonical motif. Using the example from the previous comment, it would
    # map ('RFC1', 'AAAGG') to 'GAAAG' rather than 'AAAGG' because 'GAAAG' is observed 805 times while 'AAAGG' is only
    # observed 774 times in df.
    most_common_motif_lookup = {}
    for (locus_id, motif), counter in motif_counts.items():
        key = (locus_id, compute_canonical_motif(motif))
        if key not in most_common_motif_lookup:
            most_common_motif_lookup[key] = (motif, counter)
            continue

        previous_motif, previous_counter = most_common_motif_lookup[key]
        if previous_counter < counter:
            most_common_motif_lookup[key] = (motif, counter)

    # Drop the counter from the value
    most_common_motif_lookup = {key: motif for key, (motif, _) in most_common_motif_lookup.items()}

    return most_common_motif_lookup


def add_motif_classification_field(gnomad_json, most_common_motif_lookup):
    """For repeats where the pathogenic motif differs from the reference motif, add info on which motifs are known to be
    disease-associated and which are benign.

    Args:
        gnomad_json (dict): The main .json structure being generated by this script.
        most_common_motif_lookup (dict): The dictionary generated by compute_most_common_motif_lookup_dict(..)
    """

    non_ref_pathogenic_motif_info = json.loads(pkgutil.get_data("str_analysis", "data/non_ref_motif.locus_info.json"))

    for locus_id in gnomad_json:
        if not gnomad_json[locus_id].get("GeneName"):
            continue

        gene_name = gnomad_json[locus_id]["GeneName"]
        if gene_name not in non_ref_pathogenic_motif_info:
            continue

        gnomad_json[locus_id]["RepeatUnitClassification"] = {}
        for classification, motifs in non_ref_pathogenic_motif_info[gene_name]["Motifs"].items():
            for motif in motifs:
                canonical_motif = compute_canonical_motif(motif)
                motif_key = most_common_motif_lookup.get((locus_id, canonical_motif))
                if motif_key is None:
                    # If this known-benign or known-pathogenic motif wasn't detected in any gnomAD samples, just
                    # include it as-is, the way it's recorded in data/locus_info.json
                    motif_key = motif
                gnomad_json[locus_id]["RepeatUnitClassification"][motif_key] = classification


def add_histograms_and_compute_readviz_paths(
        df, gnomad_json, most_common_motif_lookup, existing_readviz_filename_list=None, no_readviz_images=False,
        include_all_age_and_pcr_info=False,
):
    """Populate the AlleleCountHistogram, AlleleCountScatterPlot and AgeDistribution. Also, compute encrypted readviz
    paths and add these & other metadata to readviz_json.

    Args:
        df (pandas.DataFrame): Combined DataFrame generated by load_data_df(..)
        gnomad_json (dict): The main .json structure being generated by this script.
        most_common_motif_lookup (dict): The dictionary generated by compute_most_common_motif_lookup_dict(..)
        existing_readviz_filename_list (str): Path of a text file that lists all readviz .svg filenames that exist
            (one per line). These are the encrypted public filenames that don't contain sample ids -
            for example: ffa0880117e0791d51b0ef85b56f3a54216.svg. This argument is mutually-exclusive with
            no_readviz_images=True
        no_readviz_images (bool): If True, this method will assume the REViewer images were not generated for this
            dataset.
        include_all_age_and_pcr_info (bool): If True, add age and pcr info to the user-friendly genotypes table
            for all samples for which age and pcr info is shown in the gnomAD browser.
    Return:
        (list, dict): 2-tuple containing (readviz_paths_to_rename, readviz_json) where
            readviz_paths_to_rename is a list of 2-tuples that matches the original readviz svg filename with the
                corresponding encrypted filename that can be made public.
            readviz_json is the .json data structure that will be loaded into the gnomAD browser to generate the
                readviz section of the STR pages. It contains the encrypted readviz filename and associated metadata
                for each sample.
            user_friendly_genotypes_json contains a flat list of records that will be shared on
                the gnomAD Downloads page. It's main purpose is to provide a user-friendly view of the data
                underlying the browser STR pages. IT contains the subset of columns that may be useful for downstream
                analyses.
    """

    readviz_paths_to_rename = set()
    readviz_json = collections.defaultdict(list)
    user_friendly_genotypes_json = []
    user_friendly_genotypes_json_binned = collections.defaultdict(list)
    age_counter = collections.defaultdict(int)
    public_sample_id_counter = collections.defaultdict(set)

    total_readviz_counter = collections.defaultdict(int)
    missing_readviz_counter = collections.defaultdict(int)

    existing_readviz_filenames_set = None
    if existing_readviz_filename_list:
        existing_readviz_filenames_df = pd.read_table(existing_readviz_filename_list, names=["filenames"], low_memory=False)
        existing_readviz_filenames_list = existing_readviz_filenames_df["filenames"]
        existing_readviz_filenames_set = set(existing_readviz_filenames_list)

        if len(existing_readviz_filenames_list) > len(existing_readviz_filenames_set):
            raise ValueError(f"{existing_readviz_filename_list} contains duplicate entries")

    df = df.sort_values(["Num Repeats: Allele 2", "Num Repeats: Allele 1", "Motif: Allele 2", "Motif: Allele 1"], ascending=False)
    for _, row in tqdm.tqdm(df.iterrows(), unit=" rows", total=len(df)):
        locus_id = row["LocusId"]
        variant_id = row["VariantId"]

        is_adjacent_repeat = variant_id in ADJACENT_REPEAT_LABELS
        adjacent_repeat_label = ADJACENT_REPEAT_LABELS[variant_id] if is_adjacent_repeat else None
        canonical_motif1 = compute_canonical_motif(row["Motif: Allele 1"])
        canonical_motif2 = compute_canonical_motif(row["Motif: Allele 2"])
        motif1 = most_common_motif_lookup[locus_id, canonical_motif1]
        motif2 = most_common_motif_lookup[locus_id, canonical_motif2]

        # Get gnomAD fields
        sex_karyotype = row["sex_imputation.sex_karyotype"]
        population = row["population_inference.pop"]
        pcr_protocol = row["pcr_protocol"]
        read_length = int(row["read_length"])
        public_project_id = row["public_project_id"] if not pd.isna(row["public_project_id"]) else None
        public_sample_id = row["public_sample_id"] if not pd.isna(row["public_sample_id"]) else None
        is_neuro_case = row["is_neuro_case"]

        # Compute age_range
        if row["age"] == AGE_NOT_AVAILABLE:
            age_range = AGE_NOT_AVAILABLE
        else:
            age = int(row["age"])
            age_lower_bound = AGE_RANGE_SIZE * math.floor(age/AGE_RANGE_SIZE)
            age_upper_bound = AGE_RANGE_SIZE * math.ceil((age + 0.1)/AGE_RANGE_SIZE)
            assert age_lower_bound != age_upper_bound
            if age_upper_bound <= LOWER_AGE_CUTOFF:
                age_range = f"<{LOWER_AGE_CUTOFF}"
            elif age_lower_bound >= UPPER_AGE_CUTOFF:
                age_range = f">{UPPER_AGE_CUTOFF}"
            else:
                age_range = f"{age_lower_bound}-{age_upper_bound}"

        age_range_to_show_in_readviz_section = AGE_NOT_AVAILABLE
        if (population in POPULATIONS_WITH_AGE_DISPLAYED_IN_READVIZ_SECTION
                and age_counter[locus_id, sex_karyotype] < MAX_AGES_PER_BUCKET_TO_DISPLAY_IN_THE_READVIZ_SECTION):
            age_counter[locus_id, sex_karyotype] += 1
            age_range_to_show_in_readviz_section = age_range

        age_range_for_user_friendly_genotypes_file = age_range if include_all_age_and_pcr_info else age_range_to_show_in_readviz_section

        # Get num_repeats1, num_repeats2
        try:
            num_repeats1 = int(float(row["Num Repeats: Allele 1"]))
            num_repeats2 = float(row["Num Repeats: Allele 2"])
        except ValueError as e:
            print("Num Repeats parse error", e, row["Genotype"], row["GenotypeConfidenceInterval"], ". Skipping..")
            continue

        if sex_karyotype == "XY" and "X" in row["ReferenceRegion"]:
            is_hemizygous = True
            if math.isnan(num_repeats2) or num_repeats2 == num_repeats1:
                num_repeats2 = num_repeats1
            else:
                print(f"ERROR: Locus is male and on chrX, but has different values for allele1, allele2: {row.to_dict()}")
                continue
        else:
            is_hemizygous = False

        num_repeats2 = int(num_repeats2)

        # Update histogram and scatter plot counts
        histogram_key1 = f"{population}/{sex_karyotype}/{motif1}"
        histogram_key2 = f"{population}/{sex_karyotype}/{motif2}"
        scatter_plot_key = f"{population}/{sex_karyotype}/{motif1}/{motif2}"
        age_distribution_key = f"{age_range}"

        if is_adjacent_repeat:
            data_dict = gnomad_json[locus_id]["AdjacentRepeats"][adjacent_repeat_label]
        else:
            data_dict = gnomad_json[locus_id]

        for histogram_key in histogram_key1, histogram_key2:
            if histogram_key not in data_dict["AlleleCountHistogram"]:
                data_dict["AlleleCountHistogram"][histogram_key] = collections.defaultdict(int)

        data_dict["AlleleCountHistogram"][histogram_key1][f"{num_repeats1}"] += 1
        if not is_hemizygous:
            data_dict["AlleleCountHistogram"][histogram_key2][f"{num_repeats2}"] += 1

        if scatter_plot_key not in data_dict["AlleleCountScatterPlot"]:
            data_dict["AlleleCountScatterPlot"][scatter_plot_key] = collections.defaultdict(int)
        data_dict["AlleleCountScatterPlot"][scatter_plot_key][f"{num_repeats1}/{num_repeats2}"] += 1

        if age_range != AGE_NOT_AVAILABLE:
            if age_distribution_key not in data_dict["AgeDistribution"]:
                data_dict["AgeDistribution"][age_distribution_key] = collections.defaultdict(int)
            for num_repeats in num_repeats1, num_repeats2:
                data_dict["AgeDistribution"][age_distribution_key][f"{num_repeats}"] += 1

        # Update readviz metadata
        encrypted_svg_prefix = hashlib.sha512(f"{locus_id}_{row['SampleId']}_{salt}".encode("UTF-8")).hexdigest()
        # The sha digest is 128 letters long - which is too long for a filename. Use only the first 35 letters.
        encrypted_svg_filename = f"{encrypted_svg_prefix[:35]}.svg"

        # check whether this image was successfully generated by REViewer, and set to None if not
        if no_readviz_images or (existing_readviz_filenames_set is not None
                                 and encrypted_svg_filename not in existing_readviz_filenames_set):
            encrypted_svg_filename = None
        elif not is_adjacent_repeat:
            original_svg_filename = row["ReadvizFilename"]
            readviz_paths_to_rename.add((original_svg_filename, f"{locus_id}/{encrypted_svg_filename}"))

        if not is_adjacent_repeat:
            total_readviz_counter[locus_id] += 1
            if encrypted_svg_filename is None:
                missing_readviz_counter[locus_id] += 1

        combined_motif = motif1 if motif1 == motif2 else f"{motif1}/{motif2}"
        reference_region_chrom, reference_region_start_0based, reference_region_end = parse_interval(row["ReferenceRegion"])
        user_friendly_genotypes_record = {
            "Id": variant_id,  # use variant_id instead of locus_id here because variant_id differs for main vs. adjacent repeats while locus_id is the same for both
            "LocusId": locus_id,
            "Chrom": reference_region_chrom,
            "Start_0based": int(reference_region_start_0based),
            "End": int(reference_region_end),
            "Motif": combined_motif,
            "IsAdjacentRepeat": is_adjacent_repeat,
            "Population": population,
            "Sex": sex_karyotype,
            "Age": age_range_for_user_friendly_genotypes_file,
            "PcrProtocol": pcr_protocol,
            "ReadLength": read_length,
            "ReadvizFilename": encrypted_svg_filename,
            "PublicProjectId": public_project_id,
            "PublicSampleId": public_sample_id,
            "IsNeuroCase": is_neuro_case,
        }

        if public_sample_id:
            public_sample_id_counter[public_project_id].add(public_sample_id)

        # transfer additional columns from the table row. The values are: (to_column, to_column_type, warn_if_missing)
        extra_columns = {
            "ReferenceRegion": ("ReferenceRegion", str, True),
            "Genotype": ("Genotype", str, True),
            "GenotypeConfidenceInterval": ("GenotypeConfidenceInterval", str, True),
            "Num Repeats: Allele 1": ("Allele1", int, True),
            "Num Repeats: Allele 2": ("Allele2", int, False),  # expected to be missing in hemizygous genotypes
            "Genotype__UsingOfftargetRegions": ("GenotypeUsingOfftargetRegions", str, True),
            "GenotypeConfidenceInterval__UsingOfftargetRegions": ("GenotypeConfidenceIntervalUsingOfftargetRegions", str, True),
            "Num Repeats: Allele 1__UsingOfftargetRegions": ("Allele1UsingOfftargetRegions", int, True),
            "Num Repeats: Allele 2__UsingOfftargetRegions": ("Allele2UsingOfftargetRegions", int, False), # expected to be missing in hemizygous genotypes
        }
        for from_column, (to_column, to_column_type, warn_if_missing) in extra_columns.items():
            if from_column not in set(df.columns):
                print(f"WARNING: data table is missing a '{from_column}' column")
            elif row[from_column] is None or pd.isna(row[from_column]):
                if warn_if_missing:
                    print(f"WARNING: missing '{from_column}' value in row {row.to_dict()}")
                user_friendly_genotypes_record[to_column] = ""
            else:
                user_friendly_genotypes_record[to_column] = to_column_type(row[from_column])

        user_friendly_genotypes_json.append(user_friendly_genotypes_record)

        user_friendly_genotypes_json_binned[
            (variant_id, population, sex_karyotype, age_range_for_user_friendly_genotypes_file, pcr_protocol)
        ].append(user_friendly_genotypes_record)

        if not is_adjacent_repeat:
            readviz_json[locus_id].append({
                "Allele1Motif": motif1,
                "Allele2Motif": motif2,
                "Allele1HistogramKey": histogram_key1,
                "Allele2HistogramKey": histogram_key2 if not is_hemizygous else None,
                "ScatterPlotKey": scatter_plot_key,
                "ScatterPlotX": num_repeats2,
                "ScatterPlotY": num_repeats1,
                "Sex": sex_karyotype,
                "Age": age_range_to_show_in_readviz_section,
                "Population": population,
                "PcrProtocol": pcr_protocol,
                "Genotype": row["Genotype"],
                "GenotypeConfidenceInterval": row["GenotypeConfidenceInterval"],
                "ReadLength": read_length,
                "ReadvizFilename": encrypted_svg_filename,
                "PublicProjectId": public_project_id,
                "PublicSampleId": public_sample_id,
                "IsNeuroCase": is_neuro_case,
            })

    for project_id, sample_id_set in sorted(public_sample_id_counter.items(), key=lambda x: -len(x[1])):
        print(f"{len(sample_id_set):,d} public sample ids in gnomAD project '{project_id}'")

    # check whether an unexpected number of readviz images are missing
    if not no_readviz_images:
        for locus_id, missing_count in missing_readviz_counter.items():
            total = total_readviz_counter[locus_id]
            message = (f"{locus_id:20s}:  {missing_count} out of {total} ({100*missing_count/total:0.2f}%) readviz images removed "
                       f"because they are missing from {existing_readviz_filename_list}")

            if missing_count/total > MISSING_READVIZ_ERROR_THRESHOLD:
                raise ValueError(message)

            print(message)

    # Check for samples that fall into a unique metadata bin - where they are the only samples with a particular
    # combination of values for population, sex, age-bin, and pcr-protocol. Having only 1 individual in a bin like this
    # could allow someone to get all STR genotypes for that single individual by looking up its metadata signature
    # at each STR locus.
    # For added caution, and in-keeping with the principle of only sharing per-variant information, this loop discards
    # the Age and PCR-protocol values from the user_friendly_genotypes_json data structure. This should move these
    # samples into more common metadata bins that have other samples and prevent phasing of STR variants across loci.
    if not include_all_age_and_pcr_info:
        modified_records = []
        for bin_key, user_friendly_genotypes_records in user_friendly_genotypes_json_binned.items():
            if len(user_friendly_genotypes_records) == 1:
                user_friendly_genotypes_record = user_friendly_genotypes_records[0]
                (variant_id, population, sex_karyotype, age, pcr_protocol) = bin_key

                # move this sample to a more common metadata bin by dropping age
                new_bin_key = (variant_id, population, sex_karyotype, AGE_NOT_AVAILABLE, pcr_protocol)
                user_friendly_genotypes_record["Age"] = AGE_NOT_AVAILABLE
                if len(user_friendly_genotypes_json_binned.get(new_bin_key, [])) < 2:
                    # dropping age didn't move this sample to a more common bin, so also drop PCR protocol info
                    new_bin_key = (variant_id, population, sex_karyotype, AGE_NOT_AVAILABLE, PCR_INFO_NOT_AVAILABLE)
                    user_friendly_genotypes_record["PcrProtocol"] = PCR_INFO_NOT_AVAILABLE

                modified_records.append((bin_key, new_bin_key, user_friendly_genotypes_record))

        visited_individual_metadata_bins = set()
        for bin_key, new_bin_key, modified_user_friendly_genotypes_record in modified_records:
            (_, population, sex_karyotype, age, pcr_protocol) = bin_key
            individual_bin = (population, sex_karyotype, age, pcr_protocol)
            if individual_bin not in visited_individual_metadata_bins:
                print(f"Moving sample from unique metadata bin {bin_key} to {new_bin_key}")
                visited_individual_metadata_bins.add(individual_bin)

            del user_friendly_genotypes_json_binned[bin_key]
            user_friendly_genotypes_json_binned[new_bin_key].append(modified_user_friendly_genotypes_record)

        unable_to_avoid_unique_metadata_bins = False
        for bin_key, user_friendly_genotypes_records in user_friendly_genotypes_json_binned.items():
            if len(user_friendly_genotypes_records) == 1:
                print(f"Unable to avoid a metadata bin with only 1 sample. Bin key: {bin_key} ")
                print(user_friendly_genotypes_records[0])
                unable_to_avoid_unique_metadata_bins = True
        if unable_to_avoid_unique_metadata_bins:
            raise ValueError("Failed to avoid metadata bins wiht only 1 sample")

    return list(readviz_paths_to_rename), readviz_json, user_friendly_genotypes_json


def sort_keys(gnomad_json):
    """Sort keys in the output json. Python built-in dictionaries preserve key order since python3.6, so this works.
    For example, for the "AFF2" locus this sorts keys so that

    {
      "ReferenceRegion": "chrX:148500631-148500691",
      "ReferenceRepeatUnit": "GCC",
      "LocusId": "AFF2",
      "AgeDistribution": {
         "35-40": {
            "4": 3,
            "0": 10,
            "2": 2,
            ...
         },
         "20-25": {
            "5": 5,
            "6": 12,
            "0": 10,
            "2": 2,
            "4": 3,
            ...
        },
      }
    }

    is converted to

    {
      "LocusId": "AFF2",
      "ReferenceRegion": "chrX:148500631-148500691",
      "ReferenceRepeatUnit": "GCC",
      "AgeDistribution": {
         "20-25": {
            "0": 10,
            "2": 2,
            "4": 3,
            "5": 5,
            "6": 12,
            ...
        },
         "35-40": {
            "0": 10,
            "2": 2,
            "4": 3,
            ...
         },
      }
    }

    Args:
        gnomad_json (dict): The main .json structure being generated by this script.
    """

    def sort_by_key(key_type=str):
        def key_lookup(key_value):
            return key_type(key_value[0])
        return key_lookup

    def top_level_sort_order(key_value):
        # Sort top-level keys so that the histograms are last
        return key_value[0] in ("AlleleCountHistogram", "AlleleCountScatterPlot", "AgeDistribution"), key_value[0]

    for locus_id, locus_data in gnomad_json.items():
        for histogram_name, histogram_key_type in (
            ("AlleleCountHistogram", int),
            ("AlleleCountScatterPlot", str),
            ("AgeDistribution", int),
        ):
            # `histogram_key` here refers to the top level keys in the histogram
            # For example the "20-25" age range is a `histogram_key` within "AgeDistribution"
            # Each of the age ranges within "AgeDistribution" has a nested dict. 
            # e.g., "0" is a key within "20-25": 
            # gnomad_json["AFF2"]["AgeDistribution"]["20-25"]["0"] = 10
            # This first `for` loop below sorts these nested dicts (e.g., sorts the dicts within the age range "20-25")
            for histogram_key in locus_data[histogram_name]:
                locus_data[histogram_name][histogram_key] = {
                    key: value for key, value in sorted(
                        locus_data[histogram_name][histogram_key].items(), key=sort_by_key(key_type=histogram_key_type))
                }
            # This sorts the `histogram_key` values
            # e.g, this sorts "20-25", "25-30", "30-35", etc. within "AgeDistribution" 
            locus_data[histogram_name] = {
                key: value for key, value in sorted(
                    locus_data[histogram_name].items(), key=sort_by_key())
            }

        # This sorts the top level keys, which contain the histogram names above (e.g., AgeDistribution)
        # and other keys, like "Diseases", "GeneID", "GeneName", etc.
        gnomad_json[locus_id] = {
            key: value for key, value in sorted(gnomad_json[locus_id].items(), key=top_level_sort_order)
        }


def validate_json(df, gnomad_json, readviz_json, user_friendly_genotypes_json, no_readviz_images=False):
    """Perform basic checks to validate the gnomad_json and readviz_json data structure.

    Args:
        df (pandas.DataFrame): Combined DataFrame generated by load_data_df(..).
        gnomad_json (dict): The main .json structure being generated by this script.
        readviz_json (dict): The .json data structure that will be loaded into the gnomAD browser to generate the
                readviz section of the STR pages. It contains the encrypted readviz filename and associated
                metadata for each sample.
        no_readviz_images (bool): If True, this method will skip validation of the "ReadvizFilename" field. This is used
            for datasets where REViewer images were not generated.
    """

    total_samples = len(set(df["SampleId"]))
    if len(gnomad_json) != EXPECTED_N_KNOWN_PATHOGENIC_REPEATS:
        print(f"WARNING: gnomad_json contains {len(gnomad_json)} pathogenic loci. "
            f"Expected {EXPECTED_N_KNOWN_PATHOGENIC_REPEATS} loci.")

    gnomad_json_str_loci = set(gnomad_json)
    readviz_json_str_loci = set(readviz_json)
    if gnomad_json_str_loci != readviz_json_str_loci:
        raise ValueError(f"gnomad_json locus ids are different from readviz_json locus ids:\n"
                         f"{gnomad_json_str_loci} \n{readviz_json_str_loci}")

    fraction_male_samples = sum(df["sex_imputation.sex_karyotype"] == "XY")/len(df)
    for locus_id, data in gnomad_json.items():
        # Check that expected keys are present and have non-null values.
        for key in "ReferenceRepeatUnit", "LocusId", "GeneName", "GeneId", "ReferenceRegion", \
                   "AlleleCountHistogram", "AlleleCountScatterPlot", "AgeDistribution", "Diseases":
            if data.get(key) is None:
                print(f"WARNING: {locus_id} {key} is None")

        # Check that expected keys are present in the data["Diseases"] dictionary and have non-null values.
        for key in "Symbol", "Name", "Inheritance", "PathogenicMin", "OMIM":
            diseases = data.get("Diseases", [])
            if not diseases:
                print(f"WARNING: {locus_id} has no Diseases data")
                break
            for i, disease_data in enumerate(diseases):
                if disease_data.get(key) is None:
                    print(f"WARNING: {locus_id} disease #{i} {key} is None")

        # Check that total counts in the histogram and scatter plot roughly match expectation, taking into account
        # hemizygous genotypes (which only contribute 1 count) and missing genotypes due to low coverage in some samples
        if "X" in data["ReferenceRegion"]:
            expected_counts_in_histogram = total_samples * (2 - fraction_male_samples)
        else:
            expected_counts_in_histogram = total_samples * 2

        for key, expected_counts in [
            ("AlleleCountHistogram", expected_counts_in_histogram),
            ("AlleleCountScatterPlot", total_samples)
        ]:
            total_counts_in_plot = sum([sum(d.values()) for d in data[key].values()])
            if not ((1 - MISSING_GENOTYPES_ERROR_THRESHOLD) * expected_counts < total_counts_in_plot <= expected_counts):
                print(f"WARNING: {locus_id} total counts in {key} = {total_counts_in_plot} while expected counts = {expected_counts}")

        total_readviz_samples = len(readviz_json[locus_id])
        if total_readviz_samples < (1 - MISSING_GENOTYPES_ERROR_THRESHOLD) * total_samples:
            print(f"WARNING: {locus_id}: only {total_readviz_samples} readviz records. Expected {total_samples}")
        if total_readviz_samples > total_samples:
            print(f"WARNING: {locus_id}: found {total_readviz_samples} readviz records which is more than the total "
                             f"number of samples ({total_samples})")

        if not no_readviz_images:
            total_readviz_samples_with_image = sum(1 for r in readviz_json[locus_id] if r["ReadvizFilename"] is not None)
            if total_readviz_samples_with_image < (1 - MISSING_READVIZ_ERROR_THRESHOLD) * total_readviz_samples:
                raise ValueError(f"{locus_id}: found {total_readviz_samples_with_image} readviz images. Expected at "
                                 f"least {total_readviz_samples}.")

    if len(user_friendly_genotypes_json) < sum(1 for locus_id in readviz_json for _ in readviz_json[locus_id]):
        raise ValueError("user_friendly_genotypes_json contains fewer records than readviz_json")


def export_to_tsv(df, output_tsv_filename):
    """Write the given pandas DataFrame to the output_tsv_filename"""
    if not output_tsv_filename.endswith(".tsv"):
        raise ValueError(f"{output_tsv_filename} must end with '.tsv'")
    print(f"Writing {len(df)} records to {output_tsv_filename}.gz")
    df.to_csv(output_tsv_filename, sep="\t", index=False, header=True)
    os.system(f"bgzip -f {output_tsv_filename}")


def export_readviz_rename_list(readviz_paths_to_rename, readviz_rename_list_output_path):
    """Utility function for writing out the readviz_paths_to_rename data structure.

    Args:
        readviz_paths_to_rename (list): list of 2-tuples that matches the original readviz svg filename with the
            corresponding encrypted filename that can be made public.
        readviz_rename_list_output_path (str): Local output path where to write the readviz_paths_to_rename table.
    """

    if not readviz_rename_list_output_path.endswith(".gz"):
        raise ValueError(f"{readviz_rename_list_output_path} needs to end in .gz")

    print(f"Writing {readviz_rename_list_output_path}")
    with gzip.open(readviz_rename_list_output_path, "wt") as f:
        for a, b in readviz_paths_to_rename:
            f.write(f"{a}\t{b}\n")


def main():
    """Generate 3 files: the main gnomAD STR json data file which will be loaded into the gnomAD browser to populate
    the STR pages, a readviz_rename_list table which maps REViewer image filenames to the corresponding encrypted
    filenames which can be made public without revealing sample ids, and a readviz_paths_json file which contains
    metadata on samples and readviz images.
    """

    args = parse_args()

    # Generate the 3 data structures
    df = load_data_df(args)
    gnomad_json = init_gnomad_json(df)
    add_known_pathogenic_STR_annotations(args, gnomad_json)
    add_gene_ids(gnomad_json)
    most_common_motif_lookup = compute_most_common_motif_lookup_dict(df)
    add_motif_classification_field(gnomad_json, most_common_motif_lookup)
    readviz_paths_to_rename, readviz_json, user_friendly_genotypes_json = add_histograms_and_compute_readviz_paths(
        df, gnomad_json, most_common_motif_lookup,
        existing_readviz_filename_list=args.existing_readviz_filename_list,
        no_readviz_images=args.no_readviz,
        include_all_age_and_pcr_info=args.include_all_age_and_pcr_info,
    )

    sort_keys(gnomad_json)

    # Perform validity checks
    validate_json(df, gnomad_json, readviz_json, user_friendly_genotypes_json, no_readviz_images=args.no_readviz)

    # Write out the data structures
    date_stamp = datetime.now().strftime("%Y_%m_%d")
    local_output_dir = args.local_output_dir
    if local_output_dir is None:
        local_output_dir = os.path.expanduser(os.path.dirname(args.expansion_hunter_tsv))
    print(f"Output dir: {local_output_dir}")

    output_filename_label = f"__{args.output_filename_suffix}" if args.output_filename_suffix else ""

    # write out the full df for debugging and internal analyses that need sample ids
    export_to_tsv(df, f"{local_output_dir}/gnomAD_STR_calls_with_gnomAD_metadata_and_sample_ids{output_filename_label}__{date_stamp}.tsv")

    # write out readviz_records to a .tsv mostly for debugging
    readviz_metadata_df = pd.DataFrame([
        {**readviz_record, **{"LocusId": locus_id}}
        for locus_id, readviz_records in readviz_json.items() for readviz_record in readviz_records
    ])
    export_to_tsv(readviz_metadata_df, f"{local_output_dir}/gnomAD_STR_readviz_metadata{output_filename_label}__{date_stamp}.tsv")

    # write out JSON files for the gnomAD browser
    export_json(readviz_json, f"{local_output_dir}/gnomAD_STR_readviz_metadata{output_filename_label}__{date_stamp}.json.gz", args.output_dir)
    export_json(gnomad_json, f"{local_output_dir}/gnomAD_STR_distributions{output_filename_label}__{date_stamp}.json.gz", args.output_dir)

    # write out a table of original readviz image filenames together with their encrypted filenames to use for renaming the files
    export_readviz_rename_list(readviz_paths_to_rename, f"{local_output_dir}/readviz_rename_list{output_filename_label}__{date_stamp}.tsv.gz")

    # write out user_friendly_genotypes_json to a .tsv so it can be shared publicly on the gnomAD downloads page
    def user_friendly_genotypes_sort_key(record):
        # sort the output table by the locus coordinates, then by the size of the long and short alleles in desc. order
        return (
            record["Id"],
            -1*int(record["Allele2"]) if record["Allele2"] else 0,
            -1*int(record["Allele1"]) if record["Allele1"] else 0,
        )

    user_friendly_genotypes_json.sort(key=user_friendly_genotypes_sort_key)
    user_friendly_genotypes_df = pd.DataFrame(user_friendly_genotypes_json)
    user_friendly_genotypes_df = user_friendly_genotypes_df.filter(items=[
        "Id", "LocusId",
        "ReferenceRegion",
        "Chrom", "Start_0based", "End",
        "Motif",
        "IsAdjacentRepeat",
        "Population", "Sex", "Age", "PcrProtocol",
        "Genotype",
        "Allele1",
        "Allele2",
        "GenotypeConfidenceInterval",
        "GenotypeUsingOfftargetRegions",
        "Allele1UsingOfftargetRegions",
        "Allele2UsingOfftargetRegions",
        "GenotypeConfidenceIntervalUsingOfftargetRegions",
        "ReadvizFilename",
        "ReadLength",
        "IsNeuroCase",
        "PublicProjectId",
        "PublicSampleId",
    ])
    output_path = f"{local_output_dir}/gnomAD_STR_genotypes{output_filename_label}__"
    if args.include_all_age_and_pcr_info:
        output_path += "including_all_age_and_pcr_info__"
    output_path += f"{date_stamp}.tsv"
    export_to_tsv(
        user_friendly_genotypes_df,
        output_path,
    )

    print("Done")


if __name__ == "__main__":
    main()
    