"""When filter_vcf_to_tandem_repeats.py is applied to more than one single-sample VCF, this script can be used to
combine the resulting *.variants.tsv.gz tables into a single table with one record per locus. It takes the
per-sample tables as well as the DipCall VCFs and confidence regions and then outputs a single combined table by:
1. performing an outer-join on the loci in the per-sample tables
2. setting missing genotypes in each sample to homozygous reference for loci that are within the DipCall confidence
   regions and don't overlap any variants in the input VCF
"""

import argparse
import collections
import gzip
from intervaltree import IntervalTree, Interval
import numpy as np
import os
import polars as pl
import subprocess
import tqdm

ALL_CHROMS = list(map(str, range(1, 23))) + ["X", "Y"]

"""
$1                    Chrom : chr1
$2              Start1Based : 674824
$3                End1Based : 674832
$4                    Locus : 1:674824-674832
$5                  LocusId : 1-674823-674832-AGG
$6               INS_or_DEL : INS
$7                    Motif : AGG
$8   MotifInterruptionIndex : 2.0
$9           CanonicalMotif : AGG
$10               MotifSize : 3
$11   NumRepeatsInReference : 3.0
$12                  VcfPos : 674829
$13                  VcfRef : C
$14                  VcfAlt : CAGG
$15             VcfGenotype : 1|0
$16           SummaryString : 3bp:AGG:INS:3=>4:HET:not-pure
$17      IsFoundInReference : True
$18            IsPureRepeat : False
$19          IsMultiallelic : False
$20              NumRepeats : 4
$21         RepeatSize (bp) : 12
$22          NumPureRepeats : 3
$23     PureRepeatSize (bp) : 9
$24     FractionPureRepeats : 0.75
"""

PER_LOCUS_COLUMNS = [
    "Chrom",
    "Start1Based",
    "End1Based",
    "Locus",
    "LocusId",
    #"INS_or_DEL",
    "Motif",
    "CanonicalMotif",
    "MotifSize",
    "NumRepeatsInReference",
    "IsFoundInReference",
]
SAMPLE_SPECIFIC_COLUMNS = [
    #"VcfPos",
    #"VcfRef",
    #"VcfAlt",
    #"VcfGenotype",
    #"SummaryString",
    #"IsMultiallelic",
    #"HET_or_HOM_or_HEMI_or_MULTI",
    "NumRepeatsShortAllele",
    "NumRepeatsLongAllele",
    #"RepeatSizeShortAllele (bp)",
    #"RepeatSizeLongAllele (bp)",
    "IsPureRepeat",
    #"MotifInterruptionIndex",

    #"NumRepeats",
    #"RepeatSize (bp)",
    #"NumPureRepeats",
    #"PureRepeatSize (bp)",
    #"FractionPureRepeats",
]

# Define schema for polars to avoid mixed types warning
SCHEMA = {
    "Chrom": pl.Utf8,
    "Start1Based": pl.Int32,
    "End1Based": pl.Int32, 
    "Locus": pl.Utf8,
    "LocusId": pl.Utf8,
    "Motif": pl.Utf8,
    "CanonicalMotif": pl.Utf8,
    "MotifSize": pl.Int32,
    "NumRepeatsInReference": pl.Float32,
    "IsFoundInReference": pl.Boolean,
    "NumRepeatsShortAllele": pl.Float32,
    "NumRepeatsLongAllele": pl.Float32,
    "IsPureRepeat": pl.Boolean,
}

SCHEMA_WITH_INT = SCHEMA.copy()
SCHEMA_WITH_INT["NumRepeatsInReference"] = pl.Int32
SCHEMA_WITH_INT["MotifSize"] = pl.Int32
SCHEMA_WITH_INT["NumRepeatsShortAllele"] = pl.Int32
SCHEMA_WITH_INT["NumRepeatsLongAllele"] = pl.Int32



def parse_input_tsv(sample_id, input_tsv, exclude_homopolymers = False, discard_impure_genotypes = False, chroms_to_process = None):

    # Read with polars
    df = pl.read_csv(input_tsv, separator="\t", columns=PER_LOCUS_COLUMNS + SAMPLE_SPECIFIC_COLUMNS, schema_overrides=SCHEMA)
    df = df.with_columns(
        pl.col("NumRepeatsInReference").cast(pl.Int32),
        pl.col("MotifSize").cast(pl.Int32),
        pl.col("NumRepeatsShortAllele").cast(pl.Int32),
        pl.col("NumRepeatsLongAllele").cast(pl.Int32),
    )

    # normalize the chromosome name
    df = df.with_columns(
        pl.col("Chrom").str.replace("chr", "").str.to_uppercase()
    )

    if chroms_to_process:
        df = df.filter(pl.col("Chrom").is_in(chroms_to_process))

    if exclude_homopolymers:
        df = df.filter(pl.col("MotifSize") > 1)
    
    missing_columns = set(PER_LOCUS_COLUMNS + SAMPLE_SPECIFIC_COLUMNS) - set(df.columns)
    if len(missing_columns) > 0:
        raise ValueError(f"{input_tsv} is missing these columns: {missing_columns}. Its columns are: {df.columns}")

    if discard_impure_genotypes:
        df = df.filter(pl.col("IsPureRepeat"))

    # Rename columns efficiently - build rename dict once
    rename_dict = {}
    for column in SAMPLE_SPECIFIC_COLUMNS:
        rename_dict[column] = f"{column}:{sample_id}"

    df = df.rename(rename_dict)

    # sort by chrom, start, end
    df = df.sort(PER_LOCUS_COLUMNS)

    return df


def merge_two_rows(keep_row, drop_row, locus_id_to_sample_ids_with_non_missing_genotypes, merge_per_sample_columns=False, verbose=False):
    # transfer all values from the "drop_row" to the "keep_row"
    keep_row = keep_row.copy()

    # update the ShortAllelesList and LongAllelesList. The allele_size should always be non-negative anyway 
    # since keep_row is the one with the larger reference interval, but use max(0, ..) to be safe.
    num_repeats_offset = keep_row["NumRepeatsInReference"] - drop_row["NumRepeatsInReference"]

    new_short_allele_list = []
    new_long_allele_list = []
    for keep_row_short_allele, keep_row_long_allele, drop_row_short_allele, drop_row_long_allele in zip(
        keep_row["ShortAllelesList"], keep_row["LongAllelesList"], drop_row["ShortAllelesList"], drop_row["LongAllelesList"]
    ):
        drop_row_placeholders = {None, drop_row["NumRepeatsInReference"]}
        keep_row_placeholders = {None, keep_row["NumRepeatsInReference"]}

        # at each position in the lists, either keep the allele from the keep_row or the allele from the drop_row
        if (drop_row_short_allele is None and drop_row_long_allele is None) or (drop_row_short_allele == drop_row["NumRepeatsInReference"] and drop_row_long_allele == drop_row["NumRepeatsInReference"]):
            new_short_allele_list.append(keep_row_short_allele)
            new_long_allele_list.append(keep_row_long_allele)
        elif (keep_row_short_allele is None and keep_row_long_allele is None) or (keep_row_short_allele == keep_row["NumRepeatsInReference"] and keep_row_long_allele == keep_row["NumRepeatsInReference"]):
            new_short_allele_list.append(drop_row_short_allele)
            new_long_allele_list.append(drop_row_long_allele)
        elif keep_row_short_allele not in keep_row_placeholders and drop_row_short_allele not in drop_row_placeholders:
            raise ValueError(f"keep_row_short_allele {keep_row_short_allele} and drop_row_short_allele {drop_row_short_allele} are both non-placeholders")

    #if keep_row["LocusId"] == "9-133392583-133392639-TGGA":
    #    print("keep_row ShortAllelesList, LongAllelesList", keep_row["ShortAllelesList"], keep_row["LongAllelesList"])
    #    print("drop_row ShortAllelesList, LongAllelesList", drop_row["ShortAllelesList"], drop_row["LongAllelesList"])

    #    print("new_row new_short_allele_list", new_short_allele_list)
    #    print("new_row new_long_allele_list", new_long_allele_list)

    assert len(keep_row["ShortAllelesList"]) == len(keep_row["LongAllelesList"])
    assert len(new_short_allele_list) == len(keep_row["ShortAllelesList"])
    assert len(new_long_allele_list) == len(keep_row["LongAllelesList"])

    keep_row["ShortAllelesList"] = new_short_allele_list
    keep_row["LongAllelesList"] = new_long_allele_list

    # update the "locus_id_to_sample_ids_with_non_missing_genotypes" dictionary
    locus_id_to_sample_ids_with_non_missing_genotypes[keep_row["LocusId"]] |= locus_id_to_sample_ids_with_non_missing_genotypes[drop_row["LocusId"]]
    del locus_id_to_sample_ids_with_non_missing_genotypes[drop_row["LocusId"]]

    # if they are present, transfer any per-sample NumRepeatsShortAllele or NumRepeatsLongAllele values from the "drop_row" to the "keep_row"
    if merge_per_sample_columns:
        all_sample_ids = {
            c.replace("NumRepeatsShortAllele:", "").replace("NumRepeatsLongAllele:", "") for c in keep_row.keys() if c.startswith("NumRepeatsShortAllele:") or c.startswith("NumRepeatsLongAllele:")
        }

        drop_row_placeholders = {None, drop_row["NumRepeatsInReference"]}
        keep_row_placeholders = {None, keep_row["NumRepeatsInReference"]}
        for sample_id in all_sample_ids:
            if not (
                drop_row[f"NumRepeatsShortAllele:{sample_id}"] is None and drop_row[f"NumRepeatsLongAllele:{sample_id}"] is None
                or
                drop_row[f"NumRepeatsShortAllele:{sample_id}"] == drop_row["NumRepeatsInReference"] and drop_row[f"NumRepeatsLongAllele:{sample_id}"] == drop_row["NumRepeatsInReference"]
            ):
                assert keep_row[f"NumRepeatsShortAllele:{sample_id}"] in keep_row_placeholders and keep_row[f"NumRepeatsLongAllele:{sample_id}"] in keep_row_placeholders, (f"{column_name}\nkeep_row: {keep_row}\ndrop_row:{drop_row}")
                
                if keep_row[f"NumRepeatsShortAllele:{sample_id}"] is None and keep_row[f"NumRepeatsLongAllele:{sample_id}"] is None:
                    keep_row["MissingHaplotypes"] -= 2

                keep_row[f"NumRepeatsShortAllele:{sample_id}"] = drop_row[f"NumRepeatsShortAllele:{sample_id}"] + num_repeats_offset if drop_row[f"NumRepeatsShortAllele:{sample_id}"] is not None else None
                keep_row[f"NumRepeatsLongAllele:{sample_id}"] = drop_row[f"NumRepeatsLongAllele:{sample_id}"] + num_repeats_offset if drop_row[f"NumRepeatsLongAllele:{sample_id}"] is not None else None

        assert keep_row["MissingHaplotypes"] >= 0, keep_row


    if verbose:  #  and row_i["LocusId"] == "19-29537703-29537757-CACCCACCTACTCTCCCG"
        pure_i = "Pure" if keep_row["IsPureRepeat"] else "Interrupted"
        pure_i_plus_1 = "Pure" if drop_row["IsPureRepeat"] else "Interrupted"
        pure_keep_row = "Pure" if keep_row["IsPureRepeat"] else "Interrupted"
        print(f"Merging rows: {keep_row['LocusId']} ({pure_i}) and {drop_row['LocusId']} ({pure_i_plus_1}), keeping {keep_row['LocusId']} ({pure_keep_row}). Adding {num_repeats_offset}x to dropped")
        #print(locus_id_to_sample_ids_with_non_missing_genotypes[row_i["LocusId"]] & locus_id_to_sample_ids_with_non_missing_genotypes[row_i_plus_1["LocusId"]])

    return keep_row



def merge_overlapping_loci(df, locus_id_to_sample_ids_with_non_missing_genotypes, merge_per_sample_columns=False, verbose=True):

    # sort by chrom, start, end
    df = df.sort(PER_LOCUS_COLUMNS)

    # iterate over the loci, checking for overlap between each locus and the next one, and then and create a new table with the merged loci
    result_rows = []
    i = 0
    row_i = df.row(i, named=True)
    merged_row_counter = 0
    while i < df.height - 1:
        row_i_plus_1 = df.row(i + 1, named=True)
        if (
            # The two loci overlap by 1 or more motif sizes
            row_i["Chrom"] == row_i_plus_1["Chrom"] and row_i_plus_1["Start1Based"] + row_i["MotifSize"] <= row_i["End1Based"] 
            and (
                # The two loci have the same canonical motif for STRs, or the same-sized motif for VNTRs
                row_i["CanonicalMotif"] == row_i_plus_1["CanonicalMotif"] or (
                row_i["MotifSize"] == row_i_plus_1["MotifSize"] and row_i["MotifSize"] > 6)
            )
            and (
                # There is no overlap between the samples that have a genotype at one of the loci and the samples that have a genotype at the other locus
                0 == len(locus_id_to_sample_ids_with_non_missing_genotypes[row_i["LocusId"]] & locus_id_to_sample_ids_with_non_missing_genotypes[row_i_plus_1["LocusId"]])
            )
        ):
            if row_i["NumRepeatsInReference"] >= row_i_plus_1["NumRepeatsInReference"]:
                keep_row_i = i
                drop_row_i = i + 1
            else:
                keep_row_i = i + 1
                drop_row_i = i

            merged_row = merge_two_rows(
                df.row(keep_row_i, named=True), 
                df.row(drop_row_i, named=True), 
                locus_id_to_sample_ids_with_non_missing_genotypes, 
                merge_per_sample_columns=merge_per_sample_columns, 
                verbose=verbose)

            merged_row_counter += 1
            row_i = merged_row   # compare the next row to the merged row
        else:
            result_rows.append(row_i)
            row_i = df.row(i + 1, named=True)

        i += 1
        
    # Add the last row if we haven't processed it yet
    result_rows.append(row_i)

    print(f"Merged {merged_row_counter:,d} out of {df.height:,d} ({merged_row_counter/df.height:.1%}) overlapping loci")

    if result_rows:
        merged_df = pl.DataFrame(result_rows, schema=df.schema)
    else:
        merged_df = pl.DataFrame(schema=df.schema)

    return merged_df


def get_sample_id(input_tsv):
    return os.path.basename(input_tsv).split(".")[0]


def create_table_of_all_loci(input_tsvs, exclude_homopolymers = False, discard_impure_genotypes = False, output_stats_tsv = None, chroms_to_process = None):
    output_stats = []
    combined_df = None
    for table_i, input_tsv in tqdm.tqdm(enumerate(input_tsvs), total=len(input_tsvs), unit=" tables"):
        sample_id = get_sample_id(input_tsv)

        df = parse_input_tsv(sample_id, input_tsv, exclude_homopolymers=exclude_homopolymers, discard_impure_genotypes=discard_impure_genotypes, chroms_to_process=chroms_to_process)

        if combined_df is None:
            locus_ids_before_join = pure_locus_ids_before_join = pure_found_in_reference_locus_ids_before_join = pure_found_in_reference_with_3_or_more_repeats_locus_ids_before_join = 0

            combined_df = df
            combined_df = combined_df.with_columns(pl.lit(True).alias("IsPureRepeat"))
        else:
            locus_ids_before_join = combined_df.height
            pure_locus_ids_before_join = combined_df.filter(pl.col("IsPureRepeat")).height
            pure_found_in_reference_locus_ids_before_join = combined_df.filter(pl.col("IsPureRepeat") & pl.col("IsFoundInReference")).height
            pure_found_in_reference_with_3_or_more_repeats_locus_ids_before_join = combined_df.filter(pl.col("IsPureRepeat") & pl.col("IsFoundInReference") & (pl.col("NumRepeatsInReference") >= 3)).height

            combined_df = combined_df.join(df, on=PER_LOCUS_COLUMNS, how="full", coalesce=True)

        combined_df = combined_df.with_columns(
            pl.col("IsPureRepeat").fill_null(True).and_(pl.col(f"IsPureRepeat:{sample_id}").fill_null(True)).alias("IsPureRepeat")
        )

        combined_df = combined_df.drop([f"IsPureRepeat:{sample_id}", f"NumRepeatsShortAllele:{sample_id}", f"NumRepeatsLongAllele:{sample_id}"])

        # Calculate stats
        new_locus_id_count = combined_df.height - locus_ids_before_join
        new_pure_locus_id_count = combined_df.filter(pl.col("IsPureRepeat")).height - pure_locus_ids_before_join
        new_pure_found_in_reference_locus_ids_before_join = combined_df.filter(pl.col("IsPureRepeat") & pl.col("IsFoundInReference")).height - pure_found_in_reference_locus_ids_before_join
        new_pure_found_in_reference_with_3_or_more_repeats_locus_ids_before_join = combined_df.filter(pl.col("IsPureRepeat") & pl.col("IsFoundInReference") & (pl.col("NumRepeatsInReference") >= 3)).height - pure_found_in_reference_with_3_or_more_repeats_locus_ids_before_join
        print(f"#{table_i+1:<3,d}: Added {sample_id:10s} with {df.height:8,d} loci which yielded {new_locus_id_count:8,d} new loci "
            f"({new_pure_locus_id_count:6,d} pure, {new_pure_found_in_reference_locus_ids_before_join:6,d} & found in reference, "
            f"{new_pure_found_in_reference_with_3_or_more_repeats_locus_ids_before_join:6,d} & >= 3 repeats in reference) "
            f"for an overall total of {combined_df.height:9,d} loci in the combined table.")

        if output_stats_tsv:
            output_stats.append({
                "Id": table_i + 1,
                "SampleId": sample_id,
                "Loci": df.height,
                "NewLoci": new_locus_id_count,
                "FractionNewLoci": new_locus_id_count/combined_df.height,
                "CumulativeTotalLoci": combined_df.height,
            })


    if output_stats_tsv:
        pl.DataFrame(output_stats).write_csv(output_stats_tsv, separator="\t", float_precision=2)
        print(f"Wrote {len(output_stats):,d} rows to {output_stats_tsv}")

    return combined_df


def parse_dipcall_vcf_file(vcf_file_path, chroms_to_process = None):

    interval_trees = collections.defaultdict(IntervalTree)

    fopen = gzip.open if vcf_file_path.endswith(".gz") else open
    with fopen(vcf_file_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            vcf_fields = line.strip().split("\t")
            vcf_chrom = vcf_fields[0].replace("chr", "").upper()
            if chroms_to_process is not None and vcf_chrom not in chroms_to_process:
                continue

            vcf_pos = int(vcf_fields[1])
            vcf_ref = vcf_fields[3].upper()
            if "N" in vcf_ref:
                continue

            vcf_alt = vcf_fields[4].upper()
            alt_alleles = vcf_alt.split(",")

            # filter out SNV alleles, '*' alleles, and alleles with N bases
            alt_alleles = [a for a in alt_alleles if len(a) != len(vcf_ref) and a != "*" and "N" not in a]
            if len(alt_alleles) == 0:
                continue

            max_indel_size = max(abs(len(a) - len(vcf_ref)) for a in alt_alleles)

            interval_trees[vcf_chrom].add(Interval(vcf_pos - 1, vcf_pos + len(vcf_ref) - 1, data={"indel_size": max_indel_size}))
        
    return interval_trees


def parse_dipcall_confidence_regions_bed_file(bed_file_path, chroms_to_process = None):

    confidence_region_interval_trees = collections.defaultdict(IntervalTree)

    bed_schema_overrides = {
        "Chrom": pl.Utf8,
        "Start0Based": pl.Int32,
        "End1Based": pl.Int32,
    }
    
    try:
        # Read the BED file with polars
        dipcall_confidence_regions = pl.read_csv(
            bed_file_path, 
            separator="\t", 
            has_header=False, 
            new_columns=["Chrom", "Start0Based", "End1Based"],
            schema_overrides=bed_schema_overrides
        )
    except pl.exceptions.NoDataError:
        print(f"WARNING: No records found in {bed_file_path}. Returning empty interval tree.")
        return confidence_region_interval_trees

    dipcall_confidence_regions = dipcall_confidence_regions.with_columns(
        pl.col("Chrom").str.replace("chr", "").str.to_uppercase()
    )

    if chroms_to_process:
        dipcall_confidence_regions = dipcall_confidence_regions.filter(pl.col("Chrom").is_in(chroms_to_process))

    # Add intervals to the interval tree
    total_bases = 0
    for row in dipcall_confidence_regions.iter_rows(named=True):
        if "_" in row["Chrom"]:
            continue
        normalized_chrom = row["Chrom"]
        confidence_region_interval_trees[normalized_chrom].add(Interval(row["Start0Based"], row["End1Based"]))
        total_bases += row["End1Based"] - row["Start0Based"]
    
    print(f"Parsed {total_bases:,d} bases and {len(confidence_region_interval_trees):,d} confidence regions from {bed_file_path}")

    return confidence_region_interval_trees

EMPTY_INTERVAL_TREE = IntervalTree()

def fill_in_missing_genotypes(df, sample_id, dipcall_vcf_indel_interval_trees, dipcall_confidence_region_interval_trees, verbose=True):

    def should_genotype_be_set_to_hom_ref(row):
        normalized_chrom = row["Chrom"]
        start = row["Start1Based"] - 1
        end = row["End1Based"]
        motif_size = row["MotifSize"]

        # if the genotype is not missing, return False
        if row[f"NumRepeatsShortAllele:{sample_id}"] is not None or row[f"NumRepeatsLongAllele:{sample_id}"] is not None:
            return False

        # if the locus is outside the DipCall confidence regions, return False
        dipcall_confidence_region_interval_tree = dipcall_confidence_region_interval_trees.get(normalized_chrom, EMPTY_INTERVAL_TREE)
        is_in_confidence_region = dipcall_confidence_region_interval_tree.overlaps(start, max(end, start + 1))
        if not is_in_confidence_region:
            return False

        # if the locus overlaps with any indels in the DipCall VCF, return False
        dipcall_vcf_indel_interval_tree = dipcall_vcf_indel_interval_trees.get(normalized_chrom, EMPTY_INTERVAL_TREE)
        overlapping_indels = dipcall_vcf_indel_interval_tree.overlap(start, max(end, start + 1))
        for overlapping_indel in overlapping_indels:
            if overlapping_indel.data["indel_size"] > 0.5 * motif_size:
                return False
        
        return True
        

    # Add boolean column indicating if genotype should be set to homozygous reference
    df = df.with_columns(
        pl.struct(["Chrom", "Start1Based", "End1Based", f"NumRepeatsShortAllele:{sample_id}", f"NumRepeatsLongAllele:{sample_id}"])
        .map_elements(should_genotype_be_set_to_hom_ref, return_dtype=pl.Boolean)
        .alias("SetGenotypeToHomRef")
    )

    both_alleles_are_missing = pl.col(f"NumRepeatsShortAllele:{sample_id}").is_null() & pl.col(f"NumRepeatsLongAllele:{sample_id}").is_null()

    if verbose:
        both_alleles_missing_in_confidence_region_df = df.filter(pl.col('SetGenotypeToHomRef') & both_alleles_are_missing)
        message = f"Setting {both_alleles_missing_in_confidence_region_df.height:,d} out of {df.height:,d} "
        message += f"({(both_alleles_missing_in_confidence_region_df.height/df.height):.1%}) genotypes to homozygous reference for sample {sample_id}"
        #message += f", leaving {df.filter(both_alleles_are_missing).height - both_alleles_missing_in_confidence_region_df.height:,d} missing genotypes for {sample_id}"
        print(message)

    df = df.with_columns(
        pl.when(pl.col("SetGenotypeToHomRef")).then(
            pl.when(both_alleles_are_missing).then(
                pl.col(f"NumRepeatsInReference")
            ).otherwise(
                pl.col(f"NumRepeatsShortAllele:{sample_id}")
            )
        ).otherwise(None).alias(f"NumRepeatsShortAllele:{sample_id}"),
        pl.when(pl.col("SetGenotypeToHomRef")).then(
            pl.when(both_alleles_are_missing).then(
                pl.col(f"NumRepeatsInReference")
            ).otherwise(
                pl.col(f"NumRepeatsLongAllele:{sample_id}")
            )
        ).otherwise(None).alias(f"NumRepeatsLongAllele:{sample_id}"),
    )

    df = df.drop("SetGenotypeToHomRef")

    return df


def generate_combined_columns(combined_df):
    
    def compute_columns(row):
        assert len(row["ShortAllelesList"]) == len(row["LongAllelesList"]), row
        
        all_non_null_alleles = []
        all_non_null_alleles.extend([allele_size for allele_size in row["ShortAllelesList"] if allele_size is not None])
        all_non_null_alleles.extend([allele_size for allele_size in row["LongAllelesList"] if allele_size is not None])

        allele_size_histogram = collections.Counter()
        for allele_size in all_non_null_alleles:
            allele_size_histogram[allele_size] += 1

        biallelic_frequency_histogram = collections.Counter()
        for allele_size_1, allele_size_2 in zip(row["ShortAllelesList"], row["LongAllelesList"]):
            if allele_size_1 is None and allele_size_2 is None:
                continue
            elif allele_size_1 is None:
                allele_size_1 = allele_size_2
            elif allele_size_2 is None:
                allele_size_2 = allele_size_1

            biallelic_frequency_histogram[(allele_size_1, allele_size_2)] += 1

        if len(all_non_null_alleles) == 0:
            return {
                "ModeAllele": None,
                "Stdev": None,
                "Median": None,
                "99thPercentile": None,
                "AlleleSizeHistogram": None,
                "BiallelicFrequencyHistogram": None,
            }

        allele_size_histogram_str = ",".join([f"{allele_size}x:{count}" for allele_size, count in allele_size_histogram.items()])
        biallelic_frequency_histogram_str = ",".join([f"{allele_size_1}/{allele_size_2}:{count}" for (allele_size_1, allele_size_2), count in biallelic_frequency_histogram.items()])

        return {
            "ModeAllele": float(allele_size_histogram.most_common(1)[0][0]),
            "Stdev": np.std(all_non_null_alleles),
            "Median": np.median(all_non_null_alleles),
            "99thPercentile": np.percentile(all_non_null_alleles, 99),
            "AlleleSizeHistogram": allele_size_histogram_str,
            "BiallelicFrequencyHistogram": biallelic_frequency_histogram_str,
        }

    # Drop intermediate columns "ShortAllelesList" and "LongAllelesList"
    combined_df = combined_df.with_columns(
        pl.struct(["ShortAllelesList", "LongAllelesList"]).map_elements(compute_columns, return_dtype=pl.Struct([
            pl.Field("ModeAllele", pl.Float64),
            pl.Field("Stdev", pl.Float64),
            pl.Field("Median", pl.Float64),
            pl.Field("99thPercentile", pl.Float64),
            pl.Field("AlleleSizeHistogram", pl.Utf8),
            pl.Field("BiallelicFrequencyHistogram", pl.Utf8),
        ])).alias("SummaryStats")
    ).unnest("SummaryStats")

    combined_df = combined_df.drop(["ShortAllelesList", "LongAllelesList"])

    #print(combined_df.head(1000))

    return combined_df


def add_sample_columns(
    combined_df, 
    sample_id_to_input_tsv, 
    sample_id_to_dipcall_vcf_indel_interval_trees,
    sample_id_to_dipcall_confidence_region_interval_trees, 
    include_per_sample_columns=False, 
    exclude_homopolymers=False, 
    discard_impure_genotypes=False, 
    chroms_to_process=None,
    verbose=False,
):
    all_allele_size_columns = []
    locus_id_to_sample_ids_with_non_missing_genotypes = collections.defaultdict(set)  # used later for merging overlapping loci

    # Add summary stat columns
    combined_df = combined_df.with_columns(
        pl.lit([], dtype=pl.List(pl.Int32)).alias("ShortAllelesList"),
        pl.lit([], dtype=pl.List(pl.Int32)).alias("LongAllelesList"),
        pl.lit(0).alias("TotalHaplotypes"),
        pl.lit(0).alias("MissingHaplotypes"),
    )

    for table_i, (sample_id, input_tsv) in tqdm.tqdm(enumerate(sample_id_to_input_tsv.items()), total=len(sample_id_to_input_tsv), unit=" tables"):
        if sample_id in {c.replace(f"NumRepeatsShortAllele:", "").replace(f"NumRepeatsLongAllele:", "") for c in combined_df.columns}:
            print(f"Skipping {sample_id} because it already exists in the combined table")
            continue

        df = parse_input_tsv(sample_id, input_tsv, exclude_homopolymers=exclude_homopolymers, discard_impure_genotypes=discard_impure_genotypes, chroms_to_process=chroms_to_process)
        df = df.drop(f"IsPureRepeat:{sample_id}")
        df = df.join(combined_df.select(PER_LOCUS_COLUMNS), on=PER_LOCUS_COLUMNS, how="right", coalesce=True)

        dipcall_vcf_indel_interval_trees = sample_id_to_dipcall_vcf_indel_interval_trees[sample_id]
        dipcall_confidence_region_interval_trees = sample_id_to_dipcall_confidence_region_interval_trees[sample_id]
        df = fill_in_missing_genotypes(df, sample_id, dipcall_vcf_indel_interval_trees, dipcall_confidence_region_interval_trees)

        # Left join to the previously combined table of all loci
        combined_df = combined_df.join(
            df, 
            on=PER_LOCUS_COLUMNS, 
            how="left",  # left join
            coalesce=True,
        )

        # Update summary stat columns
        combined_df = combined_df.with_columns(
            pl.concat_list("ShortAllelesList", f"NumRepeatsShortAllele:{sample_id}").alias("ShortAllelesList"), 
            pl.concat_list("LongAllelesList", f"NumRepeatsLongAllele:{sample_id}").alias("LongAllelesList"), 
            pl.col("TotalHaplotypes").add(2).alias("TotalHaplotypes"),
            pl.col("MissingHaplotypes").add(
                pl.when(pl.col(f"NumRepeatsShortAllele:{sample_id}").is_null()).then(1).otherwise(0) + 
                pl.when(pl.col(f"NumRepeatsLongAllele:{sample_id}").is_null()).then(1).otherwise(0)
            ).alias("MissingHaplotypes"),
        )
        
        # Record locus ids where the genotype is not missing and is not HOM-REF (eg. was not set to HOM-REF during the previous step)
        locus_ids_with_non_missing_genotypes_in_this_sample = df.filter(
            (pl.col(f"NumRepeatsShortAllele:{sample_id}").is_not_null() & (pl.col(f"NumRepeatsShortAllele:{sample_id}") != pl.col("NumRepeatsInReference")))
            | 
            (pl.col(f"NumRepeatsLongAllele:{sample_id}").is_not_null() & (pl.col(f"NumRepeatsLongAllele:{sample_id}") != pl.col("NumRepeatsInReference")))
        ).select(["LocusId"]).to_series()

        for locus_id in locus_ids_with_non_missing_genotypes_in_this_sample:
            #if locus_id == "19-29537703-29537757-CACCCACCTACTCTCCCG":
            #    print(f"19-29537703-29537757-CACCCACCTACTCTCCCG: {sample_id}")
            locus_id_to_sample_ids_with_non_missing_genotypes[locus_id].add(sample_id)

        per_sample_columns = [f"NumRepeatsShortAllele:{sample_id}", f"NumRepeatsLongAllele:{sample_id}"]
        if include_per_sample_columns:
            all_allele_size_columns.extend(per_sample_columns)
        else:
            combined_df = combined_df.drop(per_sample_columns)

    # Drop rows where TotalHaplotypes == MissingHaplotypes, keeping 10 example LocusIds from the dropped rows
    before = combined_df.height
    all_allele_sizes_are_null_df = combined_df.select(["LocusId", "TotalHaplotypes", "MissingHaplotypes"]).filter(pl.col("TotalHaplotypes") == pl.col("MissingHaplotypes"))
    combined_df = combined_df.filter(pl.col("TotalHaplotypes") != pl.col("MissingHaplotypes"))
    if before != combined_df.height:
        print(f"Dropped {before - combined_df.height:,d} rows where all allele sizes are null after filtering by confidence regions:", 
            ", ".join(all_allele_sizes_are_null_df.select(["LocusId"]).head(10).to_series().to_list()))

    # Merge overlapping rows
    combined_df = merge_overlapping_loci(
        combined_df, locus_id_to_sample_ids_with_non_missing_genotypes, merge_per_sample_columns=include_per_sample_columns, verbose=verbose)


    print("Generating combined allele frequencies and other summary columns")
    combined_df = generate_combined_columns(combined_df)

    # Convert to integer columns
    combined_df = combined_df.with_columns([
        pl.col(c).cast(pl.Int32) for c in all_allele_size_columns + ["NumRepeatsInReference", "ModeAllele"]
    ])

    print("Done adding sample columns")

    return combined_df

def run(command):
    print(f"Running command: {command}")
    subprocess.check_output(command, shell=True)

def main():
    parser = argparse.ArgumentParser(description="Perform an outer-join on multiple per-sample tables")
    parser.add_argument("-o", "--output-tsv", help="Combined tsv file path")
    parser.add_argument("--include-per-sample-columns", action="store_true", help="If specified, per-sample genotype columns will be included in the output TSV. "
                        "Otherwise, only the summary columns will be included.")
    parser.add_argument("-n", type=int, help="Number of tables to process")
    parser.add_argument("--exclude-homopolymers", action="store_true", help="Exclude homopolymers from the output")
    parser.add_argument("--discard-impure-genotypes", action="store_true", help="Discard genotypes that are not pure repeats")
    parser.add_argument("--output-stats-tsv", action="store_true", help="If specified, will output a table with stats")
    parser.add_argument("-b", "--batches", type=int, default=1, choices=[1, 2, 3, 4, 6, 8, 12, 24], help="When the nubmer of input sample tables if large enough to "
                        "cause out-of-memory errors, split the input tables into batches by chromosome(s). ")
    parser.add_argument("input_tsv_bed_and_vcf_files", nargs="+", help="This should be a list of the input STR variant TSV files produced by the " 
                        "filter_vcf_to_STR_variants.py script, as well as the VCF files and confidence region BED files produced by DipCall for the same samples. "
                        "The filenames should have the format: '<sample_id>.*.tsv.gz', '<sample_id>.*.vcf.gz', and '<sample_id>.*.bed.gz'")
    parser.add_argument("--force", action="store_true", help="Force re-run even if the output file already exists")
    parser.add_argument("--verbose", action="store_true", help="Verbose output")
    args = parser.parse_args()


    # Parse and validate TSV and BED file paths
    sample_id_to_input_tsv = {}
    sample_id_to_dipcall_vcf = {}
    sample_id_to_dipcall_confidence_regions_bed = {}
    for input_file_path in args.input_tsv_bed_and_vcf_files:
        if not os.path.exists(input_file_path):
            parser.error(f"Input file {input_file_path} does not exist. It must be a tsv file or a bed file.")

        sample_id = get_sample_id(input_file_path)
        if input_file_path.endswith(".tsv") or input_file_path.endswith(".tsv.gz"):
            container = sample_id_to_input_tsv
        elif input_file_path.endswith(".vcf") or input_file_path.endswith(".vcf.gz"):
            container = sample_id_to_dipcall_vcf
        elif input_file_path.endswith(".bed") or input_file_path.endswith(".bed.gz"):
            container = sample_id_to_dipcall_confidence_regions_bed
        else:
            parser.error(f"Input file {input_file_path} must have a .tsv or .tsv.gz or .bed or .bed.gz suffix.")

        if sample_id in container:
            parser.error(f"Multiple input files found for sample {sample_id}: {container[sample_id]} and {input_file_path}")

        container[sample_id] = input_file_path

    # Sort input tables by file size (largest first)
    input_tsvs = list(sample_id_to_input_tsv.values())
    input_tsvs.sort(key=os.path.getsize, reverse=True)

    if args.n:
        input_tsvs = input_tsvs[:args.n]
        sample_id_to_input_tsv = {k: v for k, v in sample_id_to_input_tsv.items() if v in set(input_tsvs)}

    # Check that all samples have a VCF and a BED file
    sample_ids_without_vcf = set(sample_id_to_input_tsv) - set(sample_id_to_dipcall_vcf)
    if len(sample_id_to_dipcall_vcf) > 0 and sample_ids_without_vcf:
        parser.error(f"A DipCall VCF file was not provided for the following sample(s): {', '.join(sample_ids_without_vcf)}")

    sample_ids_without_bed = set(sample_id_to_input_tsv) - set(sample_id_to_dipcall_confidence_regions_bed)
    if len(sample_id_to_dipcall_confidence_regions_bed) > 0 and sample_ids_without_bed:
        parser.error(f"A DipCall confidence region BED file was not provided for the following sample(s): {', '.join(sample_ids_without_bed)}")


    # Compute output path
    if not args.output_tsv:
        exclude_homopolymers_suffix = ".excluding_homopolymers" if args.exclude_homopolymers else ""
        discard_impure_genotypes_suffix = ".discarding_impure_genotypes" if args.discard_impure_genotypes else ""
        args.output_tsv = f"joined.{len(input_tsvs)}_tables{exclude_homopolymers_suffix}{discard_impure_genotypes_suffix}.tsv.gz"


    filename_prefix = args.output_tsv.replace(".tsv.gz", "").replace(".tsv", "")

    batch_tsvs = []
    for batch_i in range(0, args.batches):

        if args.batches > 10:
            batch_suffix = f".batch_{batch_i+1:02d}_of_{args.batches:02d}"
        elif args.batches > 1:
            batch_suffix = f".batch_{batch_i+1}_of_{args.batches}"
        else:
            batch_suffix = ""

        chroms_to_process_in_current_batch = None  # by default, process all chromosomes
        if args.batches > 1:
            chroms_to_process_in_current_batch = ALL_CHROMS[batch_i::args.batches]
            print("-" * 200)
            print(f"Processing batch [{batch_i+1} out of {args.batches}] with chromosomes {','.join(chroms_to_process_in_current_batch)}")

        output_stats_tsv = f"{filename_prefix}{batch_suffix}.stats.tsv" if args.output_stats_tsv else None
        args.output_tsv = f"{filename_prefix}{batch_suffix}.tsv.gz"
        print(f"Output file: {args.output_tsv}")
        batch_tsvs.append(args.output_tsv)

        # Perform an outer-join on all the input tables to create a table of all loci
        table_of_all_loci_path = f"{filename_prefix}{batch_suffix}.only_loci.tsv.gz"
        if not os.path.exists(table_of_all_loci_path) or args.force:
            combined_df = create_table_of_all_loci(
                input_tsvs=input_tsvs, 
                exclude_homopolymers=args.exclude_homopolymers, 
                discard_impure_genotypes=args.discard_impure_genotypes, 
                output_stats_tsv=output_stats_tsv,
                chroms_to_process=chroms_to_process_in_current_batch,
            )
            
            with gzip.open(table_of_all_loci_path, "wb") as f:
                combined_df.write_csv(f, separator="\t", float_precision=2)

            print(f"Wrote {combined_df.height:,d} loci to {table_of_all_loci_path}")
        else:
            combined_df = pl.read_csv(table_of_all_loci_path, separator="\t", columns=PER_LOCUS_COLUMNS + ["IsPureRepeat"], schema_overrides=SCHEMA)
            combined_df = combined_df.with_columns(
                pl.col("NumRepeatsInReference").cast(pl.Int32), 
                pl.col("MotifSize").cast(pl.Int32),
            )

            print(f"Read {combined_df.height:,d} loci from {table_of_all_loci_path}")

        # Parse the DipCall confidence region BED files
        sample_id_to_dipcall_vcf_indel_interval_trees = {}
        sample_id_to_dipcall_confidence_region_interval_trees = {}
        for sample_id, vcf_file_path in sample_id_to_dipcall_vcf.items():
            sample_id_to_dipcall_vcf_indel_interval_trees[sample_id] = parse_dipcall_vcf_file(
                vcf_file_path, 
                chroms_to_process=chroms_to_process_in_current_batch)

        for sample_id, bed_file_path in sample_id_to_dipcall_confidence_regions_bed.items():
            sample_id_to_dipcall_confidence_region_interval_trees[sample_id] = parse_dipcall_confidence_regions_bed_file(
                bed_file_path, 
                chroms_to_process=chroms_to_process_in_current_batch)

        # Add sample columns and/or allele size summary statistics to the combined table
        combined_df = add_sample_columns(
            combined_df, 
            sample_id_to_input_tsv, 
            sample_id_to_dipcall_vcf_indel_interval_trees,
            sample_id_to_dipcall_confidence_region_interval_trees, 
            include_per_sample_columns=args.include_per_sample_columns,
            exclude_homopolymers=args.exclude_homopolymers, 
            discard_impure_genotypes=args.discard_impure_genotypes, 
            chroms_to_process=chroms_to_process_in_current_batch,
            verbose=args.verbose,
        )

        # Sort by PER_LOCUS_COLUMNS
        combined_df = combined_df.sort(PER_LOCUS_COLUMNS)

        # Write as gzipped file
        with gzip.open(args.output_tsv, "wb") as f:
            combined_df.write_csv(f, separator="\t", float_precision=2)

        print(f"Wrote combined table with {combined_df.height:,d} loci to {args.output_tsv}")

    # Concatenate all the batch files into a single file
    if args.batches > 1:
        run(f"gunzip -c {batch_tsvs[0]} | head -n 1 > {filename_prefix}.tsv")
        for batch_tsv in batch_tsvs:
            run(f"gunzip -c {batch_tsv} | tail -n +2 >> {filename_prefix}.tsv") 

        combined_df = pl.read_csv(f"{filename_prefix}.tsv", separator="\t", schema_overrides=SCHEMA)
        combined_df = combined_df.sort(PER_LOCUS_COLUMNS)
        with gzip.open(f"{filename_prefix}.tsv.gz", "wb") as f:
            combined_df.write_csv(f, separator="\t", float_precision=2)
        print(f"Wrote {combined_df.height:,d} loci to {filename_prefix}.tsv.gz")

        # gzip the file
        run(f"rm -f {filename_prefix}.tsv")

if __name__ == "__main__":
    main()