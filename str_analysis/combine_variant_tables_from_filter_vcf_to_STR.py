"""When filter_vcf_to_STR.py is run on multiple samples, this script can be used to combine the resulting *.variants.tsv.gz files into a single combined table 
with one record per locus. It takes the per-sample tables as well as DipCall confidence regions as input and outputs a single combined table with one record per locus by:
1. performing an outer-join on the loci in the per-sample tables
2. setting missing genotypes in each sample to homozygous reference for loci that are within the DipCall confidence regions of that sample (based on the assumption that 
if no tandem repeat variant was detected there in the original VCF, then the locus is homozygous reference).
"""

import argparse
import collections
import gzip
from intervaltree import IntervalTree, Interval
import numpy as np
import os
import polars as pl
import tqdm

BATCH_SIZE = 100

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

#ALL_LOCUS_TABLE_CHECKPOINT_FILENAME = "temp_df_all_loci.tsv.gz"
#SAMPLES_TABLE_CHECKPOINT_FILENAME = "temp_df_with_samples.tsv.gz"

def parse_input_tsv(sample_id, input_tsv, exclude_homopolymers = False, discard_impure_genotypes = False):

    # Read with polars
    df = pl.read_csv(input_tsv, separator="\t", columns=PER_LOCUS_COLUMNS + SAMPLE_SPECIFIC_COLUMNS, schema_overrides=SCHEMA)
    df = df.with_columns(
        pl.col("NumRepeatsInReference").cast(pl.Int32),
        pl.col("MotifSize").cast(pl.Int32),
        pl.col("NumRepeatsShortAllele").cast(pl.Int32),
        pl.col("NumRepeatsLongAllele").cast(pl.Int32),
    )

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

    return df


def checkpoint_combined_df(df, temp_table_path=None):
    # write to gzipped temp table and then read it back in
    if temp_table_path is None:
        temp_table_path = f"temp_df.tsv.gz"

    with gzip.open(temp_table_path, "wb") as f:
        df.write_csv(f, separator="\t")
    print(f"Wrote {df.height:,d} loci to {temp_table_path}")

    df = pl.read_csv(temp_table_path, separator="\t", schema_overrides=SCHEMA_WITH_INT)
    print(f"Read {df.height:,d} loci from {temp_table_path}")

    return df


def get_sample_id(input_tsv):
    return os.path.basename(input_tsv).split(".")[0]


def create_table_of_all_loci(input_tsvs, exclude_homopolymers = False, discard_impure_genotypes = False, output_stats_tsv = None):
    output_stats = []
    combined_df = None
    for table_i, input_tsv in tqdm.tqdm(enumerate(input_tsvs), total=len(input_tsvs), unit=" tables"):
        sample_id = get_sample_id(input_tsv)

        df = parse_input_tsv(sample_id, input_tsv, exclude_homopolymers=exclude_homopolymers, discard_impure_genotypes=discard_impure_genotypes)

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

        #if table_i % BATCH_SIZE == 0 and table_i > 0:
        #    combined_df = checkpoint_combined_df(combined_df, temp_table_path=ALL_LOCUS_TABLE_CHECKPOINT_FILENAME)


    if output_stats_tsv:
        pl.DataFrame(output_stats).write_csv(output_stats_tsv, separator="\t", float_precision=2)
        print(f"Wrote {len(output_stats):,d} rows to {output_stats_tsv}")

    return combined_df


def parse_dipcall_confidence_regions_bed_file(bed_file_path):

    bed_schema_overrides = {
        "Chrom": pl.Utf8,
        "Start0Based": pl.Int32,
        "End1Based": pl.Int32,
    }
    
    confidence_region_interval_trees = collections.defaultdict(IntervalTree)

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

    # Add intervals to the interval tree
    total_bases = 0
    for row in dipcall_confidence_regions.iter_rows(named=True):
        if "_" in row["Chrom"]:
            continue
        normalized_chrom = row["Chrom"].replace("chr", "").upper()
        confidence_region_interval_trees[normalized_chrom].add(Interval(row["Start0Based"], row["End1Based"]))
        total_bases += row["End1Based"] - row["Start0Based"]
    
    print(f"Parsed {total_bases:,d} bases and {len(confidence_region_interval_trees):,d} confidence regions from {bed_file_path}")

    return confidence_region_interval_trees

EMPTY_INTERVAL_TREE = IntervalTree()

def fill_in_missing_genotypes(df, sample_id, dipcall_confidence_region_interval_trees, verbose=True):

    def does_locus_overlap_confidence_region(row):
        normalized_chrom = row["Chrom"].replace("chr", "").upper()
        start = row["Start1Based"] - 1
        end = row["End1Based"]

        interval_tree = dipcall_confidence_region_interval_trees.get(normalized_chrom, EMPTY_INTERVAL_TREE)
        return interval_tree.overlaps(start, max(end, start + 1))

    # Add boolean column indicating if genotype should be set to homozygous reference
    df = df.with_columns(
        pl.struct(["Chrom", "Start1Based", "End1Based", f"NumRepeatsShortAllele:{sample_id}", f"NumRepeatsLongAllele:{sample_id}"])
        .map_elements(does_locus_overlap_confidence_region, return_dtype=pl.Boolean)
        .alias("DoesLocusOverlapConfidenceRegion")
    )

    both_alleles_are_missing = pl.col(f"NumRepeatsShortAllele:{sample_id}").is_null() & pl.col(f"NumRepeatsLongAllele:{sample_id}").is_null()

    if verbose:
        both_alleles_missing_in_confidence_region_df = df.filter(pl.col('DoesLocusOverlapConfidenceRegion') & both_alleles_are_missing)
        message = f"Setting {both_alleles_missing_in_confidence_region_df.height:,d} out of {df.height:,d} "
        message += f"({(both_alleles_missing_in_confidence_region_df.height/df.height):.1%}) genotypes to homozygous reference for sample {sample_id}"
        #message += f", leaving {df.filter(both_alleles_are_missing).height - both_alleles_missing_in_confidence_region_df.height:,d} missing genotypes for {sample_id}"
        print(message)

    df = df.with_columns(
        pl.when(pl.col("DoesLocusOverlapConfidenceRegion")).then(
            pl.when(both_alleles_are_missing).then(
                pl.col(f"NumRepeatsInReference")
            ).otherwise(
                pl.col(f"NumRepeatsShortAllele:{sample_id}")
            )
        ).otherwise(None).alias(f"NumRepeatsShortAllele:{sample_id}"),
        pl.when(pl.col("DoesLocusOverlapConfidenceRegion")).then(
            pl.when(both_alleles_are_missing).then(
                pl.col(f"NumRepeatsInReference")
            ).otherwise(
                pl.col(f"NumRepeatsLongAllele:{sample_id}")
            )
        ).otherwise(None).alias(f"NumRepeatsLongAllele:{sample_id}"),
    )

    df = df.drop("DoesLocusOverlapConfidenceRegion")

    return df


def generate_summary_stats(combined_df):
    
    def compute_summary_stats(row):
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
        pl.struct(["ShortAllelesList", "LongAllelesList"]).map_elements(compute_summary_stats, return_dtype=pl.Struct([
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


def add_sample_columns(combined_df, sample_id_to_input_tsv, sample_id_to_dipcall_confidence_region_interval_trees, include_per_sample_columns=False, exclude_homopolymers=False, discard_impure_genotypes=False):
    all_allele_size_columns = []

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

        df = parse_input_tsv(sample_id, input_tsv, exclude_homopolymers=exclude_homopolymers, discard_impure_genotypes=discard_impure_genotypes)
        df = df.drop(f"IsPureRepeat:{sample_id}")
        df = df.join(combined_df.select(PER_LOCUS_COLUMNS), on=PER_LOCUS_COLUMNS, how="right", coalesce=True)

        dipcall_confidence_region_interval_trees = sample_id_to_dipcall_confidence_region_interval_trees[sample_id]
        df = fill_in_missing_genotypes(df, sample_id, dipcall_confidence_region_interval_trees)

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
            pl.col("MissingHaplotypes").add(pl.col(f"NumRepeatsShortAllele:{sample_id}").is_null() + pl.col(f"NumRepeatsLongAllele:{sample_id}").is_null()).alias("MissingHaplotypes"),
        )

        per_sample_columns = [f"NumRepeatsShortAllele:{sample_id}", f"NumRepeatsLongAllele:{sample_id}"]
        if include_per_sample_columns:
            #if table_i % BATCH_SIZE == 0 and table_i > 0:
            #    combined_df = checkpoint_combined_df(combined_df, temp_table_path=SAMPLES_TABLE_CHECKPOINT_FILENAME)

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

    print("Generating summary stats")
    combined_df = generate_summary_stats(combined_df)

    # Convert to integer columns
    for c in all_allele_size_columns + ["NumRepeatsInReference", "ModeAllele"]:
        combined_df = combined_df.with_columns(
            pl.col(c).cast(pl.Int32)
        )

    return combined_df

def main():
    parser = argparse.ArgumentParser(description="Perform an outer-join on multiple per-sample tables")
    parser.add_argument("-o", "--output-tsv", help="Combined tsv file path")
    parser.add_argument("--include-per-sample-columns", action="store_true", help="If specified, per-sample genotype columns will be included in the output TSV. "
                        "Otherwise, only the summary columns will be included.")
    parser.add_argument("-n", type=int, help="Number of tables to process")
    parser.add_argument("--exclude-homopolymers", action="store_true", help="Exclude homopolymers from the output")
    parser.add_argument("--discard-impure-genotypes", action="store_true", help="Discard genotypes that are not pure repeats")
    parser.add_argument("--output-stats-tsv", action="store_true", help="If specified, will output a table with stats")
    parser.add_argument("input_tsv_and_bed_files", nargs="+", help="Input STR variant TSV files and dipcall confidence region BED files")
    parser.add_argument("--force", action="store_true", help="Force re-run even if the output file already exists")
    args = parser.parse_args()

    
    sample_id_to_input_tsv = {}
    dipcall_confidence_regions_bed_files = []
    for input_tsv_or_bed_file in args.input_tsv_and_bed_files:
        if not os.path.exists(input_tsv_or_bed_file):
            parser.error(f"Input file {input_tsv_or_bed_file} does not exist. It must be a tsv file or a bed file.")

        if input_tsv_or_bed_file.endswith(".tsv") or input_tsv_or_bed_file.endswith(".tsv.gz"):
            sample_id = get_sample_id(input_tsv_or_bed_file)
            sample_id_to_input_tsv[sample_id] = input_tsv_or_bed_file
        elif input_tsv_or_bed_file.endswith(".bed") or input_tsv_or_bed_file.endswith(".bed.gz"):
            dipcall_confidence_regions_bed_files.append(input_tsv_or_bed_file)
        else:
            parser.error(f"Input file {input_tsv_or_bed_file} must have a .tsv or .tsv.gz or .bed or .bed.gz suffix.")

    # Sort by file size (largest first)
    input_tsvs = list(sample_id_to_input_tsv.values())
    input_tsvs.sort(key=os.path.getsize, reverse=True)

    if args.n:
        input_tsvs = input_tsvs[:args.n]
        sample_id_to_input_tsv = {k: v for k, v in sample_id_to_input_tsv.items() if v in set(input_tsvs)}

    sample_id_to_bed_file_path = {}
    for sample_id in sample_id_to_input_tsv.keys():
        # find the dipcall confidence region bed file for this sample
        matching_bed_file_path = None
        for bed_file_path in dipcall_confidence_regions_bed_files:
            bed_filename = os.path.basename(bed_file_path)
            if bed_filename.startswith(f"{sample_id}."):
                if matching_bed_file_path is not None:
                    parser.error(f"Multiple DipCall confidence region bed files found for sample {sample_id}: {matching_bed_file_path} and {bed_file_path}")

                matching_bed_file_path = bed_file_path
                
        if matching_bed_file_path is None:
            parser.error(f"No dipcall confidence region bed file found for sample {sample_id}")

        if not os.path.isfile(matching_bed_file_path):
            parser.error(f"Dipcall confidence region bed file {matching_bed_file_path} does not exist")

        sample_id_to_bed_file_path[sample_id] = matching_bed_file_path


    if not args.output_tsv:
        exclude_homopolymers_suffix = ".excluding_homopolymers" if args.exclude_homopolymers else ""
        discard_impure_genotypes_suffix = ".discarding_impure_genotypes" if args.discard_impure_genotypes else ""
        args.output_tsv = f"joined.{len(input_tsvs)}_tables{exclude_homopolymers_suffix}{discard_impure_genotypes_suffix}.tsv.gz"

    if not args.output_tsv.endswith(".gz"):
        args.output_tsv += ".gz"

    filename_prefix = args.output_tsv.replace(".tsv.gz", "").replace(".tsv", "")
    if args.output_stats_tsv:
        output_stats_tsv = f"{filename_prefix}.stats.tsv"

    table_of_all_loci_path = f"{filename_prefix}.only_loci.tsv.gz"
    if not os.path.exists(table_of_all_loci_path) or args.force:
        combined_df = create_table_of_all_loci(
            input_tsvs=input_tsvs, 
            exclude_homopolymers=args.exclude_homopolymers, 
            discard_impure_genotypes=args.discard_impure_genotypes, 
            output_stats_tsv=output_stats_tsv,
        )
        
        with gzip.open(table_of_all_loci_path, "wb") as f:
            combined_df.write_csv(f, separator="\t", float_precision=2)

        print(f"Wrote {combined_df.height:,d} loci to {table_of_all_loci_path}")
    else:
        #if os.path.isfile(SAMPLES_TABLE_CHECKPOINT_FILENAME):
        #    table_of_all_loci_path = SAMPLES_TABLE_CHECKPOINT_FILENAME
        #    combined_df = pl.read_csv(table_of_all_loci_path, separator="\t", schema_overrides=SCHEMA)
        #else:
        combined_df = pl.read_csv(table_of_all_loci_path, separator="\t", columns=PER_LOCUS_COLUMNS + ["IsPureRepeat"], schema_overrides=SCHEMA)
        combined_df = combined_df.with_columns(
            pl.col("NumRepeatsInReference").cast(pl.Int32), 
            pl.col("MotifSize").cast(pl.Int32),
        )

        print(f"Read {combined_df.height:,d} loci from {table_of_all_loci_path}")

    sample_id_to_dipcall_confidence_region_interval_trees = {}
    for sample_id, bed_file_path in sample_id_to_bed_file_path.items():
        sample_id_to_dipcall_confidence_region_interval_trees[sample_id] = parse_dipcall_confidence_regions_bed_file(bed_file_path)

    combined_df = add_sample_columns(
        combined_df, 
        sample_id_to_input_tsv, 
        sample_id_to_dipcall_confidence_region_interval_trees, 
        include_per_sample_columns=args.include_per_sample_columns,
        exclude_homopolymers=args.exclude_homopolymers, 
        discard_impure_genotypes=args.discard_impure_genotypes, 
    )

    # Sort by PER_LOCUS_COLUMNS
    combined_df = combined_df.sort(PER_LOCUS_COLUMNS)

    # Write as gzipped file
    with gzip.open(args.output_tsv, "wb") as f:
        combined_df.write_csv(f, separator="\t", float_precision=2)

    print(f"Wrote combined table with {combined_df.height:,d} loci to {args.output_tsv}")


if __name__ == "__main__":
    main()