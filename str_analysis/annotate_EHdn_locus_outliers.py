"""This script takes locus outlier results table(s) from EHdn and prioritizes outliers based on the following criteria:
1. if it's a known disease-associated locus
2. genomic region: coding with high pLI / coding / 5'UTR / 3'UTR / intronic
3. motif: matches known disease-associated motif / is 3 to 6bp / is multiple of 3bp / other
4. matches reference repeats: the EHdn interval overlaps a repeat locus in h38 with the same motif

and also whether it overlaps other intervals of interest such as exome capture regions.


Case-control analysis:
$1        contig : chrUn_KI270519v1
$2         start : 136653
$3           end : 137789
$4         motif : AATGGAATGGAGTGGAGTGG
$5        pvalue : 0.5223795458583765
$6   bonf_pvalue : 1.0
$7        counts : BEG_0761-1_01:2.6302504026354256,BEG_0761-2_01:2.0236969484185106,BEG_0761-3_01:5.655289716917204,BEG_0992-1_01:1.0656652672860623,BEG_0992-2_01:1.8669790760207234,BEG_0992-3_01:0.9827417317897249,BEG_1025-2_1:1.0703462056838686,BEG_1025-3_1:1.0588707442337135,BEG_1072-2_D3970-1_D1:2.0322522539089705,BEG_1072-3_D3967-1_D1:0.9464268183109019,BEG_1078-1_01:0.8174822673656919,BEG_1078-3_01:4.078545372332545,BEG_1111-1_01:2.9641093165118155,BEG_1230-1_01:2.1798673382240388,BEG_1230-2_01:4.8926560193740265,BEG_1230-3_01:2.622230158603387,BEG_1291-1_01:1.0012043496953118,BEG_1291-3_01:0.8110253310598092,BEG_14-3_1:4.651555972877035,BEG_14-4_1:2.0959305591064936,BEG_14-5_1:1.492702415941337,BEG_1435-1_01:1.7661506831139706,BEG_1435-2_01:1.0608670367243718,BEG_1435-3_01:1.8982736456018172,BEG_21-1_1:0.9181294948126245,BEG_21-3_1:3.059935115979922,BEG_21-4_1:3.6129365138912375,BEG_851-1_1:1.9994729921251986,BEG_851-2_1:1.091843393714633,BEG_851-3_1:2.7018847187213555,BEG_887-1_1:2.983503667458916,BEG_887-3_1:4.132964622821529,BEG_894-1_1:1.8882976047093754,BEG_894-2_1:8.014903308986883,BEG_894-3_1:5.1878261220671575,BEG_916-1_1:2.4915238551627743,BEG_916-2_1:2.20325985939489


Outlier analysis:
$1             contig : chr12
$2              start : 9893619
$3                end : 9893620
$4              motif : AAGGGAGAGAGAGGGAGG
$5    top_case_zscore : 1.05
$6   high_case_counts : BEG_1230-1_01:1.09
$7             counts : 1.09
---
$1             contig : chr18
$2              start : 43105169
$3                end : 43105577
$4              motif : AAGGGAGG
$5    top_case_zscore : 2.00
$6   high_case_counts : BEG_851-1_1:9.00,BEG_0761-1_01:7.01
$7             counts : 7.01,4.05,0.95,5.42,9.00,2.18
"""
import argparse
import collections
import os
import pandas as pd
from pprint import pprint
import re
from tqdm import tqdm

from intervaltree import IntervalTree, Interval
from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif

from str_analysis.utils.gtf_utils import compute_genomic_region_of_interval
from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.file_utils import open_file, file_exists


def parse_args():
    parser = argparse.ArgumentParser(description="Annotate locus outlier results from EHdn")
    parser.add_argument("--known-disease-associated-loci", help="Path or url of ExpansionHunter catalog .json file that "
                        "contains all known disease-associated loci, such as the catalog available @ "
                        "https://github.com/broadinstitute/str-analysis/blob/main/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json",
                        default="variant_catalog_without_offtargets.GRCh38.json")
    parser.add_argument("--genes-gtf", help="Gene models gtf file path or url. This can be specified more than once "
                        "to annotate with multiple gene sources (ie. Gencode and MANE). These gtf files can be "
                        "downloaded from https://www.gencodegenes.org/human/ and https://www.ncbi.nlm.nih.gov/refseq/MANE/",
                        action="append", default=[])
    parser.add_argument("--reference-tr-bed-file", help="A catalog of all repeats in the reference genome in BED "
                        "format where the chrom, start, and end represent the repeat interval in 0-based coordinates, "
                        "and the name field (column #4) contains the repeat unit. If an EHdn call overlaps an interval "
                        "in this bed file and shares the same normalized motif, the locus id of this matching interval "
                        "will be recorded in this column. When no matching interval is found, the value with be N/A."
                        "TandemRepeatFinder reference catalogs of various sizes are publicly available at "
                        "https://console.cloud.google.com/storage/browser/str-truth-set/hg38/ref/other/ with file names "
                        "like repeat_specs*.bed.gz",
                        required=True)
    parser.add_argument("--overlap-bed-file", action="append", help="BED file containing regions of interest. "
                        "This option can be specified more than once. For every bed file, a column will be added "
                        "with the same name as the filename, and containing True if a given EHdn call overlaps "
                        "an interval in this bed file.", default=[])
    parser.add_argument("--tr-bed-file", action="append", help="BED file that contains tandem repeat loci and a name"
                        "field (column #4) that contains the repeat unit. This option can be specified more than once. "
                        "For every bed file, a column will be added with the same name as the filename. If an EHdn call "
                        "overlaps an interval in this bed file and shares the same normalized motif, the locus id of "
                        "this matching interval will be recorded in this column. When no matching interval is found, "
                        "the value with be N/A.", default=[])

    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("-o", "--output-dir", help="Output directory")

    parser.add_argument("ehdn_locus_outlier_tsv", nargs="+", help="One or more EHdn locus outlier results .tsv file(s)")
    args = parser.parse_args()

    for file_path in args.overlap_bed_file + args.tr_bed_file + args.genes_gtf + [args.known_disease_associated_loci,
                                                                                  args.reference_tr_bed_file]:
        if not file_exists(os.path.expanduser(file_path)):
            parser.error(f"File not found: {file_path}")

    return args, parser


KNOWN_DISEASE_ASSOCIATED_LOCI_COLUMNS = [
    'LocusId', 'LocusStructure', 'RepeatUnit', 'MainReferenceRegion', 'Gene', 'Inheritance', 'GeneRegion', 'GeneId']
LOCUS_OUTLIER_COLUMNS = {"contig", "start", "end", "motif", "top_case_zscore", "high_case_counts", "counts"}
CASE_CONTROL_COLUMNS = {"contig", "start", "end", "motif", "pvalue", "bonf_pvalue", "counts"}

OVERLAP_MARGIN = 600  # bp  (fragment length)


def get_overlapping_interval_generator(interval_tree, require_motif_match=False):
    """This function takes a dictionary that maps chromosome name to IntervalTree of intervals on that chromosome.
    The intervals must either have 1) data == None or 2) data == 2-tuple (motif, locus_id).
    In the case of 1) the resulting function will check for simple overlap between EHdn outlier record(s) and the
    Intervals in the IntervalTree and return True if there's an overlap. In the case of 2) they will also check
    whether the motifs match and, if they do, return the locus id of the matching Interval.

    Args:
        interval_tree (dict): a dictionary that maps chromosome names to IntervalTrees
        require_motif_match (bool): if True, only consider an interval as overlapping if its motif also matches
    """
    def get_overlapping_interval(outlier_table_row):
        chrom = outlier_table_row["contig"].replace("chr", "")
        start = int(outlier_table_row["start"])
        end = int(outlier_table_row["end"])

        for locus_interval in interval_tree[chrom].overlap(start - OVERLAP_MARGIN, end + OVERLAP_MARGIN):
            if not require_motif_match:
                return True

            interval_canonical_motif = compute_canonical_motif(locus_interval.data[0])
            if interval_canonical_motif != outlier_table_row["CanonicalMotif"]:
                continue
            matching_known_disease_associated_locus_id = locus_interval.data[1]
            return matching_known_disease_associated_locus_id

        if not require_motif_match:
            return False
        else:
            return None

    return get_overlapping_interval


def parse_known_disease_associated_loci(args, parser):
    try:
        known_disease_loci_df = pd.read_json(args.known_disease_associated_loci)
        known_disease_loci_df = known_disease_loci_df[KNOWN_DISEASE_ASSOCIATED_LOCI_COLUMNS]
    except Exception as e:
        parser.error(f"Couldn't read known disease-associated loci catalog from {args.known_disease_associated_loci}: {e}")

    canonical_motifs = set()

    # generate an IntervalTree for known disease-associated loci
    known_disease_loci_it = collections.defaultdict(IntervalTree)
    for _, row in known_disease_loci_df.iterrows():
        chrom, start, end = parse_interval(row["MainReferenceRegion"])
        chrom = chrom.replace("chr", "")
        i = Interval(start, end, data=(row["RepeatUnit"], row["LocusId"]))
        known_disease_loci_it[chrom].add(i)

        canonical_motifs.add(compute_canonical_motif(row["RepeatUnit"]))

    return known_disease_loci_it, canonical_motifs


def parse_sample_id_and_counts_column(counts_value):
    """Splits the "high_case_counts" column in the outlier analysis or the "counts" column in the Case-control analysis
    into a list of 3-tuples (sample id, z-score, rank)"""
    counts_list = counts_value.split(",")
    counts_list = [entry.split(":") for entry in counts_list]
    counts_list = [(sample_id, float(count)) for sample_id, count in counts_list]
    counts_list.sort(key=lambda x: x[1], reverse=True)
    for i, (sample_id, count) in enumerate(counts_list):
        yield {
            "SampleId": sample_id,
            "NormalizedCount": count,
            "SampleRankAtLocus": i + 1,
            "TotalSamplesAtLocus": len(counts_list),
        }


def parse_bed_to_interval_tree(bed_file_path, name_field_is_repeat_unit=False, verbose=False):
    interval_tree = collections.defaultdict(IntervalTree)

    f = open_file(bed_file_path, download_local_copy_before_opening=True)
    if verbose:
        print(f"Parsing {bed_file_path}", "assuming the name field contains the repeat unit" if name_field_is_repeat_unit else "")
        f = tqdm(f, unit=" records", unit_scale=True)

    valid_nucleotides = set("ACGTN")

    counter = 0
    for i, line in enumerate(f):
        fields = line.strip().split("\t")
        chrom = fields[0].replace("chr", "")
        start_0based = int(fields[1])
        end_1based = int(fields[2])
        if not name_field_is_repeat_unit:
            interval = Interval(start_0based + 1, end_1based, data=None)
        else:
            repeat_unit = fields[3].strip("()*")
            if set(repeat_unit) - valid_nucleotides:
                print(f"WARNING: Unexpected characters in motif '{repeat_unit}' in line #{i + 1}: {line}. Skipping...")
                continue
            locus_id = f"{chrom}-{start_0based + 1}-{end_1based}-{repeat_unit}"
            interval = Interval(start_0based + 1, end_1based, data=(repeat_unit, locus_id))

        counter += 1
        interval_tree[chrom].add(interval)
    f.close()
    print(f"Finished parsing {counter:,d} rows from {bed_file_path}")

    return interval_tree


def compute_genomic_region_of_row(row, verbose=False):
    matched_reference_TR = row["MatchedReferenceTR"]
    if matched_reference_TR and not pd.isna(matched_reference_TR):
        chrom, start, end = parse_interval(matched_reference_TR)
    else:
        chrom, start, end = row["contig"], row["start"], row["end"]
    return compute_genomic_region_of_interval(chrom, start, end, verbose=verbose)


def main():
    args, parser = parse_args()

    if args.verbose:
        print("Args:")
        for key, value in sorted(args.__dict__.items(), key=lambda x: str(x[1])):
            key += " = "
            if isinstance(value, list):
                print(f"   {key:30s}")
                for v in value:
                    print(f"       {v}")
            else:
                print(f"   {key:30s} {str(value)}")

    # parse known disease-associated loci
    known_disease_associated_loci_interval_tree, known_disease_associated_motifs = parse_known_disease_associated_loci(
        args, parser)

    known_disease_associated_loci_column_func = get_overlapping_interval_generator(
        known_disease_associated_loci_interval_tree, require_motif_match=True)

    # prepare functions that will be used to annotate the outlier tables
    path_to_column_func = {}

    current_interval_tree = parse_bed_to_interval_tree(
        args.reference_tr_bed_file, name_field_is_repeat_unit=True, verbose=args.verbose)
    path_to_column_func["MatchedReferenceTR"] = get_overlapping_interval_generator(
        current_interval_tree, require_motif_match=True)

    for bed_file_path_list, contains_TRs in [(args.tr_bed_file, True), (args.overlap_bed_file, False)]:
        for bed_file_path in bed_file_path_list:
            # compute the column name
            column_name = re.sub(".bed(.b?gz)?$", "", os.path.basename(bed_file_path)).split(".")[0]
            if contains_TRs:
                column_name = "MatchedTRFrom:" + column_name
            else:
                column_name = "OverlapsWith:" + column_name
            if column_name in path_to_column_func:
                column_name += "_"

            # generate IntervalTrees
            current_interval_tree = parse_bed_to_interval_tree(
                bed_file_path, name_field_is_repeat_unit=contains_TRs, verbose=args.verbose)

            # store the function for computing overlap with these IntervalTrees
            path_to_column_func[column_name] = get_overlapping_interval_generator(
                current_interval_tree, require_motif_match=contains_TRs)

    # parse EHdn tables
    locus_outlier_dfs = []
    case_control_dfs = []
    for locus_outlier_tsv in args.ehdn_locus_outlier_tsv:
        print(f"Processing rows from {locus_outlier_tsv}")
        try:
            df = pd.read_table(locus_outlier_tsv)
        except Exception as e:
            parser.error(f"Couldn't read locus outlier results from {locus_outlier_tsv}: {e}")

        # make sure the input table has the expected columns
        table_type = None
        counts_column = None
        if len(LOCUS_OUTLIER_COLUMNS - set(df.columns)) == 0:
            table_type = "outlier"
            counts_column = "high_case_counts"
            df_list = locus_outlier_dfs
        elif len(CASE_CONTROL_COLUMNS - set(df.columns)) == 0:
            table_type = "case_control"
            counts_column = "counts"
            df_list = case_control_dfs
        else:
            parser.error(f"{locus_outlier_tsv} has unexpected columns. Expected "
                         f"{LOCUS_OUTLIER_COLUMNS} or {CASE_CONTROL_COLUMNS}, but found {set(df.columns)}")

        # add basic columns
        df["Source"] = os.path.basename(locus_outlier_tsv)
        df["MotifSize"] = df["motif"].str.len()
        df["CanonicalMotif"] = df["motif"].apply(compute_canonical_motif)

        # compute overlap with known disease-associated loci and motifs
        df["MatchedKnownDiseaseAssociatedLocus"] = df.apply(known_disease_associated_loci_column_func, axis=1)
        df["MatchedKnownDiseaseAssociatedMotif"] = df["CanonicalMotif"].isin(known_disease_associated_motifs)

        # compute overlap with genes
        for genes_gtf in args.genes_gtf:
            label = os.path.basename(genes_gtf).split(".")[0].title()
            df[[f"{label}GeneRegion", f"{label}GeneName", f"{label}GeneId", f"{label}TranscriptId"]] = df.apply(
                lambda input_row: compute_genomic_region_of_row(input_row, verbose=args.verbose), axis=1, result_type='expand')

        # compute overlap with other bed files
        for column_name, column_func in path_to_column_func.items():
            df[column_name] = df.apply(column_func, axis=1)

        # split the counts column and output a row for each sample
        output_rows = []
        for _, row in df.iterrows():
            for sample_and_counts_columns in parse_sample_id_and_counts_column(row[counts_column]):
                output_row = row.to_dict()
                output_row.update(sample_and_counts_columns)
                del output_row[counts_column]
                if table_type == "outlier":
                    del output_row["counts"]
                output_rows.append(output_row)

        output_df = pd.DataFrame(output_rows)

        output_dir = args.output_dir or os.path.dirname(locus_outlier_tsv)
        output_path = os.path.join(output_dir, re.sub(".tsv", "", os.path.basename(locus_outlier_tsv))+".annotated.tsv")

        with open(output_path, "wt") as f:
            output_df.to_csv(f, sep="\t", index=False)

        df_list.append(output_df)

        print(f"Wrote {len(output_df):,d} rows to {output_path} ({len(output_df)/len(df):0.1f}x the number of input rows)")

    # output a combined table
    for label, df_list, sort_by in [
        ("locus_outlier", locus_outlier_dfs, "top_case_zscore"),
        ("case_control", case_control_dfs, "bonf_pvalue"),
    ]:

        if len(df_list) <= 1:
            continue
        combined_df = pd.concat(df_list)
        combined_df = combined_df.sort_values(sort_by, ascending=True)
        #combined_df = combined_df.sort_values(["Source", "MotifSize", "CanonicalMotif"], ascending=True)
        output_path = os.path.join(output_dir, f"combined.{len(df_list)}_{label}_tables.annotated.tsv")
        print("Combining", len(df_list), label, "tables")
        with open(output_path, "wt") as f:
            combined_df.to_csv(f, sep="\t", index=False)
        print(f"Wrote {len(combined_df):,d} rows to {output_path}")


if __name__ == "__main__":
    main()
