"""This tool takes one or more tables of LPS scores and combines them into per-locus allele histograms.
The input table(s) should have the following 3 columns (where "sample1_xx" and "sample1_xy" are arbitrary sample ids):

trid	motif	sample1_xx
X-95117-95209-GTCACGGCCCGAGACTCCCTCTTCCT	GTCACGGCCCGAGACTCCCTCTTCCT	2,2
X-263540-263579-TTTA	TTTA	10,11
X-264744-264774-AAAG	AAAG	5,5

or (for male samples):

trid	motif	sample1_xy
X-95117-95209-GTCACGGCCCGAGACTCCCTCTTCCT	GTCACGGCCCGAGACTCCCTCTTCCT	3
X-263540-263579-TTTA	TTTA	9
X-264744-264774-AAAG	AAAG	4

The output table will have the following columns:

LocusId	Motif	AllAlleleHistogram	ShortAlleleHistogram	HemizygousAlleleHistogram	OutlierSampleIds_AllAlleles	OutlierSampleIds_ShortAlleles	OutlierSampleIds_HemizygousAlleles
X-95117-95209-GTCACGGCCCGAGACTCCCTCTTCCT	GTCACGGCCCGAGACTCCCTCTTCCT	2x:2	2x:1		2x:sample1_xx	2x:sample1_xx
X-263540-263579-TTTA	TTTA	10x:1,11x:1	10x:1		11x:sample1_xx,10x:sample1_xx	10x:sample1_xx
X-264744-264774-AAAG	AAAG	5x:2	5x:1		5x:sample1_xx	5x:sample1_xx

or 

LocusId	Motif	AllAlleleHistogram	ShortAlleleHistogram	HemizygousAlleleHistogram	OutlierSampleIds_AllAlleles	OutlierSampleIds_ShortAlleles	OutlierSampleIds_HemizygousAlleles
X-95117-95209-GTCACGGCCCGAGACTCCCTCTTCCT	GTCACGGCCCGAGACTCCCTCTTCCT	3x:1	3x:1	3x:1	3x:sample1_xy	3x:sample1_xy	3x:sample1_xy
X-263540-263579-TTTA	TTTA	9x:1	9x:1	9x:1	9x:sample1_xy	9x:sample1_xy	9x:sample1_xy
X-264744-264774-AAAG	AAAG	4x:1	4x:1	4x:1	4x:sample1_xy	4x:sample1_xy	4x:sample1_xy
"""



import argparse
import collections
import gzip
import os
import re
import tqdm



def update_histograms(trid, motif, allele, trid_to_allele_histogram, trid_to_allele_sample_ids, sample_id, n_outlier_sample_ids):
    allele_histogram = trid_to_allele_histogram[(trid, motif)]
    allele_sample_ids = trid_to_allele_sample_ids[(trid, motif)]
    allele_histogram[allele] += 1
    if allele_histogram[allele] <= 2 * n_outlier_sample_ids:  # * 2 since each homozygous sample is counted twice
        if sample_id not in allele_sample_ids[allele]:
            allele_sample_ids[allele].append(sample_id)


def convert_counts_to_histogram_string(allele_counts):
    return ",".join(f"{allele_size}x:{count}" for allele_size, count in sorted(allele_counts.items()))


def get_lps_filename_prefix(lps_table_path):
    return re.sub("(.lps)?(.txt|.tsv)(.gz)?$", "", os.path.basename(lps_table_path))


def parse_lps_table(
    input_table_path,
    trids_to_include,
    n_outlier_sample_ids,
    trid_to_allele_histogram,
    trid_to_allele_sample_ids,
    trid_to_short_allele_histogram,
    trid_to_short_allele_sample_ids,
    trid_to_hemizygous_allele_histogram,
    trid_to_hemizygous_allele_sample_ids,
    use_sample_id_from_header=False,
):
    """Parse an LPS table and update the histogram dictionaries.

    Args:
        input_table_path: Path to the input LPS table file.
        trids_to_include: Set of TRIDs to include, or None to include all.
        n_outlier_sample_ids: Number of outlier sample IDs to track per allele.
        trid_to_allele_histogram: Dict to update with all allele counts.
        trid_to_allele_sample_ids: Dict to update with all allele sample IDs.
        trid_to_short_allele_histogram: Dict to update with short allele counts.
        trid_to_short_allele_sample_ids: Dict to update with short allele sample IDs.
        trid_to_hemizygous_allele_histogram: Dict to update with hemizygous allele counts.
        trid_to_hemizygous_allele_sample_ids: Dict to update with hemizygous allele sample IDs.
        use_sample_id_from_header: If True, use sample ID from header; otherwise derive from filename.

    Returns:
        Tuple of (included_line_count, total_line_count, sample_ids) where sample_ids is a set
        containing the single sample ID from this table.

    Raises:
        ValueError: If the header format is unexpected.
    """
    fopen = gzip.open if input_table_path.endswith("gz") else open
    with fopen(input_table_path, "rt") as f:
        header_fields = next(f).rstrip().split("\t")
        if len(header_fields) < 3 or header_fields[0] != "trid" or header_fields[1] != "motif":
            raise ValueError(f"Unexpected header in {input_table_path}: {header_fields}. Expecting trid, motif, <sample id>")

        if use_sample_id_from_header:
            sample_id = header_fields[2]
        else:
            sample_id = get_lps_filename_prefix(input_table_path)

        sample_ids = {sample_id}
        print(f"Processing {input_table_path} with sample {sample_id}")

        total_line_count = 0
        included_line_count = 0
        for line in f:
            total_line_count += 1
            fields = line.rstrip('\n\r').split("\t")
            trid = fields[0]
            if trids_to_include is not None and trid not in trids_to_include:
                continue

            included_line_count += 1
            motif = fields[1]
            lps = fields[2]
            lps_values = [int(v) for v in lps.split(",")]

            is_hemizygous = len(lps_values) == 1
            short_allele = min(lps_values)
            long_allele = max(lps_values)

            # update trid_to_allele_histogram
            update_histograms(trid, motif, short_allele, trid_to_allele_histogram, trid_to_allele_sample_ids, sample_id, n_outlier_sample_ids)
            if not is_hemizygous:
                update_histograms(trid, motif, long_allele, trid_to_allele_histogram, trid_to_allele_sample_ids, sample_id, n_outlier_sample_ids)

            # update trid_to_short_allele_histogram
            update_histograms(trid, motif, short_allele, trid_to_short_allele_histogram, trid_to_short_allele_sample_ids, sample_id, n_outlier_sample_ids)

            # update trid_to_hemizygous_allele_histogram
            if is_hemizygous:
                update_histograms(trid, motif, short_allele, trid_to_hemizygous_allele_histogram, trid_to_hemizygous_allele_sample_ids, sample_id, n_outlier_sample_ids)

    return included_line_count, total_line_count, sample_ids


def parse_histogram_string(histogram_string):
    """Parse a histogram string like '10x:1,11x:2' into a dict {10: 1, 11: 2}."""
    if not histogram_string:
        return {}
    result = {}
    for entry in histogram_string.split(","):
        allele_str, count_str = entry.split("x:")
        result[int(allele_str)] = int(count_str)
    return result


def parse_sample_ids_string(sample_ids_string):
    """Parse a sample IDs string like '10x:sample1,11x:sample2' into a dict {10: ['sample1'], 11: ['sample2']}."""
    if not sample_ids_string:
        return {}
    result = {}
    for entry in sample_ids_string.split(","):
        allele_str, sample_id = entry.split("x:", 1)
        allele = int(allele_str)
        if allele not in result:
            result[allele] = []
        result[allele].append(sample_id)
    return result


def parse_combined_lps_allele_histograms_table(
    input_table_path,
    trids_to_include,
    n_outlier_sample_ids,
    trid_to_allele_histogram,
    trid_to_allele_sample_ids,
    trid_to_short_allele_histogram,
    trid_to_short_allele_sample_ids,
    trid_to_hemizygous_allele_histogram,
    trid_to_hemizygous_allele_sample_ids,
):
    """Parse a combined LPS allele histograms table and merge into the histogram dictionaries.

    Args:
        input_table_path: Path to the input allele histograms table file.
        trids_to_include: Set of TRIDs to include, or None to include all.
        n_outlier_sample_ids: Number of outlier sample IDs to track per allele.
        trid_to_allele_histogram: Dict to update with all allele counts.
        trid_to_allele_sample_ids: Dict to update with all allele sample IDs.
        trid_to_short_allele_histogram: Dict to update with short allele counts.
        trid_to_short_allele_sample_ids: Dict to update with short allele sample IDs.
        trid_to_hemizygous_allele_histogram: Dict to update with hemizygous allele counts.
        trid_to_hemizygous_allele_sample_ids: Dict to update with hemizygous allele sample IDs.

    Returns:
        Tuple of (included_line_count, total_line_count, sample_ids) where sample_ids is a set
        containing all sample IDs found in the outlier columns.

    Raises:
        ValueError: If the header format is unexpected.
    """
    fopen = gzip.open if input_table_path.endswith("gz") else open
    with fopen(input_table_path, "rt") as f:
        header_fields = next(f).rstrip().split("\t")
        expected_header = ["LocusId", "Motif", "AllAlleleHistogram", "ShortAlleleHistogram",
                          "HemizygousAlleleHistogram", "OutlierSampleIds_AllAlleles",
                          "OutlierSampleIds_ShortAlleles", "OutlierSampleIds_HemizygousAlleles"]
        if header_fields != expected_header:
            raise ValueError(f"Unexpected header in {input_table_path}: {header_fields}. "
                           f"Expecting {expected_header}")

        print(f"Processing combined histograms table {input_table_path}")

        total_line_count = 0
        included_line_count = 0
        sample_ids = set()
        for line in f:
            total_line_count += 1
            fields = line.rstrip('\n\r').split("\t")
            trid = fields[0]
            if trids_to_include is not None and trid not in trids_to_include:
                continue

            included_line_count += 1
            motif = fields[1]

            # Parse and merge all allele histogram
            all_allele_hist = parse_histogram_string(fields[2])
            for allele, count in all_allele_hist.items():
                trid_to_allele_histogram[(trid, motif)][allele] += count

            # Parse and merge short allele histogram
            short_allele_hist = parse_histogram_string(fields[3])
            for allele, count in short_allele_hist.items():
                trid_to_short_allele_histogram[(trid, motif)][allele] += count

            # Parse and merge hemizygous allele histogram
            hemizygous_allele_hist = parse_histogram_string(fields[4])
            for allele, count in hemizygous_allele_hist.items():
                trid_to_hemizygous_allele_histogram[(trid, motif)][allele] += count

            # Parse and merge sample IDs (respecting n_outlier_sample_ids limit)
            all_sample_ids_dict = parse_sample_ids_string(fields[5] if len(fields) > 5 else "")
            for allele, sample_list in all_sample_ids_dict.items():
                sample_ids.update(sample_list)
                existing = trid_to_allele_sample_ids[(trid, motif)][allele]
                for sample_id in sample_list:
                    if len(existing) < 2 * n_outlier_sample_ids and sample_id not in existing:
                        existing.append(sample_id)

            short_sample_ids_dict = parse_sample_ids_string(fields[6] if len(fields) > 6 else "")
            for allele, sample_list in short_sample_ids_dict.items():
                sample_ids.update(sample_list)
                existing = trid_to_short_allele_sample_ids[(trid, motif)][allele]
                for sample_id in sample_list:
                    if len(existing) < 2 * n_outlier_sample_ids and sample_id not in existing:
                        existing.append(sample_id)

            hemizygous_sample_ids_dict = parse_sample_ids_string(fields[7] if len(fields) > 7 else "")
            for allele, sample_list in hemizygous_sample_ids_dict.items():
                sample_ids.update(sample_list)
                existing = trid_to_hemizygous_allele_sample_ids[(trid, motif)][allele]
                for sample_id in sample_list:
                    if len(existing) < 2 * n_outlier_sample_ids and sample_id not in existing:
                        existing.append(sample_id)

    return included_line_count, total_line_count, sample_ids


def detect_input_table_type(input_table_path):
    """Detect whether the input table is an LPS table or a combined histograms table.

    Args:
        input_table_path: Path to the input table file.

    Returns:
        'lps' if it's an LPS table, 'histograms' if it's a combined histograms table.

    Raises:
        ValueError: If the header format is not recognized.
    """
    fopen = gzip.open if input_table_path.endswith("gz") else open
    with fopen(input_table_path, "rt") as f:
        header_fields = next(f).rstrip().split("\t")

    if len(header_fields) >= 3 and header_fields[0] == "trid" and header_fields[1] == "motif":
        return "lps"
    elif header_fields[0] == "LocusId" and header_fields[1] == "Motif":
        return "histograms"
    else:
        raise ValueError(f"Unrecognized header format in {input_table_path}: {header_fields}")


def convert_sample_ids_to_string(allele_sample_ids, n_outlier_sample_ids=10):
    output_value = []
    output_sample_id_counter = 0
    for allele, sample_list in sorted(allele_sample_ids.items(), key=lambda x: -x[0]):
        if len(sample_list) >= n_outlier_sample_ids:
            break

        for sample_id in sorted(sample_list):
            output_value.append(f"{allele}x:{sample_id}")
            output_sample_id_counter += 1

        if output_sample_id_counter > n_outlier_sample_ids:
            break

    return ",".join(output_value)

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--n-outlier-sample-ids", type=int, default=10,
                        help="Number of samples with the largest alleles to record as outliers")
    parser.add_argument("-o", "--output-path", help="Output TSV file")
    parser.add_argument("--use-sample-id-from-header", action="store_true", help="Use the sample ID from "
                        "the header of each input table instead of deriving it from the filename prefix")
    parser.add_argument("--trid-list", help="Optional path of file that lists TRIDs (one per line) to "
                        "include in the output. If not specified, all TRIDs in the input tables will be included.")
    parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    parser.add_argument("--skip-loci-without-outliers", action="store_true", help="Skip loci with outliers")
    parser.add_argument("--allow-overlapping-samples", action="store_true",
                        help="Allow sample IDs to appear in multiple input tables (will merge their histograms)")
    parser.add_argument("input_tables", nargs="+", help="Input TSV files with LPS scores")

    args = parser.parse_args()

    for input_table in args.input_tables:
        if not os.path.isfile(input_table):
            parser.error(f"Input table {input_table} does not exist")



    if not args.output_path:
        if len(args.input_tables) == 1:
            args.output_path = os.path.join(
                os.path.dirname(args.input_tables[0]),
                f"{get_lps_filename_prefix(args.input_tables[0])}.lps_allele_histograms.tsv.gz",
            )
        else:
            args.output_path = f"combined.{len(args.input_tables)}_input_tables.lps_allele_histograms.tsv.gz"
    if not args.trid_list:
        trids_to_include = None
    else:
        if not os.path.isfile(args.trid_list):
            parser.error(f"TRID list file {args.trid_list} does not exist")

        fopen = gzip.open if args.trid_list.endswith("gz") else open
        with fopen(args.trid_list, "rt") as f:
            trids_to_include = set(line.strip() for line in f if line.strip())
        print(f"Will only include the {len(trids_to_include):,d} TRIDs from {args.trid_list}")

    trid_to_allele_histogram = collections.defaultdict(lambda: collections.defaultdict(int))
    trid_to_allele_sample_ids = collections.defaultdict(lambda: collections.defaultdict(list))

    trid_to_short_allele_histogram = collections.defaultdict(lambda: collections.defaultdict(int))
    trid_to_short_allele_sample_ids = collections.defaultdict(lambda: collections.defaultdict(list))

    trid_to_hemizygous_allele_histogram = collections.defaultdict(lambda: collections.defaultdict(int))
    trid_to_hemizygous_allele_sample_ids = collections.defaultdict(lambda: collections.defaultdict(list))

    input_table_iterator = args.input_tables
    if args.show_progress_bar:
        input_table_iterator = tqdm.tqdm(input_table_iterator, unit=" tables")

    all_seen_sample_ids = set()
    for input_table in input_table_iterator:
        table_type = detect_input_table_type(input_table)
        if table_type == "lps":
            included_line_count, total_line_count, table_sample_ids = parse_lps_table(
                input_table,
                trids_to_include,
                args.n_outlier_sample_ids,
                trid_to_allele_histogram,
                trid_to_allele_sample_ids,
                trid_to_short_allele_histogram,
                trid_to_short_allele_sample_ids,
                trid_to_hemizygous_allele_histogram,
                trid_to_hemizygous_allele_sample_ids,
                use_sample_id_from_header=args.use_sample_id_from_header,
            )
        else:
            included_line_count, total_line_count, table_sample_ids = parse_combined_lps_allele_histograms_table(
                input_table,
                trids_to_include,
                args.n_outlier_sample_ids,
                trid_to_allele_histogram,
                trid_to_allele_sample_ids,
                trid_to_short_allele_histogram,
                trid_to_short_allele_sample_ids,
                trid_to_hemizygous_allele_histogram,
                trid_to_hemizygous_allele_sample_ids,
            )

        # Check for overlapping sample IDs
        overlapping_sample_ids = all_seen_sample_ids & table_sample_ids
        if overlapping_sample_ids and not args.allow_overlapping_samples:
            raise ValueError(f"Input table {input_table} contains sample IDs that were already seen in "
                           f"previous input tables: {sorted(overlapping_sample_ids)}. "
                           f"Use --allow-overlapping-samples to merge histograms from overlapping samples.")
        all_seen_sample_ids.update(table_sample_ids)

        print(f"Processed {included_line_count:,d} out of {total_line_count:,d} rows from {input_table}")


    print(f"Writing output to {args.output_path}")
    header = [
        "LocusId",
        "Motif",
        "AllAlleleHistogram",
        "ShortAlleleHistogram",
        "HemizygousAlleleHistogram",

        "OutlierSampleIds_AllAlleles",
        "OutlierSampleIds_ShortAlleles",
        "OutlierSampleIds_HemizygousAlleles",
    ]

    fopen = gzip.open if args.output_path.endswith("gz") else open
    with fopen(args.output_path, "wt") as out_f:
        out_f.write("\t".join(header) + "\n")
        output_row_counter = 0
        for (trid, motif) in sorted(trid_to_short_allele_histogram.keys()):
            all_allele_outliers = convert_sample_ids_to_string(
                trid_to_allele_sample_ids[(trid, motif)], n_outlier_sample_ids=args.n_outlier_sample_ids)
            short_allele_outliers = convert_sample_ids_to_string(
                trid_to_short_allele_sample_ids[(trid, motif)], n_outlier_sample_ids=args.n_outlier_sample_ids)
            hemizygous_allele_outliers = convert_sample_ids_to_string(
                trid_to_hemizygous_allele_sample_ids[(trid, motif)], n_outlier_sample_ids=args.n_outlier_sample_ids)

            if args.skip_loci_without_outliers and not all_allele_outliers and not short_allele_outliers and not hemizygous_allele_outliers:
                continue

            output_row = {
                "LocusId": trid,
                "Motif": motif,
                "AllAlleleHistogram": convert_counts_to_histogram_string(trid_to_allele_histogram[(trid, motif)]),
                "ShortAlleleHistogram": convert_counts_to_histogram_string(trid_to_short_allele_histogram[(trid, motif)]),
                "HemizygousAlleleHistogram": convert_counts_to_histogram_string(trid_to_hemizygous_allele_histogram[(trid, motif)]),

                "OutlierSampleIds_AllAlleles": all_allele_outliers,
                "OutlierSampleIds_ShortAlleles": short_allele_outliers,
                "OutlierSampleIds_HemizygousAlleles": hemizygous_allele_outliers,
            }
            out_f.write("\t".join(str(output_row[col]) for col in header) + "\n")
            output_row_counter += 1

    print(f"Wrote {output_row_counter:,d} rows to {args.output_path}")


if __name__ == "__main__":
    main()