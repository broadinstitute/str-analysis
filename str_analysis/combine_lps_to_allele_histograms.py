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
    if allele_histogram[allele] <= n_outlier_sample_ids:
        if sample_id not in allele_sample_ids[allele]:
            allele_sample_ids[allele].append(sample_id)
    elif allele in allele_sample_ids:
        del allele_sample_ids[allele]


def convert_counts_to_histogram_string(allele_counts):
    return ",".join(f"{allele_size}x:{count}" for allele_size, count in sorted(allele_counts.items()))

def convert_sample_ids_to_string(allele_sample_ids, n_outlier_sample_ids=10):
    output_value = []
    output_sample_id_counter = 0
    for allele, sample_list in sorted(allele_sample_ids.items(), key=lambda x: -x[0]):
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
    parser.add_argument("--output-path", help="Output TSV file")
    parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    parser.add_argument("input_tables", nargs="+", help="Input TSV files with LPS scores")

    args = parser.parse_args()
    for input_table in args.input_tables:
        if not os.path.isfile(input_table):
            parser.error(f"Input table {input_table} does not exist")

    if not args.output_path:
        if len(args.input_tables) == 1:
            args.output_path = re.sub("(.lps)?(.txt|.tsv)(.gz)?$", "", args.input_tables[0]) + ".lps_allele_histograms.tsv.gz"
        else:
            args.output_path = f"combined.{len(args.input_tables)}_input_tables.lps_allele_histograms.tsv.gz"
    else:
        args.output_path = args.output_path


    trid_to_allele_histogram = collections.defaultdict(lambda: collections.defaultdict(int))
    trid_to_allele_sample_ids = collections.defaultdict(lambda: collections.defaultdict(list))

    trid_to_short_allele_histogram = collections.defaultdict(lambda: collections.defaultdict(int))
    trid_to_short_allele_sample_ids = collections.defaultdict(lambda: collections.defaultdict(list))

    trid_to_hemizygous_allele_histogram = collections.defaultdict(lambda: collections.defaultdict(int))
    trid_to_hemizygous_allele_sample_ids = collections.defaultdict(lambda: collections.defaultdict(list))

    input_table_iterator = args.input_tables
    if args.show_progress_bar:
        input_table_iterator = tqdm.tqdm(input_table_iterator, unit=" tables")

    for input_table in input_table_iterator:
        fopen = gzip.open if input_table.endswith("gz") else open
        with fopen(input_table, "rt") as f:
            header_fields = next(f).rstrip().split("\t")
            if len(header_fields) < 3 or header_fields[0] != "trid" or header_fields[1] != "motif":
                parser.error(f"Unexpected header in {input_table}: {header_fields}. Expecting trid, motif, <sample id>")

            sample_id = header_fields[2]
            print(f"Processing {input_table} with sample {sample_id}")

            line_count = 0
            for line in f:
                line_count += 1
                fields = line.rstrip().split("\t")
                trid = fields[0]
                motif = fields[1]
                lps = fields[2]
                lps_values = [int(v) for v in lps.split(",")]

                short_allele = min(lps_values)
                long_allele = max(lps_values)

                # update trid_to_allele_histogram
                update_histograms(trid, motif, short_allele, trid_to_allele_histogram, trid_to_allele_sample_ids, sample_id, args.n_outlier_sample_ids)
                if len(lps_values) > 1:
                    update_histograms(trid, motif, long_allele, trid_to_allele_histogram, trid_to_allele_sample_ids, sample_id, args.n_outlier_sample_ids)

                # update trid_to_short_allele_histogram
                update_histograms(trid, motif, short_allele, trid_to_short_allele_histogram, trid_to_short_allele_sample_ids, sample_id, args.n_outlier_sample_ids)

                #print(trid_to_allele_sample_ids[(trid, motif)])

                # update trid_to_hemizygous_allele_histogram
                if len(lps_values) == 1:
                    update_histograms(trid, motif, short_allele, trid_to_hemizygous_allele_histogram, trid_to_hemizygous_allele_sample_ids, sample_id, args.n_outlier_sample_ids)

            print(f"Processed {line_count:,d} rows from {input_table}")

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
            output_row = {
                "LocusId": trid,
                "Motif": motif,
                "AllAlleleHistogram": convert_counts_to_histogram_string(trid_to_allele_histogram[(trid, motif)]),
                "ShortAlleleHistogram": convert_counts_to_histogram_string(trid_to_short_allele_histogram[(trid, motif)]),
                "HemizygousAlleleHistogram": convert_counts_to_histogram_string(trid_to_hemizygous_allele_histogram[(trid, motif)]),

                "OutlierSampleIds_AllAlleles": convert_sample_ids_to_string(trid_to_allele_sample_ids[(trid, motif)], n_outlier_sample_ids=args.n_outlier_sample_ids),
                "OutlierSampleIds_ShortAlleles": convert_sample_ids_to_string(trid_to_short_allele_sample_ids[(trid, motif)], n_outlier_sample_ids=args.n_outlier_sample_ids),
                "OutlierSampleIds_HemizygousAlleles": convert_sample_ids_to_string(trid_to_hemizygous_allele_sample_ids[(trid, motif)], n_outlier_sample_ids=args.n_outlier_sample_ids),
            }
            out_f.write("\t".join(str(output_row[col]) for col in header) + "\n")
            output_row_counter += 1

    print(f"Wrote {output_row_counter:,d} rows to {args.output_path}")


if __name__ == "__main__":
    main()