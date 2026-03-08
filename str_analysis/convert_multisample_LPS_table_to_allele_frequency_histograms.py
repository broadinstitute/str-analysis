"""
This script takes a wide-format TSV file with columns:

trid     (example: "10-100000859-100000887-A" or "10-100000859-100000887-A,10-100001413-100001429-T")
motif    (example: "A")
<sample1> <sample2> ...  (example: "3,3" meaning two alleles of size 3; "." for no-call)

and outputs a per-locus summary table with allele frequency histograms and statistics.
"""

import argparse
import collections
import gzip
import json
import os
import numpy as np
import tqdm

HEADER_FIELDS = [
    "LocusId",
    "Motif",
    "AlleleSizeHistogram",
    "BiallelicHistogram",
    "Min",
    "Mode",
    "Mean",
    "Stdev",
    "Median",
    "99thPercentile",
    "Max",
    "ShortAllele99thPercentile",
    "ShortAlleleMax",    
    "UniqueAlleleLengths",
    "NumCalledAlleles",
]

def _format_decimal(value):
    """Format a numeric value as an integer if whole, otherwise round to 2 decimal places."""
    if value == int(value):
        return int(value)
    return round(float(value), 2)


def compute_row(locus_id, motif, allele_sizes, alleles_by_sample_id):
    """Compute statistics for a group of allele sizes.

    Args:
        locus_id (str): the locus id (can be a comma-separated list of locus ids when it's a variation cluster)
        motif (str): the motif (in variation clusters, this will correspond to the motif at the end of one of the locus ids)
        allele_sizes (list): the allele sizes for all samples
        alleles_by_sample_id (dict): the alleles for the current key by sample id

    Returns:
        dict: a dictionary mapping HEADER_FIELDS keys to values, or None if allele_sizes is empty
    """

    if not allele_sizes:
        return None

    if "," in locus_id:
        found_locus_id = None
        for specific_locus_id in locus_id.split(","):
            if specific_locus_id.endswith(f"-{motif}"):
                found_locus_id = specific_locus_id
                break
        else:
            raise ValueError(f"Couldn't resolve locus id for motif {motif} in {locus_id}")

        locus_id = found_locus_id

    allele_sizes = sorted(allele_sizes)
    allele_counts = collections.Counter(allele_sizes)

    genotype_counts = collections.defaultdict(int)
    short_alleles = []
    for sample_id, allele_list in alleles_by_sample_id.items():
        if len(allele_list) == 1:
            short_alleles.append(allele_list[0])
            allele_list = allele_list * 2
        elif len(allele_list) == 2:
            short_alleles.append(min(allele_list))
        else:
            raise ValueError(f"Found {len(allele_list)} alleles for {sample_id} in {locus_id} {motif}")

        genotype_counts[tuple(sorted(allele_list))] += 1

    return {
        "LocusId": locus_id,
        "Motif": motif,
        "AlleleSizeHistogram": ",".join(
            f"{allele_size}x:{count}" for allele_size, count in sorted(allele_counts.items())
        ),
        "BiallelicHistogram": ",".join(
            f"{genotype[0]}/{genotype[1]}:{count}" for genotype, count in sorted(genotype_counts.items(), key=lambda x: (x[0][0], x[0][1]))
        ),
        "Min": int(min(allele_sizes)),
        "Mode": min(size for size, count in allele_counts.items() if count == allele_counts.most_common(1)[0][1]),
        "Mean": f"{np.mean(allele_sizes):.2f}",
        "Stdev": f"{np.std(allele_sizes):.2f}",
        "Median": _format_decimal(np.median(allele_sizes)),
        "99thPercentile": _format_decimal(np.percentile(allele_sizes, 99)),
        "Max": int(max(allele_sizes)),
        "ShortAllele99thPercentile": _format_decimal(np.percentile(short_alleles, 99)),
        "ShortAlleleMax": int(max(short_alleles)),
        "UniqueAlleleLengths": len(set(allele_sizes)),
        "NumCalledAlleles": len(allele_sizes),
    }


"""
Example row in 1kgp metadata TSV file:

SampleId                                                            HG00459
Sex                                                                    male
BiosampleId                                                      SAME125269
Subpopulation                                                           CHS
SubpopulationName                                      Southern Han Chinese
Population                                                              EAS
PopulationName                                          East Asian Ancestry
SubpopulationElasticId                                                  CHS
DataSource                1000 Genomes 30x on GRCh38,1000 Genomes phase ...

Population counts:
   1427 AFR
    535 AMR
    623 EAS
    672 EUR
      1 EUR,AFR
    661 SAS


 2343 female
   2653 male

Input table example row:

trid       1-19175-19184-TCC
motif                    TCC
HG00096                  3,3
HG00097                  3,3
HG00099                  3,3
...


"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-metadata-tsv", default="https://storage.googleapis.com/tandem-repeat-catalog/1kGP_metadata.tsv", help="Sample ancestry metadata TSV file (local path or URL)")
    parser.add_argument("--input-table", default="hprc-lps_2025-12-06/hprc-lps.txt.gz", help="Combine HPRC LPS dataset")
    parser.add_argument("--no-header", action="store_true", help="If set, assume the first row is data (not a header) and generate synthetic sample names (_s1, _s2, ...)")
    parser.add_argument("--population", choices=["AFR", "AMR", "EAS", "EUR", "SAS"], help="If specified, only process samples from this population")
    parser.add_argument("--sex", choices=["male", "female"], help="If specified, only process samples from this sex")
    parser.add_argument("--output-format", choices=["TSV", "JSON"], default="TSV", help="Output format (default: TSV)")
    parser.add_argument("-n", "--num-samples", type=int, default=None, help="Number of samples to process")
    parser.add_argument("-l", "--num-loci", type=int, default=None, help="Number of loci to process")
    args = parser.parse_args()

    if not os.path.isfile(args.input_table):
        parser.error(f"Input file {args.input_table} does not exist")

    fopen = gzip.open if args.input_table.endswith("gz") else open

    if args.no_header:
        # Read first line to determine number of columns, then generate synthetic header
        with fopen(args.input_table, "rt") as infile:
            first_line = next(infile).strip().split("\t")
        num_sample_columns = len(first_line) - 2
        header_fields = ["trid", "motif"] + [f"_s{i+1}" for i in range(num_sample_columns)]
        sample_ids_in_input_table = header_fields[2:]
    else:
        with fopen(args.input_table, "rt") as infile:
            header_fields = next(infile).strip().split("\t")
            if header_fields[0:2] != ["trid", "motif"]:
                parser.error(f"Expected header to start with 'trid', 'motif' columns in {args.input_table}. Found: {header_fields}")
        sample_ids_in_input_table = header_fields[2:]
        print(f"Parsed {len(sample_ids_in_input_table):,d} sample ids from {args.input_table}")

    if args.no_header:
        if args.population or args.sex:
            parser.error("--population and --sex filters require a header row with real sample ids and a --sample-metadata-tsv")
        sample_ids_to_include_list = sample_ids_in_input_table
        if args.num_samples is not None and len(sample_ids_to_include_list) > args.num_samples:
            sample_ids_to_include_list = sample_ids_to_include_list[:args.num_samples]
        sample_ids_to_include = set(sample_ids_to_include_list)
    else:
        import pandas as pd
        df_metadata = pd.read_table(args.sample_metadata_tsv)
        all_sample_ids = set(df_metadata.SampleId)
        unexpected_sample_ids = set(sample_ids_in_input_table) - all_sample_ids
        if unexpected_sample_ids:
            parser.error(f"{len(unexpected_sample_ids):,d} sample id(s) not found in {args.sample_metadata_tsv}: {', '.join(sorted(unexpected_sample_ids))}")

        df_metadata = df_metadata[df_metadata.SampleId.isin(set(sample_ids_in_input_table))]
        if args.population:
            df_metadata = df_metadata[df_metadata.Population == args.population]
            print(f"Kept {len(df_metadata):,d} samples from population {args.population}")
        if args.sex:
            df_metadata = df_metadata[df_metadata.Sex == args.sex]
            print(f"Kept {len(df_metadata):,d} samples from sex {args.sex}")

        valid_ids = set(df_metadata.SampleId)
        sample_ids_to_include_list = [s for s in sample_ids_in_input_table if s in valid_ids]
        if args.num_samples is not None and len(sample_ids_to_include_list) > args.num_samples:
            sample_ids_to_include_list = sample_ids_to_include_list[:args.num_samples]
        sample_ids_to_include = set(sample_ids_to_include_list)

        print(f"Included sample ids: {', '.join(sample_ids_to_include_list)}")

    output_dir = os.path.dirname(args.input_table)
    output_path = os.path.join(output_dir, os.path.basename(args.input_table).replace(".tsv", "").replace(".txt", "").replace(".gz", ""))
    output_path += ".per_locus_and_motif"
    if args.population:
        output_path += f".only_{args.population}"
    if args.sex:
        output_path += f".only_{args.sex}"

    output_path += f".{len(sample_ids_to_include_list)}_samples"
    if args.output_format == "JSON":
        output_path += ".json.gz"
    else:
        output_path += ".tsv.gz"
    print(f"Writing data from {len(sample_ids_to_include_list):,d} samples to {output_path}")
    with fopen(args.input_table, "rt") as infile, gzip.open(output_path, "wt") as outfile:
        if not args.no_header:
            next(infile)  # skip header

        if args.output_format == "TSV":
            outfile.write("\t".join(HEADER_FIELDS) + "\n")
        else:
            outfile.write("[\n")
            json_first_row = True

        line_count = 0
        for line_number, line in tqdm.tqdm(enumerate(infile), unit=" lines", unit_scale=True):
            fields = line.strip().split("\t")
            if len(fields) < 2:
                continue  # skip blank or malformed lines
            if fields[0] == "trid" and fields[1] == "motif":
                continue  # skip internal header rows in concatenated files

            if args.num_loci is not None and line_count >= args.num_loci:
                break
            line_count += 1

            if len(fields) != len(header_fields):
                raise ValueError(f"Line {line_number + 1} has {len(fields)} fields, expected {len(header_fields)}: {line.strip()[:200]}")

            locus_id = fields[0]
            motif = fields[1]

            alleles = []
            alleles_by_sample_id = collections.defaultdict(list)
            for sample_id, allele_sizes in zip(header_fields[2:], fields[2:]):
                if sample_id not in sample_ids_to_include:
                    continue
                if allele_sizes == ".":
                    continue
                for allele_size in allele_sizes.split(","):
                    try:
                        allele_size = int(allele_size)
                    except ValueError:
                        raise ValueError(f"Expected integer allele size, got {allele_size} in line #{line_number + 1}: {line}")
                    alleles.append(allele_size)
                    alleles_by_sample_id[sample_id].append(allele_size)

            row = compute_row(locus_id, motif, alleles, alleles_by_sample_id)
            if row:
                if args.output_format == "JSON":
                    if not json_first_row:
                        outfile.write(",\n")
                    json_first_row = False
                    outfile.write("  " + json.dumps(row, indent=2).replace("\n", "\n  "))
                else:
                    outfile.write("\t".join(str(row[field]) for field in HEADER_FIELDS) + "\n")

        if args.output_format == "JSON":
            outfile.write("\n]\n")

    print(f"Wrote {line_count:9,d} lines to {output_path}")

if __name__ == "__main__":
    main()
