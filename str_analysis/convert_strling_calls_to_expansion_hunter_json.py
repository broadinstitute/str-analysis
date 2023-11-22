"""This script converts a STRling output .txt to the .json format ExpansionHunter uses to output results.
This makes it easier to pass STRling results to downstream scripts.
"""

import argparse
import collections
import simplejson as json
import os
import re
import pandas as pd


"""
ExpansionHunter output format:

  "LocusResults": {
        "chr12-57610122-57610131-GCA": {
          "AlleleCount": 2,
          "Coverage": 50.469442942130875,
          "FragmentLength": 433,
          "LocusId": "chr12-57610122-57610131-GCA",
          "ReadLength": 151,
          "Variants": {
            "chr12-57610122-57610131-GCA": {
              "CountsOfFlankingReads": "(1, 1), (2, 4)",
              "CountsOfInrepeatReads": "()",
              "CountsOfSpanningReads": "(2, 1), (3, 48), (6, 1)",
              "Genotype": "3/3",
              "GenotypeConfidenceInterval": "3-3/3-3",
              "ReferenceRegion": "chr12:57610122-57610131",
              "RepeatUnit": "GCA",
              "VariantId": "chr12-57610122-57610131-GCA",
              "VariantType": "Repeat"
            }
          }
        },

  "SampleParameters": {
        "SampleId": "NA19239",
        "Sex": "Female"
  }

"""

"""
STRling txt calls format:

$1                    #chrom : chr11
$2                      left : 119648275
$3                     right : 119648276
$4                repeatunit : GGT
$5               allele1_est : 2.00
$6               allele2_est : 22.03
$7            anchored_reads : 3
$8            spanning_reads : 14
$9            spanning_pairs : 8
$10  expected_spanning_pairs : 68.42
$11      spanning_pairs_pctl : 0.00
$12               left_clips : 0
$13              right_clips : 0
$14           unplaced_pairs : 0
$15                    depth : 34.0
$16           sum_str_counts : 126
"""


def main():
    p = argparse.ArgumentParser()
    p.add_argument("strling_txt_path", nargs="+", help="STRling calls .txt path(s)")
    args = p.parse_args()

    for strling_txt_path in args.strling_txt_path:
        print(f"Processing {strling_txt_path}")
        locus_results = process_strling_txt(strling_txt_path)

        output_json_path = re.sub(".txt$", "", strling_txt_path) + ".json"
        print(f"Writing results for {len(locus_results['LocusResults'])} loci to {output_json_path}")
        with open(output_json_path, "wt") as f:
            json.dump(locus_results, f, indent=4)


def compute_CI(allele_size, repeat_unit_length):
    """Compute a pseudo "confidence interval" based on the criteria using for mendelian violation tests in the
    STRling pre-print - where a genotype was considered not to be a mendelian violation if it was within
    25% or 10bp of the expected genotype.

    Args:
        allele_size (int): The number of repeats.
        repeat_unit_length (str): Length of the repeat unit in base pairs.

    Returns:
         str: "{lower}-{upper}" representing the lower/upper bound in numbers of repeats.
    """
    max_repeats_distance = 10.0  # based on this utility script in the STRling repo: https://github.com/laurelhiatt/strling-MV/blob/main/denovo.py#L170
    fraction = 0.25
    lower = int(round(min(allele_size - max_repeats_distance, (1-fraction) * allele_size)))
    upper = int(round(max(allele_size + max_repeats_distance, (1+fraction) * allele_size)))

    lower = max(0, lower)

    return f"{lower}-{upper}"


def process_strling_txt(strling_txt_path):
    sample_id = re.sub("(-genotype)?.txt$", "", os.path.basename(strling_txt_path))

    locus_results = {
        "LocusResults": {},
        "SampleParameters": {
            "SampleId": sample_id,
            "Sex": None,
        },
    }

    df = pd.read_table(strling_txt_path)
    df.rename(columns={"#chrom": "chrom"}, inplace=True)

    counter = collections.defaultdict(int)
    for _, row in df.iterrows():
        counter["total loci"] += 1
        start_1based = int(row.left) + 1
        end_1based = int(row.right)

        if end_1based - start_1based == 0:
            counter["0-sized locus"] += 1

        locus_id = f"{row.chrom}-{start_1based - 1}-{end_1based}-{row.repeatunit}"
        if pd.isna(row.allele1_est) and pd.isna(row.allele2_est):
            counter["allele1 & 2 are NaN"] += 1
            continue

        elif pd.isna(row.allele1_est):
            allele2 = max(int(round(row.allele2_est)), 0)
            counter["allele1 is NaN"] += 1
            genotype = f"{allele2}"
            genotype_ci = f"{compute_CI(allele2, len(row.repeatunit))}"
        elif pd.isna(row.allele2_est):
            allele1 = max(int(round(row.allele1_est)), 0)
            genotype = f"{allele1}"
            genotype_ci = f"{compute_CI(allele1, len(row.repeatunit))}"
            counter["allele2 is NaN"] += 1
        else:
            allele1 = max(int(round(row.allele1_est)), 0)
            allele2 = max(int(round(row.allele2_est)), 0)
            genotype = f"{allele1}/{allele2}"
            genotype_ci = f"{compute_CI(allele1, len(row.repeatunit))}/{compute_CI(allele2, len(row.repeatunit))}"

        locus_results["LocusResults"][locus_id] = {
            "AlleleCount": 2,
            "LocusId": locus_id,
            "Coverage": row.depth,
            "ReadLength": None,
            "FragmentLength": None,
            "Variants": {
                locus_id: {
                    "Genotype": genotype,
                    "GenotypeConfidenceInterval": genotype_ci,
                    "ReferenceRegion": f"{row.chrom}:{start_1based - 1}-{end_1based}",
                    "RepeatUnit": row.repeatunit,
                    "VariantId": locus_id,
                    "VariantType": "Repeat",
                    #"CountsOfFlankingReads": "()",
                    #"CountsOfInrepeatReads": "()",
                    #"CountsOfSpanningReads": "()",
                }
            }
        }

    print(f"Finished processing {strling_txt_path}")
    for key, value in sorted(counter.items(), key=lambda x: -x[1]):
        print(f"{value:10d}: {key:20s} ({100*value/counter['total loci']:0.1f}%)")

    return locus_results


if __name__ == "__main__":
    main()
