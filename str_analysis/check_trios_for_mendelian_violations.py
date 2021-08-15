#!/usr/bin/env python3

import argparse
import collections
from intervaltree import Interval
import logging
import pandas as pd
import pathlib
import re

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")


def parse_args(args_list=None):
    p = argparse.ArgumentParser()
    p.add_argument(
        "-f",
        "--fam-file",
        help=".fam file that describes parent-child relationships of individual IDs in "
        "the combined_str_calls_tsv",
        type=pathlib.Path,
        required=True,
    )
    p.add_argument(
        "combined_str_calls_tsv",
        help=".tsv table that conatains str calls for all individuals, created by the "
        "combine_expansion_hunter_json_results_to_tsv script",
        type=pathlib.Path,
    )
    args = p.parse_args(args=args_list)

    for path in [args.combined_str_calls_tsv, args.fam_file]:
        if not path.is_file():
            p.error(f"{path} not found")

    return args


def check_for_duplicate_keys(df, file_path):
    num_duplicate_keys = sum(df.index.duplicated())
    if num_duplicate_keys > 0:
        for idx, row in df[df.index.duplicated()].iterrows():
            logging.info("="*100)
            logging.info(f"ERROR: duplicate key: {idx}")
            logging.info(row)

        raise ValueError(f"Found {num_duplicate_keys} duplicate keys in {file_path}")


def parse_combined_str_calls_tsv_path(combined_str_calls_tsv_path):
    combined_str_calls_df = pd.read_table(combined_str_calls_tsv_path)
    combined_str_calls_df.set_index(["SampleId", "LocusId", "VariantId"], inplace=True)
    check_for_duplicate_keys(combined_str_calls_df, combined_str_calls_tsv_path)
    return combined_str_calls_df.reset_index()


def parse_fam_file(fam_file_path):
    fam_file_df = pd.read_table(
        fam_file_path,
        names=[
            "family_id",
            "individual_id",
            "father_id",
            "mother_id",
            "sex",
            "phenotype",
        ],
    )

    fam_file_df = fam_file_df[["individual_id", "father_id", "mother_id"]]
    fam_file_df.set_index("individual_id", inplace=True)
    check_for_duplicate_keys(fam_file_df, fam_file_path)
    return fam_file_df.reset_index()


def group_rows_by_trio(combined_str_calls_df):
    """Returns a list of row tuples"""
    print(f"{len(combined_str_calls_df)} rows before")

    all_rows = {}
    for _, row in combined_str_calls_df.iterrows():
        all_rows[(row.SampleId, row.LocusId, row.VariantId)] = row
        all_rows[(row.Filename, row.LocusId, row.VariantId)] = row

    #all_sample_ids = set(combined_str_calls_df.SampleId)

    print(f"{len(combined_str_calls_df)} rows after")

    #print(f"{len(all_sample_ids)} all_sample_ids")
    #print( all_sample_ids)
    calls_counter = 0
    trio_ids = set()
    trios = []
    other_rows = []
    combined_str_calls_df = combined_str_calls_df[~combined_str_calls_df.father_id.isna() & ~combined_str_calls_df.mother_id.isna()]
    for _, row in combined_str_calls_df.iterrows():
        #if row.SampleId in all_sample_ids and row.father_id in all_sample_ids and row.mother_id in all_sample_ids:
        try:
            father_row = all_rows[(row.father_id, row.LocusId, row.VariantId)]
            mother_row = all_rows[(row.mother_id, row.LocusId, row.VariantId)]
            trio_ids.add((row.SampleId, row.father_id, row.mother_id))
            trios.append((row, father_row, mother_row))
            calls_counter += 1
        except KeyError as e:
            print(f"WARNING: skipping {row.SampleId} because table doesn't contain {e}")
            other_rows.append(row)

    print(f"Processed {calls_counter} calls in {len(trio_ids)} trios")

    return trios, other_rows


def intervals_overlap(i1, i2):
    return not (i1.end < i2.begin or i2.end < i1.begin)


def compute_min_distance_mendelian(proband_allele, parent_alleles):
    return min([abs(int(proband_allele) - int(pa)) for pa in parent_alleles])


def compute_min_distance_mendelian_ci(proband_CI, parent_CIs):
    return min([abs(proband_CI.distance_to(parent_CI)) for parent_CI in parent_CIs])


def compute_mendelian_violations(trios):
    """Compute mendelian violations using both exact genotypes and intervals"""
    counters = collections.defaultdict(int)
    counters_ci = collections.defaultdict(int)
    results = []
    results_ci = []

    results_rows = []
    for proband_row, father_row, mother_row in trios:
        locus_id = proband_row.LocusId
        proband_alleles = proband_row.Genotype.split("/")   # Num Repeats: Allele 1
        father_alleles = father_row.Genotype.split("/")
        mother_alleles = mother_row.Genotype.split("/")

        assert len(proband_alleles) in [1, 2], proband_alleles
        assert len(father_alleles) in [1, 2], (father_alleles, father_row)
        assert len(mother_alleles) in [1, 2], (mother_alleles, mother_row)

        proband_CIs = [Interval(*[int(i) for i in ci.split("-")]) for ci in proband_row.GenotypeConfidenceInterval.split("/")]
        father_CIs = [Interval(*[int(i) for i in ci.split("-")]) for ci in father_row.GenotypeConfidenceInterval.split("/")]
        mother_CIs = [Interval(*[int(i) for i in ci.split("-")]) for ci in mother_row.GenotypeConfidenceInterval.split("/")]

        if len(proband_alleles) == 1:
            # AFAIK ExpansionHunter only outputs len(proband_alleles) == 1 if the proband is male, and the locus is on the X chromosome. 
            # Strictly speaking, this means the X-chrom allele should always be inheritted from the mother. However, since there might be more complicated
            # scenarios (such as in the PAR region), call it ok (not a mendelian violation) even if it matches the father's allele
            ok_mendelian = proband_alleles[0] in father_alleles or proband_alleles[0] in mother_alleles
            distance_mendelian = compute_min_distance_mendelian(proband_alleles[0], father_alleles) + compute_min_distance_mendelian(proband_alleles[0], mother_alleles)

            ok_mendelian_ci = any([intervals_overlap(proband_CIs[0], i) for i in father_CIs]) or any([intervals_overlap(proband_CIs[0], i) for i in mother_CIs])
            distance_mendelian_ci = compute_min_distance_mendelian_ci(proband_CIs[0], father_CIs) + compute_min_distance_mendelian_ci(proband_CIs[0], mother_CIs)

        elif len(proband_alleles) == 2:
            ok_mendelian = (
                    (proband_alleles[0] in father_alleles and proband_alleles[1] in mother_alleles) or
                    (proband_alleles[1] in father_alleles and proband_alleles[0] in mother_alleles))
            distance_mendelian = min(
                compute_min_distance_mendelian(proband_alleles[0], father_alleles) + compute_min_distance_mendelian(proband_alleles[1], mother_alleles),
                compute_min_distance_mendelian(proband_alleles[1], father_alleles) + compute_min_distance_mendelian(proband_alleles[0], mother_alleles),
                )

            ok_mendelian_ci = (
                    (any([intervals_overlap(proband_CIs[0], i) for i in father_CIs]) and any([intervals_overlap(proband_CIs[1], i)for i in mother_CIs])) or
                    (any([intervals_overlap(proband_CIs[1], i) for i in father_CIs]) and any([intervals_overlap(proband_CIs[0], i) for i in mother_CIs])))
            distance_mendelian_ci = min(
                compute_min_distance_mendelian_ci(proband_CIs[0], father_CIs) + compute_min_distance_mendelian_ci(proband_CIs[1], mother_CIs),
                compute_min_distance_mendelian_ci(proband_CIs[1], father_CIs) + compute_min_distance_mendelian_ci(proband_CIs[0], mother_CIs),
                )
        else:
            raise ValueError(f"Unexpected proband_alleles value: {proband_alleles}")

        if ok_mendelian:
            assert distance_mendelian == 0, f"{locus_id}  d:{distance_mendelian}  ({father_row.Genotype} + {mother_row.Genotype} => {proband_row.Genotype})"

        mendelian_results_string = f"{locus_id}  d:{distance_mendelian}  ({father_row.Genotype} + {mother_row.Genotype} => {proband_row.Genotype}) {proband_row.SampleId}  {father_row.SampleId}  {mother_row.SampleId}    {proband_row.SampleId}*_ExpansionHunter4/*{locus_id}*.svg {father_row.SampleId}*_ExpansionHunter4/*{locus_id}*.svg {mother_row.SampleId}*_ExpansionHunter4/*{locus_id}*.svg"
        if not ok_mendelian:
            results += [mendelian_results_string]
            counters[f"{locus_id} ({proband_row.RepeatUnit})"] += 1

        mendelian_ci_results_string = f"{locus_id}   ({father_row.GenotypeConfidenceInterval} + {mother_row.GenotypeConfidenceInterval} => {proband_row.GenotypeConfidenceInterval})  {proband_row.SampleId}  {father_row.SampleId}  {mother_row.SampleId}  {proband_row.SampleId}*_ExpansionHunter4/*{locus_id}*.svg {father_row.SampleId}*_ExpansionHunter4/*{locus_id}*.svg {mother_row.SampleId}*_ExpansionHunter4/*{locus_id}*.svg"
        if not ok_mendelian_ci:
            results_ci += [mendelian_ci_results_string]
            counters_ci[f"{locus_id} ({proband_row.RepeatUnit})"] += 1

        #assert not (ok_mendelian and not ok_mendelian_ci)  # it should never be the case that mendelian inheritance is consistent for exact genotypes, and not consistent for CI interval-overlap.

        results_rows.append({
            'LocusId': f"{locus_id} ({proband_row.VariantId})",
            'ReferenceRegion': proband_row.ReferenceRegion,
            'VariantId': proband_row.VariantId,
            'RepeatUnit': proband_row.RepeatUnit,
            'IsMendelianViolation': not ok_mendelian,
            'IsMendelianViolationCI': not ok_mendelian_ci,

            'MendelianViolationDistance': distance_mendelian,
            'MendelianViolationDistanceCI': distance_mendelian_ci,

            'ProbandGenotype': proband_row.Genotype,
            'ProbandGenotypeCI': proband_row.GenotypeConfidenceInterval,

            'FatherGenotype': father_row.Genotype,
            'FatherGenotypeCI': father_row.GenotypeConfidenceInterval,

            'MotherGenotype': mother_row.Genotype,
            'MotherGenotypeCI': mother_row.GenotypeConfidenceInterval,

            'ProbandSampleId': proband_row.SampleId,
            'FatherSampleId': father_row.SampleId,
            'MotherSampleId': mother_row.SampleId,
            'ProbandSex': proband_row.Sex,

            'ProbandNumSpanningReads': proband_row.NumSpanningReads,
            'ProbandNumFlankingReads': proband_row.NumFlankingReads,
            'ProbandNumInrepeatReads': proband_row.NumInrepeatReads,
            'FatherNumSpanningReads': father_row.NumSpanningReads,
            'FatherNumFlankingReads': father_row.NumFlankingReads,
            'FatherNumInrepeatReads': father_row.NumInrepeatReads,
            'MotherNumSpanningReads': mother_row.NumSpanningReads,
            'MotherNumFlankingReads': mother_row.NumFlankingReads,
            'MotherNumInrepeatReads': mother_row.NumInrepeatReads,

            'ProbandNumAllelesSupportedBySpanningReads': int(proband_row.NumAllelesSupportedBySpanningReads) + int(proband_row.NumAllelesSupportedByFlankingReads) + int(proband_row.NumAllelesSupportedByInrepeatReads),
            'FatherNumAllelesSupportedByFlankingReads': int(father_row.NumAllelesSupportedBySpanningReads) + int(father_row.NumAllelesSupportedByFlankingReads) + int(father_row.NumAllelesSupportedByInrepeatReads),
            'MotherNumAllelesSupportedByInrepeatReads': int(mother_row.NumAllelesSupportedBySpanningReads) + int(mother_row.NumAllelesSupportedByFlankingReads) + int(mother_row.NumAllelesSupportedByInrepeatReads),

            'ProbandNumRepeatsAllele2': proband_row['Num Repeats: Allele 2'],
            'FatherNumRepeatsAllele2 ': father_row['Num Repeats: Allele 2'],
            'MotherNumRepeatsAllele2': mother_row['Num Repeats: Allele 2'],

            'ProbandCoverage': proband_row.Coverage,
            'FatherCoverage': father_row.Coverage,
            'MotherCoverage': mother_row.Coverage,
            'MinCoverage': min(float(proband_row.Coverage), float(father_row.Coverage), float(mother_row.Coverage)),

            'mendelian_results_string': mendelian_results_string,
            'mendelian_ci_results_string': mendelian_ci_results_string,
        })

        #if ok_mendelian_ci and not ok_mendelian:
        #    print("    " + mendelian_results_string)
        #    print("ci: " + mendelian_ci_results_string)

        #print(results_ci[-1])

    for name, c in sorted(counters.items(), key=lambda t: -t[1]):
        print(f"{c:5d}  ({100*float(c)/len(trios):0.0f}%) {name}")

    print("\n".join(sorted(results)))

    for name, c in sorted(counters_ci.items(), key=lambda t: -t[1]):
        print(f"{c:5d}  ({100*float(c)/len(trios):0.0f}%) {name}")

    print("\n".join(sorted(results_ci)))

    return pd.DataFrame(results_rows)


def main():
    args = parse_args()

    combined_str_calls_df = parse_combined_str_calls_tsv_path(args.combined_str_calls_tsv)
    fam_file_df = parse_fam_file(args.fam_file)

    combined_str_calls_df.set_index("SampleId", inplace=True)
    fam_file_df.set_index("individual_id", inplace=True)

    combined_str_calls_df = combined_str_calls_df.join(fam_file_df, how="left").reset_index().rename(columns={"index": "SampleId"})

    trios, other_rows = group_rows_by_trio(combined_str_calls_df)

    mendelian_violations_df = compute_mendelian_violations(trios)
    other_rows_df = pd.DataFrame((r.to_dict() for r in other_rows))

    output_tsv_path = re.sub(".tsv$", "", str(args.combined_str_calls_tsv)) + ".mendelian_violations.tsv"
    mendelian_violations_df.to_csv(output_tsv_path, index=False, header=True, sep="\t")
    print(f"Wrote {len(mendelian_violations_df)} rows to {output_tsv_path}")

    output_tsv_path = re.sub(".tsv$", "", str(args.combined_str_calls_tsv)) + ".non_trio_rows.tsv"
    other_rows_df.to_csv(output_tsv_path, index=False, header=True, sep="\t")
    print(f"Wrote {len(other_rows_df)} rows to {output_tsv_path}")


if __name__ == "__main__":
    main()
