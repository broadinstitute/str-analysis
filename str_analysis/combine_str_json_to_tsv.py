#!/usr/bin/env python3

import argparse
import collections
import simplejson as json
import logging
import os
import pathlib
import re
import sys

import numpy as np
import pandas as pd
from pprint import pformat

from str_analysis.combine_json_to_tsv import get_sample_id_column_index
from str_analysis.utils.misc_utils import parse_interval

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
ALREADY_WARNED_ABOUT = set()  # used for logging


def parse_args(args_list=None):
    """Parse command line args and return the argparse args object"""
    p = argparse.ArgumentParser()
    p.add_argument(
        "-c",
        "--variant-catalog",
        help="path of json variant catalog file. If specified, annotations for each locus in this variant catalog will "
             "be added to the output table as additional columns. This option can be specified more than once to "
             "use multiple variant catalogs as an annotation source.",
        action="append",
        default=[],
    )
    p.add_argument(
        "-m",
        "--sample-metadata",
        help="Table of sample annotations. If specified, all columns from this table will be added to the output table."
    )
    p.add_argument(
        "--sample-metadata-key",
        help="The column name in the --sample-metdata table that contains a sample id. "
             "If not specified, 'sample_id' and other variations like 'SampleId', 'sample', and 'ParticipantId' "
             "will be tried."
    )
    p.add_argument(
        "--include-extra-expansion-hunter-fields",
        action="store_true",
        help="If specified, additional fields from the ExpansionHunter will be added."
    )
    p.add_argument(
        "--include-extra-gangstr-fields",
        action="store_true",
        help="If specified, additional fields from GangSTR will be added. The input json files are expected to be the "
             "result of running convert_gangstr_vcf_to_expansion_hunter_json."
    )
    p.add_argument(
        "--include-extra-hipstr-fields",
        action="store_true",
        help="If specified, additional fields from HipSTR will be added. The input json files are expected to be the "
             "result of running convert_hipstr_vcf_to_expansion_hunter_json."
    )
    p.add_argument(
        "--include-extra-trgt-fields",
        action="store_true",
        help="If specified, additional fields from TRGT will be added. The input json files are expected to be the "
             "result of running convert_trgt_vcf_to_expansion_hunter_json."
    )

    p.add_argument(
        "-o",
        "--output-prefix",
        help="Combined table output filename prefix",
    )
    p.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print additional logging messages",
    )
    p.add_argument(
        "json_paths",
        help="EpxansionHunter output json file and/or directory path(s) to search. If not specified, this script will search for json files in "
             "the current directory and subdirectories.",
        type=pathlib.Path,
        nargs="*"
    )

    args = p.parse_args(args=args_list)

    for variant_catalog in args.variant_catalog:
        if not os.path.isfile(variant_catalog):
            p.error(f"{variant_catalog} not found")

    json_file_paths = []
    directories_to_search = []
    if not args.json_paths:
        directories_to_search.append(".")
    else:
        # check if files exist
        for json_path in args.json_paths:
            if os.path.isdir(json_path):
                directories_to_search.append(json_path)
            elif os.path.isfile(json_path):
                json_file_paths.append(json_path)
            else:
                p.error(f"File not found: {json_path}")

    for dir_path in directories_to_search:
        # find all .json files in this directory
        json_paths = [p for p in pathlib.Path(dir_path).glob("**/*.json")]
        json_file_paths.extend(json_paths)
        print(f"Found {len(json_paths)} .json files in directory: {dir_path}")

    args.json_paths = json_file_paths
    if len(args.json_paths) == 0:
        sys.exit(0)

    return args


def compute_bed_file_record(variant_record):
    chrom, start_0based, end_1based = parse_interval(variant_record["ReferenceRegion"])

    return [chrom, start_0based, end_1based, variant_record["SummaryString"], "."]


def main():
    """Main"""

    args = parse_args()

    output_prefix = args.output_prefix or f"combined_expansion_hunter"
    output_prefix += f".{len(args.json_paths)}_json_files"

    combined_variant_catalog_contents = []
    if args.variant_catalog:
        for variant_catalog in args.variant_catalog:
            variant_catalog_contents = parse_json_file(variant_catalog)
            logging.info(f"Parsed {len(variant_catalog_contents)} records from {variant_catalog}")
            combined_variant_catalog_contents.extend(variant_catalog_contents)
        #for d in variant_catalog_contents:
        #    logging.info("    " + d["LocusId"])

    sample_metadata_lookup = collections.defaultdict(list)
    if args.sample_metadata:
        sample_metadata_df = pd.read_table(args.sample_metadata)
        sample_id_column_idx = get_sample_id_column_index(sample_metadata_df, column_name=args.sample_metadata_key)
        if sample_id_column_idx == -1:
            raise ValueError(f"'sample_id' column not found in sample metadata table. The columns found were: {sample_metadata_df.columns}")
        sample_id_column = sample_metadata_df.columns[sample_id_column_idx]
        for _, row in sample_metadata_df.iterrows():
            if args.verbose and len(sample_metadata_lookup[row[sample_id_column]]) > 0:
                logging.info(f"  {row[sample_id_column]} is a duplicate sample id in {args.sample_metadata}")
            sample_metadata_lookup[row[sample_id_column]].append(row)
        logging.info(f"Parsed {len(sample_metadata_df)} rows from {args.sample_metadata}")

    variant_table_columns = []
    allele_table_columns = []
    wrote_variant_table_header = wrote_allele_table_header = False
    variant_records_counter = allele_records_counter = 0
    variant_output_file = open(f"{output_prefix}.variants.tsv", "wt")
    allele_output_file = open(f"{output_prefix}.alleles.tsv", "wt")
    bed_file_records = []
    sample_metadata_lookup_counters = {}

    for just_get_header in True, False:
        if just_get_header:
            json_paths = args.json_paths[:20]
        else:
            json_paths = args.json_paths

        for json_path in json_paths:
            json_path = str(json_path)  # convert from pathlib.Path to str

            try:
                json_contents = parse_json_file(json_path)
            except Exception as e:
                logging.info(f"Skipping {json_path}... Unable to parse json: {e}")
                continue

            if not isinstance(json_contents, dict) or "SampleParameters" not in json_contents:
                logging.info(f"Skipping {json_path}... Expected key 'SampleParameters' not found.")
                continue

            for variant_record in convert_expansion_hunter_json_to_tsv_columns(
                json_contents,
                variant_catalog_contents=combined_variant_catalog_contents,
                sample_metadata_lookup=sample_metadata_lookup,
                sample_metadata_lookup_counters=sample_metadata_lookup_counters,
                json_file_path=json_path,
                yield_allele_records=False,
                include_extra_expansion_hunter_fields=args.include_extra_expansion_hunter_fields,
                include_extra_gangstr_fields=args.include_extra_gangstr_fields,
                include_extra_hipstr_fields=args.include_extra_hipstr_fields,
                include_extra_trgt_fields=args.include_extra_trgt_fields,
            ):
                if just_get_header:
                    variant_table_columns.extend([k for k in variant_record.keys() if k not in variant_table_columns])
                else:
                    if not wrote_variant_table_header:
                        variant_output_file.write("\t".join(variant_table_columns) + "\n")
                        wrote_variant_table_header = True
                    variant_records_counter += 1
                    variant_output_file.write(
                        "\t".join([str(variant_record.get(c, "")) for c in variant_table_columns]) + "\n")

                    bed_file_records.append(compute_bed_file_record(variant_record))

            for allele_record in convert_expansion_hunter_json_to_tsv_columns(
                json_contents,
                variant_catalog_contents=combined_variant_catalog_contents,
                sample_metadata_lookup=sample_metadata_lookup,
                json_file_path=json_path,
                yield_allele_records=True,
                include_extra_expansion_hunter_fields=args.include_extra_expansion_hunter_fields,
                include_extra_gangstr_fields=args.include_extra_gangstr_fields,
                include_extra_hipstr_fields=args.include_extra_hipstr_fields,
            ):
                if just_get_header:
                    allele_table_columns.extend([k for k in allele_record.keys() if k not in allele_table_columns])
                else:
                    if not wrote_allele_table_header:
                        allele_output_file.write("\t".join(allele_table_columns) + "\n")
                        wrote_allele_table_header = True
                    allele_records_counter += 1
                    allele_output_file.write(
                        "\t".join([str(allele_record.get(c, "")) for c in allele_table_columns]) + "\n")

    if sample_metadata_lookup_counters:
        logging.info(f"Found matches for {sample_metadata_lookup_counters['sample_ids_found']} out of "
              f"{sample_metadata_lookup_counters['total_sample_ids']} "
              f"({100*sample_metadata_lookup_counters['sample_ids_found']/sample_metadata_lookup_counters['total_sample_ids']:0.1f}%) "
              f"ExpansionHunter json sample ids in {args.sample_metadata}")

        if len(sample_metadata_lookup) != len(sample_metadata_df):
            logging.info(f"{len(sample_metadata_df) - len(sample_metadata_lookup)} duplicate sample ids in {args.sample_metadata}")

    logging.info(f"Wrote {variant_records_counter:,d} records to {output_prefix}.variants.tsv")
    logging.info(f"Wrote {allele_records_counter:,d} records to {output_prefix}.alleles.tsv")

    with open(f"{output_prefix}.bed", "wt") as f:
        for bed_record in sorted(bed_file_records, key=lambda r: tuple(r[0:2])):
            f.write("\t".join(map(str, bed_record)) + "\n")

    logging.info(f"Wrote {len(bed_file_records):,d} records to {output_prefix}.bed")


class ParseError(Exception):
    """Represents an error that occurs while parsing"""
    pass


def parse_json_file(json_path):
    """Takes a json file path and loads its contents into a dictionary or list"""

    if not os.path.isfile(json_path):
        raise ValueError(f"{json_path} not found")

    with open(json_path, "rt", encoding="UTF-8") as f:
        try:
            return json.load(f)
        except Exception as e:
            raise ParseError(f"Unable to parse {json_path}: {e}")


def parse_read_count_tuples(tuple_string):
    """Parse the string value of one of the CountsOf* ExpansionHunter output fields
    (CountsOfSpanningReads, CountsOfFlankingReads, or CountsOfInrepeatReads). For example:

        "CountsOfSpanningReads": "(6, 1), (17, 10), (19, 16)"

    Return: a list of 2-tuples where t[0] = number of repeats, t[1] = number of supporting reads for this repeat size
    """
    json_string = tuple_string.replace("(", "[").replace(")", "]")

    return [tuple(t) for t in json.loads(f"[{json_string}]") if t]


def weighted_avg_and_std(values, weights):
    """Return the weighted average and standard deviation for a given array of values, and weights.
    Copied from: https://stackoverflow.com/questions/2413522/weighted-standard-deviation-in-numpy
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return average, np.sqrt(variance)


def compute_variant_summary_string(variant_record):
    """Returns a short easy-to-read string summarizing the tool's STR genotype.

    Args:
        variant_record (dict): list of 1 or more allele spec dictionaries
    Return:
        str: short summary of the variant
    """

    repeat_unit = variant_record["RepeatUnit"]
    chrom, start_0based, end_1based = parse_interval(variant_record["ReferenceRegion"])
    reference_locus_size = end_1based - start_0based
    num_repeats_ref = int(reference_locus_size/len(repeat_unit))

    if f"Repeat Size (bp): Allele 2" not in variant_record:
        het_or_hom = "HEMI"
        allele_numbers = [1]
    elif variant_record[f"Repeat Size (bp): Allele 1"] == variant_record[f"Repeat Size (bp): Allele 2"]:
        het_or_hom = "HOM"
        allele_numbers = [1]
    else:
        het_or_hom = "HET"
        allele_numbers = [1, 2]

    ins_or_del_or_ref = []
    for i in allele_numbers:
        allele_size = int(variant_record[f"Repeat Size (bp): Allele {i}"])
        if allele_size == reference_locus_size:
            ins_or_del_or_ref.append("REF")
        elif allele_size > reference_locus_size:
            ins_or_del_or_ref.append("INS")
        elif allele_size < reference_locus_size:
            ins_or_del_or_ref.append("DEL")

    summary_string = str(num_repeats_ref) + "=>" + variant_record["Genotype"]
    summary_string += "[" + variant_record["GenotypeConfidenceInterval"] + "]"
    summary_string += ":" + f"{repeat_unit}[{len(repeat_unit)}bp]"
    summary_string += ":" + het_or_hom
    summary_string += ":" + ",".join(ins_or_del_or_ref)

    return summary_string


def convert_expansion_hunter_json_to_tsv_columns(
    json_contents,
    variant_catalog_contents=None,
    sample_metadata_lookup=None,
    sample_metadata_lookup_counters=None,
    variant_info=None,
    json_file_path="",
    yield_allele_records=True,
    include_extra_expansion_hunter_fields=False,
    include_extra_gangstr_fields=False,
    include_extra_hipstr_fields=False,
    include_extra_trgt_fields=False,
):
    """Converts a dictionary that represents the contents of an ExpansionHunter v3 or v4 json output file to
    a dictionary of tsv column values.

    Args:
        json_contents (dict): a dict with the contents of an ExpansionHunter v3 or v4 json output file
        variant_catalog_contents (dict): optional dict with the contents of the variant catalog used when running
            expansion hunter. If provided, the fields will be added to the output table.
        sample_metadata_lookup (dict): maps sample id to row of sample metadata
        sample_metadata_lookup_counters (dict): dictionary for accumulating counters and stats about the
            ability to find samples in sample_metadata_lookup
        variant_info (dict): if provided, results will be added to this dict. Otherwise, a new dict will be created.
        json_file_path (str): if provided, it will be added as a field to the output table, and also used for logging
        yield_allele_records (bool): if True, the returned list will have one record per allele rather than per variant
        include_extra_expansion_hunter_fields (bool): if True, include additional fields provided by ExpansionHunter.
        include_extra_gangstr_fields (bool): if True, include additional fields provided by GangSTR.
        include_extra_hipstr_fields (bool): if True, include additional fields provided by HipSTR.
        include_extra_trgt_fields (bool): if True, include additional fields provided by TRGT.
    Yields:
        dict: dictionary representing the output tsv row
    """

    if variant_info is None:
        variant_info = collections.OrderedDict()

    if variant_catalog_contents:
        variant_catalog_contents = {c['LocusId']: c for c in variant_catalog_contents}

    if json_file_path:
        variant_info["Filename"] = json_file_path

    try:
        variant_info["Sex"] = json_contents["SampleParameters"]["Sex"]
        variant_info["SampleId"] = json_contents["SampleParameters"]["SampleId"]
        locus_results_list = json_contents["LocusResults"].values()
    except KeyError as e:
        raise ValueError(f"Key {e} not found in json file {json_file_path}")

    if sample_metadata_lookup:
        sample_metadata_rows = sample_metadata_lookup.get(variant_info["SampleId"])
        if sample_metadata_rows is None:
            sample_id_from_filename = re.sub("([._]expansion_hunter[0-9]?)?.json.*$", "", variant_info["Filename"])
            sample_metadata_rows = sample_metadata_lookup.get(sample_id_from_filename)

        if sample_metadata_rows is None:
            if sample_metadata_lookup_counters is not None:
                sample_metadata_lookup_counters["total_sample_ids"] = sample_metadata_lookup_counters.get(
                    "total_sample_ids", 0) + 1
                sample_metadata_lookup_counters["sample_ids_found"] = sample_metadata_lookup_counters.get(
                    "sample_ids_found", 0)
        else:
            if sample_metadata_lookup_counters is not None:
                sample_metadata_lookup_counters["total_sample_ids"] = sample_metadata_lookup_counters.get(
                    "total_sample_ids", 0) + len(sample_metadata_rows)
                sample_metadata_lookup_counters["sample_ids_found"] = sample_metadata_lookup_counters.get(
                    "sample_ids_found", 0) + len(sample_metadata_rows)

            for sample_metadata_row in sample_metadata_rows:
                row_dict = sample_metadata_row.to_dict()
                for key, value in row_dict.items():
                    key = f"Sample_{key}"

                    if value is not None and not pd.isna(value):
                        output_value = str(value)
                    else:
                        output_value = ""

                    if key not in variant_info:
                        variant_info[key] = output_value
                    else:
                        variant_info[key] += f"; {output_value}"

    for locus_json in locus_results_list:
        locus_record = collections.OrderedDict(variant_info)
        for key in "LocusId", "ReadLength", "Coverage", "AlleleCount":
            locus_record[key] = locus_json[key]

        # if variant catalog specified, transfer info from it
        locus_id = locus_record["LocusId"]
        if variant_catalog_contents:
            locus_record["InVariantCatalog"] = locus_id in variant_catalog_contents
            if locus_record["InVariantCatalog"]:
                excluded_keys = "LocusId", "ReferenceRegion", "VariantId", "VariantType", "OfftargetRegions"
                locus_record.update({f"VariantCatalog_{k}": v for k, v in variant_catalog_contents[locus_id].items() if k not in excluded_keys})
                locus_record["VariantCatalog_NumOfftargetRegions"] = len(variant_catalog_contents[locus_id].get("OfftargetRegions", []))
            else:
                if locus_id not in ALREADY_WARNED_ABOUT:
                    ALREADY_WARNED_ABOUT.add(locus_id)
                    logging.warning(f"{json_file_path} locus id {locus_id} not found in variant catalog.")

        for variant_json in locus_json.get("Variants", {}).values():
            if "Genotype" not in variant_json:
                continue

            # docs @ https://github.com/Illumina/ExpansionHunter/blob/master/docs/05_OutputJsonFiles.md
            variant_record = collections.OrderedDict(locus_record)
            variant_record["VariantId"] = variant_json.get("VariantId", "")
            variant_record["RepeatUnit"] = variant_json["RepeatUnit"]
            variant_record["RepeatUnitLength"] = len(variant_json["RepeatUnit"])
            variant_record["ReferenceRegion"] = variant_json["ReferenceRegion"]

            if include_extra_expansion_hunter_fields:
                variant_record["CountsOfSpanningReads"] = variant_json["CountsOfSpanningReads"]
                variant_record["CountsOfFlankingReads"] = variant_json["CountsOfFlankingReads"]
                variant_record["CountsOfInrepeatReads"] = variant_json["CountsOfInrepeatReads"]

                spanning_read_tuples = parse_read_count_tuples(variant_json["CountsOfSpanningReads"])
                flanking_read_tuples = parse_read_count_tuples(variant_json["CountsOfFlankingReads"])
                inrepeat_read_tuples = parse_read_count_tuples(variant_json["CountsOfInrepeatReads"])

                variant_record["NumSpanningReads"] = sum(t[1] for t in spanning_read_tuples)
                variant_record["NumFlankingReads"] = sum(t[1] for t in flanking_read_tuples)
                variant_record["NumInrepeatReads"] = sum(t[1] for t in inrepeat_read_tuples)
                variant_record["NumReadsTotal"] = sum([variant_record[k] for k in (
                    "NumSpanningReads", "NumFlankingReads", "NumInrepeatReads")])

                variant_record["NumAllelesSupportedBySpanningReads"] = len(spanning_read_tuples)
                variant_record["NumAllelesSupportedByFlankingReads"] = len(flanking_read_tuples)
                variant_record["NumAllelesSupportedByInrepeatReads"] = len(inrepeat_read_tuples)
                variant_record["NumAllelesSupportedTotal"] = sum([variant_record[k] for k in (
                    "NumAllelesSupportedBySpanningReads", "NumAllelesSupportedByFlankingReads", "NumAllelesSupportedByInrepeatReads")])

            if include_extra_gangstr_fields:
                variant_record["RC"] = variant_json["RC"]
                variant_record["ENCLREADS"] = variant_json["ENCLREADS"]
                variant_record["FLNKREADS"] = variant_json["FLNKREADS"]
                enclosing, spanning, FRR, flanking = variant_json["RC"].split(",")
                variant_record["NumSpanningReads"] = int(enclosing)
                variant_record["NumFlankingReads"] = int(flanking)  # + int(spanning)
                variant_record["NumInrepeatReads"] = int(FRR)
                variant_record["NumReadsTotal"] = sum([variant_record[k] for k in (
                    "NumSpanningReads", "NumFlankingReads", "NumInrepeatReads")])
                variant_record["Q"] = float(variant_json["Q"])

            if include_extra_hipstr_fields:
                variant_record["Q"] = float(variant_json["Q"])
                variant_record["DP"] = float(variant_json["DP"])
                variant_record["AB"] = float(variant_json["AB"])
                variant_record["FS"] = float(variant_json["FS"])
                variant_record["DFLANKINDEL"] = float(variant_json["DFLANKINDEL"])
                variant_record["DSTUTTER"] = float(variant_json["DSTUTTER"])

            if include_extra_trgt_fields:
                variant_record["AllelePurity"] = variant_json["AP"]
                variant_record["MeanMethylation"] = variant_json["AM"]
                variant_record["SpanningReadsPerAllele"] = variant_json["SD"]

            variant_record["Genotype"] = variant_json["Genotype"]
            variant_record["GenotypeConfidenceInterval"] = variant_json["GenotypeConfidenceInterval"]
            genotype_tuples = list(zip(
                variant_json["Genotype"].split("/"),
                variant_json["GenotypeConfidenceInterval"].split("/"),
            ))

            for i, (genotype, genotypeCI) in enumerate(genotype_tuples):
                if yield_allele_records:
                    suffix = ""
                    output_record = collections.OrderedDict(variant_record)
                else:
                    suffix = f": Allele {i+1}"
                    output_record = variant_record

                try:
                    confidence_interval_start, confidence_interval_end = genotypeCI.split("-")
                except ValueError as e:
                    raise ValueError(f"Unexpected confidence interval format: {genotypeCI} in: {pformat(variant_json)}")
                output_record[f"Allele Number{suffix}"] = i + 1
                output_record[f"Num Repeats{suffix}"] = int(genotype)
                output_record[f"Repeat Size (bp){suffix}"] = int(genotype) * len(variant_json.get("RepeatUnit", ""))
                output_record[f"CI start{suffix}"] = int(confidence_interval_start)
                output_record[f"CI end{suffix}"] = int(confidence_interval_end)
                output_record[f"CI size{suffix}"] = int(confidence_interval_end) - int(confidence_interval_start)
                output_record[f"CI ratio{suffix}"] = output_record[f"CI size{suffix}"]/(output_record[f"Num Repeats{suffix}"] or 1)

                if include_extra_expansion_hunter_fields:
                    is_homozygous = len(genotype_tuples) > 1 and genotype_tuples[0][0] == genotype_tuples[1][0]
                    divisor = 2 if is_homozygous else 1
                    read_length = locus_record.get("ReadLength", 150)
                    ru_length = variant_record["RepeatUnitLength"]
                    for label, read_tuples in [("Spanning", spanning_read_tuples), ("Flanking", flanking_read_tuples), ("Inrepeat", inrepeat_read_tuples)]:
                        output_record[f"Num{label}ReadsThatSupportGenotype{suffix}"] = sum(
                            t[1] for t in read_tuples
                            if (t[0] == int(genotype)) or (int(genotype) > t[0] and t[0]*ru_length > 0.8*read_length)
                        ) / divisor

                    output_record[f"NumReadsTotalThatSupportGenotype{suffix}"] = (
                        output_record[f"NumSpanningReadsThatSupportGenotype{suffix}"] +
                        output_record[f"NumFlankingReadsThatSupportGenotype{suffix}"] +
                        output_record[f"NumInrepeatReadsThatSupportGenotype{suffix}"]
                    )
                    output_record[f"FractionOfReadsThatSupportsGenotype{suffix}"] = (
                        output_record[f"NumReadsTotalThatSupportGenotype{suffix}"] / float(output_record["NumReadsTotal"]) if int(output_record["NumReadsTotal"]) > 0 else 0
                    )

                    # ExpansionHunter Q score based on EnsemblTR https://github.com/gymrek-lab/EnsembleTR/blob/main/ensembletr/utils.py#L53-L59
                    output_record[f"Q{suffix}"] = 1/np.exp(4*output_record[f"CI ratio{suffix}"])

                if include_extra_gangstr_fields:
                    # ex. "2,10|3,7|4,14"
                    for source_field, dest_field in [("ENCLREADS", f"NumSpanningReadsThatSupportGenotype{suffix}"),
                                                     ("FLNKREADS", f"NumFlankingReadsThatSupportGenotype{suffix}")]:
                        if variant_json[source_field] != "NULL":
                            read_support = dict([key_value.split(",") for key_value in variant_json[source_field].split("|")])
                            output_record[dest_field] = int(read_support.get(genotype, 0))
                        else:
                            output_record[dest_field] = 0

                    output_record[f"NumReadsTotalThatSupportGenotype{suffix}"] = (
                            output_record.get(f"NumSpanningReadsThatSupportGenotype{suffix}", 0) +
                            output_record.get(f"NumFlankingReadsThatSupportGenotype{suffix}", 0)
                    )
                    output_record[f"FractionOfReadsThatSupportsGenotype{suffix}"] = (
                        output_record[f"NumReadsTotalThatSupportGenotype{suffix}"] / float(output_record["NumReadsTotal"]) if int(output_record["NumReadsTotal"]) > 0 else 0
                    )
                if yield_allele_records:
                    yield output_record

            if not yield_allele_records:
                output_record["SummaryString"] = compute_variant_summary_string(output_record)
                yield output_record


if __name__ == "__main__":
    main()
