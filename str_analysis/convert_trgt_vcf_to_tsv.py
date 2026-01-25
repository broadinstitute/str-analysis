#!/usr/bin/env python3

"""Converts a TRGT output VCF to TSV format"""

import argparse
import gzip
import re

import numpy as np
import pandas as pd
from tqdm import tqdm

from str_analysis.utils.misc_utils import parse_interval


def parse_args(args_list=None):
    """Parse command line arguments."""
    p = argparse.ArgumentParser(
        description="Convert TRGT VCF to TSV format",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument("--discard-hom-ref", action="store_true",
                   help="Discard hom-ref calls")
    p.add_argument("--dont-output-REF-ALT-fields", action="store_true",
                   help="Exclude the VCF REF and ALT fields from the output.")
    p.add_argument("--parse-genotype-from-AL-field", action="store_true",
                   help="Parse genotype from AL field instead of MC field for single-motif loci.")
    grp = p.add_mutually_exclusive_group()
    grp.add_argument("--parse-reference-region-from-locus-id", action="store_true",
                     help="Parse reference region from locus ID instead of VCF position.")
    grp.add_argument("--set-locus-id", action="store_true",
                     help="Set locus ID to '{chrom}-{start_0based}-{end}-{motif}'.")

    p.add_argument("-m", "--sample-metadata",
                   help="Table of sample annotations to add to output.")
    p.add_argument("--sample-metadata-key", default="sample_id",
                   help="Column name in sample metadata containing sample ID.")

    p.add_argument("--dont-output-allele-table", action="store_true",
                   help="Don't output the per-allele table.")
    p.add_argument("--dont-output-variant-table", action="store_true",
                   help="Don't output the per-variant table.")
    p.add_argument("--dont-output-bed-file", action="store_true",
                   help="Don't output the BED file.")

    p.add_argument("-n", "--n-loci", type=int,
                   help="Only process the first N loci from the VCF.")
    p.add_argument("-o", "--output-prefix",
                   help="Output filename prefix. Defaults to input filename without .vcf extension.")
    p.add_argument("--show-progress-bar", action="store_true",
                   help="Show a progress bar during processing.")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Print additional logging messages.")
    p.add_argument("--sample-id",
                   help="Sample ID. If not specified, parsed from VCF header.")
    p.add_argument("vcf_path", help="TRGT VCF file path")

    args = p.parse_args(args=args_list)

    return args


def get_variant_columns(include_ref_alt, sample_metadata_columns=None):
    """Return the list of variant table columns."""
    cols = [
        "SampleId", "LocusId", "AlleleCount", "VariantId", "RepeatUnit",
        "RepeatUnitLength", "ReferenceRegion", "Genotype", "GenotypeConfidenceInterval",
        "AP", "AM", "SummaryString",
    ]
    if include_ref_alt:
        cols.extend(["Ref", "Alt"])

    for allele_i in 1, 2:
        cols.extend([
            f"Num Repeats: Allele {allele_i}", f"CI start: Allele {allele_i}", f"CI end: Allele {allele_i}",
            f"CI size: Allele {allele_i}", f"Q: Allele {allele_i}",
            f"AllelePurity: Allele {allele_i}", f"MeanMethylation: Allele {allele_i}", f"SpanningReadsPerAllele: Allele {allele_i}",
        ])

    if sample_metadata_columns:
        cols.extend([f"Sample_{c}" for c in sample_metadata_columns])

    return cols


def get_allele_columns(include_ref_alt, sample_metadata_columns=None):
    """Return the list of allele table columns."""
    cols = [
        "SampleId", "LocusId", "AlleleCount", "VariantId", "RepeatUnit",
        "RepeatUnitLength", "ReferenceRegion", "Genotype", "GenotypeConfidenceInterval",
        "AP", "AM",
        "AllelePurity", "MeanMethylation", "SpanningReadsPerAllele",
        "Num Repeats", "CI start", "CI end", "CI size", "Q",
    ]
    if include_ref_alt:
        cols.extend(["Ref", "Alt"])

    if sample_metadata_columns:
        cols.extend([f"Sample_{c}" for c in sample_metadata_columns])

    return cols


def compute_summary_string(repeat_unit, ref_region, genotype, genotype_ci, num_alleles):
    """Compute variant summary string efficiently."""
    chrom, start_0based, end_1based = parse_interval(ref_region)
    reference_locus_size = end_1based - start_0based
    num_repeats_ref = int(reference_locus_size / len(repeat_unit))

    genotype_parts = genotype.split("/")
    allele1 = int(genotype_parts[0])
    allele2 = int(genotype_parts[1]) if len(genotype_parts) > 1 else None

    if allele2 is None:
        het_or_hom = "HEMI"
        allele_numbers = [allele1]
    elif allele1 == allele2:
        het_or_hom = "HOM"
        allele_numbers = [allele1]
    else:
        het_or_hom = "HET"
        allele_numbers = [allele1, allele2]

    ins_or_del_or_ref = []
    for allele_repeats in allele_numbers:
        allele_size = allele_repeats * len(repeat_unit)
        if allele_size == reference_locus_size:
            ins_or_del_or_ref.append("REF")
        elif allele_size > reference_locus_size:
            ins_or_del_or_ref.append("INS")
        else:
            ins_or_del_or_ref.append("DEL")

    return (f"{num_repeats_ref}=>{genotype}[{genotype_ci or 'noCIs'}]"
            f":{repeat_unit}[{len(repeat_unit)}bp]"
            f":{het_or_hom}"
            f":{','.join(ins_or_del_or_ref)}")


def process_vcf_line(line, sample_id, args, sample_metadata_row, variant_columns, allele_columns):
    """Process a single VCF line and return variant/allele rows.

    Returns:
        tuple: (variant_rows, allele_rows)
    """
    fields = line.rstrip('\n').split('\t')
    chrom = fields[0]
    start_1based = int(fields[1])
    ref_seq = fields[3]
    alt_seq = fields[4]
    info = fields[7]

    if len(fields) < 10 or not fields[9] or fields[9] == ".":
        return None, None

    # Parse INFO and FORMAT
    info_dict = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            info_dict[k] = v

    genotype_fields = fields[8].split(":")
    genotype_values = fields[9].split(":")
    genotype_dict = dict(zip(genotype_fields, genotype_values))

    if args.discard_hom_ref and genotype_dict.get("GT") == "0/0":
        return None, None

    if genotype_dict.get("AL", ".") == ".":
        return None, None

    end_1based = int(info_dict["END"])
    motifs = info_dict["MOTIFS"].split(",")
    locus_id = info_dict["TRID"]

    # Sort alleles by size (smaller first)
    allele_sizes_bp = [int(s) for s in genotype_dict["AL"].split(",")]
    flip_alleles = len(allele_sizes_bp) == 2 and allele_sizes_bp[0] > allele_sizes_bp[1]
    if flip_alleles:
        for key in ("AL", "ALLR", "SD", "MC", "MS", "AP", "AM"):
            if key in genotype_dict and genotype_dict[key]:
                genotype_dict[key] = ",".join(genotype_dict[key].split(",")[::-1])

    if args.parse_reference_region_from_locus_id:
        locus_id_fields = re.split("[_-]", locus_id)
        if len(locus_id_fields) >= 3:
            chrom = locus_id_fields[0]
            start_1based = int(locus_id_fields[1]) + 1
            end_1based = int(locus_id_fields[2])

    variant_rows = []
    allele_rows = []

    # Generate variant records for each motif
    motif_records = []
    if args.parse_genotype_from_AL_field and len(motifs) == 1:
        repeat_unit = motifs[0]
        if args.set_locus_id:
            locus_id = f"{chrom}-{start_1based - 1}-{end_1based}-{repeat_unit}"

        genotype_list = [str(int(s) // len(repeat_unit)) for s in genotype_dict["AL"].split(",")]
        genotype = "/".join(genotype_list)

        ci_list = []
        for ci in genotype_dict["ALLR"].split(","):
            ci_parts = ci.split("-")
            ci_parts = [str(int(v) // len(repeat_unit)) for v in ci_parts]
            ci_list.append("-".join(ci_parts))
        genotype_ci = "/".join(ci_list)

        motif_records.append((locus_id, locus_id, repeat_unit, genotype, genotype_ci))
    else:
        # Use MC field
        motif_counts = [[int(c) for c in mc.split("_")] for mc in genotype_dict["MC"].split(",")]
        for motif_i, motif in enumerate(motifs):
            genotype_list = [str(mc[motif_i]) for mc in motif_counts]
            genotype = "/".join(genotype_list)
            genotype_ci = "/".join([f"{g}-{g}" for g in genotype_list])
            variant_id = f"{locus_id}-m{motif_i}-{motif}" if len(motifs) > 1 else locus_id
            motif_records.append((locus_id, variant_id, motif, genotype, genotype_ci))

    # Parse purity values once
    allele_purities = genotype_dict.get("AP", "").split(",") if genotype_dict.get("AP") else ["", ""]
    mean_methylation = genotype_dict.get("AM", "").split(",") if genotype_dict.get("AM") else ["", ""]
    spanning_reads = genotype_dict.get("SD", "").split(",") if genotype_dict.get("SD") else ["", ""]

    ref_region = f"{chrom}:{start_1based - 1}-{end_1based}"
    allele_count = len(allele_sizes_bp)

    for locus_id_rec, variant_id, repeat_unit, genotype, genotype_ci in motif_records:
        genotype_parts = genotype.split("/")
        ci_parts = genotype_ci.split("/") if genotype_ci else ["", ""]

        # Build variant row as a dict for column lookup
        row = {
            "SampleId": sample_id,
            "LocusId": locus_id_rec,
            "AlleleCount": allele_count,
            "VariantId": variant_id,
            "RepeatUnit": repeat_unit,
            "RepeatUnitLength": len(repeat_unit),
            "ReferenceRegion": ref_region,
            "Genotype": genotype,
            "GenotypeConfidenceInterval": genotype_ci,
        }

        if not args.dont_output_REF_ALT_fields:
            row["Ref"] = ref_seq
            row["Alt"] = alt_seq

        row["AP"] = genotype_dict.get("AP", "")
        row["AM"] = genotype_dict.get("AM", "")

        # Process each allele
        for i, (gt, ci) in enumerate(zip(genotype_parts, ci_parts)):
            suffix = f": Allele {i + 1}"
            num_repeats = int(gt)
            row[f"Num Repeats{suffix}"] = num_repeats

            if ci and "-" in ci:
                ci_start, ci_end = map(int, ci.split("-"))
                row[f"CI start{suffix}"] = ci_start
                row[f"CI end{suffix}"] = ci_end
                ci_size = ci_end - ci_start
                row[f"CI size{suffix}"] = ci_size
                ci_ratio = ci_size / (num_repeats or 1)
                row[f"Q{suffix}"] = round(1 / np.exp(4 * ci_ratio), 3)
            else:
                row[f"CI start{suffix}"] = ""
                row[f"CI end{suffix}"] = ""
                row[f"CI size{suffix}"] = ""
                row[f"Q{suffix}"] = ""

            purity = allele_purities[i] if i < len(allele_purities) else ""
            row[f"AllelePurity{suffix}"] = purity if purity != "." else ""
            meth = mean_methylation[i] if i < len(mean_methylation) else ""
            row[f"MeanMethylation{suffix}"] = meth if meth != "." else ""
            sd = spanning_reads[i] if i < len(spanning_reads) else ""
            row[f"SpanningReadsPerAllele{suffix}"] = sd if sd != "." else ""

            # Build allele row
            if not args.dont_output_allele_table:
                allele_row = {
                    "SampleId": sample_id,
                    "LocusId": locus_id_rec,
                    "AlleleCount": allele_count,
                    "VariantId": variant_id,
                    "RepeatUnit": repeat_unit,
                    "RepeatUnitLength": len(repeat_unit),
                    "ReferenceRegion": ref_region,
                    "Genotype": genotype,
                    "GenotypeConfidenceInterval": genotype_ci,
                }
                if not args.dont_output_REF_ALT_fields:
                    allele_row["Ref"] = ref_seq
                    allele_row["Alt"] = alt_seq
                allele_row["AP"] = genotype_dict.get("AP", "")
                allele_row["AM"] = genotype_dict.get("AM", "")
                allele_row["Num Repeats"] = num_repeats
                allele_row["CI start"] = row.get(f"CI start{suffix}", "")
                allele_row["CI end"] = row.get(f"CI end{suffix}", "")
                allele_row["CI size"] = row.get(f"CI size{suffix}", "")
                allele_row["Q"] = row.get(f"Q{suffix}", "")
                allele_row["AllelePurity"] = row.get(f"AllelePurity{suffix}", "")
                allele_row["MeanMethylation"] = row.get(f"MeanMethylation{suffix}", "")
                allele_row["SpanningReadsPerAllele"] = row.get(f"SpanningReadsPerAllele{suffix}", "")
                if sample_metadata_row is not None:
                    for k, v in sample_metadata_row.items():
                        allele_row[f"Sample_{k}"] = v if pd.notna(v) else ""
                allele_rows.append(allele_row)

        # Compute summary string
        row["SummaryString"] = compute_summary_string(repeat_unit, ref_region, genotype, genotype_ci, allele_count)

        # Add sample metadata
        if sample_metadata_row is not None:
            for k, v in sample_metadata_row.items():
                row[f"Sample_{k}"] = v if pd.notna(v) else ""

        variant_rows.append(row)

    return variant_rows, allele_rows


def main():
    """Main entry point."""
    args = parse_args()

    print(f"Processing {args.vcf_path}")

    # Determine output prefix
    output_prefix = args.output_prefix
    if output_prefix is None:
        output_prefix = re.sub(r"\.vcf(\.gz)?$", "", args.vcf_path)

    # Load sample metadata if provided
    sample_metadata_row = None
    sample_metadata_columns = None
    if args.sample_metadata:
        sample_metadata_df = pd.read_table(args.sample_metadata)
        sample_id_col = args.sample_metadata_key
        if sample_id_col not in sample_metadata_df.columns:
            for col in ["sample_id", "SampleId", "sample", "ParticipantId"]:
                if col in sample_metadata_df.columns:
                    sample_id_col = col
                    break
        sample_metadata_columns = list(sample_metadata_df.columns)
        sample_metadata_lookup = {row[sample_id_col]: row for _, row in sample_metadata_df.iterrows()}
        print(f"Loaded {len(sample_metadata_df)} rows from {args.sample_metadata}")

    # Parse sample ID from VCF header
    sample_id = args.sample_id
    fopen = gzip.open if args.vcf_path.endswith("gz") else open
    with fopen(args.vcf_path, "rt") as vcf:
        for line in vcf:
            if line.startswith("#CHROM"):
                header_fields = line.strip().split("\t")
                if sample_id is None and len(header_fields) >= 10:
                    sample_id = header_fields[9]
                    print(f"Got sample ID '{sample_id}' from VCF header")
                break

    # Look up sample metadata
    if args.sample_metadata and sample_id:
        sample_metadata_row = sample_metadata_lookup.get(sample_id)

    # Define columns upfront
    include_ref_alt = not args.dont_output_REF_ALT_fields
    variant_columns = get_variant_columns(include_ref_alt, sample_metadata_columns)
    allele_columns = get_allele_columns(include_ref_alt, sample_metadata_columns)

    # Open output files
    variant_output_file = None
    allele_output_file = None
    if not args.dont_output_variant_table:
        variant_output_file = gzip.open(f"{output_prefix}.variants.tsv.gz", "wt")
        variant_output_file.write("\t".join(variant_columns) + "\n")
    if not args.dont_output_allele_table:
        allele_output_file = gzip.open(f"{output_prefix}.alleles.tsv.gz", "wt")
        allele_output_file.write("\t".join(allele_columns) + "\n")

    variant_records_counter = 0
    allele_records_counter = 0
    loci_processed = 0
    bed_file_records = []

    with fopen(args.vcf_path, "rt") as vcf:
        line_iter = vcf
        if args.show_progress_bar:
            line_iter = tqdm(vcf, unit=" lines", unit_scale=True)

        for line in line_iter:
            if args.n_loci is not None and loci_processed >= args.n_loci:
                break
            if line.startswith("#"):
                continue

            try:
                variant_rows, allele_rows = process_vcf_line(
                    line, sample_id, args, sample_metadata_row, variant_columns, allele_columns
                )
            except Exception as e:
                if args.verbose:
                    print(f"Error parsing line: {e}")
                continue

            if variant_rows is None:
                continue

            loci_processed += 1

            # Write variant rows
            if not args.dont_output_variant_table:
                for row in variant_rows:
                    values = [str(row.get(c, "")) for c in variant_columns]
                    variant_output_file.write("\t".join(values) + "\n")
                    variant_records_counter += 1

                    # Collect BED records
                    if not args.dont_output_bed_file:
                        chrom, start, end = parse_interval(row["ReferenceRegion"])
                        bed_file_records.append([chrom, start, end, row["SummaryString"], "."])

            # Write allele rows
            if not args.dont_output_allele_table and allele_rows:
                for row in allele_rows:
                    values = [str(row.get(c, "")) for c in allele_columns]
                    allele_output_file.write("\t".join(values) + "\n")
                    allele_records_counter += 1

    # Close output files
    if variant_output_file:
        variant_output_file.close()
    if allele_output_file:
        allele_output_file.close()

    # Write BED file
    if not args.dont_output_bed_file:
        with open(f"{output_prefix}.bed", "wt") as f:
            for bed_record in sorted(bed_file_records, key=lambda r: (r[0], r[1])):
                f.write("\t".join(map(str, bed_record)) + "\n")
        print(f"Wrote {len(bed_file_records):,d} records to {output_prefix}.bed")

    if not args.dont_output_variant_table:
        print(f"Wrote {variant_records_counter:,d} records to {output_prefix}.variants.tsv.gz")

    if not args.dont_output_allele_table:
        print(f"Wrote {allele_records_counter:,d} records to {output_prefix}.alleles.tsv.gz")


if __name__ == "__main__":
    main()
