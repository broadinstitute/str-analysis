"""This script reads a variant catalog in BED or JSON format and writes a new catalog with motifs adjusted
to maximize purity for each interval. The output is written in the same format as the input (or in the format
specified by --output-format), preserving input order. BED output is bgzipped and tabix-indexed if already sorted.
"""

import argparse
import gzip
import os
import re
import simplejson as json
import pyfaidx
import tqdm

from str_analysis.utils.eh_catalog_utils import (
    get_variant_catalog_iterator,
    parse_motifs_from_locus_structure,
)
from str_analysis.utils.find_motif_utils import adjust_motif_to_maximize_purity_in_interval
from str_analysis.utils.misc_utils import parse_interval


def parse_args():
    parser = argparse.ArgumentParser(
        description="Reads a variant catalog and writes a new catalog with motifs adjusted to maximize purity "
                    "for each interval.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-R", "--reference-fasta", required=True, help="Reference FASTA file path")
    parser.add_argument("-o", "--output-file", help="Output file path. If not specified, will be derived from "
                        "the input file name.")
    parser.add_argument("--output-format", choices=("JSON", "BED"), help="Output file format. If not specified, "
                        "the format will match the input file.")
    parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose output")
    parser.add_argument("variant_catalog", help="Path to a variant catalog in JSON or BED format")

    args = parser.parse_args()

    if not os.path.isfile(args.variant_catalog):
        parser.error(f"File not found: {args.variant_catalog}")

    if not os.path.isfile(args.reference_fasta):
        parser.error(f"Reference FASTA not found: {args.reference_fasta}")

    return args


def determine_file_format(path):
    """Determine whether a file is JSON or BED based on its extension."""
    if any(path.endswith(suffix) for suffix in [".json", ".json.gz"]):
        return "JSON"
    elif any(path.endswith(suffix) for suffix in [".bed", ".bed.gz", ".bed.bgz"]):
        return "BED"
    else:
        raise ValueError(f"Unrecognized file extension: {path}")


def process_record(record, fasta_obj, counters, verbose=False):
    """Process a single record, adjusting its motif to maximize purity.

    Args:
        record (dict): The input record
        fasta_obj: pyfaidx.Fasta object for the reference genome
        counters (dict): Counter dictionary to update
        verbose (bool): Whether to print verbose output

    Returns:
        dict: The processed record (possibly with adjusted motif)
    """
    counters["total"] += 1

    locus_structure = record["LocusStructure"]
    motifs = parse_motifs_from_locus_structure(locus_structure)

    if len(motifs) != 1:
        if verbose:
            print(f"Skipping locus with multiple motifs: {record['LocusId']}")
        counters["skipped"] += 1
        return record

    original_motif = motifs[0]
    reference_region = record["ReferenceRegion"]

    if isinstance(reference_region, list):
        if verbose:
            print(f"Skipping locus with multiple reference regions: {record['LocusId']}")
        counters["skipped"] += 1
        return record

    chrom, start_0based, end_1based = parse_interval(reference_region)

    # Normalize chromosome name for pyfaidx lookup
    if chrom.startswith("chr"):
        fasta_chrom = chrom
    else:
        fasta_chrom = f"chr{chrom}"

    if fasta_chrom not in fasta_obj:
        fasta_chrom = chrom.replace("chr", "")
        if fasta_chrom not in fasta_obj:
            if verbose:
                print(f"Chromosome {chrom} not found in reference FASTA, skipping {record['LocusId']}")
            counters["skipped"] += 1
            return record

    adjusted_motif, purity = adjust_motif_to_maximize_purity_in_interval(
        fasta_obj, fasta_chrom, start_0based, end_1based, original_motif
    )

    if adjusted_motif != original_motif:
        counters["modified"] += 1
        if verbose:
            print(f"{record['LocusId']}: {original_motif} -> {adjusted_motif} (purity: {purity:.3f})")

        new_record = dict(record)
        new_record["LocusStructure"] = f"({adjusted_motif})*"

        if original_motif in new_record["LocusId"]:
            new_record["LocusId"] = new_record["LocusId"].replace(original_motif, adjusted_motif)

        return new_record

    return record


def record_to_bed_tuple(record):
    """Convert a single JSON record to a BED tuple.

    Args:
        record (dict): The record to convert

    Returns:
        tuple: (chrom, start_0based, end_1based, motif, motif_length)
    """
    locus_structure = record["LocusStructure"]
    reference_region = record["ReferenceRegion"]
    motifs = parse_motifs_from_locus_structure(locus_structure)
    motif = motifs[0]
    chrom, start_0based, end_1based = parse_interval(reference_region)
    return (chrom, start_0based, end_1based, motif, len(motif))


def process_catalog(input_path, reference_fasta_path, output_path, output_format, show_progress_bar=False, verbose=False):
    """Process a variant catalog, adjusting motifs to maximize purity.

    Args:
        input_path (str): Path to the input variant catalog
        reference_fasta_path (str): Path to the reference FASTA file
        output_path (str): Path for the output catalog
        output_format (str): Output format ("JSON" or "BED")
        show_progress_bar (bool): Whether to show a progress bar
        verbose (bool): Whether to print verbose output
    """
    fasta_obj = pyfaidx.Fasta(reference_fasta_path, one_based_attributes=False, as_raw=True)

    print(f"Processing {input_path}")

    counters = {"total": 0, "modified": 0, "skipped": 0}
    iterator = get_variant_catalog_iterator(input_path, show_progress_bar=show_progress_bar)

    if output_format == "JSON":
        fopen = gzip.open if output_path.endswith("gz") else open
        with fopen(output_path, "wt") as f:
            f.write("[\n")
            first_record = True
            for record in iterator:
                processed = process_record(record, fasta_obj, counters, verbose)
                if not first_record:
                    f.write(",\n")
                f.write("    " + json.dumps(processed, ignore_nan=True))
                first_record = False
            f.write("\n]\n")

    elif output_format == "BED":
        output_path_uncompressed = re.sub(r"\.b?gz$", "", output_path)
        is_sorted = True
        prev_chrom = None
        prev_start = None
        seen_chroms = set()

        with open(output_path_uncompressed, "wt") as f:
            for record in iterator:
                processed = process_record(record, fasta_obj, counters, verbose)
                bed_tuple = record_to_bed_tuple(processed)
                chrom, start = bed_tuple[0], bed_tuple[1]

                if is_sorted:
                    if chrom == prev_chrom:
                        # Within same chromosome, check start coordinate order
                        if start < prev_start:
                            is_sorted = False
                    elif chrom in seen_chroms:
                        # Chromosome appeared before but not contiguously
                        is_sorted = False

                seen_chroms.add(chrom)
                prev_chrom = chrom
                prev_start = start

                f.write("\t".join(map(str, bed_tuple)) + "\n")

        if is_sorted:
            os.system(f"bgzip -f {output_path_uncompressed}")
            os.system(f"tabix -f -p bed {output_path_uncompressed}.gz")
            output_path = f"{output_path_uncompressed}.gz"
        else:
            print("Note: Output is not sorted, skipping bgzip/tabix indexing")
            output_path = output_path_uncompressed

    print(f"Processed {counters['total']:,d} loci")
    print(f"  Modified: {counters['modified']:,d} ({counters['modified']/max(1, counters['total']):.1%})")
    print(f"  Skipped: {counters['skipped']:,d} ({counters['skipped']/max(1, counters['total']):.1%})")
    print(f"Wrote {counters['total']:,d} records to {output_path}")


def main():
    args = parse_args()

    input_format = determine_file_format(args.variant_catalog)
    output_format = args.output_format or input_format

    if args.output_file:
        output_path = args.output_file
    else:
        # Generate output filename based on input
        base_name = re.sub(r"\.(json|bed)(\.gz)?$", "", args.variant_catalog, flags=re.IGNORECASE)
        if output_format == "JSON":
            output_path = f"{base_name}.maximized_purity.json.gz"
        else:
            output_path = f"{base_name}.maximized_purity.bed"

    process_catalog(
        args.variant_catalog,
        args.reference_fasta,
        output_path,
        output_format,
        show_progress_bar=args.show_progress_bar,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
