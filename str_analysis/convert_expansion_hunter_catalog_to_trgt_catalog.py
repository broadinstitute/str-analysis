"""This script converts an ExpansionHunter variant catalog to a TRGT catalog, including proper handling of
 variant catalog entries that include adjacent repeats (eg. with LocusStructure like "(A)*(ACG)*").
"""

import argparse
import gzip
import ijson
import os
import pysam
import re
import tqdm

from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.fasta_utils import create_normalize_chrom_function


def main():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-R", "--reference-fasta", required=True, help="Reference genome FASTA file")
    p.add_argument("-s", "--split-adjacent-repeats", action="store_true", help="If a locus is defined "
                   "as having several adjacent repeats in the ExpansionHunter catalog, split it into separate "
                   "entries in the TRGT catalog. This can simplify downstream analysis.")
    p.add_argument("--keep-wide-boundaries", action="store_true", help="When the --split-adjacent-repeats option is used, "
                   "use the outer boundaries of the adjacent repeats list as the locus boundaries.")
    p.add_argument("--set-locus-id", action="store_true", help=f"Set the ID field to 'chom-start_0based-end-motif'"
                                                               f"instead of using the LocusId field from the "
                                                               f"ExpansionHunter catalog.")
    p.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    p.add_argument("-o", "--output-file", help="BED file output path")
    p.add_argument("expansion_hunter_catalog", help="ExpansionHunter variant catalog in JSON format")
    args = p.parse_args()

    if not args.output_file:
        args.output_file = re.sub(".json(.gz)?$", "", args.expansion_hunter_catalog) + ".trgt.bed"
    elif not args.output_file.endswith(".bed") and not args.output_file.endswith(".bed.gz"):
        p.error(f"Output file does not have a .bed extension: {args.output_file}")
    elif args.output_file.endswith(".gz"):
        args.output_file = args.output_file[:-3]

    if args.keep_wide_boundaries and not args.split_adjacent_repeats:
        p.error("The --wide-boundaries option can only be used in combination with the --split-adjacent-repeats option")


    if not os.path.isfile(args.reference_fasta):
        p.error(f"Reference genome FASTA file not found: {args.reference_fasta}")

    process_expansion_hunter_catalog(args.reference_fasta, args.expansion_hunter_catalog, args.output_file,
                                     split_adjacent_repeats=args.split_adjacent_repeats,
                                     keep_wide_boundaries=args.keep_wide_boundaries,
                                     set_locus_id=args.set_locus_id,
                                     show_progress_bar=args.show_progress_bar)


def convert_expansion_hunter_record_to_trgt_row(i, record, split_adjacent_repeats=False,
                                                keep_wide_boundaries=False, set_locus_id=False):
    locus_id = record["LocusId"]
    locus_structure = record["LocusStructure"]
    motifs = re.findall("[(]([A-Z]+)[)]", locus_structure)
    if not motifs:
        raise ValueError(f"Unable to parse LocusStructure '{locus_structure}' in variant catalog "
                         f"record #{i+1}: {record}")

    reference_regions = record["ReferenceRegion"]
    if not isinstance(reference_regions, list):
        reference_regions = [reference_regions]

    if len(motifs) != len(reference_regions):
        raise ValueError(f"LocusStructure elements != # of entries in the list of ReferenceRegions in "
                         f"variant catalog record #{i+1}: {record}")

    chrom = locus_start_0based = locus_end_1based = None
    for motif, reference_region in zip(motifs, reference_regions):
        chrom, start_0based, end_1based = parse_interval(reference_region)

        locus_start_0based = min(locus_start_0based, start_0based) if locus_start_0based is not None else start_0based
        locus_end_1based = max(locus_end_1based, end_1based) if locus_end_1based is not None else end_1based

    if locus_start_0based >= locus_end_1based:
        print(f"WARNING: Skipping locus {locus_id} because its ReferenceRegion "
              f"{chrom}:{locus_start_0based}-{locus_end_1based} has a width = "
              f"{locus_end_1based - locus_start_0based}bp")
        return []

    if "|" in locus_structure:
        print(f"WARNING: Skipping locus {locus_id} @ {chrom}:{locus_start_0based+1}-{locus_end_1based} because "
              f"its LocusStructure {locus_structure} contains a sequence swap operation '|' which is not "
              f"supported by TRGT.")
        return []

    results = []
    if split_adjacent_repeats:
        for motif, reference_region in zip(motifs, reference_regions):
            chrom, start_0based, end_1based = parse_interval(reference_region)

            if set_locus_id:
                locus_label = f"{chrom.replace('chr', '')}-{start_0based}-{end_1based}-{motif}"
            else:
                locus_label = f"{locus_id}-{motif}" if len(motifs) > 1 else locus_id
            struc = f"({motif})n"
            results.append([
                chrom,
                start_0based if not keep_wide_boundaries else locus_start_0based,
                end_1based if not keep_wide_boundaries else locus_end_1based,
                f"ID={locus_label};MOTIFS={motif};STRUC={struc}",
            ])
    else:
        motif_string = ",".join(motifs)
        # Handle the different LocusStructure regular expression operations described in
        #   https://github.com/Illumina/ExpansionHunter/blob/master/docs/04_VariantCatalogFiles.md#using-regular-expressions-to-define-locus-structure
        struc = locus_structure.replace("*", "n").replace("+", "n").replace("?", "n")

        results.append([
            chrom,
            locus_start_0based,
            locus_end_1based,
            f"ID={locus_id};MOTIFS={motif_string};STRUC={struc}",
        ])

    return results

def process_expansion_hunter_catalog(reference_fasta_path, expansion_hunter_catalog_path, output_file_path,
                                     split_adjacent_repeats=False, keep_wide_boundaries=False, set_locus_id=False,
                                     show_progress_bar=False):

    ref_fasta = pysam.FastaFile(reference_fasta_path)
    does_chrom_start_with_chr = ref_fasta.references[0].startswith("chr")
    normalize_chrom = create_normalize_chrom_function(does_chrom_start_with_chr)

    print(f"Parsing {expansion_hunter_catalog_path}")
    counter = 0
    total = 0
    fopen = gzip.open if expansion_hunter_catalog_path.endswith("gz") else open
    with fopen(expansion_hunter_catalog_path, "rt") as f:
        iterator = ijson.items(f, "item", use_float=True)
        if show_progress_bar:
            iterator = tqdm.tqdm(iterator, unit=" variant catalog records", unit_scale=True)

        fopen2 = gzip.open if output_file_path.endswith("gz") else open
        with fopen2(output_file_path, "wt") as f2:
            for i, record in enumerate(iterator):
                total += 1
                for output_row in convert_expansion_hunter_record_to_trgt_row(
                    i, record, split_adjacent_repeats=split_adjacent_repeats,
                    keep_wide_boundaries=keep_wide_boundaries, set_locus_id=set_locus_id):
                    counter += 1
                    output_row[0] = normalize_chrom(output_row[0])
                    f2.write("\t".join(map(str, output_row)) + "\n")

    bgzip_step = "| bgzip" if output_file_path.endswith("gz") else ""
    os.system(f"bedtools sort -i {output_file_path} {bgzip_step} > {output_file_path}.sorted")
    os.system(f"mv {output_file_path}.sorted {output_file_path}")
    #os.system(f"bgzip -f {output_file_path}")
    #os.system(f"tabix -f {output_file_path}.gz")
    print(f"Wrote {counter:,d} out of {total:,d} rows to {output_file_path}")


if __name__ == "__main__":
    main()
