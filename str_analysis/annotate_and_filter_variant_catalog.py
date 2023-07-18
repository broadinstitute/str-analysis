"""This script takes a variant catalog either in .json or .bed format and filters it.
"""


import argparse
import collections
import ijson
import json
import os
import pandas as pd
import re
from tqdm import tqdm

from intervaltree import IntervalTree, Interval
from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.file_utils import open_file, file_exists
from str_analysis.utils.gtf_utils import compute_genomic_region_of_interval
from str_analysis.utils.misc_utils import parse_interval


VALID_GENE_REGIONS = {"CDS", "UTR", "5UTR", "3UTR", "promoter", "exon", "intron", "intergenic"}


def parse_args():
    parser = argparse.ArgumentParser(description="Annotate and filter a variant catalog.")
    parser.add_argument("--known-disease-associated-loci",
                        help="ExpansionHunter catalog .json file with all known disease-associated loci",
                        default="~/code/str-analysis/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json")
    #default="https://raw.githubusercontent.com/broadinstitute/str-analysis/main/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json")
    parser.add_argument("--genes-gtf", help="Gene models gtf file path or url.",
                        default="~/code/str-truth-set/ref/other/MANE.v1.0.ensembl_genomic.sorted.gtf.gz")
    parser.add_argument("--gene-models-source", help="Source of the genes-gtf file. If not specified, it will be "
                        "computed based on the filename", choices=["gencode", "mane", "refseq"])

    parser.add_argument("-o", "--output-path", help="Output path")
    parser.add_argument("--verbose", action="store_true", help="Print verbose output.")

    annotations_group = parser.add_argument_group("Annotations")
    annotations_group.add_argument("--skip-gene-annotations", action="store_true", help="Don't addd gene annotations to "
                                                                                        "the output catalog")
    annotations_group.add_argument("--skip-disease-loci-annotations", action="store_true", help="Don't annotate known "
                                   "disease associated loci in the output catalog")

    filter_group = parser.add_argument_group("Filters")
    filter_group.add_argument("--min-motif-size", type=int, help="Minimum motif size to include in the output catalog")
    filter_group.add_argument("--max-motif-size", type=int, help="Maximum motif size to include in the output catalog")
    filter_group.add_argument("-ms", "--motif-size", type=int, action="append", help="Only include loci with these motif sizes")
    filter_group.add_argument("-m", "--motif", action="append", help="Only include loci whose canonical motif matched this motif")

    filter_group.add_argument("-t", "--region-type", action="append", help="Gene region(s) to include in the output "
                              "catalog", choices=VALID_GENE_REGIONS)
    filter_group.add_argument("-xt", "--exclude-region-type", action="append", choices=VALID_GENE_REGIONS,
                              help="Gene region(s) to exclude from the output catalog",)

    filter_group.add_argument("-g", "--gene-name", action="append", help="Only include loci in this gene")
    filter_group.add_argument("-xg", "--exclude-gene-name", action="append", help="Exclude loci in this gene")

    filter_group.add_argument("--gene-id", action="append", help="Only include loci with this gene id")
    filter_group.add_argument("--exclude-gene-id", action="append", help="Exclude loci with this gene id")

    filter_group.add_argument("-ok", "--only-known-disease-associated-loci", action="store_true",
                        help="Only include known disease-associated loci")
    filter_group.add_argument("-xk", "--exclude-known-disease-associated-loci", action="store_true",
                        help="Exclude known disease-associated loci")

    filter_group.add_argument("-om", "--only-known-disease-associated-motifs", action="store_true",
                        help="Only include loci that have known disease-associated motifs")

    filter_group.add_argument("-l", "--locus-id", action="append", help="Only include the locus with this locus id")
    filter_group.add_argument("-xl", "--exclude-locus-id", action="append", help="Exclude the locus with this locus id")

    parser.add_argument("variant_catalog_json_or_bed", help="A catalog of repeats to annotate and filter, either "
                        "in JSON or BED format. For BED format, the chrom, start, and end should represent the repeat "
                        "interval in 0-based coordinates, and the name field (column #4) should be the repeat unit.")

    args = parser.parse_args()

    for file_path in [args.known_disease_associated_loci, args.genes_gtf, args.variant_catalog_json_or_bed]:
        if not file_exists(os.path.expanduser(file_path)):
            parser.error(f"File not found: {file_path}")

    return args, parser


KNOWN_DISEASE_ASSOCIATED_LOCI_COLUMNS = [
    'LocusId', 'LocusStructure', 'RepeatUnit', 'MainReferenceRegion', 'Gene', 'Inheritance', 'GeneRegion', 'GeneId']


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
        canonical_motif = compute_canonical_motif(row["RepeatUnit"])
        i = Interval(start, end, data=argparse.Namespace(
            canonical_motif=canonical_motif,
            motif=row["RepeatUnit"],
            locus_id=row["LocusId"],
        ))
        known_disease_loci_it[chrom].add(i)

        canonical_motifs.add(canonical_motif)

    return known_disease_loci_it, canonical_motifs


def get_overlap(interval_tree, chrom, start_0based, end, canonical_motif=None):
    """Returns overlapping interval(s) from the interval tree, if any.

    Args:
        interval_tree (dict): a dictionary that maps chromosome names to IntervalTrees
        chrom (str): chromosome name
        start_0based (int): start position of the interval to check for overlap
        end (int): end position of the interval to check for overlap
        canonical_motif (str): the canonical motif to match. If specified, only consider an interval as overlapping if
            its motif also matches

    Returns:
        list of strings: locus ids of entries in the IntervalTree that overlap the given interval chrom:start_0based-end
    """
    chrom = chrom.replace("chr", "")
    for locus_interval in interval_tree[chrom].overlap(start_0based, end):
        if not canonical_motif:
            return True

        if locus_interval.data.canonical_motif != canonical_motif:
            continue
        matching_locus_id = locus_interval.data.locus_id
        return matching_locus_id

    if not canonical_motif:
        return False
    else:
        return None


def get_variant_catalog_iterator(variant_catalog_json_or_bed):
    """Takes the path of a JSON or BED file and returns an iterator over variant catalog records parsed from that file.

    Args:
        variant_catalog_json_or_bed (str): path to a JSON or BED file containing variant catalog records

    Yields:
        dict: a variant catalog record parsed from the file
    """

    if ".json" in variant_catalog_json_or_bed:
        with open_file(variant_catalog_json_or_bed, "rt") as f:
            for record in ijson.items(f, "item"):
                yield record
    else:        
        with open_file(variant_catalog_json_or_bed, "rt") as input_variant_catalog:
            for line in input_variant_catalog:
                fields = line.strip().split("\t")
                chrom = fields[0].replace("chr", "")
                start_0based = int(fields[1])
                end_1based = int(fields[2])
                motif = fields[3].strip("()*+")
                record = {
                    "LocusId": f"{chrom}-{start_0based + 1}-{end_1based}-{motif}",
                    "ReferenceRegion": f"{chrom}:{start_0based + 1}-{end_1based}",
                    "LocusStructure": f"({motif})*",
                    "VariantType": "Repeat",
                }
                yield record


def main():
    args, parser = parse_args()

    if args.verbose:
        print("Args:")
        for key, value in sorted(args.__dict__.items(), key=lambda x: str(x[1])):
            if value is None:
                continue
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

    if not args.output_path:
        filename_prefix = re.sub("(.json|.bed)(.b?gz)?$", "", os.path.basename(args.variant_catalog_json_or_bed))
        args.output_path = f"{filename_prefix}.filtered.json"

    output_records = []
    input_variant_catalog_iterator = get_variant_catalog_iterator(args.variant_catalog_json_or_bed)
    if args.verbose:
        input_variant_catalog_iterator = tqdm(input_variant_catalog_iterator, unit=" variant catalog records")

    for i, input_variant_catalog_record in enumerate(input_variant_catalog_iterator):
        # validate input_variant_catalog_record
        error = False
        for key in ["LocusId", "ReferenceRegion", "LocusStructure", "VariantType"]:
            if key not in input_variant_catalog_record:
                print(f"ERROR: {key} not found in variant catalog record #{i+1}: {input_variant_catalog_record}. Skipping...")
                error = True
                break
        if error:
            continue

        # parse input_variant_catalog_record based on whether it has adjacent repeats
        locus_id = input_variant_catalog_record["LocusId"]
        if args.locus_id and locus_id not in args.locus_id:
            continue
        if args.exclude_locus_id and locus_id in args.exclude_locus_id:
            continue

        if isinstance(input_variant_catalog_record["ReferenceRegion"], list):
            motifs = []
            for m in input_variant_catalog_record["LocusStructure"].strip("()*+").split(")"):
                if "(" in m:
                    m = m.split("(")[1]
                motifs.append(m.strip("()*+"))
            reference_regions = input_variant_catalog_record["ReferenceRegion"]
            variant_types = input_variant_catalog_record["VariantType"]
        else:
            motifs = [input_variant_catalog_record["LocusStructure"].strip("()*+")]
            reference_regions = [input_variant_catalog_record["ReferenceRegion"]]
            variant_types = [input_variant_catalog_record["VariantType"]]

        # validate input_variant_catalog_record with adjacent repeats
        if len(motifs) != len(reference_regions):
            print(f"ERROR: number of motifs != number of reference regions in variant catalog record #{i+1}: {input_variant_catalog_record}. Skipping...")
            continue
        if len(motifs) != len(variant_types):
            print(f"ERROR: number of motifs != number of variant types in variant catalog record #{i+1}: {input_variant_catalog_record}. Skipping...")
            continue

        # parse intervals
        chroms_start_0based_ends = [parse_interval(reference_region) for reference_region in reference_regions]
        canonical_motifs = [compute_canonical_motif(motif) for motif in motifs]

        # apply motif size and motif filters
        if args.min_motif_size and all(len(motif) < args.min_motif_size for motif in motifs):
            continue
        if args.max_motif_size and all(len(motif) > args.max_motif_size for motif in motifs):
            continue
        if args.motif_size and all(len(motif) not in args.motif_size for motif in motifs):
            continue
        if args.motif and all(motif not in args.motif for motif in motifs):
            continue

        if len(canonical_motifs) == 1:
            input_variant_catalog_record["CanonicalMotif"] = canonical_motifs[0]
            input_variant_catalog_record["MotifSize"] = len(canonical_motifs[0])

        # parse input variant catalog record and check for overlap with known disease-associated loci
        matched_known_disease_associated_locus = False
        matched_known_disease_associated_motif = False
        for canonical_motif, reference_region, chrom_start_0based_end, variant_type in zip(
                canonical_motifs, reference_regions, chroms_start_0based_ends, variant_types):
            chrom, start_0based, end = chrom_start_0based_end
            matched_known_disease_associated_locus |= bool(get_overlap(
                known_disease_associated_loci_interval_tree,  chrom, start_0based, end, canonical_motif=canonical_motif))
            matched_known_disease_associated_motif |= canonical_motif in known_disease_associated_motifs

        input_variant_catalog_record["KnownDiseaseAssociatedLocus"] = matched_known_disease_associated_locus
        input_variant_catalog_record["KnownDiseaseAssociatedMotif"] = matched_known_disease_associated_motif

        if args.only_known_disease_associated_loci and not matched_known_disease_associated_locus:
            continue
        if args.exclude_known_disease_associated_loci and matched_known_disease_associated_locus:
            continue
        if args.only_known_disease_associated_motifs and not matched_known_disease_associated_motif:
            continue

        # filter by gene region
        if len(chroms_start_0based_ends) == 1:
            spanning_interval_chrom, spanning_interval_start0_based, spanning_interval_end = chroms_start_0based_ends[0]
        else:
            spanning_interval_chrom = spanning_interval_start0_based = spanning_interval_end = None
            for chrom, start_0based, end in chroms_start_0based_ends:
                spanning_interval_chrom = chrom
                spanning_interval_start0_based = min(start_0based, spanning_interval_start0_based) if spanning_interval_start0_based is not None else start_0based
                spanning_interval_end = max(end, spanning_interval_end) if spanning_interval_end is not None else end

        if not args.gene_models_source:
            args.gene_models_source = os.path.basename(args.genes_gtf).split(".")[0].title()
        (
            input_variant_catalog_record[f"{args.gene_models_source}GeneRegion"],
            input_variant_catalog_record[f"{args.gene_models_source}GeneName"],
            input_variant_catalog_record[f"{args.gene_models_source}GeneId"],
            input_variant_catalog_record[f"{args.gene_models_source}TranscriptId"],
        ) = compute_genomic_region_of_interval(
            spanning_interval_chrom,
            spanning_interval_start0_based + 1,
            spanning_interval_end,
            args.genes_gtf,
            verbose=args.verbose)

        if args.region_type and input_variant_catalog_record[f"{args.gene_models_source}GeneRegion"] not in args.region_type:
            continue
        if args.exclude_region_type and input_variant_catalog_record[f"{args.gene_models_source}GeneRegion"] in args.exclude_region_type:
            continue
        if args.gene_name and input_variant_catalog_record[f"{args.gene_models_source}GeneName"] not in args.gene_name:
            continue
        if args.exclude_gene_name and input_variant_catalog_record[f"{args.gene_models_source}GeneName"] in args.exclude_gene_name:
            continue
        if args.gene_id and input_variant_catalog_record[f"{args.gene_models_source}GeneId"] not in args.gene_id:
            continue
        if args.exclude_gene_id and input_variant_catalog_record[f"{args.gene_models_source}GeneId"] in args.exclude_gene_id:
            continue

        # add this record to the output variant catalog
        output_records.append(input_variant_catalog_record)

    with open(args.output_path, "wt") as f:
        json.dump(output_records, f, indent=4)

    print(f"Wrote {len(output_records):,d} records to {args.output_path}")


if __name__ == "__main__":
    main()
