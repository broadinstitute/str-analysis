"""This script takes a variant catalog either in .json or .bed format and filters it."""

import argparse
import collections
import gzip
import ijson
from intervaltree import IntervalTree, Interval
import simplejson as json
import os
import pandas as pd
from pprint import pformat
import pyBigWig
import pysam
import re
from tqdm import tqdm

from str_analysis.compute_catalog_stats import compute_catalog_stats
from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.eh_catalog_utils import (parse_motifs_from_locus_structure,
                                                 convert_json_records_to_bed_format_tuples, get_variant_catalog_iterator)
from str_analysis.utils.file_utils import open_file, file_exists, download_local_copy
from str_analysis.utils.find_repeat_unit import find_repeat_unit_without_allowing_interruptions
from str_analysis.utils.gtf_utils import compute_genomic_region_of_interval, GENE_MODELS
from str_analysis.utils.misc_utils import parse_interval

VALID_GENE_REGIONS = {"CDS", "UTR", "5UTR", "3UTR", "promoter", "exon", "intron", "intergenic"}

# The GEM mappability tracks  are based on 36bp, 50bp, 75bp, or 100bp kmers.
# These tracks can be viewed @ https://tgg-viewer.broadinstitute.org
MAPPABILITY_TRACK_KMER_SIZE = 36
MAPPABILITY_TRACK_BIGWIG_URL = f"gs://tgg-viewer/ref/GRCh38/mappability/" \
                               f"GRCh38_no_alt_analysis_set_GCA_000001405.15-k{MAPPABILITY_TRACK_KMER_SIZE}_m2.bw"
FLANK_MAPPABILITY_WINDOW_SIZE = 150

ACGT_REGEX = re.compile("^[ACGT]+$", re.IGNORECASE)
ACGTN_REGEX = re.compile("^[ACGTN]+$", re.IGNORECASE)

def parse_args():
    parser = argparse.ArgumentParser(description="Annotate and filter an STR variant catalog.")
    parser.add_argument("-r", "--reference-fasta", required=True, help="Reference fasta file path or url.")
    parser.add_argument("-o", "--output-path", help="Output json path. Any additional output paths such as tsv or bed "
                                                    "will be computed based on this path.")
    parser.add_argument("--output-tsv", action="store_true", help="Output the catalog as a tab-separated table "
                        "in addition to the standard ExpansionHunter JSON file format")
    parser.add_argument("--output-bed", action="store_true", help="Output the catalog as a BED file "
                        "in addition to the standard ExpansionHunter JSON file format")
    parser.add_argument("--output-stats", action="store_true", help="Output detailed stats about loci in the catalog")
    parser.add_argument("--verbose", action="store_true", help="Print verbose output.")
    parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")

    annotations_group = parser.add_argument_group("annotations")
    annotations_group.add_argument("--known-disease-associated-loci",
        help="ExpansionHunter catalog .json file with all known disease-associated loci",
        default="~/code/str-analysis/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json")
    #default="https://raw.githubusercontent.com/broadinstitute/str-analysis/main/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json")
    annotations_group.add_argument("--genes-gtf", help="Gene models gtf file path or url.")
    annotations_group.add_argument("--gene-models-source", action="append", help="Source of the genes-gtf file. "
        "Can be specified more than once to annotate with multiple gene models.", choices=["gencode", "mane", "refseq"])
    annotations_group.add_argument("--mappability-track-bigwig", default=MAPPABILITY_TRACK_BIGWIG_URL,
        help="Path or URL of the bigWig file containing mappability scores for each base in the reference genome")
    annotations_group.add_argument("--skip-gene-annotations", action="store_true", help="Don't addd gene annotations to "
        "the output catalog")
    annotations_group.add_argument("--skip-disease-loci-annotations", action="store_true", help="Don't annotate known "
        "disease associated loci in the output catalog")
    annotations_group.add_argument("--skip-mappability-annotations", action="store_true", help="Don't annotate known "
        "loci with mappability")

    annotations_group.add_argument("--add-gene-region-to-locus-id", action="store_true", help="Append the gene region "
        "to each locus id")
    annotations_group.add_argument("--add-canonical-motif-to-locus-id", action="store_true", help="Append the motif "
        "to each locus id")

    filter_group = parser.add_argument_group("filters")
    filter_group.add_argument("--min-motif-size", type=int, help="Minimum motif size to include in the output catalog")
    filter_group.add_argument("--max-motif-size", type=int, help="Maximum motif size to include in the output catalog")
    filter_group.add_argument("-ms", "--motif-size", type=int, action="append", help="Only include loci with these motif sizes")
    filter_group.add_argument("-m", "--motif", action="append", help="Only include loci whose canonical motif matched "
                                                                     "the canonical motif version of the given motif")

    filter_group.add_argument("--min-interval-size-bp", type=int, help="Filter out any loci with a reference interval "
                                                                       "that's smaller than this many base pairs")
    filter_group.add_argument("--max-interval-size-bp", type=int, help="Filter out any loci with a reference interval "
                                                                       "that's larger than this many base pairs")
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

    #filter_group.add_argument("--max-interruptions", type=int, help="Maximum number of interruptions allowed in the "
    #                                                                "reference repeat sequence")
    filter_group.add_argument("--min-reference-repeat-purity", type=float, help="Minimum purity of the "
                              "reference repeat sequence, computed as the fraction of bases that deviate from pure "
                              "repeats of the given motif")
    #filter_group.add_argument("--min-repeat-purity", type=float, help="Minimum purity of the reference repeat sequence, "
    #                          "computed as the fraction of repeats that deviate from pure repeats of the given motif")
    filter_group.add_argument("--min-mappability", type=float, help="Minimum overall mappability of the repeat region")
    filter_group.add_argument("--discard-overlapping-intervals-with-similar-motifs", action="store_true",
                              help="Discard intervals if they overlap another interval with a similar motif")
    filter_group.add_argument("--discard-loci-with-non-ACGT-bases-in-motif", action="store_true", help="Discard loci "
                              "that include non-ACGT bases within their motif")
    filter_group.add_argument("--discard-loci-with-non-ACGTN-bases-in-motif", action="store_true", help="Discard loci "
                              "that include non-ACGTN bases within their motif")
    filter_group.add_argument("--discard-loci-with-non-ACGT-bases-in-reference", action="store_true", help="Discard "
                              "loci if their reference sequence contains non-ACGT bases")
    filter_group.add_argument("--discard-loci-with-non-ACGTN-bases-in-reference", action="store_true", help="Discard "
                              "loci if their reference sequence contains non-ACGTN bases")

    modification_group = parser.add_argument_group("modifications")
    modification_group.add_argument("--dont-simplify-motifs", action="store_true", help="By default, repeat motifs within the "
                              "LocusStructure will be simplified to their most compact form. For example, a "
                              "LocusStructure defined as '(TTT)*CGAG(CAGCAG)*' will be replaced with '(T)*CGAG(CAG)*'. "
                              "This option disables this simplification step.")
    modification_group.add_argument("--trim-loci", action="store_true", help="Trim locus end coordinate so that the"
                                                                             "interval width is an exact multiple of the "
                                                                             "motif size")
    parser.add_argument("variant_catalog_json_or_bed", help="A catalog of repeats to annotate and filter, either "
                        "in JSON or BED format. For BED format, the chrom, start, and end should represent the repeat "
                        "interval in 0-based coordinates, and the name field (column #4) should be the repeat unit.")

    args = parser.parse_args()

    if not args.skip_gene_annotations and args.genes_gtf and not file_exists(os.path.expanduser(args.genes_gtf)):
        parser.error(f"File not found: {args.genes_gtf}")
    if not args.skip_disease_loci_annotations and not file_exists(os.path.expanduser(args.known_disease_associated_loci)):
        parser.error(f"File not found: {args.known_disease_associated_loci}")
    if not file_exists(os.path.expanduser(args.variant_catalog_json_or_bed)):
        parser.error(f"File not found: {args.variant_catalog_json_or_bed}")


    return args, parser


KNOWN_DISEASE_ASSOCIATED_LOCI_COLUMNS = [
    'LocusId', 'LocusStructure', 'RepeatUnit', 'MainReferenceRegion', 'Gene', 'GeneRegion', 'GeneId']


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



def output_tsv(output_path, output_records):
    output_tsv_rows = []
    for record in output_records:
        if not isinstance(record["ReferenceRegion"], list):
            output_tsv_rows.append(record)
            continue

        fields_that_are_lists = [
            "ReferenceRegion", "VariantType", "ReferenceRepeatPurity", #"FractionPureRepeats", "InterruptionBaseCount", 
        ]
        for i in range(len(record["ReferenceRegion"])):
            output_row = dict(record)
            for field in fields_that_are_lists:
                output_row[field] = record[field][i]
            output_tsv_rows.append(output_row)

    pd.DataFrame(output_records).to_csv(output_path, sep="\t", index=False, header=True)


def compute_sequence_purity_stats(nucleotide_sequence, motif):
    """Compute the number of interruptions, base-level purity, and repeat-level purity of a nucleotide sequence,
    compared to an equal-length pure repeat sequence with the given motif. If the nucleotide sequence length is not an
    exact multiple of the motif length, the sequence will be trimmed to the nearest multiple of the motif length before
    computing stats.

    Args:
        nucleotide_sequence (str): a nucleotide sequence
        motif (str): a repeat motif

    Return:
        tuple: (interruption_base_count, fraction_pure_bases, fraction_pure_repeats)
    """
    num_repeats = int(len(nucleotide_sequence)/len(motif))
    pure_sequence = motif.upper() * num_repeats
    nucleotide_sequence_trimmed = nucleotide_sequence[:len(pure_sequence)]
    num_matching_bases = sum(1 for nuc1, nuc2 in zip(nucleotide_sequence_trimmed, pure_sequence) if nuc1 == nuc2)

    if len(nucleotide_sequence) < len(motif):
        return 0, 0, 0

    interruption_base_count = len(nucleotide_sequence_trimmed) - num_matching_bases
    fraction_pure_bases = round(num_matching_bases / len(pure_sequence), 2)
    fraction_pure_repeats = round(nucleotide_sequence.count(motif) / num_repeats, 2)

    return interruption_base_count, fraction_pure_bases, fraction_pure_repeats


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

    # get fasta chrom sizes
    ref_fasta = pysam.FastaFile(args.reference_fasta)

    ref_fasta_chromosome_sizes = dict(zip(ref_fasta.references, ref_fasta.lengths))

    # download and open the mappability bigWig file
    if not args.skip_mappability_annotations:
        if args.mappability_track_bigwig and any(args.mappability_track_bigwig.startswith(prefix) for prefix in ("http", "gs://")):
            mappability_bigwig_path = download_local_copy(args.mappability_track_bigwig)
        else:
            mappability_bigwig_path = args.mappability_track_bigwig
        mappability_bigwig = pyBigWig.open(mappability_bigwig_path)

    # parse known disease-associated loci
    if not args.skip_disease_loci_annotations:
        known_disease_associated_loci_interval_tree, known_disease_associated_motifs = parse_known_disease_associated_loci(
            args, parser)
    else:
        known_disease_associated_motifs = set()
        known_disease_associated_loci_interval_tree = collections.defaultdict(IntervalTree)

    if not args.output_path:
        filename_prefix = re.sub("(.json|.bed)(.b?gz)?$", "", os.path.basename(args.variant_catalog_json_or_bed))
        args.output_path = f"{filename_prefix}.annotated_and_filtered.json.gz"

    output_records = []
    if args.verbose:
        print(f"Parsing {args.variant_catalog_json_or_bed}")

    input_variant_catalog_iterator = get_variant_catalog_iterator(args.variant_catalog_json_or_bed)
    if args.show_progress_bar:
        input_variant_catalog_iterator = tqdm(input_variant_catalog_iterator, unit=" variant catalog records", unit_scale=True)

    filter_counters = collections.Counter()
    modification_counters = collections.Counter()
    interval_trees = collections.defaultdict(IntervalTree)  # used to check for overlap between records in the catalog
    warning_counter = 0
    for i, variant_catalog_record in enumerate(input_variant_catalog_iterator):
        filter_counters["total"] += 1

        # validate variant_catalog_record
        error = None
        for key in ["LocusId", "ReferenceRegion", "LocusStructure", "VariantType"]:
            if key not in variant_catalog_record:
                print(f"ERROR: {key} not found in variant catalog record #{i+1}: {variant_catalog_record}. Skipping...")
                error = f"missing {key} key"
                break
        if error:
            filter_counters[f"row {error}"] += 1
            continue

        # parse variant_catalog_record based on whether it has adjacent repeats
        locus_id = variant_catalog_record["LocusId"]
        if args.locus_id and locus_id not in args.locus_id:
            filter_counters[f"row LocusId != {args.locus_id}"] += 1
            continue
        if args.exclude_locus_id and locus_id in args.exclude_locus_id:
            filter_counters[f"row LocusId == {args.exclude_locus_id}"] += 1
            continue

        motifs = parse_motifs_from_locus_structure(variant_catalog_record["LocusStructure"])
        if isinstance(variant_catalog_record["ReferenceRegion"], list):
            reference_regions = variant_catalog_record["ReferenceRegion"]
            variant_types = variant_catalog_record["VariantType"]
        else:
            reference_regions = [variant_catalog_record["ReferenceRegion"]]
            variant_types = [variant_catalog_record["VariantType"]]

        # validate variant_catalog_record with adjacent repeats
        if len(motifs) != len(reference_regions):
            print(f"ERROR: number of motifs != number of reference regions in variant catalog record #{i+1}: {variant_catalog_record}. Skipping...")
            filter_counters[f"row LocusStructure motif count != number of reference regions"] += 1
            continue
        if len(motifs) != len(variant_types):
            print(f"ERROR: number of motifs != number of variant types in variant catalog record #{i+1}: {variant_catalog_record}. Skipping...")
            filter_counters[f"row LocusStructure motif count != number of VariantTypes"] += 1
            continue

        # Check if the motif is actually composed of multiple repeats of a simpler motif.
        # For example, a "TTT" motif can be simplified to just a "T" homopolymer.
        if not args.dont_simplify_motifs:
            for i, (reference_region, motif) in enumerate(zip(reference_regions, motifs)):
                simplified_motif = None
    
                # the above elif clauses are an optimization to speed up the
                simplified_motif, num_repeats, _ = find_repeat_unit_without_allowing_interruptions(
                    motif, allow_partial_repeats=False)
                if num_repeats > 1:
                    variant_catalog_record["LocusStructure"] = variant_catalog_record["LocusStructure"].replace(
                        f"({motif})", f"({simplified_motif})")
    
                    motifs[i] = simplified_motif
                    counter_key = f"replaced {len(motif)}bp motif with a simplified {len(simplified_motif)}bp motif"
                    modification_counters[counter_key] += 1

                    #if args.verbose:
                    #    warning_counter += 1
                    #    print(f"WARNING #{warning_counter}: collapsing "
                    #          f"{locus_id} motif from {motif} to just {simplified_motif}:",
                    #          variant_catalog_record["LocusStructure"])

                if args.trim_loci:
                    # trim locus end coordinate so that the interval width is an exact multiple of the motif size
                    chrom, start_0based, end_1based = parse_interval(reference_region)
                    motif_size = len(simplified_motif)
                    interval_width = end_1based - start_0based
                    if interval_width % motif_size != 0:
                        new_end_1based = end_1based - interval_width % motif_size
                        if new_end_1based == start_0based:
                            new_end_1based += motif_size
                        variant_catalog_record["ReferenceRegion"] = f"{chrom}:{start_0based}-{new_end_1based}"
                        reference_regions[i] = variant_catalog_record["ReferenceRegion"]
                        counter_key = f"trimmed locus end coordinate to an exact multiple of the motif size"
                        modification_counters[counter_key] += 1

        # parse intervals
        chroms_start_0based_ends = [parse_interval(reference_region) for reference_region in reference_regions]
        canonical_motifs = [compute_canonical_motif(motif) for motif in motifs]

        # apply motif size and motif filters
        if args.min_motif_size and all(len(motif) < args.min_motif_size for motif in motifs):
            filter_counters[f"motif size(s) < {args.min_motif_size}"] += 1
            continue
        if args.max_motif_size and all(len(motif) > args.max_motif_size for motif in motifs):
            filter_counters[f"motif size(s) > {args.max_motif_size}"] += 1
            continue
        if args.motif_size and all(len(motif) not in args.motif_size for motif in motifs):
            filter_counters[f"motif size(s) not in {args.motif_size}"] += 1
            continue
        if args.motif:
            args_canonical_motif_set = {compute_canonical_motif(motif) for motif in args.motif}
            if all(canonical_motif not in args_canonical_motif_set for canonical_motif in canonical_motifs):
                filter_counters[f"row canonical motif(s) not among {args_canonical_motif_set}"] += 1
                continue

        if len(canonical_motifs) == 1:
            variant_catalog_record["CanonicalMotif"] = canonical_motifs[0]
            #variant_catalog_record["MotifSize"] = len(canonical_motifs[0])

        # apply interval size filters
        if args.min_interval_size_bp and all(end - start_0based < args.min_interval_size_bp for chrom, start_0based, end in chroms_start_0based_ends):
            filter_counters[f"row interval size < {args.min_interval_size_bp}bp"] += 1
            continue

        # parse input variant catalog record and check for overlap with known disease-associated loci
        matched_known_disease_associated_locus = False
        matched_known_disease_associated_motif = False
        for canonical_motif, reference_region, chrom_start_0based_end, variant_type in zip(
                canonical_motifs, reference_regions, chroms_start_0based_ends, variant_types):
            chrom, start_0based, end = chrom_start_0based_end
            matched_known_disease_associated_locus |= bool(get_overlap(
                known_disease_associated_loci_interval_tree,  chrom, start_0based, end, canonical_motif=canonical_motif))
            matched_known_disease_associated_motif |= canonical_motif in known_disease_associated_motifs

        variant_catalog_record["KnownDiseaseAssociatedLocus"] = matched_known_disease_associated_locus
        variant_catalog_record["KnownDiseaseAssociatedMotif"] = matched_known_disease_associated_motif

        if args.only_known_disease_associated_loci and not matched_known_disease_associated_locus:
            filter_counters[f"row isn't a known disease-associated locus"] += 1
            continue
        if args.exclude_known_disease_associated_loci and matched_known_disease_associated_locus:
            filter_counters[f"row is a known disease-associated locus"] += 1
            continue
        if args.only_known_disease_associated_motifs and not matched_known_disease_associated_motif:
            filter_counters[f"row motif doesn't match the motif of a known disease-associated locus"] += 1
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

        if not args.skip_gene_annotations:
            if args.genes_gtf:
                gene_models_source = os.path.basename(args.genes_gtf).split(".")[0].lower()
                args.gene_models_source = [gene_models_source]
                GENE_MODELS[gene_models_source] = args.genes_gtf

            for gene_models_source in args.gene_models_source:
                genes_gtf = GENE_MODELS[gene_models_source]

                gene_models_source = gene_models_source.title()

                (
                    variant_catalog_record[f"{gene_models_source}GeneRegion"],
                    variant_catalog_record[f"{gene_models_source}GeneName"],
                    variant_catalog_record[f"{gene_models_source}GeneId"],
                    variant_catalog_record[f"{gene_models_source}TranscriptId"],
                ) = compute_genomic_region_of_interval(
                    spanning_interval_chrom,
                    spanning_interval_start0_based + 1,
                    spanning_interval_end,
                    genes_gtf,
                    verbose=args.verbose,
                    show_progress_bar=args.show_progress_bar)

                if variant_catalog_record[f"{gene_models_source}GeneId"] == variant_catalog_record[f"{gene_models_source}GeneName"]:
                    del variant_catalog_record[f"{gene_models_source}GeneName"]

        if args.region_type and variant_catalog_record.get(f"{args.gene_models_source}GeneRegion") not in args.region_type:
            filter_counters[f"row {args.gene_models_source}GeneRegion isn't among {args.region_type}"] += 1
            continue
        if args.exclude_region_type and variant_catalog_record.get(f"{args.gene_models_source}GeneRegion") in args.exclude_region_type:
            filter_counters[f"row {args.gene_models_source}GeneRegion is one of {args.exclude_region_type}"] += 1
            continue
        if args.gene_name and variant_catalog_record.get(f"{args.gene_models_source}GeneName") not in args.gene_name:
            filter_counters[f"row {args.gene_models_source}GeneName isn't among {args.gene_name}"] += 1
            continue
        if args.exclude_gene_name and variant_catalog_record.get(f"{args.gene_models_source}GeneName") in args.exclude_gene_name:
            filter_counters[f"row {args.gene_models_source}GeneName is one of {args.exclude_gene_name}"] += 1
            continue
        if args.gene_id and variant_catalog_record.get(f"{args.gene_models_source}GeneId") not in args.gene_id:
            filter_counters[f"row {args.gene_models_source}GeneId isn't among {args.gene_id}"] += 1
            continue
        if args.exclude_gene_id and variant_catalog_record.get(f"{args.gene_models_source}GeneId") in args.exclude_gene_id:
            filter_counters[f"row {args.gene_models_source}GeneId is one of {args.exclude_gene_id}"] += 1
            continue

        # annotate repeat purity in the reference genome
        reference_repeats_list = []
        fraction_pure_bases_list = []
        overlaps_other_interval = False
        overlaps_other_interval_with_similar_motif = False
        has_invalid_bases = False
        for (chrom, start_0based, end), motif in zip(chroms_start_0based_ends, motifs):
            ref_fasta_sequence = ref_fasta.fetch(chrom, start_0based, end)  # fetch uses 0-based coords
            ref_fasta_sequence = ref_fasta_sequence.upper()

            if ((args.discard_loci_with_non_ACGT_bases_in_motif and not ACGT_REGEX.match(motif)) or
                (args.discard_loci_with_non_ACGTN_bases_in_motif and not ACGTN_REGEX.match(motif))):
                has_invalid_bases = True
                filter_counters["row motif has invalid bases"] += 1
                break
            if len(ref_fasta_sequence) > 0:
                if ((args.discard_loci_with_non_ACGT_bases_in_reference and not ACGT_REGEX.match(ref_fasta_sequence)) or
                    (args.discard_loci_with_non_ACGTN_bases_in_reference and not ACGTN_REGEX.match(ref_fasta_sequence))):
                    has_invalid_bases = True
                    filter_counters["row reference sequence has invalid bases"] += 1
                    break
            

            reference_repeats_list.append(len(ref_fasta_sequence) // len(motif))

            interrupted_bases_count, fraction_pure_bases, fraction_pure_repeats = compute_sequence_purity_stats(
                ref_fasta_sequence, motif)
            
            fraction_pure_bases_list.append( fraction_pure_bases )

            # check for overlap
            canonical_motif = compute_canonical_motif(motif, include_reverse_complement=False)
            for overlapping_interval in interval_trees[chrom].overlap(start_0based, end):
                overlapping_interval_motif = overlapping_interval.data
                larger_motif_size = max(len(canonical_motif), len(overlapping_interval_motif))
                if overlapping_interval.overlap_size(start_0based, end) >= 2*larger_motif_size:
                    overlaps_other_interval_with_similar_motif = canonical_motif == overlapping_interval_motif
                    #if args.verbose and overlaps_other_interval_with_similar_motif:
                    #    print(f"    {locus_id} overlaps another interval with a similar motif: "
                    #          f"{chrom}-{overlapping_interval.begin}-{overlapping_interval.end}-{overlapping_interval_motif} "
                    #          f"by 2 or more repeats of the larger motif. Source: {variant_catalog_record.get('Source')}")
                    break

            interval_trees[chrom].add(Interval(start_0based, end, data=canonical_motif))

        if has_invalid_bases:
            continue

        variant_catalog_record["ReferenceRepeats"] = reference_repeats_list[0] if len(chroms_start_0based_ends) == 1 else reference_repeats_list
        variant_catalog_record["ReferenceRepeatPurity"] = fraction_pure_bases_list[0] if len(chroms_start_0based_ends) == 1 else fraction_pure_bases_list

        if args.discard_overlapping_intervals_with_similar_motifs and overlaps_other_interval_with_similar_motif:
            filter_counters[f"overlapping interval"] += 1
            continue

        if args.min_reference_repeat_purity is not None and variant_catalog_record["ReferenceRepeatPurity"] < args.min_reference_repeat_purity:
            filter_counters[f"row ReferenceRepeatPurity < {args.min_reference_repeat_purity}"] += 1
            continue

        # compute mappability of left and right flanking sequence
        if not args.skip_mappability_annotations:
            chrom, left_flank_end, _ = chroms_start_0based_ends[0]
            _, _, right_flank_start = chroms_start_0based_ends[-1]

            mappability_left_flank = None
            left_flank_mappability_interval_start = max(1, left_flank_end - FLANK_MAPPABILITY_WINDOW_SIZE - MAPPABILITY_TRACK_KMER_SIZE)
            # subtract the kmer size so that mappability scores are retrieved from an interval that doesn't overlap the repeat
            # region and so only measure mappability of the flank itself
            left_flank_mappability_interval_end = max(2, left_flank_end - MAPPABILITY_TRACK_KMER_SIZE)
            try:

                (mappability_left_flank, ) = mappability_bigwig.stats(chrom,
                                                                      left_flank_mappability_interval_start,
                                                                      left_flank_mappability_interval_end,
                                                                      type="mean")
            except Exception as e:
                warning_counter += 1
                print(f"WARNING #{warning_counter}: Couldn't compute mappability of left flank interval: {chrom}:{left_flank_mappability_interval_start}-{left_flank_mappability_interval_end}: {e}")

            mappability_right_flank = None
            right_flank_mappability_interval_start = right_flank_start
            right_flank_mappability_interval_end = min(ref_fasta_chromosome_sizes[chrom], right_flank_start + FLANK_MAPPABILITY_WINDOW_SIZE)
            try:
                (mappability_right_flank, ) = mappability_bigwig.stats(chrom,
                                                                       right_flank_start,
                                                                       right_flank_mappability_interval_end,
                                                                       type="mean")
            except Exception as e:
                warning_counter += 1
                print(f"WARNING #{warning_counter}: Couldn't compute mappability of right flank interval: {chrom}:{right_flank_start}-{right_flank_mappability_interval_end}: {e}")

            mappability_overall = None
            try:
                (mappability_overall, ) = mappability_bigwig.stats(chrom,
                                                                   left_flank_mappability_interval_start,
                                                                   right_flank_mappability_interval_end,
                                                                   type="mean")
            except Exception as e:
                warning_counter += 1
                print(f"WARNING #{warning_counter}: Couldn't compute mappability of overall interval: {chrom}:{left_flank_mappability_interval_start}-{right_flank_mappability_interval_end}: {e}")

            variant_catalog_record[f"LeftFlankMappability"] = round(mappability_left_flank, 2)
            variant_catalog_record[f"FlanksAndLocusMappability"] = round(mappability_overall, 2)
            variant_catalog_record[f"RightFlankMappability"] = round(mappability_right_flank, 2)

            if args.min_mappability is not None and variant_catalog_record["FlanksAndLocusMappability"] < args.min_mappability:
                filter_counters[f"row FlanksAndLocusMappability < {args.min_mappability}"] += 1
                continue

        if args.add_gene_region_to_locus_id:
            gene_region = variant_catalog_record.get(f"{args.gene_models_source}GeneRegion")
            if gene_region:
                variant_catalog_record["LocusId"] += "-" + gene_region.replace(" ", "").replace("'", "")

        if args.add_canonical_motif_to_locus_id:
            variant_catalog_record["LocusId"] += "-" + "-".join(canonical_motifs)

        for key in list(variant_catalog_record.keys()):
            value = variant_catalog_record[key]
            if value is None or value is False:
                del variant_catalog_record[key]

        # add this record to the output variant catalog
        filter_counters["passed all filters"] += 1
        output_records.append(variant_catalog_record)

    fopen = gzip.open if args.output_path.endswith(".gz") else open
    with fopen(args.output_path, "wt") as f:
        json.dump(output_records, f, indent=4)
    print(f"Wrote {len(output_records):,d} records to {args.output_path}")

    output_path_prefix = re.sub("(.json|.bed)(.b?gz)?$", "", args.output_path)
    if args.output_stats:
        if args.verbose: print("Calculating catalog stats..")
        compute_catalog_stats(args.output_path, output_records, verbose=args.verbose)

    if args.verbose:
        print("\nFilter stats:")
        for key, value in sorted(filter_counters.items(), key=lambda x: -x[1]):
            if key == "total":
                print(f" {value:9,d} total input rows")
            else:
                print(f" {value:9,d} out of {filter_counters['total']:,d} ({value/filter_counters['total']:3.0%})  {key}")
        if sum(modification_counters.values()) > 0:
            print("\nModification stats:")
            for key, value in sorted(modification_counters.items(), key=lambda x: -x[1]):
                print(f" {value:9,d} out of {filter_counters['total']:,d} ({value/filter_counters['total']:3.0%})  {key}")


    if args.output_bed:
        output_path = f"{output_path_prefix}.bed"
        if args.verbose: print(f"Writing to {output_path}")
        with open(output_path, "wt") as output_bed:
            total = 0
            for bed_record in sorted(convert_json_records_to_bed_format_tuples(output_records)):
                total += 1
                output_bed.write("\t".join(map(str, bed_record)) + "\n")

        os.system(f"bgzip -f {output_path}")
        os.system(f"tabix -f -p bed {output_path}.gz")
        print(f"Done writing {total:,d} output records to {output_path}.gz")

    if args.output_tsv:
        if args.verbose: print(f"Writing to {output_path_prefix}.tsv.gz")
        output_tsv(f"{output_path_prefix}.tsv.gz", output_records)



if __name__ == "__main__":
    main()
