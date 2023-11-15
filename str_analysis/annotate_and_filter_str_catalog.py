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
from str_analysis.utils.eh_catalog_utils import parse_motifs_from_locus_structure, convert_json_records_to_bed_format_tuples
from str_analysis.utils.file_utils import open_file, file_exists, download_local_copy
from str_analysis.utils.gtf_utils import compute_genomic_region_of_interval
from str_analysis.utils.misc_utils import parse_interval

VALID_GENE_REGIONS = {"CDS", "UTR", "5UTR", "3UTR", "promoter", "exon", "intron", "intergenic"}

# The GEM mappability tracks  are based on 36bp, 50bp, 75bp, or 100bp kmers.
# These tracks can be viewed @ https://tgg-viewer.broadinstitute.org
MAPPABILITY_TRACK_KMER_SIZE = 100
MAPPABILITY_TRACK_BIGWIG_URL = f"gs://tgg-viewer/ref/GRCh38/mappability/" \
                               f"GRCh38_no_alt_analysis_set_GCA_000001405.15-k{MAPPABILITY_TRACK_KMER_SIZE}_m2.bw"
FLANK_MAPPABILITY_WINDOW_SIZE = 100

ACGT_REGEX = re.compile("^[ACGT]+$", re.IGNORECASE)

def parse_args():
    parser = argparse.ArgumentParser(description="Annotate and filter an STR variant catalog.")
    parser.add_argument("-r", "--reference-fasta", required=True, help="Reference fasta file path or url.")
    parser.add_argument("-o", "--output-path", help="Output path")
    parser.add_argument("--output-tsv", action="store_true", help="Output the catalog as a tab-separated table "
                        "in addition to the standard ExpansionHunter JSON file format")
    parser.add_argument("--output-bed", action="store_true", help="Output the catalog as a BED file "
                        "in addition to the standard ExpansionHunter JSON file format")
    parser.add_argument("--output-stats", action="store_true", help="Output detailed stats about loci in the catalog")
    parser.add_argument("--verbose", action="store_true", help="Print verbose output.")

    annotations_group = parser.add_argument_group("annotations")
    annotations_group.add_argument("--known-disease-associated-loci",
        help="ExpansionHunter catalog .json file with all known disease-associated loci",
        default="~/code/str-analysis/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json")
    #default="https://raw.githubusercontent.com/broadinstitute/str-analysis/main/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json")
    annotations_group.add_argument("--genes-gtf", help="Gene models gtf file path or url.",
        default="~/code/str-truth-set/ref/other/MANE.v1.0.ensembl_genomic.sorted.gtf.gz")
    annotations_group.add_argument("--gene-models-source", help="Source of the genes-gtf file. If not specified, "
        "it will be computed based on the filename", choices=["gencode", "mane", "refseq"])
    annotations_group.add_argument("--skip-gene-annotations", action="store_true", help="Don't addd gene annotations to "
        "the output catalog")
    annotations_group.add_argument("--skip-disease-loci-annotations", action="store_true", help="Don't annotate known "
        "disease associated loci in the output catalog")

    filter_group = parser.add_argument_group("filters")
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

    filter_group.add_argument("--max-interruptions", type=int, help="Maximum number of interruptions allowed in the "
                                                                    "reference repeat sequence")
    filter_group.add_argument("--min-base-purity", type=float, help="Minimum purity of the reference repeat sequence, "
                              "computed as the fraction of bases that deviate from pure repeats of the given motif")
    filter_group.add_argument("--min-repeat-purity", type=float, help="Minimum purity of the reference repeat sequence, "
                              "computed as the fraction of repeats that deviate from pure repeats of the given motif")
    filter_group.add_argument("--min-mappability", type=float, help="Minimum overall mappability of the repeat region")
    filter_group.add_argument("--discard-overlapping-intervals-with-similar-motifs", action="store_true",
                              help="Discard intervals if they overlap another interval with a similar motif")
    filter_group.add_argument("--discard-loci-with-non-acgt-bases", action="store_true", help="Discard loci that have "
                              "non-ACGT bases within their motif and/or their reference repeat sequence")
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
                unmodified_chrom = fields[0]
                chrom = unmodified_chrom.replace("chr", "")
                start_0based = int(fields[1])
                end_1based = int(fields[2])
                motif = fields[3].strip("()*+").upper()
                if not ACGT_REGEX.match(motif):
                    print(f"WARNING: skipping line with invalid motif: {line.strip()}")
                    continue

                record = {
                    "LocusId": f"{chrom}-{start_0based + 1}-{end_1based}-{motif}",
                    "ReferenceRegion": f"{unmodified_chrom}:{start_0based}-{end_1based}",
                    "LocusStructure": f"({motif})*",
                    "VariantType": "Repeat",
                }
                yield record


def output_tsv(output_path, output_records):
    output_tsv_rows = []
    for record in output_records:
        if not isinstance(record["ReferenceRegion"], list):
            output_tsv_rows.append(record)
            continue

        fields_that_are_lists = [
            "ReferenceRegion", "VariantType", "InterruptionBaseCount", "FractionPureBases", "FractionPureRepeats",
        ]
        for i in range(len(record["ReferenceRegion"])):
            output_row = dict(record)
            for field in fields_that_are_lists:
                output_row[field] = record[field][i]
            output_tsv_rows.append(output_row)

    pd.DataFrame(output_records).to_csv(output_path, sep="\t", index=False, header=True)


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

    ref_fasta = pysam.FastaFile(args.reference_fasta)

    # download and open the mappability bigWig file
    mappability_bigwig_path = download_local_copy(MAPPABILITY_TRACK_BIGWIG_URL)
    mappability_bigwig = pyBigWig.open(mappability_bigwig_path)

    # parse known disease-associated loci
    known_disease_associated_loci_interval_tree, known_disease_associated_motifs = parse_known_disease_associated_loci(
        args, parser)

    if not args.output_path:
        filename_prefix = re.sub("(.json|.bed)(.b?gz)?$", "", os.path.basename(args.variant_catalog_json_or_bed))
        args.output_path = f"{filename_prefix}.annotated_and_filtered.json.gz"

    output_records = []
    if args.verbose:
        print(f"Parsing {args.variant_catalog_json_or_bed}")

    input_variant_catalog_iterator = get_variant_catalog_iterator(args.variant_catalog_json_or_bed)
    if args.verbose:
        input_variant_catalog_iterator = tqdm(input_variant_catalog_iterator, unit=" variant catalog records", unit_scale=True)

    interval_trees = collections.defaultdict(IntervalTree)  # used to check for overlap between records in the catalog
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

        motifs = parse_motifs_from_locus_structure(input_variant_catalog_record["LocusStructure"])
        if isinstance(input_variant_catalog_record["ReferenceRegion"], list):
            reference_regions = input_variant_catalog_record["ReferenceRegion"]
            variant_types = input_variant_catalog_record["VariantType"]
        else:
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

        # annotate repeat purity in the reference genome
        interruption_base_count = []
        fraction_pure_bases = []
        fraction_pure_repeats = []
        overlaps_other_interval = False
        overlaps_other_interval_with_similar_motif = False
        has_non_acgt_bases = False
        for (chrom, start_0based, end), motif in zip(chroms_start_0based_ends, motifs):
            trimmed_end = end - (end - start_0based) % len(motif)
            ref_fasta_sequence = ref_fasta.fetch(chrom, start_0based, trimmed_end)  # fetch uses 0-based coords
            ref_fasta_sequence = ref_fasta_sequence.upper()

            if args.discard_loci_with_non_acgt_bases and not ACGT_REGEX.match(motif + ref_fasta_sequence):
                has_non_acgt_bases = True
                break

            pure_sequence = motif.upper() * int(len(ref_fasta_sequence)/len(motif))
            if len(pure_sequence) != len(ref_fasta_sequence):
                raise ValueError(f"Invalid input locus spec: {pformat(input_variant_catalog_record)}. Reference "
                                 f"interval is not a multiple of the motif size")

            num_matching_bases = sum(1 for nuc1, nuc2 in zip(ref_fasta_sequence, pure_sequence) if nuc1 == nuc2)

            interruption_base_count.append( len(ref_fasta_sequence) - num_matching_bases )
            fraction_pure_bases.append( round(num_matching_bases / len(ref_fasta_sequence), 2) )
            fraction_pure_repeats.append( round(ref_fasta_sequence.count(motif) / int(len(ref_fasta_sequence) / len(motif)), 2) )

            # check for overlap
            for overlapping_interval in interval_trees[chrom].overlap(start_0based, end):
                overlapping_interval_motif_size = overlapping_interval.data
                larger_motif_size = max(len(motif), overlapping_interval_motif_size)
                if overlapping_interval.overlap_size(start_0based, end) >= 2*larger_motif_size:
                    overlaps_other_interval = True
                    overlaps_other_interval_with_similar_motif = len(motif) == overlapping_interval_motif_size
                    break

            interval_trees[chrom].add(Interval(start_0based, end, data=len(motif)))

        if args.discard_loci_with_non_acgt_bases and has_non_acgt_bases:
            continue

        input_variant_catalog_record["InterruptionBaseCount"] = interruption_base_count[0] if len(chroms_start_0based_ends) == 1 else interruption_base_count
        input_variant_catalog_record["FractionPureBases"] = fraction_pure_bases[0] if len(chroms_start_0based_ends) == 1 else fraction_pure_bases
        input_variant_catalog_record["FractionPureRepeats"] = fraction_pure_repeats[0] if len(chroms_start_0based_ends) == 1 else fraction_pure_repeats
        input_variant_catalog_record["OverlapsOtherInterval"] = overlaps_other_interval

        if args.discard_overlapping_intervals_with_similar_motifs and overlaps_other_interval_with_similar_motif:
            continue

        if args.max_interruptions is not None and input_variant_catalog_record["InterruptionBaseCount"] > args.max_interruptions:
            continue
        if args.min_base_purity is not None and input_variant_catalog_record["FractionPureBases"] < args.min_base_purity:
            continue
        if args.min_repeat_purity is not None and input_variant_catalog_record["FractionPureRepeats"] < args.min_repeat_purity:
            continue

        # compute mappability of left and right flanking sequence
        chrom, left_flank_end, _ = chroms_start_0based_ends[0]
        _, right_flank_start, _ = chroms_start_0based_ends[-1]
        (mappability_left_flank, ) = mappability_bigwig.stats(
            chrom,
            left_flank_end - FLANK_MAPPABILITY_WINDOW_SIZE - MAPPABILITY_TRACK_KMER_SIZE,
            left_flank_end - MAPPABILITY_TRACK_KMER_SIZE)  # subtract the kmer size so that mappability scores are
                                                           # retrieved from an interval that doesn't overlap the repeat
                                                           # region and so only measure mappability of the flank itself
        (mappability_right_flank, ) = mappability_bigwig.stats(
            chrom,
            right_flank_start,
            right_flank_start + FLANK_MAPPABILITY_WINDOW_SIZE)

        (mappability_overall, ) = mappability_bigwig.stats(
            chrom,
            left_flank_end - FLANK_MAPPABILITY_WINDOW_SIZE,
            right_flank_start + FLANK_MAPPABILITY_WINDOW_SIZE)

        input_variant_catalog_record["LeftFlankMappability"] = round(mappability_left_flank, 2)
        input_variant_catalog_record["RightFlankMappability"] = round(mappability_right_flank, 2)
        input_variant_catalog_record["OverallMappability"] = round(mappability_overall, 2)

        if args.min_mappability is not None and input_variant_catalog_record["OverallMappability"] < args.min_mappability:
            continue

        # add this record to the output variant catalog
        output_records.append(input_variant_catalog_record)

    fopen = gzip.open if args.output_path.endswith(".gz") else open
    with fopen(args.output_path, "wt") as f:
        json.dump(output_records, f, indent=4)
    print(f"Wrote {len(output_records):,d} records to {args.output_path}")

    output_path_prefix = re.sub("(.json|.bed)(.b?gz)?$", "", args.output_path)
    if args.output_stats:
        if args.verbose: print("Calculating catalog stats..")
        compute_catalog_stats(args.output_path, output_records, verbose=args.verbose)

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
