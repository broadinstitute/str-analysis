#!/usr/bin/env python3

"""
This script takes a VCF (either single-sample or multi-sample) and filters it to the subset of insertions and deletions
that represent tandem repeat (TR) expansions or contractions. It does this by checking each indel to see if the inserted
or deleted sequence consists entirely of repeats of some motif, and if yes, whether these repeats can be extended into
the flanking reference sequences immediately to the left or right of the variant. The output is a set of tandem repeat
loci (including their motifs and reference start and end coordinates) that can then be used for downstream analyses,
such as genotyping.

This script is the next iteration of the filter_vcf_to_STR_variants.py script. It implements multiple approaches to
detecting repeat sequences within each variant - first doing a fast, brute-force scan for perfect (or nearly perfect)
repeats. If no repeats are detected in this first step, it runs TandemRepeatFinder to discover more imperfect repeats
(particularly VNTRs). The script then merges overlapping tandem repeat alleles that have very similar motifs and writes
the results to output files. Unlike the original filter_vcf_to_STR_variants.py script, it separates tandem repeat
locus discovery (the 'catalog' subcommand) from genotyping (the optional downstream 'genotype' subcommand, which
genotypes the loci of a catalog from the indel genotypes in a single-sample VCF).

---

Pseudocode: 

 1. for each allele, check if the allele is a tandem repeat using a simple brute force scan for perfect (or nearly perfect) repeats
     - each allele will either be: 
           a. a tandem repeat  (add it to results)
           b. not a tandem repeat (go to #2) 
           c. a tandem repeat too long for the flanking sequence (increase flanking sequence size and redo #1)

  2. write each allele + flanking sequence to a FASTA file and run TRF on all of them using multiple threads, then parse the TRF output
     - each allele will either be: 
           a. tandem repeat (add it to results)
           b. not a tandem repeat 
           c. a tandem repeat too long for the flanking sequence (increase flanking sequence size and redo #2)

  3. write results to output files

  NOTE: Merging of overlapping TR loci is handled by the separate 'merge' subcommand, which should be run
  after the 'catalog' subcommand when combining results from multiple VCFs or when deduplication is desired.
"""

import argparse
import collections
import configargparse
import datetime
import functools
import gzip
import importlib.util
import itertools
import json
import math
import multiprocessing
import os
import shlex
import tempfile

import intervaltree
import pyfaidx
import pysam
import re
import shutil
import tqdm

from concurrent.futures import ThreadPoolExecutor
from pprint import pformat

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.fasta_utils import create_normalize_chrom_function, normalize_chromosome_name
from str_analysis.utils.find_repeat_unit import find_repeat_unit_allowing_interruptions
from str_analysis.utils.find_repeat_unit import find_repeat_unit_without_allowing_interruptions
from str_analysis.utils.find_repeat_unit import extend_repeat_into_sequence_allowing_interruptions
from str_analysis.utils.find_repeat_unit import extend_repeat_into_sequence_without_allowing_interruptions
from str_analysis.utils.file_utils import open_file, file_exists
from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.trf_runner import TRFRunner
from str_analysis.utils.find_motif_utils import compute_repeat_purity, compute_most_common_motif
from str_analysis.utils.find_motif_utils import compute_best_phase_repeat_purity, compute_sequence_periodicity
from str_analysis.utils.find_motif_utils import compute_partial_copy_purity, EDIT_DISTANCE_METRIC
from str_analysis.utils.find_motif_utils import split_sequence_into_motifs, format_motifs_as_sequence_string

DETECTION_MODE_PURE_REPEATS = "pure"
DETECTION_MODE_ALLOW_INTERRUPTIONS = "interrupted"
DETECTION_MODE_TRF = "trf"

# Methods used to parse an allele sequence into an ordered list of motifs for --add-motif-composition
MOTIF_DETECTION_METHOD_TRF = "trf"
MOTIF_DETECTION_METHOD_TRVIZ = "trviz"
MOTIF_DETECTION_METHOD_BASIC_SPLIT = "basic-split"

# Valid IUPAC DNA bases for allele validation
DNA_BASES = set("ACGTNRYSWKMBDHV")


class NonIUPACAlleleError(ValueError):
    """Raised when a REF or ALT allele contains characters that aren't IUPAC nucleotide codes.

    The usual cause is a symbolic allele such as <DEL> or a breakend such as N[chr2:123[. It is a ValueError so
    that callers catching ValueError from convert_variants_to_haplotype_sequence() still catch it.
    """


class RefAlleleMismatchError(ValueError):
    """Raised when a variant's REF allele doesn't match the reference sequence at its position.

    The usual cause is a VCF called against a different reference than the FASTA given with -R. It is a
    ValueError so that callers catching ValueError from convert_variants_to_haplotype_sequence() still catch it.
    """

CURRENT_TIMESTAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S.%f")
TRF_WORKING_DIR = f"trf_working_dir"


MAX_INDEL_SIZE = 100_000  # bp
MAX_FLANKING_SEQUENCE_SIZE = 1_000_000  # bp


FILTER_ALLELE_WITH_N_BASES = "contains Ns in the variant sequence"
FILTER_ALLELE_INDEL_WITHOUT_REPEATS = "INDEL without repeats"
FILTER_ALLELE_TOO_BIG = f"INDEL > {MAX_INDEL_SIZE:,d}bp"
FILTER_TR_ALLELE_NOT_ENOUGH_REPEATS = "contains < {:,d} full repeats"
FILTER_TR_ALLELE_NOT_ENOUGH_REPEATS_IN_REFERENCE = "contains < {:,d} full repeats in reference"
FILTER_TR_ALLELE_TOO_MANY_REPEATS = "contains > {:,d} repeats"
FILTER_TR_ALLELE_DOESNT_SPAN_ENOUGH_BASE_PAIRS = "spans < {:,d} bp"
FILTER_TR_ALLELE_SPANS_TOO_MANY_BASE_PAIRS = "spans > {:,d} bp"
FILTER_TR_ALLELE_PURITY_IS_TOO_LOW = "purity < {:.2f}"

TRF_MAX_REPEATS_IN_REFERENCE_THRESHOLD = 3_500  
TRF_MAX_SPAN_IN_REFERENCE_THRESHOLD = 10_000      # 10Kb

# Thresholds for deciding whether inserted bases belong to a locus tandem repeat, set from a held-out measurement on
# HG002. They are named rather than written inline because parse_args(), build_insertion_filter()'s fallbacks, and the
# tests all need the same values, and retuning one copy while missing another would leave the tests green against a
# threshold the CLI no longer uses.
DEFAULT_MIN_INSERTION_SIZE_TO_CHECK = 20
DEFAULT_MIN_INSERTION_PURITY = 0.9
DEFAULT_MIN_INSERTION_PERIODICITY = 0.55

# The edit distance between an allele and a pure repeat of the same length costs O(length^2), and the phase
# search computes it once per motif rotation, so a 19kb allele at a 170bp-motif locus took several seconds. The
# resulting purity is only an output annotation, so above this length it is left missing instead.
DEFAULT_MAX_ALLELE_LENGTH_FOR_EDIT_DISTANCE_PURITY = 10_000

#FILTER_TR_ALLELE_PARTIAL_REPEAT = "ends in partial repeat"

# Reasons why an inserted sequence at a tandem repeat locus was judged not to be part of the repeat
INSERTION_FILTER_REASON_CONTAINS_NS = "inserted sequence contains Ns"
INSERTION_FILTER_REASON_NOT_REPEAT_LIKE = "inserted sequence is not a tandem repeat"

# Reasons a whole locus is set to no call before its haplotypes are built. Each one describes a locus the VCF
# cannot answer for, as opposed to one the VCF says is homozygous reference.
# Each value is also the counter key used to summarize a run, so it stays free of per-locus detail: the
# specific bases or coordinates go in the parenthesized suffix that describe_no_call() appends.
NO_CALL_REASON_HET_AT_HAPLOID_LOCUS = "HET variant at a haploid locus"
NO_CALL_REASON_HAPLOID_GENOTYPE_AT_DIPLOID_LOCUS = "haploid genotype at a diploid locus"
NO_CALL_REASON_AMBIGUOUS_PHASING = "contains more than one HET variant with unclear phasing between them"
NO_CALL_REASON_NON_REPEAT_INSERTION = "inserted sequence is not sufficiently repetitive"
NO_CALL_REASON_HAPLOTYPE_BUILD_ERROR = "the haplotype sequence couldn't be built"
NO_CALL_REASON_NON_IUPAC_ALLELE = "a REF or ALT allele contains characters that aren't IUPAC nucleotide codes"
NO_CALL_REASON_REF_ALLELE_MISMATCH = "a REF allele doesn't match the reference FASTA"
NO_CALL_REASON_CONTIG_NOT_IN_REFERENCE = "contig not found in the reference FASTA"
NO_CALL_REASON_LOCUS_PAST_CONTIG_END = "locus extends past the end of the contig"
NO_CALL_REASON_MISSING_GENOTYPE = "genotype includes uncalled allele"
NO_CALL_REASON_CHROMOSOME_ABSENT = "sample lacks this chromosome"

# Haplotype build failures caused by the catalog or reference rather than by the sample's variants. They are left
# out of the count of alleles whose variants couldn't be applied (GenotypedTandemRepeat.num_alleles_with_build_errors).
REFERENCE_BUILD_ERROR_REASONS = (NO_CALL_REASON_CONTIG_NOT_IN_REFERENCE, NO_CALL_REASON_LOCUS_PAST_CONTIG_END)

# Pseudoautosomal regions (PARs) as 0-based half-open intervals, keyed by the length of chrX in the reference,
# which is what identifies the genome build. Every GRCh38 flavor (UCSC hg38, the NCBI analysis sets, the Broad
# Homo_sapiens_assembly38) shares the same primary chromosome lengths and so the same PARs, and likewise for
# GRCh37 and hg19. Outside the PARs, chrX and chrY are haploid in a male, so on a chromosome detected as
# haploid (see detect_sex_chromosome_ploidy) a genotype that names only one allele (1, 1|., .|1, ./1) is a
# valid haploid call rather than an uncalled haplotype.
PAR_REGIONS_BY_CHRX_LENGTH = {
    156_040_895: ("GRCh38", {
        "X": [(10_000, 2_781_479), (155_701_382, 156_030_895)],
        "Y": [(10_000, 2_781_479), (56_887_902, 57_217_415)],
    }),
    155_270_560: ("GRCh37", {
        "X": [(60_000, 2_699_520), (154_931_043, 155_260_560)],
        "Y": [(10_000, 2_649_520), (59_034_049, 59_363_566)],
    }),
    154_259_566: ("T2T-CHM13", {
        "X": [(0, 2_394_410), (153_925_834, 154_259_566)],
        "Y": [(0, 2_458_320), (62_122_809, 62_460_029)],
    }),
}

# Thresholds for detect_sex_chromosome_ploidy's first two chrX tests. Among non-PAR chrX records that call
# both alleles and carry an alt, 4-19% are heterozygous in males whose assembly haplotypes aren't split by
# parent (HGSVC h1/h2), versus 50-71% in females, so 35% sits in the middle of the gap. Each test is a fraction
# of the records it looks at, so each needs enough of them for the fraction to mean something: a whole-genome
# VCF has 100,000+ non-PAR chrX records, so 1,000 only guards against deciding from a VCF that barely covers
# chrX, where a handful of ".|1" records would otherwise mark the whole sample as haploid.
MAX_CHRX_HET_FRACTION_FOR_HAPLOID_X = 0.35
MIN_CHRX_RECORDS_FOR_PLOIDY_DETECTION = 1000

# Called non-PAR chrY records needed to count chrY as present. In a survey of 331 DipCall VCFs every sample
# had either none or at least 2,161.
MIN_CHRY_RECORDS_FOR_HAPLOID_Y = 1000


def describe_no_call(no_call_reason, detail=None):
    """Combine a NO_CALL_REASON_* value with the per-locus detail behind it, for the output columns.

    Args:
        no_call_reason (str): one of the NO_CALL_REASON_* values
        detail (str): what specifically went wrong at this locus, or None

    Returns:
        str: the reason, with the detail in parentheses when there is one
    """
    return f"{no_call_reason} ({detail})" if detail else no_call_reason

FILTER_TR_ALLELE_REPEAT_UNIT_TOO_SHORT = "repeat unit < {:,d} bp"
FILTER_TR_ALLELE_REPEAT_UNIT_TOO_LONG = "repeat unit > {:,d} bp"

# Output columns for genotype subcommand TSV output
GENOTYPE_TSV_OUTPUT_COLUMNS = [
    "Chrom",
    "Start0Based",
    "End",
    "Locus",           # Format: "chr1:1000-1050" (0-based start)
    "LocusId",         # Format: "chr1-1000-1050-CAG" (0-based start)
    "Motif",
    "CanonicalMotif",
    "MotifSize",
    "NumRepeatsInReference",
    "NumRepeatsShortAllele",
    "NumRepeatsLongAllele",
    "RepeatSizeShortAlleleBp",
    "RepeatSizeLongAlleleBp",
    "Zygosity",        # HOM, HET, or HEMI
    "IsPureRepeat",
    "RepeatPurity",
    "RepeatPurityShortAllele",
    "RepeatPurityLongAllele",
    "RepeatPurityViaEditDistance",
    "RepeatPurityViaEditDistanceShortAllele",
    "RepeatPurityViaEditDistanceLongAllele",
    "Allele1Sequence",
    "Allele2Sequence",
    "Allele1MotifSequence",
    "Allele1SequenceMotifSplittingMethod",
    "Allele2MotifSequence",
    "Allele2SequenceMotifSplittingMethod",
    "NumOverlappingVariants",
    "VariantPositions",
    "NoCallReason",
]


def parse_args():
    """Parse command-line arguments."""

    p = configargparse.ArgumentParser(formatter_class=configargparse.ArgumentDefaultsHelpFormatter)
    subparsers = p.add_subparsers(dest="subcommand", required=True, title="Subcommand")
    catalog_p = subparsers.add_parser("catalog", help="step1: Discover tandem repeat loci in a VCF file.",
                                      formatter_class=configargparse.ArgumentDefaultsHelpFormatter)
    catalog_p.add_argument("-c", "--config-file", is_config_file=True, help="Optional config file path")

    # catalog subcommand
    catalog_p.add_argument("-R", "--reference-fasta-path", help="Reference genome fasta path.", required=True)
    catalog_p.add_argument("--dont-allow-interruptions", action="store_true", help="Only detect perfect repeats. This implicitly "
                   "also enables --dont-run-trf since detection of pure repeats does not require running TandemRepeatFinder (TRF).")
    catalog_p.add_argument("--dont-run-trf", action="store_true", help="Don't use TandemRepeatFinder (TRF) to help detect imperfect TRs. Instead, only "
                   "use the simpler algorithm that allows one position within the repeat unit to vary across repeats.")
    catalog_p.add_argument("--trf-executable-path", help="Path to the TandemRepeatFinder (TRF) executable. This is "
                   "required unless --dont-run-trf is specified.")
    catalog_p.add_argument("-t", "--trf-threads", default=max(1, multiprocessing.cpu_count() - 2), type=int, help="Number of TandemRepeatFinder (TRF) "
                   "instances to run in parallel.")
    catalog_p.add_argument("--allow-multiple-trf-results-per-locus", action="store_true",
                           help="At some loci, TRF returns multiple valid results that have different locus boundaries and motif sizes, "
                           "with no obvious best choice. By default, the result with the shortest motif (>= 3bp) is selected. "
                           "This option changes the behavior to return separate entries for all results that pass specified thresholds.")
    catalog_p.add_argument("--min-indel-size-to-run-trf", default=7, type=int, help="Only run TandemRepeatFinder (TRF) "
        "on insertions and deletions that are at least this many base pairs.")

    catalog_p.add_argument("--trf-min-repeats-in-reference", type=int, default=2, help="For TRF results, require a locus to span "
        "at least this many repeats in the reference genome. This helps filter out noisy, low-quality TRF results.")
    catalog_p.add_argument("--trf-min-purity", type=float, default=0.2, help="For TRF results, filter out locus definitions where the "
        "repeat purity is below this threshold (defined as the fraction of bases that correspond to perfect repeats "
        "of the locus motif).")
    catalog_p.add_argument("--trf-mismatch-penalty", default=7, type=int, help="TandemRepeatFinder (TRF) mismatch penalty.")
    catalog_p.add_argument("--trf-indel-penalty", default=7, type=int, help="TandemRepeatFinder (TRF) indel penalty.")
    catalog_p.add_argument("--trf-min-score", default=20, type=int, help="TandemRepeatFinder (TRF) minimum alignment score.")
    catalog_p.add_argument("--trf-working-dir", default=TRF_WORKING_DIR, help="Directory to store intermediate files "
        "for TandemRepeatFinder (TRF).")
    catalog_p.add_argument("--min-tandem-repeat-length", type=int, default=9, help="Only detect tandem repeat variants "
        "that are at least this long (in base pairs). This threshold will be applied to the total repeat sequence "
        "including any repeats in the flanking sequence to the left and right of the variant in addition to the "
        "inserted or deleted bases.")
    catalog_p.add_argument("--min-repeats", type=int, default=3, help="Only detect tandem repeat loci that consist of at least this many repeats. "
                   "This threshold will be applied to the total repeat sequence including any repeats in the flanking sequence to "
                   "the left and right of the variant in addition to the inserted or deleted bases")
    catalog_p.add_argument("--min-repeat-unit-length", type=int, default=1, help="Minimum repeat unit length in base pairs.")
    catalog_p.add_argument("--max-repeat-unit-length", type=int, default=10**9, help="Max repeat unit length in base pairs.")
    catalog_p.add_argument("--show-progress-bar", help="Show a progress bar in the terminal when processing variants.",
                   action="store_true")
    catalog_p.add_argument("-v", "--verbose", help="Print detailed logs.", action="store_true")
    catalog_p.add_argument( "--debug", help="Print any debugging info and don't delete intermediate files.", action="store_true")

    catalog_p.add_argument("--offset", default=0, type=int, help="Skip the first N variants in the VCF file. This is useful for testing ")
    catalog_p.add_argument("-n", type=int, help="Only process N rows from the VCF (after applying --offset). Useful for testing.")

    catalog_p.add_argument("-o", "--output-prefix", help="Output file prefix. If not specified, it will be computed based on "
                   "the input vcf filename")

    catalog_p.add_argument("--write-detailed-bed", help="Output a second BED file (in addition to the main output BED file of all TR loci) "
                           "where the name field (ie. column 4) contains additional info besides the repeat unit.", action="store_true")
    catalog_p.add_argument("--write-vcf", help="Output a VCF file with all variants that were found to be TRs.", action="store_true")
    catalog_p.add_argument("--write-filtered-variants-to-vcf", help="Output a VCF file with indel variants that were not identified as TRs. "
                   "The FILTER column will contain the reason each variant was filtered.", action="store_true")
    catalog_p.add_argument("--write-fasta", help="Output a FASTA file containing all TR alleles", action="store_true")
    catalog_p.add_argument("--write-tsv", help="Output a TSV file containing all TR alleles", action="store_true")
    catalog_p.add_argument("-ik", "--copy-info-field-keys-to-tsv", help="Copy the values of these INFO field keys from the input "
                                                                 "VCF to the output TSV files.", action="append")
    catalog_p.add_argument("-L", "--interval", help="Only process variants in this genomic interval (format: chrom:start-end) or BED file",
                   action="append")
    
    catalog_p.add_argument("input_vcf_path", help="Input single-sample VCF file path. This script was designed and tested on VCFs produced by DipCall, "
                   "but should work with any single-sample VCF.")

    # merge subcommand
    merge_p = subparsers.add_parser("merge", help="step2: Optionally merge catalogs produced in step1 from two or more VCF files.",
                                    formatter_class=configargparse.ArgumentDefaultsHelpFormatter)
    merge_p.add_argument("-c", "--config-file", is_config_file=True, help="Optional config file path")

    merge_p.add_argument("-R", "--reference-fasta-path", help="Reference genome fasta path.", required=True)
    merge_p.add_argument("--write-detailed-bed", help="Output a second BED file (in addition to the main output BED file of all TR loci after merging) "
                         "where the name field (ie. column 4) contains additional info besides the repeat unit.", action="store_true")    
    merge_p.add_argument("-L", "--interval", help="Only process loci in this genomic interval (format: chrom:start-end)", action="append")
    merge_p.add_argument("-o", "--output-prefix", help="Output file prefix. If not specified, it will be computed based on "
                   "the input vcf filename")
    merge_p.add_argument("-v", "--verbose", help="Print detailed logs.", action="store_true")
    merge_p.add_argument("--show-progress-bar", help="Show a progress bar in the terminal when processing variants.", action="store_true")
    merge_p.add_argument("--batch-size", type=int, default=10_000_000, help="Merge tandem repeat loci in batches of this size")
    merge_p.add_argument("input_bed_paths", help="Input BED files generated by the 'catalog' subcommand.", nargs="+")

    # genotype subcommand
    genotype_p = subparsers.add_parser("genotype", help="step3: Genotype tandem repeat loci by looking at the genotypes of indels in the input single-sample VCF file.",
                                       formatter_class=configargparse.ArgumentDefaultsHelpFormatter)
    genotype_p.add_argument("-c", "--config-file", is_config_file=True, help="Optional config file path")

    genotype_p.add_argument("-R", "--reference-fasta-path", help="Reference genome fasta path.", required=True)
    genotype_p.add_argument("--catalog-bed", help="Input BED file containing a tandem repeat catalog where the name field contains "
                            "(or starts with) the locus repeat unit. This catalog can be generated by the 'catalog' subcommand, or can be from other sources.", required=True)
    genotype_p.add_argument("-o", "--output-prefix", help="Output file prefix. If not specified, it will be computed based on "
                            "the input vcf filename")
    genotype_p.add_argument("-L", "--interval", help="Only genotype loci in this genomic interval (format: chrom:start-end). "
                            "Requires the catalog BED to be bgzip compressed and tabix indexed.", action="append")
    genotype_p.add_argument("-v", "--verbose", help="Print detailed logs.", action="store_true")
    genotype_p.add_argument("--show-progress-bar", help="Show a progress bar in the terminal when processing variants.", action="store_true")
    genotype_p.add_argument("--write-vcf", help="Output a VCF file with the subset of variants that contributed to TR genotyping.", action="store_true")
    genotype_p.add_argument("--write-json", help="Output a JSON file containing all genotyped TR loci.", action="store_true")
    genotype_p.add_argument("--add-motif-composition", choices=["trviz", "trf", "basic"],
                            help="Add the parsed motif sequence for each allele to the output (eg. '[CAG][CAG][CCG][CAG]'). "
                            "'trviz' (recommended) aligns the allele sequence to the annotated motif using the "
                            "decomposition algorithm from the trviz python library (which must be installed separately "
                            "with 'pip3 install trviz'), and so allows for insertions or deletions within the repeat "
                            "sequence. 'trf' uses TandemRepeatsFinder for more flexible detection that also allows for "
                            "insertions or deletions within the repeat sequence, while 'basic' naively splits the allele "
                            "sequence into subsequences the size of the annotated motif. The motif sequence is added to "
                            "both the TSV and JSON outputs, while per-allele motif counts are added to the JSON output only.")
    genotype_p.add_argument("--trf-executable-path", help="Path to the TandemRepeatsFinder (TRF) executable. "
                            "Required if --add-motif-composition trf is specified.")
    genotype_p.add_argument("--min-allele-length-for-trf-motif-splitting", type=int, default=12,
                            help="Only used with --add-motif-composition trf. Allele sequences shorter than "
                            "max(this value, 2 * motif_size) are split into motifs using the basic splitting method "
                            "rather than running TandemRepeatsFinder (TRF).")
    genotype_p.add_argument("--max-allele-length-for-edit-distance-purity", type=int,
                            default=DEFAULT_MAX_ALLELE_LENGTH_FOR_EDIT_DISTANCE_PURITY,
                            help="Alleles longer than this many base pairs get an empty RepeatPurityViaEditDistance "
                                 "value, since the edit distance to a pure repeat is expensive to compute for long "
                                 "sequences. The position-by-position RepeatPurity is still computed for them.")
    genotype_p.add_argument("-t", "--trf-threads", default=max(1, multiprocessing.cpu_count() - 2), type=int,
                            help="Only used with --add-motif-composition trf. Number of TandemRepeatsFinder (TRF) "
                            "instances to run in parallel (one per thread) when splitting allele sequences into motifs.")
    genotype_p.add_argument("--threads", default=1, type=int,
                            help="Number of worker processes to use for genotyping loci. The catalog is split into "
                                 "chunks of loci that are genotyped in parallel, and the output is the same as with "
                                 "a single thread.")
    genotype_p.add_argument("--skip-hom-ref-loci", action="store_true",
                            help="Skip loci that were genotyped as homozygous reference, meaning no variant "
                                 "overlapped them, and don't include them in the output. Loci that got no "
                                 "call are always kept, since they are the opposite of a confirmed reference "
                                 "call rather than a case of it.")

    genotype_p.add_argument("--min-insertion-size-to-check-for-repetitiveness", type=int,
                            default=DEFAULT_MIN_INSERTION_SIZE_TO_CHECK,
                            help="A locus where either allele contains an insertion that doesn't look like part of "
                            "the tandem repeat (eg. an Alu element inserted into a poly-A tract, or an assembly "
                            "error) gets no call, since there is no way to say how many repeats it carries. "
                            "Insertions shorter than this many base pairs are always accepted as part of the "
                            "repeat. Short insertions rarely inflate a repeat count much, and there aren't enough "
                            "bases in them to tell a tandem repeat apart from random sequence. Set to a very large "
                            "value to turn this check off.")
    genotype_p.add_argument("--min-insertion-purity", type=float, default=DEFAULT_MIN_INSERTION_PURITY,
                            help="Accept an insertion as part of the repeat if at least this fraction of its bases "
                            "match a pure repeat of the locus motif (trying every starting offset within the motif).")
    genotype_p.add_argument("--min-insertion-periodicity", type=float, default=DEFAULT_MIN_INSERTION_PERIODICITY,
                            help="Also accept an insertion if at least this fraction of its bases match a copy of "
                            "itself shifted by the best-fitting number of bases. This keeps expansions made up of a "
                            "different motif than the one annotated for the locus (eg. an AAGGG expansion at the "
                            "AAAAG reference repeat in RFC1). Random sequence scores around 0.35 by this measure, "
                            "while a tandem repeat of any motif scores close to 1.")
    genotype_p.add_argument("input_vcf_path", help="Input VCF single-sample VCF file containing variant genotypes from "
                                                   "which to compute the TR genotypes.")

    args = p.parse_args()

    if args.subcommand == "catalog":
        if args.dont_allow_interruptions:
            args.dont_run_trf = True

        if not args.dont_run_trf and not args.trf_executable_path:
            p.error(f"Must specify --trf-executable-path or --dont-run-trf")

        if args.copy_info_field_keys_to_tsv:
            args.copy_info_field_keys_to_tsv = {key: 0 for key in args.copy_info_field_keys_to_tsv}
        
    if args.subcommand == "genotype":
        if args.add_motif_composition == "trf" and not args.trf_executable_path:
            p.error("--add-motif-composition trf requires --trf-executable-path to be specified")
        if args.add_motif_composition == "trviz" and importlib.util.find_spec("trviz") is None:
            p.error("--add-motif-composition trviz requires the trviz python library. Install it with "
                    "'pip3 install trviz'")
        if args.threads < 1:
            p.error(f"--threads must be at least 1, not {args.threads}")

    if args.subcommand == "catalog" or args.subcommand == "genotype":
        args.input_vcf_prefix = re.sub(".vcf(.gz|.bgz)?$", "", os.path.basename(args.input_vcf_path))

    return args


class Allele:
    """Represents a single VCF allele."""

    def __init__(self, chrom, pos, ref, alt, fasta_obj, order=-1, info_field_dict=None):
        """Initialize an Allele object

        Args:
            chrom (str): Chromosome name
            pos (int): Position (1-based)
            ref (str): Reference allele sequence
            alt (str): Alt allele sequence
            fasta_obj (pyfaidx.Fasta): Reference fasta object initialized with one_based_attributes = False
            order (int): Optional order of this allele in the VCF file
            info_field_dict (dict): Optional VCF info fields dict
        """
        self._chrom = chrom
        self._pos = pos
        self._ref = ref
        self._alt = alt
        self._fasta_obj = fasta_obj
        self._order = order
        self._info_field_dict = info_field_dict
        
        if fasta_obj.faidx.one_based_attributes:
            raise ValueError("Fasta object should be created with one_based_attributes set to False")
        
        self._left_flank_start_0based = None
        self._left_flank_end = None
        self._right_flank_start_0based = None
        self._right_flank_end = None

        self._left_flanking_reference_sequence = None
        self._right_flanking_reference_sequence = None
        self._left_flank_stops_at_N = False
        self._right_flank_stops_at_N = False

        if len(self._ref) == len(self._alt):
            raise ValueError(f"Logic error: variant {self._chrom}:{self._pos}:{self._ref}:{self._alt} is a SNV/MNV")
        elif len(self._ref) < len(self._alt) and self._alt.startswith(self._ref):
            self._ins_or_del = "INS"
            self._variant_bases = self._alt[len(self._ref):]
        elif len(self._alt) < len(self._ref) and self._ref.startswith(self._alt):
            self._ins_or_del = "DEL"
            self._variant_bases = self._ref[len(self._alt):]
        else:
            raise ValueError(f"Logic error: variant {self._chrom}:{self._pos}:{self._ref}:{self._alt} is a complex MNV insertion/deletion")

        self._k = {"left": 0, "right": 0}  # controls the size of the left and right flanking sequences
        self._already_retrieved_k = {"left": None, "right": None}  # tracks the current size of the left and right flanking sequences
        self._previously_increased_flanking_sequence_size = False  # tracks whether either the left or right flanking sequence size was increased from its starting size
        self._shortened_variant_id = None

    def _get_flanking_sequence_size(self, left_or_right):
        """Computes the size of the left or right flanking sequence (in base pairs) based on the current value of k.
        
        Args:
            left_or_right (str): "left" or "right"

        Returns:
            int: the current size of the left or right flanking sequence (in base pairs)
        """

        # if k = 0, this formula sets the size multiplier to 3, if k = 1 ==> 10, k = 2 ==> 30, k = 3 ==> 100, etc.
        exponent = self._k[left_or_right] // 2
        if self._k[left_or_right] % 2 == 0:
            size_multiplier = 3 * 10**exponent
        else:
            size_multiplier = 10**(exponent+1)

        num_flanking_bases = size_multiplier * max(len(self._variant_bases), 100)

        return num_flanking_bases

    def _retrieve_flanking_sequence(self, left_or_right):
        """Loads the left or right flanking sequence from the reference genome into memory, and updates stored coordinates.
        
        Args:
            left_or_right (str): "left" or "right"
        """

        if self._k[left_or_right] == self._already_retrieved_k[left_or_right]:
            return
        
        self._already_retrieved_k[left_or_right] = self._k[left_or_right]

        if self._ins_or_del == "INS":
            self._left_flank_end = self._pos + len(self._ref) - 1
            self._right_flank_start_0based = self._left_flank_end
        elif self._ins_or_del == "DEL":
            self._left_flank_end = self._pos + len(self._alt) - 1
            self._right_flank_start_0based = self._pos + len(self._ref) - 1
        else:
            raise ValueError(f"Logic error: variant {self._chrom}:{self._pos}:{self._ref}:{self._alt} is a complex MNV insertion/deletion")

        num_flanking_bases = self._get_flanking_sequence_size(left_or_right)
        if left_or_right == "left":
            self._left_flank_start_0based = max(self._left_flank_end - num_flanking_bases, 0)
            self._left_flanking_reference_sequence = str(self._fasta_obj[self._chrom][self._left_flank_start_0based : self._left_flank_end]).upper()

            # Stop at N's in the sequence (search in reverse from variant towards left)
            n_pos = self._left_flanking_reference_sequence[::-1].find('N')
            if n_pos != -1:
                # N found - keep only the portion to the right of the N (closer to variant)
                self._left_flank_stops_at_N = True
                self._left_flanking_reference_sequence = self._left_flanking_reference_sequence[len(self._left_flanking_reference_sequence) - n_pos:]
                self._left_flank_start_0based = self._left_flank_end - len(self._left_flanking_reference_sequence)
        elif left_or_right == "right":
            chrom_size = len(self._fasta_obj[self._chrom])
            self._right_flank_end = min(self._right_flank_start_0based + num_flanking_bases, chrom_size)
            self._right_flanking_reference_sequence = str(self._fasta_obj[self._chrom][self._right_flank_start_0based : self._right_flank_end]).upper()

            # Stop at N's in the sequence (search from variant towards right)
            n_pos = self._right_flanking_reference_sequence.find('N')
            if n_pos != -1:
                # N found - keep only the portion before the N
                self._right_flank_stops_at_N = True
                self._right_flanking_reference_sequence = self._right_flanking_reference_sequence[:n_pos]
                self._right_flank_end = self._right_flank_start_0based + n_pos
        else:
            raise ValueError(f"Logic error: left_or_right must be 'left' or 'right', not {left_or_right}")

    def get_left_flanking_sequence(self):
        self._retrieve_flanking_sequence("left")
        return self._left_flanking_reference_sequence

    def get_right_flanking_sequence(self):
        self._retrieve_flanking_sequence("right")
        return self._right_flanking_reference_sequence


    def get_left_flank_start_0based(self):
        self._retrieve_flanking_sequence("left")
        return self._left_flank_start_0based
    
    def get_left_flank_end(self):
        self._retrieve_flanking_sequence("left")
        return self._left_flank_end


    def get_right_flank_start_0based(self):
        self._retrieve_flanking_sequence("right")
        return self._right_flank_start_0based
    
    def get_right_flank_end(self):
        self._retrieve_flanking_sequence("right")
        return self._right_flank_end

    def get_left_flank_stops_at_N(self):
        self._retrieve_flanking_sequence("left")
        return self._left_flank_stops_at_N

    def get_right_flank_stops_at_N(self):
        self._retrieve_flanking_sequence("right")
        return self._right_flank_stops_at_N

    def increase_left_flanking_sequence_size(self):
        """Increases the size of the left flanking sequence without reading it in yet."""
        self._previously_increased_flanking_sequence_size = True
        self._k["left"] += 1

    def increase_right_flanking_sequence_size(self):
        """Increases the size of the right flanking sequence without reading it in yet."""
        self._previously_increased_flanking_sequence_size = True
        self._k["right"] += 1

    def get_expected_left_flanking_sequence_size(self):
        """Returns the expected size of the left flanking sequence (in base pairs) based on the current value of k."""
        return self._get_flanking_sequence_size("left")
    
    def get_expected_right_flanking_sequence_size(self):
        """Returns the expected size of the right flanking sequence (in base pairs) based on the current value of k."""
        return self._get_flanking_sequence_size("right")

    def increase_flanking_sequence_size(self):
        """Increases the size of both the left and right flanking sequences without reading them in yet."""
        self.increase_left_flanking_sequence_size()
        self.increase_right_flanking_sequence_size()

    def __str__(self):
        return f"{self._chrom}:{self._pos:} {self._ref}>{self._alt} ({self._ins_or_del})"
    
    def __repr__(self):
        return self.__str__()
    
    @property
    def chrom(self):
        return self._chrom
    
    @property
    def pos(self):
        return self._pos
    
    @property
    def ref(self):
        return self._ref
    
    @property
    def alt(self):
        return self._alt
    
    @property
    def ins_or_del(self):
        return self._ins_or_del
    
    @property
    def variant_bases(self):
        return self._variant_bases
    
    @property
    def order(self):
        return self._order

    @property
    def previously_increased_flanking_sequence_size(self):
        return self._previously_increased_flanking_sequence_size

    @property
    def number_of_times_flanking_sequence_size_was_increased(self):
        return self._k["left"] + self._k["right"]
    
    @property
    def variant_id(self):
        return f"{self._chrom}-{self._pos}-{self._ref}-{self._alt}"

    @property
    def shortened_variant_id(self):
        if self._shortened_variant_id is None:
            self._shortened_variant_id = f"{self._chrom}-{self._pos}-"
            self._shortened_variant_id += self._ref if len(self._ref) < 24 else (f"{self._ref[0]}..{len(self._ref)}..{self._ref[-1]}")
            self._shortened_variant_id += "-"
            self._shortened_variant_id += self._alt if len(self._alt) < 24 else (f"{self._alt[0]}..{len(self._alt)}..{self._alt[-1]}")
            self._shortened_variant_id += f"-h{abs(hash(self.variant_id))}"  # add a hash to make the ID unique

        return self._shortened_variant_id

    @property
    def info_field_dict(self):
        return self._info_field_dict


class TandemRepeatAllele:
    """Stores additional information about a VCF insertion or deletion allele that
    represents a tandem repeat expansion or contraction.
    """

    def __init__(
            self, 
            allele,
            repeat_unit,
            adjust_repeat_unit,
            num_repeat_bases_in_left_flank,
            num_repeat_bases_in_variant, 
            num_repeat_bases_in_right_flank, 
            detection_mode,
    ):
        """Initialize a TandemRepeatAllele object.

        Args:
            allele (Allele): the allele record that this TandemRepeatAllele object is based on
            repeat_unit (str): the repeat unit of the tandem repeat allele
            adjust_repeat_unit (bool): whether to set the repeat unit to the most common motif in the variant sequence of the same length as the given repeat unit
            num_repeat_bases_in_left_flank (int): the number of repeat bases in the left flanking sequence
            num_repeat_bases_in_variant (int): the number of repeat bases in the variant
            num_repeat_bases_in_right_flank (int): the number of repeat bases in the right flanking sequence
            detection_mode (str): the detection mode used to find the tandem repeat allele
        """

        self._allele = allele
        self._repeat_unit = repeat_unit
        self._num_repeat_bases_in_left_flank = num_repeat_bases_in_left_flank
        self._num_repeat_bases_in_variant = num_repeat_bases_in_variant
        self._num_repeat_bases_in_right_flank = num_repeat_bases_in_right_flank

        # Set start/end before _adjust_repeat_unit_to_maximize_purity() because its ValueError handler
        # references self.start_0based/end_1based, which would otherwise raise a secondary AttributeError
        # that masks the original error.
        self._start_0based = self._allele.get_left_flank_end() - self._num_repeat_bases_in_left_flank
        self._end_1based = self._allele.get_right_flank_start_0based() + self._num_repeat_bases_in_right_flank

        if adjust_repeat_unit:
            self._adjust_repeat_unit_to_maximize_purity()

        self._repeat_unit_length = len(self._repeat_unit)
        self._detection_mode = detection_mode
        self._summary_string = None

        self._canonical_repeat_unit = None
        self._repeat_purity = None

        if self._start_0based > self._end_1based:
            raise ValueError(f"Logic error: start_0based ({self._start_0based}) > end_1based ({self._end_1based})")

    def _adjust_repeat_unit_to_maximize_purity(self):
        if self._num_repeat_bases_in_left_flank + self._num_repeat_bases_in_variant + self._num_repeat_bases_in_right_flank < len(self._repeat_unit):
            return

        try:
            most_common_motif = compute_most_common_motif(self.variant_and_flanks_repeat_sequence, len(self._repeat_unit))
        except ValueError as e:
            raise ValueError(
                f"Error computing most common motif for allele at {self.chrom}:{self.start_0based}-{self.end_1based} "
                f"with repeat unit '{self._repeat_unit}': {str(e)}"
            ) from e

        simplified_motif, _, _ = find_repeat_unit_without_allowing_interruptions(most_common_motif, allow_partial_repeats=False)
        if simplified_motif != self._repeat_unit:
            self.repeat_unit_adjusted = True  # for debugging
            self._repeat_unit = simplified_motif

    @property
    def chrom(self):
        return self._allele.chrom

    @property
    def start_0based(self):
        return self._start_0based

    @property
    def end_1based(self):
        return self._end_1based

    @property
    def ref_interval_size(self):
        return self._end_1based - self._start_0based

    @property
    def num_repeats_ref(self):
        num_repeat_bases = self._num_repeat_bases_in_left_flank + self._num_repeat_bases_in_right_flank
        if self._allele.ins_or_del == "DEL":
            num_repeat_bases += self._num_repeat_bases_in_variant
        return num_repeat_bases // len(self._repeat_unit)

    @property
    def num_repeats_alt(self):
        num_repeat_bases = self._num_repeat_bases_in_left_flank + self._num_repeat_bases_in_right_flank
        if self._allele.ins_or_del == "INS":
            num_repeat_bases += self._num_repeat_bases_in_variant
        return num_repeat_bases // len(self._repeat_unit)

    @property
    def ref_allele_repeat_sequence(self):
        ref_allele_repeat_sequence = ""
        if self._num_repeat_bases_in_left_flank:
            left_flanking_sequence = self._allele.get_left_flanking_sequence()[-self._num_repeat_bases_in_left_flank:]
            ref_allele_repeat_sequence += left_flanking_sequence

        if self._allele.ins_or_del == "DEL":
            ref_allele_repeat_sequence += self._allele.variant_bases

        if self._num_repeat_bases_in_right_flank:
            right_flanking_sequence = self._allele.get_right_flanking_sequence()[:self._num_repeat_bases_in_right_flank]
            ref_allele_repeat_sequence += right_flanking_sequence

        return ref_allele_repeat_sequence

    @property
    def alt_allele_repeat_sequence(self):
        alt_allele_repeat_sequence = ""
        if self._num_repeat_bases_in_left_flank:
            left_flanking_sequence = self._allele.get_left_flanking_sequence()[-self._num_repeat_bases_in_left_flank:]
            alt_allele_repeat_sequence += left_flanking_sequence

        if self._allele.ins_or_del == "INS":
            alt_allele_repeat_sequence += self._allele.variant_bases

        if self._num_repeat_bases_in_right_flank:
            right_flanking_sequence = self._allele.get_right_flanking_sequence()[:self._num_repeat_bases_in_right_flank]
            alt_allele_repeat_sequence += right_flanking_sequence

        return alt_allele_repeat_sequence

    @property
    def variant_and_flanks_repeat_sequence(self):
        variant_and_flanks_repeat_sequence = ""
        if self._num_repeat_bases_in_left_flank:
            left_flanking_sequence = self._allele.get_left_flanking_sequence()[-self._num_repeat_bases_in_left_flank:]
            variant_and_flanks_repeat_sequence += left_flanking_sequence

        variant_and_flanks_repeat_sequence += self._allele.variant_bases

        if self._num_repeat_bases_in_right_flank:
            right_flanking_sequence = self._allele.get_right_flanking_sequence()[:self._num_repeat_bases_in_right_flank]
            variant_and_flanks_repeat_sequence += right_flanking_sequence

        return variant_and_flanks_repeat_sequence

    @property
    def num_repeats_in_variant_and_flanks(self):
        return self.num_repeats_alt if self.ins_or_del == "INS" else self.num_repeats_ref
    
    @property
    def allele(self):
        return self._allele

    @property
    def repeat_unit_length(self):
        return self._repeat_unit_length

    @property
    def repeat_unit(self):
        return self._repeat_unit

    @property
    def canonical_repeat_unit(self):
        if self._canonical_repeat_unit is None:
            self._canonical_repeat_unit = compute_canonical_motif(self.repeat_unit, include_reverse_complement=True)
        return self._canonical_repeat_unit

    @property
    def detection_mode(self):
        return self._detection_mode

    @property
    def locus_id(self):
        return f"{self._allele.chrom}-{self._start_0based}-{self._end_1based}-{self.repeat_unit}"

    @property
    def repeat_purity(self):
        if self._repeat_purity is None:
            self._repeat_purity, _ = compute_repeat_purity(
                self.variant_and_flanks_repeat_sequence, self.repeat_unit, include_partial_repeats=True)
        return self._repeat_purity

    @property
    def summary_string(self):
        if self._summary_string is None:                
            self._summary_string = f"{self.repeat_unit_length}bp:"
            self._summary_string += f"{self.num_repeats_in_variant_and_flanks:0.1f}x:"
            if self.repeat_unit_length > 30:
                self._summary_string += f"{self.repeat_unit[:30]}...:"
            else:
                self._summary_string += f"{self.repeat_unit}:"
            self._summary_string += f"{self.detection_mode}"
            self._summary_string += f":p{self.repeat_purity:0.2f}"

        return self._summary_string

    @property
    def num_repeat_bases_in_left_flank(self):
        return self._num_repeat_bases_in_left_flank
    
    @property
    def num_repeat_bases_in_variant(self):
        return self._num_repeat_bases_in_variant
    
    @property
    def num_repeat_bases_in_right_flank(self):
        return self._num_repeat_bases_in_right_flank
    

    @property
    def num_repeats_in_left_flank(self):
        return self._num_repeat_bases_in_left_flank // self._repeat_unit_length

    @property
    def num_repeats_in_variant(self):
        return self._num_repeat_bases_in_variant // self._repeat_unit_length

    @property
    def num_repeats_in_right_flank(self):
        return self._num_repeat_bases_in_right_flank // self._repeat_unit_length
    
    @property
    def ins_or_del(self):
        return self._allele.ins_or_del

    @property
    def is_pure_repeat(self):
        return self.repeat_purity > 0.99999
    
    @property
    def order(self):
        return self._allele.order

    @property
    def variant_id(self):
        return f"{self._allele.chrom}-{self._allele.pos}-{self._allele.ref}-{self._allele.alt}"

    @property
    def info_field_dict(self):
        return self._allele.info_field_dict

    def do_repeats_cover_entire_left_flanking_sequence(self):
        return self._num_repeat_bases_in_left_flank > len(self._allele.get_left_flanking_sequence()) - self._repeat_unit_length

    def do_repeats_cover_entire_right_flanking_sequence(self):
        return self._num_repeat_bases_in_right_flank > len(self._allele.get_right_flanking_sequence()) - self._repeat_unit_length

    def do_repeats_cover_entire_flanking_sequence(self):
        return self.do_repeats_cover_entire_left_flanking_sequence() or self.do_repeats_cover_entire_right_flanking_sequence()
    

    def __str__(self):
        return (f"{self._allele.chrom}:{self._start_0based}-{self._end_1based} "
                f"{self.num_repeats_ref}x{self.repeat_unit} ({self.repeat_unit_length}bp) [{self._detection_mode}]")
    
    def __repr__(self):
        return self.__str__()
    


class ReferenceTandemRepeat:
    """Represents a tandem repeat locus in the reference genome.

    This class stores the coordinates and motif of a TR locus from a catalog.
    It uses 0-based half-open coordinates (BED format) where start_0based is
    inclusive and end_1based is exclusive.

    Attributes:
        chrom (str): Chromosome name
        start_0based (int): Start position (0-based, inclusive)
        end_1based (int): End position (1-based inclusive / 0-based exclusive)
        repeat_unit (str): The repeat motif sequence
        repeat_unit_length (int): Length of the repeat unit in bp
        canonical_repeat_unit (str): Canonical form of the motif
        locus_id (str): Unique identifier in format 'chr-start-end-motif'
        num_repeats_ref (int): Number of full repeats in the reference
    """

    def __init__(
            self,
            chrom,
            start_0based,
            end_1based,
            repeat_unit,
            detection_mode=None,
        ):
        """Initialize a ReferenceTandemRepeat object.

        Args:
            chrom (str): Chromosome name (e.g., "chr1" or "1")
            start_0based (int): Start position (0-based, inclusive)
            end_1based (int): End position (1-based, inclusive; same as 0-based exclusive)
            repeat_unit (str): The repeat motif sequence (e.g., "CAG", "AAGGG")
            detection_mode (str): Detection mode used to identify this TR locus
                (e.g., "pure_repeats", "allow_interruptions")

        Raises:
            ValueError: If start_0based > end_1based
        """

        if start_0based > end_1based:
            raise ValueError(f"start_0based ({start_0based}) > end_1based ({end_1based})")

        self._chrom = chrom
        self._start_0based = start_0based
        self._end_1based = end_1based

        self._repeat_unit = repeat_unit
        self._repeat_unit_length = len(repeat_unit)
        self._detection_mode = detection_mode

        self._canonical_repeat_unit = None
        self._summary_string = None

    @property
    def chrom(self):
        """Chromosome name."""
        return self._chrom

    @property
    def start_0based(self):
        """Start position (0-based, inclusive)."""
        return self._start_0based

    @property
    def end_1based(self):
        """End position (1-based inclusive, same as 0-based exclusive)."""
        return self._end_1based

    @property
    def repeat_unit_length(self):
        """Length of the repeat unit in base pairs."""
        return self._repeat_unit_length

    @property
    def repeat_unit(self):
        """The repeat motif sequence (e.g., 'CAG')."""
        return self._repeat_unit

    @property
    def canonical_repeat_unit(self):
        """Canonical form of the repeat unit.

        Computed lazily by normalizing for rotations and reverse complement.
        """
        if self._canonical_repeat_unit is None:
            self._canonical_repeat_unit = compute_canonical_motif(self.repeat_unit, include_reverse_complement=True)
        return self._canonical_repeat_unit

    @property
    def detection_mode(self):
        """Detection mode used to identify this TR locus."""
        return self._detection_mode

    @property
    def ref_interval_size(self):
        """Size of the reference interval in base pairs."""
        return self._end_1based - self._start_0based

    @property
    def locus_id(self):
        """Unique identifier in format 'chr-start-end-motif' (0-based start)."""
        return f"{self._chrom}-{self._start_0based}-{self._end_1based}-{self.repeat_unit}"

    @property
    def num_repeats_ref(self):
        """Number of full repeats in the reference interval."""
        return self.ref_interval_size // self.repeat_unit_length

    @property
    def summary_string(self):
        """Human-readable summary string with motif size, repeat count, motif, and detection mode."""
        if self._summary_string is None:
            self._summary_string = f"{self.repeat_unit_length}bp:"
            self._summary_string += f"{self.ref_interval_size/self.repeat_unit_length:0.1f}x:"
            if self.repeat_unit_length > 30:
                self._summary_string += f"{self.repeat_unit[:30]}...:"
            else:
                self._summary_string += f"{self.repeat_unit}:"
            self._summary_string += f"{self.detection_mode}"

        return self._summary_string

    def __str__(self):
        return self.locus_id

    def __repr__(self):
        return self.__str__()


def min_with_None_check(value1, value2):
    """Return the smaller of two values, ignoring a None, or None if both are None.

    Args:
        value1 (float): a value, or None
        value2 (float): a value, or None

    Returns:
        float: min(value1, value2), the one that isn't None, or None
    """
    if value1 is None:
        return value2
    if value2 is None:
        return value1
    return min(value1, value2)


class GenotypedTandemRepeat:
    """Represents a genotyped tandem repeat locus with computed allele information.

    This class stores the TR locus information from a catalog along with the VCF
    variant(s) that overlap the locus and the computed genotype information including
    repeat counts and full allele sequences for both haplotypes.

    Attributes:
        chrom (str): Chromosome name
        start_0based (int): Start position (0-based, inclusive)
        end (int): End position (0-based exclusive / 1-based inclusive)
        motif (str): The repeat motif
        locus_id (str): Unique identifier in format 'chr-start-end-motif'
        zygosity (str): 'HOM', 'HET', 'HEMI', or None if missing
        num_repeats_short_allele (int): Repeat count in shorter allele
        num_repeats_long_allele (int): Repeat count in longer allele
        allele1_sequence (str): Full sequence for haplotype 0, or None if that allele has no call
        allele2_sequence (str): Full sequence for haplotype 1, or None if that allele has no call

    Example:
        >>> tr_locus = ReferenceTandemRepeat("chr4", 3074876, 3074933, "CAG")
        >>> genotyped = genotype_single_locus(tr_locus, vcf_file, fasta_obj)
        >>> print(f"{genotyped.locus_id}: {genotyped.zygosity} "
        ...       f"{genotyped.num_repeats_short_allele}/{genotyped.num_repeats_long_allele}")
        chr4-3074876-3074933-CAG: HET 19/22
    """

    def __init__(
            self,
            tr_locus,
            overlapping_variants=None,
            allele1_sequence=None,
            allele2_sequence=None,
            num_repeats_allele1=None,
            num_repeats_allele2=None,
            allele1_purity=None,
            allele2_purity=None,
            num_alleles_with_non_repeat_insertions=0,
            num_alleles_with_build_errors=0,
            no_call_reason=None,
            no_call_detail=None,
            allele1_purity_via_edit_distance=None,
            allele2_purity_via_edit_distance=None,
        ):
        """Initialize a GenotypedTandemRepeat object.

        Args:
            tr_locus (ReferenceTandemRepeat): The tandem repeat locus from the catalog
            overlapping_variants (list): OverlappingVariant tuples (chrom, pos, ref, alts), one per VCF record
                that changes at least one base inside this locus (see get_overlapping_vcf_variants)
            allele1_sequence (str): Full haplotype sequence for allele 1 (haplotype 0)
            allele2_sequence (str): Full haplotype sequence for allele 2 (haplotype 1)
            num_repeats_allele1 (int): Number of repeats in allele 1
            num_repeats_allele2 (int): Number of repeats in allele 2
            allele1_purity (float): Repeat purity for allele 1 (0.0-1.0), from comparing it position by position
                to a pure repeat of the same length
            allele2_purity (float): Repeat purity for allele 2 (0.0-1.0), computed the same way
            num_alleles_with_non_repeat_insertions (int): How many of the two alleles contained an insertion
                that isn't part of the tandem repeat. Any nonzero value means the whole locus has no call.
            num_alleles_with_build_errors (int): How many of the two alleles had variants that couldn't be
                applied to the reference. Any nonzero value means the whole locus has no call.
            no_call_reason (str): Why the locus was set to no call (one of the NO_CALL_REASON_* values),
                or None if it wasn't. This is the value runs are summarized by, so it carries no per-locus
                detail.
            no_call_detail (str): What specifically went wrong at this locus, or None. Reported alongside
                no_call_reason in the output columns.
            allele1_purity_via_edit_distance (float): Repeat purity for allele 1 (0.0-1.0), from the edit
                distance to a pure repeat of the same length, so an indel isn't charged for every base it shifts
            allele2_purity_via_edit_distance (float): Repeat purity for allele 2 (0.0-1.0), computed the same way
        """
        self._tr_locus = tr_locus
        self._overlapping_variants = overlapping_variants if overlapping_variants else []
        self._allele1_sequence = allele1_sequence
        self._allele2_sequence = allele2_sequence
        self._num_repeats_allele1 = num_repeats_allele1
        self._num_repeats_allele2 = num_repeats_allele2
        self._allele1_purity = allele1_purity
        self._allele2_purity = allele2_purity
        self._num_alleles_with_non_repeat_insertions = num_alleles_with_non_repeat_insertions
        self._num_alleles_with_build_errors = num_alleles_with_build_errors
        self._no_call_reason = no_call_reason
        self._no_call_detail = no_call_detail
        self._allele1_purity_via_edit_distance = allele1_purity_via_edit_distance
        self._allele2_purity_via_edit_distance = allele2_purity_via_edit_distance

    # Properties from the underlying TR locus
    @property
    def chrom(self):
        return self._tr_locus.chrom

    @property
    def start_0based(self):
        return self._tr_locus.start_0based

    @property
    def end(self):
        """End coordinate (0-based half-open, same as end_1based for BED format)."""
        return self._tr_locus.end_1based

    @property
    def motif(self):
        return self._tr_locus.repeat_unit

    @property
    def canonical_motif(self):
        return self._tr_locus.canonical_repeat_unit

    @property
    def motif_size(self):
        return self._tr_locus.repeat_unit_length

    @property
    def locus(self):
        """Locus string in format chr:start-end (0-based start)."""
        return f"{self.chrom}:{self.start_0based}-{self.end}"

    @property
    def locus_id(self):
        """Locus ID in format chr-start-end-motif (0-based start)."""
        return self._tr_locus.locus_id

    @property
    def num_repeats_in_reference(self):
        return self._tr_locus.num_repeats_ref

    # Genotype-specific properties
    @property
    def overlapping_variants(self):
        return self._overlapping_variants

    @property
    def num_overlapping_variants(self):
        return len(self._overlapping_variants)

    @property
    def variant_positions(self):
        """List of 1-based positions of overlapping variants."""
        return [v.pos for v in self._overlapping_variants]

    @property
    def allele1_sequence(self):
        return self._allele1_sequence

    @property
    def allele2_sequence(self):
        return self._allele2_sequence

    @property
    def num_repeats_allele1(self):
        return self._num_repeats_allele1

    @property
    def num_repeats_allele2(self):
        return self._num_repeats_allele2

    @property
    def allele1_purity(self):
        return self._allele1_purity

    @property
    def allele2_purity(self):
        return self._allele2_purity

    @property
    def allele1_purity_via_edit_distance(self):
        return self._allele1_purity_via_edit_distance

    @property
    def allele2_purity_via_edit_distance(self):
        return self._allele2_purity_via_edit_distance

    @property
    def num_alleles_with_non_repeat_insertions(self):
        """How many of the two alleles contained an insertion that isn't part of the tandem repeat. Any nonzero
        value means the whole locus was set to no-call."""
        return self._num_alleles_with_non_repeat_insertions

    @property
    def num_alleles_with_build_errors(self):
        """How many of the two alleles had variants that couldn't be applied to the reference sequence. Any
        nonzero value means the whole locus was set to no-call."""
        return self._num_alleles_with_build_errors

    @property
    def no_call_reason(self):
        """Why the locus was set to no call, or None if it wasn't. One of the NO_CALL_REASON_* values, with
        no per-locus detail, so a run can be summarized by it. This distinguishes a locus the VCF cannot
        answer for from one it reports as homozygous reference."""
        return self._no_call_reason

    @property
    def no_call_detail(self):
        """What specifically went wrong at this locus, or None."""
        return self._no_call_detail

    @property
    def no_call_description(self):
        """The no-call reason together with its per-locus detail, as written to the output columns."""
        if self._no_call_reason is None:
            return None
        return describe_no_call(self._no_call_reason, self._no_call_detail)

    def _purity_in_short_long_order(self, p1, p2):
        """Return (short_allele_purity, long_allele_purity) ordered so the values line up with the other
        short/long columns.

        Repeat counts are truncated to whole motifs, so two alleles that differ by less than one motif tie on
        count while still differing in length. Ordering on count alone would then leave the purity columns
        paired with the opposite alleles from repeat_size_short_allele_bp / repeat_size_long_allele_bp, so
        sequence length breaks the tie here the same way it decides the size columns.

        For HEMI loci (only one allele present), the present allele's purity is used for both. Returns
        (None, None) if the genotype is missing.

        Args:
            p1 (float): allele 1's purity (either kind), or None
            p2 (float): allele 2's purity, computed the same way as p1, or None
        """
        n1, n2 = self._num_repeats_allele1, self._num_repeats_allele2
        if n1 is None and n2 is None:
            return None, None
        if n1 is None:
            return p2, p2
        if n2 is None:
            return p1, p1
        return (p1, p2) if (n1, self.allele1_size_bp) <= (n2, self.allele2_size_bp) else (p2, p1)

    @property
    def repeat_purity_short_allele(self):
        """Repeat purity of the shorter allele (matches the short-allele repeat-count and size columns)."""
        return self._purity_in_short_long_order(self._allele1_purity, self._allele2_purity)[0]

    @property
    def repeat_purity_long_allele(self):
        """Repeat purity of the longer allele (matches the long-allele repeat-count and size columns)."""
        return self._purity_in_short_long_order(self._allele1_purity, self._allele2_purity)[1]

    @property
    def repeat_purity_via_edit_distance_short_allele(self):
        """Edit-distance repeat purity of the shorter allele (matches the short-allele columns)."""
        return self._purity_in_short_long_order(
            self._allele1_purity_via_edit_distance, self._allele2_purity_via_edit_distance)[0]

    @property
    def repeat_purity_via_edit_distance_long_allele(self):
        """Edit-distance repeat purity of the longer allele (matches the long-allele columns)."""
        return self._purity_in_short_long_order(
            self._allele1_purity_via_edit_distance, self._allele2_purity_via_edit_distance)[1]

    @property
    def num_repeats_short_allele(self):
        """Number of repeats in the shorter allele."""
        if self._num_repeats_allele1 is None and self._num_repeats_allele2 is None:
            return None
        if self._num_repeats_allele1 is None:
            return self._num_repeats_allele2
        if self._num_repeats_allele2 is None:
            return self._num_repeats_allele1
        return min(self._num_repeats_allele1, self._num_repeats_allele2)

    @property
    def num_repeats_long_allele(self):
        """Number of repeats in the longer allele."""
        if self._num_repeats_allele1 is None and self._num_repeats_allele2 is None:
            return None
        if self._num_repeats_allele1 is None:
            return self._num_repeats_allele2
        if self._num_repeats_allele2 is None:
            return self._num_repeats_allele1
        return max(self._num_repeats_allele1, self._num_repeats_allele2)

    @property
    def allele1_size_bp(self):
        """Size in bp of allele 1."""
        if self._allele1_sequence is None:
            return None
        return len(self._allele1_sequence)

    @property
    def allele2_size_bp(self):
        """Size in bp of allele 2."""
        if self._allele2_sequence is None:
            return None
        return len(self._allele2_sequence)

    @property
    def repeat_size_short_allele_bp(self):
        """Size in bp of the shorter allele (based on actual sequence length)."""
        size1 = self.allele1_size_bp
        size2 = self.allele2_size_bp
        if size1 is None and size2 is None:
            return None
        if size1 is None:
            return size2
        if size2 is None:
            return size1
        return min(size1, size2)

    @property
    def repeat_size_long_allele_bp(self):
        """Size in bp of the longer allele (based on actual sequence length)."""
        size1 = self.allele1_size_bp
        size2 = self.allele2_size_bp
        if size1 is None and size2 is None:
            return None
        if size1 is None:
            return size2
        if size2 is None:
            return size1
        return max(size1, size2)

    @property
    def zygosity(self):
        """Zygosity based on allele length: HOM, HET, or HEMI.

        Repeat counts are truncated to whole motifs, so two alleles that differ by less than one motif copy
        tie on count while still being different lengths. Deciding on count alone would report such a locus
        as HOM while repeat_size_short_allele_bp and repeat_size_long_allele_bp disagree, which hides real
        heterozygous indels: at a VNTR locus with a long motif, any indel smaller than one copy lands in the
        same count bin. Length is what the size columns report, so length decides here too, the same way
        _purity_in_short_long_order() breaks count ties.

        Returns:
            str: 'HOM' if both alleles are the same length,
                 'HET' if the alleles are different lengths,
                 'HEMI' if only one allele is present (hemizygous),
                 None if genotype is missing (both alleles are None)
        """
        if self._num_repeats_allele1 is None and self._num_repeats_allele2 is None:
            return None
        if self._num_repeats_allele1 is None or self._num_repeats_allele2 is None:
            return "HEMI"
        if self.allele1_size_bp == self.allele2_size_bp:
            return "HOM"
        return "HET"

    @property
    def is_pure_repeat(self):
        """Whether every allele with a defined purity is a pure repeat (purity > 0.99).

        An allele counts here only if its purity is defined. An allele deleted down to an empty sequence, or one
        shorter than a single motif copy, has no purity to judge, so it neither passes nor fails: a locus where
        no allele has a defined purity returns None (unknown) rather than False.
        """
        purity_threshold = 0.99

        # Both alleles present with a defined purity: both must be pure
        if self._allele1_purity is not None and self._allele2_purity is not None:
            return self._allele1_purity > purity_threshold and self._allele2_purity > purity_threshold
        # Only one allele has a defined purity: judge on that one
        if self._allele1_purity is not None:
            return self._allele1_purity > purity_threshold
        if self._allele2_purity is not None:
            return self._allele2_purity > purity_threshold
        # Missing genotype, or no allele long enough to have a purity
        return None

    @property
    def repeat_purity(self):
        """Overall repeat purity (minimum of both alleles, or the one present)."""
        return min_with_None_check(self._allele1_purity, self._allele2_purity)

    @property
    def repeat_purity_via_edit_distance(self):
        """Overall edit-distance repeat purity (minimum of both alleles, or the one present), or None if it
        wasn't computed for an allele that is present (see --max-allele-length-for-edit-distance-purity)."""
        # An allele that has a position-by-position purity but no edit-distance purity was skipped for being
        # too long, not absent, so the other allele's value can't stand in as the minimum over both.
        if ((self._allele1_purity is not None and self._allele1_purity_via_edit_distance is None)
                or (self._allele2_purity is not None and self._allele2_purity_via_edit_distance is None)):
            return None
        return min_with_None_check(self._allele1_purity_via_edit_distance, self._allele2_purity_via_edit_distance)

    def to_tsv_dict(self, motif_lists=None):
        """Convert this genotyped locus to a dictionary for TSV output.

        Args:
            motif_lists (dict): Optional dict with keys 'allele1' and 'allele2', each a parsed-motif entry
                {"motifs": [...], "prefix": str, "suffix": str} (or None), plus 'allele1_method' and
                'allele2_method' naming the method used ("trf" or "basic-split"). If provided, the parsed
                motif sequence (eg. "CA[GCA][GCA][GCC]G") and detection method are added for each allele.

        Returns:
            dict: Dictionary with keys matching GENOTYPE_TSV_OUTPUT_COLUMNS.
                Values are formatted appropriately for TSV output:
                - None values become empty strings
                - Lists are comma-joined
                - Floats are formatted to 4 decimal places
        """
        # Format purity: None -> empty string, otherwise 4 decimal places
        def format_purity(purity):
            return f"{purity:.4f}" if purity is not None else ""

        # Format is_pure_repeat: None -> empty, bool -> True/False
        is_pure_str = ""
        if self.is_pure_repeat is not None:
            is_pure_str = str(self.is_pure_repeat)

        # Format variant positions as comma-separated list
        variant_positions_str = ",".join(str(p) for p in self.variant_positions)

        return {
            "Chrom": self.chrom,
            "Start0Based": self.start_0based,
            "End": self.end,
            "Locus": self.locus,
            "LocusId": self.locus_id,
            "Motif": self.motif,
            "CanonicalMotif": self.canonical_motif,
            "MotifSize": self.motif_size,
            "NumRepeatsInReference": self.num_repeats_in_reference if self.num_repeats_in_reference is not None else "",
            "NumRepeatsShortAllele": self.num_repeats_short_allele if self.num_repeats_short_allele is not None else "",
            "NumRepeatsLongAllele": self.num_repeats_long_allele if self.num_repeats_long_allele is not None else "",
            "RepeatSizeShortAlleleBp": self.repeat_size_short_allele_bp if self.repeat_size_short_allele_bp is not None else "",
            "RepeatSizeLongAlleleBp": self.repeat_size_long_allele_bp if self.repeat_size_long_allele_bp is not None else "",
            "Zygosity": self.zygosity if self.zygosity is not None else "",
            "IsPureRepeat": is_pure_str,
            "RepeatPurity": format_purity(self.repeat_purity),
            "RepeatPurityShortAllele": format_purity(self.repeat_purity_short_allele),
            "RepeatPurityLongAllele": format_purity(self.repeat_purity_long_allele),
            "RepeatPurityViaEditDistance": format_purity(self.repeat_purity_via_edit_distance),
            "RepeatPurityViaEditDistanceShortAllele": format_purity(self.repeat_purity_via_edit_distance_short_allele),
            "RepeatPurityViaEditDistanceLongAllele": format_purity(self.repeat_purity_via_edit_distance_long_allele),
            "Allele1Sequence": self.allele1_sequence if self.allele1_sequence is not None else "",
            "Allele2Sequence": self.allele2_sequence if self.allele2_sequence is not None else "",
            "Allele1MotifSequence": (format_motif_entry_as_sequence_string(motif_lists.get("allele1")) or "") if motif_lists else "",
            "Allele1SequenceMotifSplittingMethod": (motif_lists.get("allele1_method") or "") if motif_lists else "",
            "Allele2MotifSequence": (format_motif_entry_as_sequence_string(motif_lists.get("allele2")) or "") if motif_lists else "",
            "Allele2SequenceMotifSplittingMethod": (motif_lists.get("allele2_method") or "") if motif_lists else "",
            "NumOverlappingVariants": self.num_overlapping_variants,
            "VariantPositions": variant_positions_str,
            "NoCallReason": self.no_call_description if self.no_call_reason is not None else "",
        }

    def to_json_dict(self, motif_lists=None):
        """Convert this genotyped locus to a dictionary for JSON output.

        Args:
            motif_lists (dict): Optional dict with keys 'allele1' and 'allele2', each a parsed-motif entry
                {"motifs": [...], "prefix": str, "suffix": str} (or None), plus 'allele1_method' and
                'allele2_method' naming the method used ("trf" or "basic-split"). If provided, the parsed
                motif sequence (eg. "CA[GCA][GCA][GCC]G"), per-allele motif counts (full motifs only,
                excluding prefix/suffix), and detection method are added for each allele.

        Returns:
            dict: Dictionary with all fields for JSON serialization. Unlike to_tsv_dict(),
                this preserves native types (int, float, bool, list) rather than converting
                to strings. None values are preserved as None (which becomes null in JSON).
        """
        result = {
            "Chrom": self.chrom,
            "Start0Based": self.start_0based,
            "End": self.end,
            "Locus": self.locus,
            "LocusId": self.locus_id,
            "Motif": self.motif,
            "CanonicalMotif": self.canonical_motif,
            "MotifSize": self.motif_size,
            "NumRepeatsInReference": self.num_repeats_in_reference,
            "NumRepeatsShortAllele": self.num_repeats_short_allele,
            "NumRepeatsLongAllele": self.num_repeats_long_allele,
            "RepeatSizeShortAlleleBp": self.repeat_size_short_allele_bp,
            "RepeatSizeLongAlleleBp": self.repeat_size_long_allele_bp,
            "Zygosity": self.zygosity,
            "IsPureRepeat": self.is_pure_repeat,
            "RepeatPurity": round(self.repeat_purity, 3) if self.repeat_purity is not None else None,
            "RepeatPurityShortAllele": round(self.repeat_purity_short_allele, 3) if self.repeat_purity_short_allele is not None else None,
            "RepeatPurityLongAllele": round(self.repeat_purity_long_allele, 3) if self.repeat_purity_long_allele is not None else None,
            "RepeatPurityViaEditDistance": (round(self.repeat_purity_via_edit_distance, 3)
                                            if self.repeat_purity_via_edit_distance is not None else None),
            "RepeatPurityViaEditDistanceShortAllele": (
                round(self.repeat_purity_via_edit_distance_short_allele, 3)
                if self.repeat_purity_via_edit_distance_short_allele is not None else None),
            "RepeatPurityViaEditDistanceLongAllele": (
                round(self.repeat_purity_via_edit_distance_long_allele, 3)
                if self.repeat_purity_via_edit_distance_long_allele is not None else None),
            "Allele1Sequence": self.allele1_sequence,
            "Allele2Sequence": self.allele2_sequence,
            "NumOverlappingVariants": self.num_overlapping_variants,
            "VariantPositions": self.variant_positions,
            "NoCallReason": self.no_call_description,
        }

        if motif_lists:
            allele1_entry = motif_lists.get("allele1")
            allele2_entry = motif_lists.get("allele2")
            result["Allele1MotifCounts"] = dict(collections.Counter(allele1_entry["motifs"])) if allele1_entry and allele1_entry["motifs"] else None
            result["Allele2MotifCounts"] = dict(collections.Counter(allele2_entry["motifs"])) if allele2_entry and allele2_entry["motifs"] else None
            result["Allele1MotifSequence"] = format_motif_entry_as_sequence_string(allele1_entry)
            result["Allele2MotifSequence"] = format_motif_entry_as_sequence_string(allele2_entry)
            result["Allele1SequenceMotifSplittingMethod"] = motif_lists.get("allele1_method")
            result["Allele2SequenceMotifSplittingMethod"] = motif_lists.get("allele2_method")

        return result

    def __str__(self):
        return f"{self.locus_id}:{self.zygosity}:{self.num_repeats_short_allele}/{self.num_repeats_long_allele}"

    def __repr__(self):
        return self.__str__()


def get_PAR_region_coordinates(fasta_obj, fasta_contig_lookup):
    """Look up the pseudoautosomal regions of the reference, identifying the genome build by the length of chrX.

    Args:
        fasta_obj (pyfaidx.Fasta): the reference genome
        fasta_contig_lookup (dict): the reference's contig names indexed by build_contig_name_lookup

    Returns:
        dict: "X" and "Y" to lists of (start_0based, end) PAR intervals, or an empty dict when the reference
            has no chrX or is not one of the builds in PAR_REGIONS_BY_CHRX_LENGTH. With no PARs known, all of
            chrX and chrY is treated as non-PAR.
    """
    fasta_chrx = fasta_contig_lookup.get("X")
    if fasta_chrx is None:
        return {}

    chrx_length = len(fasta_obj[fasta_chrx])
    if chrx_length not in PAR_REGIONS_BY_CHRX_LENGTH:
        print(f"WARNING: chrX length {chrx_length} doesn't match GRCh37, GRCh38 or T2T-CHM13, so the "
              f"pseudoautosomal regions are unknown. All of chrX and chrY will be treated as non-PAR.")
        return {}

    genome_version, par_regions = PAR_REGIONS_BY_CHRX_LENGTH[chrx_length]
    print(f"Reference chrX length matches {genome_version}; using its pseudoautosomal regions")
    return par_regions


def overlaps_par(chrom, start_0based, end, par_regions):
    """Check whether an interval overlaps a pseudoautosomal region, even partly.

    Args:
        chrom (str): the chromosome, in any naming convention
        start_0based (int): interval start position (0-based, inclusive)
        end (int): interval end position (0-based, exclusive)
        par_regions (dict): "X" and "Y" to lists of (start_0based, end) PAR intervals (see get_PAR_region_coordinates),
            or None when no PARs are known

    Returns:
        bool: True if the interval overlaps a PAR
    """
    return any(par_start < end and start_0based < par_end
               for par_start, par_end in (par_regions or {}).get(normalize_chromosome_name(chrom), []))


def get_locus_ploidy(chrom, start_0based, end, par_regions, sex_chromosome_ploidy):
    """Look up how many copies of a locus this sample carries.

    Autosomes and the PARs are diploid. On chrX and chrY outside the PARs the ploidy is whatever
    detect_sex_chromosome_ploidy found for that chromosome. A locus that overlaps a PAR even partly counts as
    diploid, since a "." allele there may be an uncalled haplotype. The exception is a chromosome the sample
    lacks altogether (chrY in an XX sample): its PAR loci are as absent as the rest of it, so they stay at
    ploidy 0 rather than being reported as diploid reference.

    Args:
        chrom (str): the locus chromosome, in any naming convention
        start_0based (int): locus start position (0-based, inclusive)
        end (int): locus end position (0-based, exclusive)
        par_regions (dict): "X" and "Y" to lists of (start_0based, end) PAR intervals (see get_PAR_region_coordinates),
            or None when no PARs are known
        sex_chromosome_ploidy (dict): "X" and "Y" to the sample's ploidy of each outside the PARs (see
            detect_sex_chromosome_ploidy), or None to treat both as diploid

    Returns:
        int: the ploidy of the locus. 0 means the sample lacks the chromosome.
    """
    normalized_chrom = normalize_chromosome_name(chrom)
    if normalized_chrom not in ("X", "Y") or sex_chromosome_ploidy is None:
        return 2

    chromosome_ploidy = sex_chromosome_ploidy[normalized_chrom]
    if chromosome_ploidy == 0:
        return 0

    return 2 if overlaps_par(chrom, start_0based, end, par_regions) else chromosome_ploidy


def get_called_non_par_genotypes(vcf_file, vcf_chrom, par_regions):
    """Yield the GT of each record on a sex chromosome outside the PARs that calls at least one allele.

    Args:
        vcf_file (pysam.VariantFile): the open single-sample VCF
        vcf_chrom (str): the chromosome, spelled the way the VCF spells it, or None if the VCF lacks it
        par_regions (dict): "X" and "Y" to lists of (start_0based, end) PAR intervals (see get_PAR_region_coordinates)

    Yields:
        tuple: the GT tuple of each such record
    """
    if vcf_chrom is None:
        return

    try:
        records = vcf_file.fetch(vcf_chrom)
    except ValueError:
        return

    for variant in records:
        if overlaps_par(vcf_chrom, variant.start, variant.stop, par_regions):
            continue
        gt = variant.samples[0].get("GT")
        if gt and any(allele is not None for allele in gt):
            yield gt


def detect_sex_chromosome_ploidy(vcf_file, vcf_contig_lookup, par_regions):
    """Detect the sample's ploidy of chrX and chrY outside the PARs from its own records.

    Cutoffs come from a survey of 331 DipCall VCFs. chrX is haploid if any of three signals says so:

    1. When a male's assembly haplotypes are split by parent (as in HPRC), the paternal haplotype has no chrX
       sequence and the maternal one no chrY, so DipCall writes non-PAR chrX calls as ".|1" (and chrY as
       "1|."), and nearly every called record there has one uncalled allele (99.9-100%). In a female, one uncalled allele
       instead means one assembly didn't cover the site, which affects only 1-18% of chrX records. More than
       half marks chrX as haploid, once there are at least MIN_CHRX_RECORDS_FOR_PLOIDY_DETECTION called
       records for the fraction to mean something. A single-allele GT such as "1", which callers run with
       ploidy 1 on chrX write, calls one allele just as ".|1" does and counts the same way.
    2. When a male's assembly haplotypes aren't split by parent (HGSVC h1/h2), both carry chrX and chrY
       sequence, so DipCall writes non-PAR chrX as diploid, mostly "1|1". Its male option doesn't change these
       genotypes; it only adds a DIPX or DIPY filter. The first signal misses such samples. Among records
       that call both alleles and carry an alt, only 4-19% are
       heterozygous in such males versus 50-71% in females. Below MAX_CHRX_HET_FRACTION_FOR_HAPLOID_X marks
       chrX as haploid, once there are enough such records for the fraction to mean something.
    3. When more than half of the called records have a single-allele GT such as "1", the caller itself
       treated chrX as haploid. Unlike ".|1", which can also mean one assembly didn't cover the site, a
       single-allele GT has no other explanation, so this needs no minimum record count. Without it, a
       VCF from such a caller that covers little of chrX would fall back to diploid, and every chrX
       locus it calls would become a no call.

    chrY is haploid if it has at least MIN_CHRY_RECORDS_FOR_HAPLOID_Y called records, and absent otherwise.
    Detecting the two separately is what lets an XXY sample keep a diploid chrX and a haploid chrY.

    Cases these signals can't see: an XXY sample whose two X copies are identical looks like XY, and extra
    copies beyond two (XXX, XYY) don't show up in a two-haplotype VCF at all.

    Args:
        vcf_file (pysam.VariantFile): the open single-sample VCF
        vcf_contig_lookup (dict): the VCF's contig names indexed by build_contig_name_lookup
        par_regions (dict): "X" and "Y" to lists of (start_0based, end) PAR intervals (see get_PAR_region_coordinates)

    Returns:
        dict: "X" to 1 or 2, and "Y" to 0 or 1. chrX is 2 when no test marks it haploid, including when the
            VCF has too few called non-PAR chrX records for the first two.
    """
    num_called = 0
    num_with_one_called_allele = 0
    num_single_allele_gt = 0
    num_het = 0
    num_hom_alt = 0
    for gt in get_called_non_par_genotypes(vcf_file, vcf_contig_lookup.get("X"), par_regions):
        num_called += 1
        if len(gt) == 1:
            num_single_allele_gt += 1
        if sum(allele is not None for allele in gt) == 1:
            num_with_one_called_allele += 1
        elif len(gt) == 2 and gt[0] != gt[1]:
            num_het += 1
        elif len(gt) == 2 and gt[0] != 0:
            num_hom_alt += 1

    num_diploid_alt = num_het + num_hom_alt
    is_haploid_x_by_uncalled_alleles = (
        num_called >= MIN_CHRX_RECORDS_FOR_PLOIDY_DETECTION
        and num_with_one_called_allele > num_called / 2)
    is_haploid_x_by_het_fraction = (
        num_diploid_alt >= MIN_CHRX_RECORDS_FOR_PLOIDY_DETECTION
        and num_het / num_diploid_alt < MAX_CHRX_HET_FRACTION_FOR_HAPLOID_X)
    is_haploid_x_by_single_allele_genotypes = num_single_allele_gt > num_called / 2
    chrx_ploidy = 1 if (is_haploid_x_by_uncalled_alleles or is_haploid_x_by_het_fraction
                        or is_haploid_x_by_single_allele_genotypes) else 2

    num_called_chry = sum(1 for _ in get_called_non_par_genotypes(
        vcf_file, vcf_contig_lookup.get("Y"), par_regions))
    chry_ploidy = 1 if num_called_chry >= MIN_CHRY_RECORDS_FOR_HAPLOID_Y else 0

    karyotype = "X" * chrx_ploidy + ("Y" * chry_ploidy if chry_ploidy else ("0" if chrx_ploidy == 1 else ""))
    uncalled_fraction = f" ({num_with_one_called_allele / num_called:.1%})" if num_called else ""
    het_fraction = f" ({num_het / num_diploid_alt:.1%})" if num_diploid_alt else ""
    print(f"Non-PAR chrX: {num_with_one_called_allele:,d} of {num_called:,d} called records call only one "
          f"allele{uncalled_fraction} ({num_single_allele_gt:,d} with a single-allele GT), and {num_het:,d} of {num_diploid_alt:,d} records that call both alleles "
          f"with an alt are heterozygous{het_fraction}. Non-PAR chrY: {num_called_chry:,d} called records. "
          f"Detected sex chromosome karyotype: {karyotype} (chrX ploidy {chrx_ploidy}, chrY ploidy "
          f"{chry_ploidy})")
    return {"X": chrx_ploidy, "Y": chry_ploidy}


def match_interval_to_catalog_contig(interval, catalog_contigs):
    """Rewrite an interval's contig to the spelling the catalog uses, or return None if it has no match.

    A catalog written with "chr1" and an interval given as "1:1-1000000" name the same region, but tabix
    matches contig names literally. Without this, such an interval looks exactly like a contig the catalog has
    no records for, and the run would quietly genotype nothing. Names are matched the same way as for the VCF
    and reference (see build_contig_name_lookup), so "chrMT" also finds a "chrM" catalog.

    Args:
        interval (str): a region string ("chr1:1-1000000") or a bare contig name
        catalog_contigs (collection): the contig names present in the catalog's tabix index

    Returns:
        str: the interval with its contig rewritten to the catalog's spelling, or None if no spelling of it is
            in the catalog
    """
    chrom, separator, span = interval.partition(":")
    if chrom in catalog_contigs:
        return interval

    catalog_chrom = build_contig_name_lookup(catalog_contigs).get(normalize_chromosome_name(chrom))
    if catalog_chrom is None:
        return None

    return f"{catalog_chrom}{separator}{span}"


def fetch_catalog_records_within_intervals(tabix_file, intervals, catalog_bed_path):
    """Yield the catalog BED lines that overlap any of the given intervals.

    The interval's contig is first matched against the catalog's own spelling, so a catalog written with
    "chr1" still answers an interval given as "1:1-1000000". Without that, tabix's literal name matching makes
    a naming mismatch look exactly like a contig the catalog has no loci on, and the run genotypes nothing
    while reporting only that the interval was empty.

    An interval on a contig the catalog genuinely has no loci on is reported and skipped rather than aborting
    the run, since that is routine when a genome-wide run is sharded by chromosome: a catalog built from a
    female sample has no chrY loci, one built on the primary assembly has no decoy contigs.

    Args:
        tabix_file (pysam.TabixFile): the open, indexed catalog
        intervals (list): genomic intervals to fetch (eg. ["chr1:1-100000"])
        catalog_bed_path (str): path to the catalog, used in the warning messages

    Yields:
        str: each BED line that overlaps any of the intervals
    """
    catalog_contigs = set(tabix_file.contigs)
    for interval in intervals:
        matched_interval = match_interval_to_catalog_contig(interval, catalog_contigs)
        if matched_interval is None:
            print(f"WARNING: the contig in interval {interval} is not present in {catalog_bed_path} under "
                  f"any of its usual names. Skipping that interval.")
            continue

        # The contig is known to be in the index at this point, so tabix raising means the interval itself is
        # malformed (eg. end before start, or non-numeric coordinates). An empty region just yields nothing.
        try:
            interval_iterator = tabix_file.fetch(matched_interval)
        except ValueError as e:
            raise ValueError(f"Invalid interval '{interval}': {e}")

        yield from interval_iterator


def parse_catalog_bed_line(line, line_num, catalog_bed_path):
    """Split one catalog BED line into its fields and its validated motif.

    Comment, header, UCSC track/browser and blank lines are skipped. tabix already drops '#' lines when a
    catalog is read by interval, so without this the same catalog would parse with -L and fail without it.

    Args:
        line (str): the BED line
        line_num (int): its 1-based line number, for error messages
        catalog_bed_path (str): path to the catalog, for error messages

    Returns:
        2-tuple (list, str): the tab-separated fields and the upper-cased motif, which is the name field's first
            ":"-separated token (eg. "CAG" from "CAG:3bp:19.0x:pure_repeats"), or None for a skipped line

    Raises:
        ValueError: if the line has fewer than 4 columns, or its motif is empty or not made of DNA bases
    """
    if not line.strip() or line.startswith("#") or line.startswith(("track ", "track\t", "browser ")):
        return None

    fields = line.strip().split("\t")
    if len(fields) < 4:
        raise ValueError(f"Invalid BED file format in {catalog_bed_path} on line {line_num}: "
                       f"expected at least 4 columns, got {len(fields)}: {line.strip()}")

    repeat_unit = fields[3].split(":")[0].upper()

    # An empty name field, or one that starts with the ":" separator, yields an empty motif that would otherwise
    # pass the DNA-base check below (the set difference of an empty string is empty) and then divide by zero
    # downstream.
    if not repeat_unit:
        raise ValueError(f"Missing repeat unit in {catalog_bed_path} on line {line_num}: the name field "
                       f"'{fields[3]}' has no motif before the first ':' separator. "
                       f"Line contents: {line.strip()}")

    invalid_bases = set(repeat_unit) - DNA_BASES
    if invalid_bases:
        raise ValueError(f"Invalid repeat unit in {catalog_bed_path} on line {line_num}: "
                       f"'{repeat_unit}' contains non-DNA characters {invalid_bases}. "
                       f"Line contents: {line.strip()}")

    return fields, repeat_unit


def parse_catalog_bed_file(catalog_bed_path, intervals=None, verbose=False):
    """Parse a catalog BED file into a list of ReferenceTandemRepeat objects.

    This function opens a BED catalog file and parses each locus into a
    ReferenceTandemRepeat object. The BED file must have at least 4 columns
    where the name field (column 4) contains (or starts with) the repeat unit.

    Args:
        catalog_bed_path (str): Path to the catalog BED file (can be .gz/.bgz compressed)
        intervals (list): Optional list of genomic intervals to filter to (e.g., ["chr1:100-200"]).
            If specified, the BED file MUST be bgzip compressed and tabix indexed.
        verbose (bool): If True, print progress information

    Returns:
        list: List of ReferenceTandemRepeat objects representing TR loci

    Raises:
        ValueError: If intervals are specified but the BED file is not tabix-indexed
        ValueError: If the BED file format is invalid

    Note:
        Loading all loci into memory may require significant RAM for large catalogs
        (e.g., ~1-2 GB for catalogs with millions of loci). Consider using interval
        filtering with -L to reduce memory usage when processing specific regions.
    """
    tr_loci = []
    input_files_to_close = []

    if intervals:
        # Check if file is tabix-indexed
        tabix_index_path = catalog_bed_path + ".tbi"
        if not file_exists(tabix_index_path):
            raise ValueError(
                f"Cannot filter by intervals: {catalog_bed_path} is not tabix-indexed.\n"
                f"Please index the file using:\n"
                f"  bgzip -c {catalog_bed_path} > {catalog_bed_path}.gz\n"
                f"  tabix -p bed {catalog_bed_path}.gz\n"
                f"Then use the .gz file as input."
            )

        if verbose:
            print(f"Parsing {', '.join(intervals)} from {catalog_bed_path}")

        tabix_file = pysam.TabixFile(catalog_bed_path)
        bed_iterator = fetch_catalog_records_within_intervals(tabix_file, intervals, catalog_bed_path)
        input_files_to_close.append(tabix_file)
    else:
        if verbose:
            print(f"Parsing catalog: {catalog_bed_path}")
        bed_iterator = open_file(catalog_bed_path, is_text_file=True)
        input_files_to_close.append(bed_iterator)

    # Parse the BED file into a list of ReferenceTandemRepeat objects.
    # A locus is emitted once even if several -L intervals cover it: tabix returns every record that overlaps
    # a region, so a locus spanning the boundary between two adjacent intervals is fetched by both, and
    # without this it would be genotyped twice and written twice to every output.
    seen_loci = set()
    for line_num, line in enumerate(bed_iterator, start=1):
        parsed_line = parse_catalog_bed_line(line, line_num, catalog_bed_path)
        if parsed_line is None:
            continue
        fields, repeat_unit = parsed_line

        locus_key = (fields[0], int(fields[1]), int(fields[2]), repeat_unit)
        if locus_key in seen_loci:
            continue
        seen_loci.add(locus_key)

        # Use the motif exactly as specified in the catalog - do NOT simplify
        # This preserves the original annotation and avoids confusion
        tr_loci.append(ReferenceTandemRepeat(
            chrom=fields[0],
            start_0based=int(fields[1]),
            end_1based=int(fields[2]),
            repeat_unit=repeat_unit,
        ))

    # Close input files
    for input_file in input_files_to_close:
        input_file.close()

    if verbose:
        print(f"Parsed {len(tr_loci):,d} TR loci from {catalog_bed_path}")

    return tr_loci


def open_vcf_for_genotyping(vcf_path):
    """Open a VCF file for genotyping and validate it's single-sample.

    This function opens the VCF file, validates that it contains exactly one
    sample (as required for genotyping), and indexes the contig names in its header so that catalog loci
    can be matched to the spelling this VCF uses.

    Args:
        vcf_path (str): Path to the VCF file (can be .gz/.bgz compressed)

    Returns:
        tuple: (pysam.VariantFile, str, dict) containing:
            - The opened VCF file object
            - The sample name
            - The contig name lookup for this VCF (see build_contig_name_lookup)

    Raises:
        ValueError: If the VCF is not a single-sample VCF
        FileNotFoundError: If the VCF file doesn't exist
    """
    if not file_exists(vcf_path):
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    vcf_file = pysam.VariantFile(vcf_path)

    # Genotyping fetches variants per-locus, which requires a tabix/csi index. An unindexed
    # VCF makes every fetch raise "fetch requires an index", which get_overlapping_vcf_variants
    # would otherwise silently swallow and report every locus as homozygous-reference. Fail loudly.
    if vcf_file.index is None:
        raise ValueError(
            f"VCF file must be bgzip-compressed and tabix/csi-indexed for genotyping, but no index "
            f"was found for: {vcf_path}. Run e.g. `bgzip {vcf_path} && tabix -p vcf {vcf_path}.gz`")

    # Validate single-sample VCF
    sample_names = list(vcf_file.header.samples)
    if len(sample_names) == 0:
        raise ValueError(f"VCF file has no samples: {vcf_path}")
    if len(sample_names) > 1:
        raise ValueError(
            f"VCF file must be single-sample for genotyping, but found {len(sample_names)} samples: "
            f"{', '.join(sample_names[:5])}{'...' if len(sample_names) > 5 else ''}"
        )

    sample_name = sample_names[0]

    return vcf_file, sample_name, build_contig_name_lookup(vcf_file.header.contigs)


def get_called_alt_alleles(variant):
    """Return the alt alleles the sample actually carries, or None if its genotype leaves that unknown.

    Only the alleles named by GT matter. A record's other ALT alleles belong to other samples or to alleles
    this sample does not carry, so letting them decide whether the record affects a locus would keep records
    whose called allele changes nothing there. The star allele is skipped because it stands for a deletion
    described by its own separate record, and the reference allele changes nothing by definition.

    Args:
        variant (pysam.VariantRecord): the record to inspect

    Returns:
        list: the alt allele strings this sample carries, or None if the record has no GT or a missing
            allele, which leaves the haplotype unknown rather than unchanged
    """
    gt = variant.samples[0].get("GT")
    if gt is None:
        return None

    called_alts = []
    for allele_index in gt:
        if allele_index is None:
            return None
        if allele_index == 0 or allele_index >= len(variant.alleles):
            continue
        allele = variant.alleles[allele_index]
        if allele is None or allele == "*":
            continue
        called_alts.append(allele)

    return called_alts


def is_sequence_made_of_motif(sequence, repeat_unit):
    """Whether a sequence is one or more copies of a repeat unit, starting at any rotation of it.

    The last copy may be partial, so "AGCAGCA" is made of "CAG", while "CA" is not (less than one full copy).

    Args:
        sequence (str): the sequence to test, upper case
        repeat_unit (str): the repeat unit, upper case

    Returns:
        bool: True if the sequence is a prefix of some rotation of the repeat unit repeated
    """
    if not repeat_unit or len(sequence) < len(repeat_unit):
        return False

    num_copies = len(sequence) // len(repeat_unit) + 1
    return any(((repeat_unit[i:] + repeat_unit[:i]) * num_copies).startswith(sequence)
               for i in range(len(repeat_unit)))


def is_position_of_inserted_bases_inside_locus(position_of_inserted_bases, start_0based, end, inserted_bases=None,
                                               repeat_unit=None):
    """Whether bases inserted at the given position land inside the locus.

    The position of inserted bases is the end of the variant's suffix-trimmed reference span, so a left-aligned
    repeat-unit insertion just before the locus is positioned at start_0based and lands inside it, while one
    positioned at end sits in the right flank. Two kinds of insertion positioned at end belong to the locus all
    the same:

    - At a zero-width locus (start_0based == end), a repeat that catalog discovery found only in an insertion
      with no repeat bases in the reference, the insertion positioned at its single coordinate is the locus.
    - An insertion whose bases are copies of the locus motif. Catalog discovery ends an interrupted tract such
      as CAGCAGCAT at its last base, and an inserted CAG after the CAT cannot be left-aligned into the tract,
      so it is positioned exactly at end even though it is the expansion the locus was discovered from. Requiring
      the motif to match keeps an insertion of a different motif with the neighbouring locus it belongs to when
      two loci are adjacent, and leaves a non-repeat insertion right after the locus in the flank.

    Args:
        position_of_inserted_bases (int): 0-based genomic position of the inserted bases (the end of the
            variant's suffix-trimmed reference span)
        start_0based (int): locus start position (0-based, inclusive)
        end (int): locus end position (0-based, exclusive)
        inserted_bases (str): the inserted bases, upper case, or None to skip the motif test
        repeat_unit (str): the locus motif, upper case, or None to skip the motif test

    Returns:
        bool: True if the inserted bases fall inside the locus
    """
    if start_0based <= position_of_inserted_bases < end:
        return True

    if position_of_inserted_bases != end:
        return False

    return start_0based == end or (
        inserted_bases is not None and repeat_unit is not None
        and is_sequence_made_of_motif(inserted_bases, repeat_unit))


def does_alt_allele_change_bases_inside_locus(variant_pos_1based, ref, alt, start_0based, end, repeat_unit=None):
    """Whether applying one alt allele changes any base between the locus boundaries.

    A record can overlap a locus through its reference span while every base it actually changes lies outside.
    The common shape is an indel that left-alignment has anchored on the locus's last base: its anchor base is
    unchanged and all of its inserted or deleted bases sit in the flank. Counting such a record as overlapping
    inflates NumOverlappingVariants, keeps the locus out of --skip-hom-ref-loci, writes the flank record into
    the contributing-variants VCF, and can turn a genotypable locus into an ambiguous-phasing no call.

    Substituted and deleted bases sit at their own reference positions; inserted bases are positioned at the end
    of the suffix-trimmed reference span, the same convention compute_locus_start_and_end_offsets_in_haplotype()
    uses.

    Args:
        variant_pos_1based (int): the variant's 1-based position
        ref (str): the variant's reference allele
        alt (str): the alt allele to test
        start_0based (int): locus start position (0-based, inclusive)
        end (int): locus end position (0-based, exclusive)
        repeat_unit (str): the locus motif, which lets an insertion of that motif positioned exactly at end
            count as inside the locus (see is_position_of_inserted_bases_inside_locus), or None

    Returns:
        bool: True if this alt changes at least one base inside [start_0based, end)
    """
    trimmed_ref, trimmed_alt = trim_shared_suffix(ref.upper(), alt.upper())
    variant_start_0based = variant_pos_1based - 1

    for i, reference_base in enumerate(trimmed_ref):
        replacement = trimmed_alt[i] if i < len(trimmed_alt) else None
        if replacement != reference_base and start_0based <= variant_start_0based + i < end:
            return True

    if len(trimmed_alt) > len(trimmed_ref):
        if is_position_of_inserted_bases_inside_locus(variant_start_0based + len(trimmed_ref), start_0based, end,
                                                      trimmed_alt[len(trimmed_ref):], repeat_unit):
            return True

    return False


def does_variant_change_the_locus(variant, start_0based, end, repeat_unit=None):
    """Whether a record changes the locus, judged on the alleles the sample actually carries.

    Args:
        variant (pysam.VariantRecord): the record to test
        start_0based (int): locus start position (0-based, inclusive)
        end (int): locus end position (0-based, exclusive)
        repeat_unit (str): the locus motif (see does_alt_allele_change_bases_inside_locus), or None

    Returns:
        bool: True if the record changes a base inside the locus, or leaves a haplotype unknown there
    """
    called_alts = get_called_alt_alleles(variant)
    if called_alts is None:
        # The genotype does not say which allele this sample carries, so the record matters only if some
        # allele it offers could change the locus at all. A flank-only indel with a "." allele leaves the
        # locus untouched whichever allele turns out to be real, so it must not force a no call there.
        called_alts = [alt for alt in variant.alleles[1:] if alt is not None and alt != "*"]

    return any(does_alt_allele_change_bases_inside_locus(variant.pos, variant.ref, alt, start_0based, end, repeat_unit)
               for alt in called_alts)


def get_overlapping_vcf_variants(vcf_file, chrom, start_0based, end, repeat_unit=None):
    """Fetch VCF variants that overlap a genomic interval and can change a haplotype there.

    Args:
        vcf_file (pysam.VariantFile): An open pysam VariantFile object
        chrom (str): Chromosome name, spelled the way this VCF spells it (see build_contig_name_lookup)
        start_0based (int): Start position (0-based, inclusive)
        end (int): End position (0-based, exclusive / 1-based inclusive)
        repeat_unit (str): the locus motif (see does_alt_allele_change_bases_inside_locus), or None

    Returns:
        list: List of pysam.VariantRecord objects sorted by position

    Note:
        The tabix index keys each record by its REF span [POS, POS+len(REF)), so a plain
        fetch(start_0based, end) returns variants whose REF extends into the locus (e.g.
        deletions anchored in the left flank) but NOT insertions. A normalized (left-aligned)
        repeat-unit insertion anchors to the base immediately before the tract, i.e. its REF
        span is [start_0based - 1, start_0based), which does not overlap [start_0based, end)
        even though the inserted bases belong to the locus. To catch these, the fetch window
        is widened one base to the left and the returned records are filtered back down to
        those that actually affect the locus.
    """
    try:
        # pysam fetch uses 0-based half-open coordinates. Widen the window one base to the
        # left so left-anchored insertions (REF span ending exactly at start_0based) are seen.
        candidates = list(vcf_file.fetch(chrom, max(0, start_0based - 1), end))
    except ValueError:
        # Chromosome not found in VCF - return empty list
        return []

    # Keep only variants that actually affect the locus: those whose REF span overlaps
    # [start_0based, end), plus left-anchored insertions that put bases inside the locus even
    # though their REF span ends exactly at start_0based. A variant whose REF span ends at
    # start_0based without inserting anything inside the locus (a SNV on the last flank base, or
    # an insertion whose inserted bases fall back into the flank once the shared suffix is trimmed)
    # is dropped so it cannot perturb the locus, inflate the overlapping-variant count, or spuriously trigger the
    # multi-variant phasing-ambiguity check. Records genotyped as homozygous reference are dropped for the
    # same reasons, as are records whose called allele changes nothing inside the locus even though its
    # reference span reaches into it (see does_variant_change_the_locus).
    variants = [variant for variant in candidates
                if does_variant_change_the_locus(variant, start_0based, end, repeat_unit)]

    # Sort by position
    variants.sort(key=lambda v: v.pos)

    return variants


def trim_shared_suffix(ref, alt):
    """Trim the shared trailing bases of a variant's ref and alt alleles.

    For example, ref="ATCAG" alt="TTTCAG" share the suffix "CAG", so this returns
    ("AT", "TTT"). Trimming is only applied when both alleles are longer than 1 bp,
    matching the convention used when applying variants to a haplotype sequence.

    Args:
        ref (str): Reference allele
        alt (str): Alternate allele

    Returns:
        tuple: (trimmed_ref, trimmed_alt)
    """
    if len(ref) > 1 and len(alt) > 1:
        common_suffix_length = 0
        for j in range(1, min(len(ref), len(alt)) + 1):
            if ref[-j] == alt[-j]:
                common_suffix_length += 1
            else:
                break
        if common_suffix_length > 0:
            return ref[:-common_suffix_length], alt[:-common_suffix_length]
    return ref, alt


def get_inserted_bases_inside_locus(variant_pos_1based, ref, alt, start_0based, end, repeat_unit=None):
    """Return the bases an alt allele inserts between the locus boundaries, or None if it inserts none there.

    Inserted bases are positioned at the end of the variant's suffix-trimmed reference span, the same convention
    compute_locus_start_and_end_offsets_in_haplotype() uses to decide which bases fall inside the locus. A
    variant can overlap the locus through its reference span while inserting its bases outside it, and trimming
    a shared suffix can pull the position of the inserted bases back into the left flank, so that position
    rather than the variant position decides.

    Args:
        variant_pos_1based (int): the variant's 1-based position
        ref (str): the variant's reference allele
        alt (str): the alt allele to check
        start_0based (int): locus start position (0-based, inclusive)
        end (int): locus end position (0-based, exclusive / 1-based inclusive)
        repeat_unit (str): the locus motif (see is_position_of_inserted_bases_inside_locus), or None

    Returns:
        str: the inserted bases that land inside the locus, or None if this alt inserts none there
    """
    trimmed_ref, trimmed_alt = trim_shared_suffix(ref.upper(), alt.upper())
    inserted_base_count = len(trimmed_alt) - len(trimmed_ref)
    if inserted_base_count <= 0:
        return None

    inserted_bases = trimmed_alt[-inserted_base_count:]
    if not is_position_of_inserted_bases_inside_locus(variant_pos_1based - 1 + len(trimmed_ref), start_0based,
                                                      end, inserted_bases, repeat_unit):
        return None

    return inserted_bases


def convert_variants_to_haplotype_sequence(pos_1based, reference_sequence, variants, verbose=False):
    """Convert a reference sequence and overlapping variants into an alternate haplotype sequence.

    This function takes a reference sequence at a genomic locus and a list of one or more
    variants that overlap this sequence, and returns the alternate haplotype sequence with
    all variants applied.

    For example, if the reference sequence at position 100 is "ACAGCAG" and the variants are:
        [(103, "G", "A"), (106, "G", "AT")]
    then the output sequence would be "ACAACAAT".

    Args:
        pos_1based (int): The 1-based genomic position of the start of the reference sequence
        reference_sequence (str): The reference sequence that these variants overlap
        variants (list): List of (pos_1based, ref, alt) tuples representing the variants
            to apply. Variants must be sorted by position and not overlap each other.
        verbose (bool): If True, print detailed output for debugging

    Returns:
        str: The alternate haplotype sequence with all variants applied

    Raises:
        NonIUPACAlleleError: If a ref or alt allele contains characters that aren't IUPAC nucleotide codes
            (a subclass of ValueError)
        RefAlleleMismatchError: If a variant's ref allele doesn't match the reference sequence at the expected
            position (a subclass of ValueError)
        ValueError: If variants are out of order or overlap, or if a variant extends beyond the reference
            sequence
    """
    # Calculate expected output length for validation
    expected_output_length = len(reference_sequence) + sum(len(alt) - len(ref) for _, ref, alt in variants)

    reference_sequence = reference_sequence.upper()
    output_sequence = ""

    if verbose:
        output_sequence_with_separators = ""
        print(f"Processing sequence at position {pos_1based:,d}: {reference_sequence} with {len(variants)} variants")

    # Iterate through variants, building the output sequence
    offset = 0
    next_offset = 0

    for i, (variant_pos_1based, variant_ref, variant_alt) in enumerate(variants):
        variant_ref = variant_ref.upper()
        variant_alt = variant_alt.upper()

        # Validate alleles contain only IUPAC nucleotide codes. The offending characters are sorted so the message,
        # which ends up in the NoCallReason output column, is the same from run to run.
        if set(variant_ref) - DNA_BASES:
            raise NonIUPACAlleleError(
                f"Invalid ref allele '{variant_ref}' at position {variant_pos_1based:,d}: contains non-IUPAC "
                f"characters {''.join(sorted(set(variant_ref) - DNA_BASES))}")
        if set(variant_alt) - DNA_BASES:
            raise NonIUPACAlleleError(
                f"Invalid alt allele '{variant_alt}' at position {variant_pos_1based:,d}: contains non-IUPAC "
                f"characters {''.join(sorted(set(variant_alt) - DNA_BASES))}")

        # Check that variant position is not before current position (variants must be in order)
        if pos_1based + next_offset > variant_pos_1based:
            raise ValueError(f"Variant at position {variant_pos_1based:,d} is before or overlaps the current "
                           f"position {pos_1based + next_offset:,d}. Variants must be sorted and non-overlapping.")

        # Add reference sequence up to the current variant
        next_offset = variant_pos_1based - pos_1based
        if next_offset > len(reference_sequence):
            raise ValueError(f"Variant at position {variant_pos_1based:,d} is beyond the end of the reference "
                           f"sequence interval ({pos_1based:,d}-{pos_1based - 1 + len(reference_sequence):,d})")

        output_sequence += reference_sequence[offset:next_offset]
        if verbose:
            print(f"Adding ref sequence [{offset}:{next_offset}] between variants: "
                  f"'{reference_sequence[offset:next_offset]}'")
            output_sequence_with_separators += reference_sequence[offset:next_offset]

        # Validate that ref allele matches reference sequence at this position
        if not reference_sequence[next_offset:].startswith(variant_ref):
            actual_ref = reference_sequence[next_offset:next_offset + len(variant_ref)]
            raise RefAlleleMismatchError(
                f"Reference sequence at position {variant_pos_1based:,d} does not match variant ref: "
                f"found '{actual_ref}' but variant has '{variant_ref}'")

        # Handle shared suffix trimming between ref and alt alleles
        # This handles cases like ref="CAGCAG" alt="CAG" where they share a "CAG" suffix
        trimmed_ref, trimmed_alt = trim_shared_suffix(variant_ref, variant_alt)
        if verbose and (trimmed_ref != variant_ref or trimmed_alt != variant_alt):
            print(f"Removing shared suffix from ref and alt alleles:")
            print(f"  ref: {variant_ref} -> {trimmed_ref}")
            print(f"  alt: {variant_alt} -> {trimmed_alt}")
        variant_ref, variant_alt = trimmed_ref, trimmed_alt

        # Apply the variant: add alt allele to output
        output_sequence += variant_alt
        next_offset = next_offset + len(variant_ref)

        if verbose:
            output_sequence_with_separators += f"|var{i + 1}:{variant_alt}|"
            print(f"Variant #{i + 1}: pos={variant_pos_1based:,d} {variant_ref} -> {variant_alt} : "
                  f"{output_sequence_with_separators}  (offset: {offset}, next_offset: {next_offset})")

        offset = next_offset

    # Add remaining reference sequence after the last variant
    output_sequence += reference_sequence[offset:]

    # Validate output length matches expected
    if len(output_sequence) != expected_output_length:
        raise ValueError(f"Output sequence length ({len(output_sequence)}) does not match expected length "
                        f"({expected_output_length}). This indicates a bug in the variant processing logic.")

    return output_sequence


# Locus motif and thresholds used to decide whether inserted bases represent tandem repeats.
InsertionFilter = collections.namedtuple("InsertionFilter", [
    "motif",
    "min_insertion_size_to_check",
    "min_insertion_purity",
    "min_insertion_periodicity",
])


def build_insertion_filter(motif, args):
    """Build an InsertionFilter for a locus.

    Args:
        motif (str): the locus repeat motif
        args (argparse.Namespace): command-line arguments parsed by parse_args()

    Returns:
        InsertionFilter: the thresholds to apply at this locus
    """
    return InsertionFilter(
        motif=motif,
        min_insertion_size_to_check=getattr(args, "min_insertion_size_to_check_for_repetitiveness",
                                            DEFAULT_MIN_INSERTION_SIZE_TO_CHECK),
        min_insertion_purity=getattr(args, "min_insertion_purity", DEFAULT_MIN_INSERTION_PURITY),
        min_insertion_periodicity=getattr(args, "min_insertion_periodicity", DEFAULT_MIN_INSERTION_PERIODICITY),
    )


def check_if_inserted_sequence_is_repetitive(inserted_sequence, insertion_filter):
    """Decide whether inserted bases at a tandem repeat locus are part of that locus tandem repeat.

    An inserted sequence is accepted if any of the following is true:
      1. it is shorter than insertion_filter.min_insertion_size_to_check, since short insertions can't inflate a
         repeat count by much and are too short to tell apart from random sequence
      2. enough of its bases match a pure repeat of the locus motif, at the best starting offset within the motif.
         The threshold is high (see DEFAULT_MIN_INSERTION_PURITY) because purity is maximized over every rotation of
         the motif, so a sequence barely one copy long finds some rotation that fits it. Requiring the insertion to
         span two motif copies instead was measured and rejected: it lowers the non-repeat base pairs accepted by
         less than the measurement's own uncertainty, and it discards single-copy VNTR insertions outright, since
         one copy of a non-repetitive unit has near-perfect purity but no internal periodicity to fall back on.
         An insertion shorter than one motif copy is scored against the best-matching stretch of the motif
         instead, since comparing it to whole rotations yields nan and would leave only the periodicity test,
         which a fragment of a non-repetitive VNTR unit cannot pass
      3. it looks like a tandem repeat of some other motif, which keeps expansions where the expanded motif
         differs from the one annotated for the locus (eg. an AAGGG expansion at the AAAAG repeat in RFC1) as well
         as expansions of large, highly degenerate VNTR units

    Args:
        inserted_sequence (str): the inserted bases
        insertion_filter (InsertionFilter): locus motif and thresholds

    Returns:
        2-tuple (bool, str):
            bool: True if the inserted bases are part of the tandem repeat
            str: the reason the bases were rejected, or None if they were accepted
    """
    if len(inserted_sequence) < insertion_filter.min_insertion_size_to_check:
        return True, None

    if "N" in inserted_sequence.upper():
        return False, INSERTION_FILTER_REASON_CONTAINS_NS

    purity, _, _ = compute_best_phase_repeat_purity(inserted_sequence, insertion_filter.motif)
    if purity != purity:  # nan, ie. the insertion is shorter than one copy of the motif
        purity = compute_partial_copy_purity(inserted_sequence, insertion_filter.motif)
    if purity == purity and purity >= insertion_filter.min_insertion_purity:
        return True, None

    max_period = min(500, max(100, 3 * len(insertion_filter.motif)))
    _, periodicity = compute_sequence_periodicity(inserted_sequence, max_period=max_period)
    if periodicity >= insertion_filter.min_insertion_periodicity:
        return True, None

    return False, INSERTION_FILTER_REASON_NOT_REPEAT_LIKE


def find_insufficiently_repetitive_insertions_inside_locus(variant_list, start_0based, end, insertion_filter,
                                                           verbose=False):
    """Find the insertions in a haplotype's variants that aren't part of the locus tandem repeat.

    An insertion that isn't a tandem repeat (an Alu element dropped into a poly-A tract, an assembly error)
    makes the whole haplotype unusable as a repeat measurement, since there is no way to say how many repeats
    the sample really carries. Callers set the affected allele to a no-call rather than reporting a count.

    Only insertions whose bases land between the locus boundaries are judged. A variant can overlap the locus
    through its reference span while inserting its bases outside it, and such an insertion says nothing about
    the repeat. Inserted bases are positioned at the end of the variant's reference span, the same convention
    compute_locus_start_and_end_offsets_in_haplotype() uses to decide which bases fall inside the locus.

    Args:
        variant_list (list): list of (pos_1based, ref, alt) tuples for one haplotype
        start_0based (int): locus start position (0-based, inclusive)
        end (int): locus end position (0-based, exclusive / 1-based inclusive)
        insertion_filter (InsertionFilter): locus motif and thresholds
        verbose (bool): if True, print each rejected insertion

    Returns:
        list: (inserted_sequence, reason) for each insertion inside the locus that isn't part of the tandem
            repeat. Empty if every such insertion belongs to the repeat.
    """
    rejected_insertions = []
    for variant_pos_1based, variant_ref, variant_alt in variant_list:
        inserted_sequence = get_inserted_bases_inside_locus(
            variant_pos_1based, variant_ref, variant_alt, start_0based, end, insertion_filter.motif)
        if inserted_sequence is None:
            continue

        is_insertion_sufficiently_repetitive, reason = check_if_inserted_sequence_is_repetitive(
            inserted_sequence, insertion_filter)
        if not is_insertion_sufficiently_repetitive:
            rejected_insertions.append((inserted_sequence, reason))
            if verbose:
                print(f"  The {len(inserted_sequence):,d} bases inserted at position {variant_pos_1based:,d} "
                      f"are not sufficiently repetitive ({reason}), so this allele has no call")

    return rejected_insertions


# Per-haplotype result of building the sequence of a tandem repeat locus from VCF variants.
HaplotypeSequence = collections.namedtuple("HaplotypeSequence", [
    "sequence",             # locus sequence with all variants applied, or None if this allele has no call
    "rejected_insertions",  # (inserted_sequence, reason) for each insertion that isn't part of the repeat
    "build_error",          # why the sequence couldn't be built, or None if it was
    "build_error_reason",   # which NO_CALL_REASON_* the build_error belongs to, or None for the generic one
], defaults=(None, None))

MISSING_HAPLOTYPE_SEQUENCES = HaplotypeSequence(None, ())


def compute_locus_start_and_end_offsets_in_haplotype(locus_start_0based, locus_end_0based,
                                                     fetched_reference_sequence_start, variant_list,
                                                     repeat_unit=None):
    """Find where the locus starts and ends within a haplotype sequence built from the reference plus variants.

    The haplotype sequence is the one built by convert_variants_to_haplotype_sequence() from the reference
    starting at fetched_reference_sequence_start plus the given variants. Each variant's (suffix-trimmed) alt
    bases are assigned genomic positions: alt base i aligns to variant_start + i for i < len(ref) (matched,
    substituted or deleted bases) and to variant_start + len(ref) for any extra inserted bases (i >= len(ref)).
    An alt base falls on the near side of a boundary iff its assigned position is strictly before it. This
    places a left-anchored insertion's inserted bases (positioned at the variant's ref-span end, == the locus
    start) inside the locus, while keeping a boundary-spanning deletion's surviving bases on the correct side.

    At the locus end, inserted bases positioned exactly at the end also count as inside the locus when
    is_position_of_inserted_bases_inside_locus() says they belong to it (a zero-width locus, or an insertion of
    the locus motif), matching the overlap test that admitted the variant.

    Args:
        locus_start_0based (int): locus start position (0-based, inclusive)
        locus_end_0based (int): locus end position (0-based, exclusive)
        fetched_reference_sequence_start (int): 0-based genomic start of the reference sequence the haplotype
            was built from
        variant_list (list): list of (pos_1based, ref, alt) tuples that were applied, sorted by position
        repeat_unit (str): the locus motif, used for the locus-end insertion test, or None

    Returns:
        2-tuple (int, int): the offsets of the locus start and the locus end within the haplotype sequence, so
            that haplotype_sequence[start_offset:end_offset] is the locus
    """
    def compute_offset(locus_boundary_0based, is_locus_end):
        output_pos = 0
        genomic_pos = fetched_reference_sequence_start  # next reference position not yet consumed
        for variant_pos_1based, variant_ref, variant_alt in variant_list:
            # Upper-case before trimming, the same way convert_variants_to_haplotype_sequence() does, so that the
            # offsets computed here line up with the sequence it built. VCFs written against a soft-masked
            # reference (eg. DipCall output) give the ref allele in lower case and the alt allele in upper case,
            # so a case-sensitive comparison would find no shared suffix where there is one.
            variant_ref, variant_alt = trim_shared_suffix(variant_ref.upper(), variant_alt.upper())
            variant_start_0based = variant_pos_1based - 1
            ref_len = len(variant_ref)
            counts_insertion_at_boundary = (
                is_locus_end and len(variant_alt) > ref_len
                and variant_start_0based + ref_len == locus_boundary_0based
                and is_position_of_inserted_bases_inside_locus(locus_boundary_0based, locus_start_0based,
                                                               locus_end_0based, variant_alt[ref_len:], repeat_unit))
            if locus_boundary_0based <= genomic_pos and not counts_insertion_at_boundary:
                return output_pos
            # Reference bases between the current position and this variant, up to the boundary
            ref_run_end = min(variant_start_0based, locus_boundary_0based)
            if ref_run_end > genomic_pos:
                output_pos += ref_run_end - genomic_pos
            if locus_boundary_0based <= variant_start_0based and not counts_insertion_at_boundary:
                return output_pos
            # The variant starts before the boundary: count only the alt bases whose assigned genomic position
            # lies before it
            for i in range(len(variant_alt)):
                if i >= ref_len and counts_insertion_at_boundary:
                    output_pos += 1
                elif variant_start_0based + (i if i < ref_len else ref_len) < locus_boundary_0based:
                    output_pos += 1
            genomic_pos = max(genomic_pos, variant_start_0based + ref_len)
        # Trailing reference bases after the last variant
        if locus_boundary_0based > genomic_pos:
            output_pos += locus_boundary_0based - genomic_pos
        return output_pos

    return compute_offset(locus_start_0based, is_locus_end=False), compute_offset(locus_end_0based, is_locus_end=True)


def get_haplotype_variant_list(vcf_variants, haplotype, verbose=False):
    """Collect the variants that apply to one haplotype of a single-sample VCF.

    Args:
        vcf_variants (list): list of pysam.VariantRecord objects that overlap the locus
        haplotype (int): 0 for the first haplotype, 1 for the second
        verbose (bool): if True, print detailed output for debugging

    Returns:
        list: (pos_1based, ref, alt) tuples for the alt alleles this haplotype carries, upper-cased so that
            downstream comparisons don't depend on the soft-masking of the reference the VCF was written
            against, or None if the genotype is missing for this haplotype. Star alleles are skipped, since the
            deletion they refer to is represented by its own VCF record.
    """
    variant_list = []
    for variant in vcf_variants:
        gt = variant.samples[0].get("GT")
        if gt is None:
            if verbose:
                print(f"Haplotype {haplotype}: Missing GT for variant at {variant.pos}")
            return None

        # Haploid genotypes (e.g. a chrX/chrY call with GT == (1,)) have no second haplotype, so treat the
        # absent haplotype as missing and let the existing HEMI logic handle it instead of indexing out of range.
        if haplotype >= len(gt):
            if verbose:
                print(f"Haplotype {haplotype}: Haploid genotype at {variant.pos}; no allele for this haplotype")
            return None

        gt_value = gt[haplotype]

        # Handle missing genotype (.)
        if gt_value is None:
            if verbose:
                print(f"Haplotype {haplotype}: Missing genotype (.) at {variant.pos}")
            return None

        # Skip reference allele (0) - no change to reference sequence
        if gt_value == 0:
            continue

        alt_allele = variant.alleles[gt_value]

        # Skip star alleles - they indicate overlap with a previously specified deletion. The actual deletion is
        # represented by a separate VCF record which pysam will fetch if it overlaps this locus.
        if alt_allele == "*":
            if verbose:
                print(f"Haplotype {haplotype}: Skipping star allele at {variant.pos}")
            continue

        if verbose:
            print(f"Haplotype {haplotype}: Variant at {variant.pos} {variant.ref} -> {alt_allele}")

        variant_list.append((variant.pos, variant.ref.upper(), alt_allele.upper()))

    return variant_list


def is_heterozygous_genotype(gt):
    """Check whether a genotype carries two different alleles, and so needs phasing to be resolved.

    A genotype that names only one distinct allele puts that allele on every haplotype no matter how the record
    is phased. That covers homozygous calls (0/0, 1/1), haploid calls (1), and calls where one allele is missing
    (./1), none of which need a phase.

    Args:
        gt (tuple): the GT tuple from a pysam sample record, or None if the record has no GT

    Returns:
        bool: True if the genotype names more than one distinct allele
    """
    if gt is None:
        return False

    return len({allele for allele in gt if allele is not None}) > 1


def build_contig_name_lookup(contig_names):
    """Index a file's contig names by their normalized form so a catalog contig can be resolved to whatever
    spelling that file uses.

    A "chr1" catalog has to line up with a VCF or reference that calls the same contig "1", and a chrM locus
    with an MT-named one. normalize_chromosome_name() collapses both of those differences, so indexing by it
    turns per-locus resolution into a single dict lookup. Building the index once matters because a catalog
    holds millions of loci but only a couple of dozen distinct contigs.

    Args:
        contig_names (iterable): the contig names a VCF header or reference FASTA declares

    Returns:
        dict: normalized contig name (see normalize_chromosome_name) to the spelling the file uses
    """
    return {normalize_chromosome_name(contig_name): contig_name for contig_name in contig_names}


def find_records_with_a_missing_genotype_allele(vcf_variants, require_every_allele_missing=False):
    """Return the 1-based positions of records that leave a haplotype uncalled at a diploid locus.

    A record with no GT field at all, or a GT carrying a "." allele (".|1", "1|.", "./.", or a haploid "."),
    does not say what one or both haplotypes are. DipCall emits these in bulk, including ".|." with FILTER
    GAP1;GAP2 for sites uncalled on both haplotypes.

    Args:
        vcf_variants (list): pysam.VariantRecord objects that overlap the locus
        require_every_allele_missing (bool): when True, only report records where no allele is called. This is
            what a haploid locus needs: there ".|1" is the expected shape of a haploid call and the called
            allele is used, while ".|." still says nothing at all.

    Returns:
        list: the 1-based positions of the records with a missing allele, in the order given
    """
    positions = []
    for variant in vcf_variants:
        gt = variant.samples[0].get("GT")
        if gt is None or len(gt) == 0:
            positions.append(variant.pos)
            continue
        missing_alleles = [allele for allele in gt if allele is None]
        if not missing_alleles:
            continue
        if require_every_allele_missing and len(missing_alleles) < len(gt):
            continue
        positions.append(variant.pos)

    return positions


def find_records_with_a_haploid_genotype(vcf_variants):
    """Return the 1-based positions of records whose GT names a single allele, such as "1" or "0".

    Callers run with ploidy 1 write these on a male's non-PAR chrX and chrY. At a diploid locus such a record
    says nothing about the second haplotype, so it is as incomplete as "1|.", even though no allele is ".".

    Args:
        vcf_variants (list): pysam.VariantRecord objects that overlap the locus

    Returns:
        list: the 1-based positions of the records with a single-allele GT, in the order given
    """
    return [variant.pos for variant in vcf_variants if len(variant.samples[0].get("GT") or ()) == 1]


def are_variants_unambiguously_phased(vcf_variants, verbose=False):
    """Check whether a set of variants overlapping one locus can be split into two haplotypes without guessing.

    Only heterozygous variants need a phase: a homozygous or hom-ref record contributes the same allele to both
    haplotypes however it is written, and a single heterozygous variant can go on either haplotype, which yields
    the same pair of allele sequences either way. So the genotype is ambiguous only when two or more
    heterozygous variants overlap the locus and the VCF does not say which haplotype each one is on.

    Two heterozygous variants can be placed on haplotypes relative to each other only if both are phased AND
    both belong to the same phase set. Phasing tools such as whatshap and HiPhase split a chromosome into independent phase
    blocks and tag each with a PS value, and "0|1" in one block says nothing about which physical haplotype
    "0|1" refers to in another. Assembly-based callers such as dipcall phase the whole chromosome at once and
    emit no PS at all, which reads here as a single shared phase set.

    Records that came from splitting one multiallelic site are the exception. `bcftools norm -m -any` rewrites
    a single "1/2" record as two records at the same position, genotyped "1/0" and "0/1". Those two alt alleles
    came from one diploid genotype, so they are necessarily on opposite haplotypes no matter that each record
    is written unphased. Counting them as two independently unphased heterozygous variants would make a
    normalized callset a no call at every multiallelic tandem repeat, while the same sample genotypes fine
    before normalization, so records sharing a position are counted once here.

    The position alone is the key, not the position and the REF allele. Run with -f, which is the normal
    usage, bcftools minimizes each split allele separately, so a deletion pair from one "1/2" record comes out
    with different REF strings at the same POS. Keying on REF as well would miss exactly the shapes
    normalization produces at repeat loci. In a single-sample VCF two heterozygous records at one position can
    only have come from splitting one multiallelic record, so the position is enough.

    Args:
        vcf_variants (list): list of pysam.VariantRecord objects that overlap the locus
        verbose (bool): if True, print the reason the genotype is ambiguous

    Returns:
        bool: True if the variants can be assigned to haplotypes unambiguously
    """
    heterozygous_variants = [
        variant for variant in vcf_variants if is_heterozygous_genotype(variant.samples[0].get("GT"))
    ]
    if len({(variant.chrom, variant.pos) for variant in heterozygous_variants}) < 2:
        return True

    for variant in heterozygous_variants:
        # pysam returns phased=False if GT uses the '/' separator
        if not variant.samples[0].phased:
            if verbose:
                print(f"Multiple heterozygous variants overlap locus and the GT at {variant.pos} is unphased: "
                      f"returning missing genotype")
            return False

    phase_sets = {variant.samples[0].get("PS") for variant in heterozygous_variants}
    if len(phase_sets) > 1:
        if verbose:
            print(f"Multiple heterozygous variants overlap locus but belong to different phase sets "
                  f"({', '.join(str(phase_set) for phase_set in sorted(phase_sets, key=str))}): "
                  f"returning missing genotype")
        return False

    return True


def extract_haplotype_sequences_and_insertions_from_vcf(chrom, start_0based, end, fasta_obj, vcf_variants,
                                                        verbose=False, insertion_filter=None,
                                                        fasta_contig_lookup=None, repeat_unit=None):
    """Extract diploid haplotype sequences for a genomic locus from VCF variants, optionally rejecting a
    haplotype that contains an insertion which isn't part of the locus tandem repeat.

    This is the implementation behind extract_haplotype_sequences_from_vcf(). It additionally reports, for each
    haplotype, which insertions were judged not to be part of the repeat. A haplotype with any such insertion
    has no sequence, since there is no way to say how many repeats it carries.

    Args:
        chrom (str): Chromosome name
        start_0based (int): Locus start position (0-based, inclusive)
        end (int): Locus end position (0-based, exclusive / 1-based inclusive)
        fasta_obj (pyfaidx.Fasta): Reference genome fasta object (with one_based_attributes=False)
        vcf_variants (list): List of pysam.VariantRecord objects that overlap this locus
        verbose (bool): If True, print detailed output for debugging
        insertion_filter (InsertionFilter): thresholds for deciding whether inserted bases belong to the repeat,
            or None to accept every inserted base. The genotype subcommand always passes a filter (see
            build_insertion_filter); None is only for direct callers such as extract_haplotype_sequences_from_vcf.
        fasta_contig_lookup (dict): the reference's contig names indexed by build_contig_name_lookup, so a
            "chr1" catalog lines up with a "1"-named reference and a chrM locus with an MT-named one. When
            None, chrom is used as given.
        repeat_unit (str): the locus motif, which lets an insertion of that motif positioned exactly at end
            count as part of the locus (see is_position_of_inserted_bases_inside_locus), or None

    Returns:
        tuple: (HaplotypeSequence, HaplotypeSequence) for haplotype 0 and haplotype 1
    """
    # Check for phasing ambiguity first
    if not are_variants_unambiguously_phased(vcf_variants, verbose=verbose):
        return (MISSING_HAPLOTYPE_SEQUENCES, MISSING_HAPLOTYPE_SEQUENCES)

    # Determine the actual region we need to fetch from the reference
    # This may be larger than the locus if variants extend beyond its boundaries
    fetched_reference_sequence_start = start_0based
    fetched_reference_sequence_end = end

    for variant in vcf_variants:
        # variant.pos is 1-based, convert to 0-based
        variant_start_0based = variant.pos - 1
        variant_end_0based = variant_start_0based + len(variant.ref)

        # Expand fetch region if variant extends beyond locus
        if variant_start_0based < fetched_reference_sequence_start:
            fetched_reference_sequence_start = variant_start_0based
        if variant_end_0based > fetched_reference_sequence_end:
            fetched_reference_sequence_end = variant_end_0based

    # Fetch the (potentially expanded) reference sequence under whatever name this reference gives the
    # contig. A contig the fasta simply does not have is a build error rather than a silent missing
    # genotype, so the locus says why it got no call.
    fasta_chrom = fasta_contig_lookup.get(normalize_chromosome_name(chrom), chrom) if fasta_contig_lookup \
        else chrom
    try:
        reference_sequence = str(
            fasta_obj[fasta_chrom][fetched_reference_sequence_start:fetched_reference_sequence_end]).upper()
    except KeyError:
        reference_sequence = None
    if reference_sequence is None:
        build_error = f"{chrom} is not present under any of its usual names"
        if verbose:
            print(f"  {NO_CALL_REASON_CONTIG_NOT_IN_REFERENCE}: {build_error}")
        missing = HaplotypeSequence(None, (), build_error, NO_CALL_REASON_CONTIG_NOT_IN_REFERENCE)
        return (missing, missing)

    # pyfaidx clips a slice that runs past the end of a contig instead of raising, so a catalog locus whose
    # end exceeds the contig length would otherwise yield an allele shorter than the locus span, with no
    # warning. Its repeat count would come from the truncated sequence while NumRepeatsInReference came from
    # the full catalog span, which reads as a real contraction.
    if len(reference_sequence) != fetched_reference_sequence_end - fetched_reference_sequence_start:
        build_error = (f"{chrom}:{fetched_reference_sequence_start}-{fetched_reference_sequence_end} returned "
                       f"only {len(reference_sequence):,d}bp of reference")
        if verbose:
            print(f"  {NO_CALL_REASON_LOCUS_PAST_CONTIG_END}: {build_error}")
        clipped = HaplotypeSequence(None, (), build_error, NO_CALL_REASON_LOCUS_PAST_CONTIG_END)
        return (clipped, clipped)

    if verbose:
        print(f"Extracting haplotypes for {chrom}:{start_0based}-{end}")
        print(f"Fetched reference region {chrom}:{fetched_reference_sequence_start}-"
              f"{fetched_reference_sequence_end}: {reference_sequence}")
        print(f"Processing {len(vcf_variants)} overlapping variants")

    # Build haplotype sequences for both haplotypes (0 and 1)
    haplotype_results = []

    for haplotype in (0, 1):
        variant_list = get_haplotype_variant_list(vcf_variants, haplotype, verbose=verbose)
        if variant_list is None:
            haplotype_results.append(MISSING_HAPLOTYPE_SEQUENCES)
            continue

        # If no variants affect this haplotype, use reference sequence
        if not variant_list:
            # Return just the locus portion, not the expanded fetch region
            haplotype_seq = reference_sequence[
                start_0based - fetched_reference_sequence_start:end - fetched_reference_sequence_start]
            haplotype_results.append(HaplotypeSequence(haplotype_seq, ()))
            continue

        if insertion_filter is not None:
            rejected_insertions = find_insufficiently_repetitive_insertions_inside_locus(
                variant_list, start_0based, end, insertion_filter, verbose=verbose)
            if rejected_insertions:
                haplotype_results.append(HaplotypeSequence(None, tuple(rejected_insertions)))
                continue

        haplotype_seq, build_error, build_error_reason = build_locus_sequence_from_variants(
            fetched_reference_sequence_start, reference_sequence, variant_list, start_0based, end, verbose=verbose,
            repeat_unit=repeat_unit)
        if haplotype_seq is None:
            haplotype_results.append(HaplotypeSequence(None, (), build_error, build_error_reason))
            continue

        haplotype_results.append(HaplotypeSequence(haplotype_seq, ()))

    return tuple(haplotype_results)


def build_locus_sequence_from_variants(fetched_reference_sequence_start, reference_sequence, variant_list,
                                       start_0based, end, verbose=False, repeat_unit=None):
    """Apply a haplotype's variants to the fetched reference sequence and trim the result to the locus.

    Args:
        fetched_reference_sequence_start (int): 0-based genomic start of reference_sequence
        reference_sequence (str): reference sequence covering the locus and any variants that extend beyond it
        variant_list (list): (pos_1based, ref, alt) tuples for this haplotype, sorted by position
        start_0based (int): locus start position (0-based, inclusive)
        end (int): locus end position (0-based, exclusive)
        verbose (bool): if True, print detailed output for debugging
        repeat_unit (str): the locus motif (see is_position_of_inserted_bases_inside_locus), or None

    Returns:
        3-tuple (str, str, str):
            str: the locus sequence for this haplotype, or None if the variants could not be applied
            str: why the variants could not be applied, or None if they applied cleanly
            str: the NO_CALL_REASON_* for that failure when it has a specific one (a non-IUPAC allele, or a REF
                allele that doesn't match the reference), or None for the generic
                NO_CALL_REASON_HAPLOTYPE_BUILD_ERROR or when there was no failure
    """
    try:
        # pos_1based for convert_variants_to_haplotype_sequence is fetched_reference_sequence_start + 1
        full_haplotype_seq = convert_variants_to_haplotype_sequence(
            fetched_reference_sequence_start + 1, reference_sequence, variant_list, verbose=verbose
        )
    except ValueError as e:
        if verbose:
            print(f"Error converting variants: {e}")
        if isinstance(e, NonIUPACAlleleError):
            build_error_reason = NO_CALL_REASON_NON_IUPAC_ALLELE
        elif isinstance(e, RefAlleleMismatchError):
            build_error_reason = NO_CALL_REASON_REF_ALLELE_MISMATCH
        else:
            build_error_reason = None
        return None, str(e), build_error_reason

    # Trim the haplotype sequence back to the original locus boundaries. We may have built the haplotype over an
    # expanded fetch region to cover variants that extend beyond the locus, so map the genomic locus interval
    # [start_0based, end) to offsets in the output sequence. An insertion positioned exactly at end that
    # is_position_of_inserted_bases_inside_locus() assigns to the locus (a zero-width locus, or an insertion of the
    # locus motif) must fall inside the slice, the same way the overlap test that admitted it treats it.
    locus_start_offset, locus_end_offset = compute_locus_start_and_end_offsets_in_haplotype(
        start_0based, end, fetched_reference_sequence_start, variant_list, repeat_unit=repeat_unit)
    return full_haplotype_seq[locus_start_offset:locus_end_offset], None, None


def extract_haplotype_sequences_from_vcf(chrom, start_0based, end, fasta_obj, vcf_variants, verbose=False):
    """Extract diploid haplotype sequences for a genomic locus from VCF variants.

    This function takes the coordinates of a TR locus, a reference genome, and a list of
    VCF variant records that overlap the locus. It returns two haplotype sequences derived
    by applying the variants to the reference sequence based on their diploid genotypes.

    For example, if the reference sequence at chr1:100 is "ACAGCAG" and the VCF variants are:
        CHROM   POS    REF   ALT   FORMAT    sample
        chr1    103    G     A     GT        0|1
        chr1    106    G     AT    GT        1|0

    then the two output haplotype sequences would be:
        ("ACAGCAAT", "ACAACAG")

    Args:
        chrom (str): Chromosome name
        start_0based (int): Locus start position (0-based, inclusive)
        end (int): Locus end position (0-based, exclusive / 1-based inclusive)
        fasta_obj (pyfaidx.Fasta): Reference genome fasta object (with one_based_attributes=False)
        vcf_variants (list): List of pysam.VariantRecord objects that overlap this locus
        verbose (bool): If True, print detailed output for debugging

    Returns:
        tuple: (haplotype0_seq, haplotype1_seq) where each is:
            - A string containing the haplotype sequence, or
            - None if the genotype is missing (e.g., due to unresolved phasing between
              multiple heterozygous variants, or a missing GT field)

    Notes:
        - If two or more heterozygous variants overlap and the VCF does not say which haplotype
          each one is on (any of them unphased, or they span more than one phase set), returns
          (None, None) to indicate missing genotype. Homozygous records never need a phase, and a single
          heterozygous variant gives the same pair of alleles whichever haplotype it goes on.
        - Variants whose REF extends beyond locus boundaries are handled by fetching
          an expanded reference region and trimming the result
        - Star alleles ('*') are skipped as they indicate overlap with a deletion
    """
    haplotype_results = extract_haplotype_sequences_and_insertions_from_vcf(
        chrom, start_0based, end, fasta_obj, vcf_variants, verbose=verbose)

    return tuple(haplotype_result.sequence for haplotype_result in haplotype_results)


def compute_repeat_counts_from_sequence(sequence, repeat_unit,
                                        max_allele_length_for_edit_distance_purity=DEFAULT_MAX_ALLELE_LENGTH_FOR_EDIT_DISTANCE_PURITY):
    """Compute repeat count and purity metrics from a haplotype sequence.

    Given a haplotype sequence and a repeat unit (motif), this function calculates
    the number of full repeats, the total size in base pairs, and the purity of the
    sequence relative to the expected repeat pattern.

    The repeat count uses simple integer division: len(sequence) // len(repeat_unit).
    This is an approximation that assumes the sequence starts at a motif boundary.
    The purity fields indicate how well the sequence matches a pure repeat pattern,
    revealing any interruptions or motif changes. Both are computed at whichever starting
    offset within the motif fits the sequence best, so an allele that starts in the
    middle of a motif is not penalized. purity compares the sequence position by position, so an indel inside
    the tract counts every base it shifts out of phase. purity_via_edit_distance uses the edit distance
    instead, which charges the indel itself plus the length difference it leaves at the end.

    Args:
        sequence (str): The haplotype sequence to analyze. Can be None for missing
            genotypes.
        repeat_unit (str): The expected repeat unit/motif (e.g., "CAG", "AAGGG").
        max_allele_length_for_edit_distance_purity (int): purity_via_edit_distance is left as None for
            sequences longer than this, since the edit distance costs O(length^2) per motif rotation.

    Returns:
        dict: A dictionary with the following keys:
            - num_repeats (int or None): Number of full repeats of the motif
            - repeat_size_bp (int or None): Total length of the sequence in base pairs
            - purity (float or None): Fraction of bases matching pure repeat pattern
                (0.0 to 1.0)
            - purity_via_edit_distance (float or None): (length - edit distance to a pure repeat of the
                same length) / length (0.0 to 1.0), or None if the sequence is longer than
                max_allele_length_for_edit_distance_purity
            - is_pure (bool or None): True if purity > 0.99, indicating essentially
                no interruptions

        Returns a dict with all None values if sequence is None (missing genotype).

    Example:
        >>> compute_repeat_counts_from_sequence("CAGCAGCAGCAG", "CAG")
        {'num_repeats': 4, 'repeat_size_bp': 12, 'purity': 1.0, 'purity_via_edit_distance': 1.0, 'is_pure': True}

        >>> compute_repeat_counts_from_sequence("CAGCAACAGCAG", "CAG")
        {'num_repeats': 4, 'repeat_size_bp': 12, 'purity': 0.917, 'purity_via_edit_distance': 0.917,
         'is_pure': False}

        >>> compute_repeat_counts_from_sequence(None, "CAG")
        {'num_repeats': None, 'repeat_size_bp': None, 'purity': None, 'purity_via_edit_distance': None,
         'is_pure': None}
    """
    # Handle missing genotype
    if sequence is None:
        return {
            "num_repeats": None,
            "repeat_size_bp": None,
            "purity": None,
            "purity_via_edit_distance": None,
            "is_pure": None,
        }

    # Handle empty sequence
    if not sequence:
        return {
            "num_repeats": 0,
            "repeat_size_bp": 0,
            "purity": None,
            "purity_via_edit_distance": None,
            "is_pure": None,
        }

    # Calculate number of full repeats using simple integer division
    # Repeat count formula
    num_repeats = len(sequence) // len(repeat_unit)
    repeat_size_bp = len(sequence)

    # Calculate purity - how well the sequence matches a pure repeat.
    # Include partial repeats since the sequence may not be an exact multiple of motif length, and try every
    # starting offset within the motif since an allele often starts in the middle of a motif (eg. after a 1bp
    # deletion just upstream of the locus, "CAGCAGCAG" becomes "AGCAGCAG").
    purity, _, _ = compute_best_phase_repeat_purity(sequence, repeat_unit, include_partial_repeats=True)
    if len(sequence) > max_allele_length_for_edit_distance_purity:
        purity_via_edit_distance = None
    else:
        purity_via_edit_distance, _, _ = compute_best_phase_repeat_purity(
            sequence, repeat_unit, include_partial_repeats=True, distance_metric=EDIT_DISTANCE_METRIC)

    # Handle NaN purity (can occur if sequence is shorter than motif)
    if purity != purity:  # NaN check
        purity = None
        purity_via_edit_distance = None
        is_pure = None
    else:
        # Threshold for considering a repeat "pure" (no significant interruptions)
        is_pure = purity > 0.99

    return {
        "num_repeats": num_repeats,
        "repeat_size_bp": repeat_size_bp,
        "purity": purity,
        "purity_via_edit_distance": purity_via_edit_distance,
        "is_pure": is_pure,
    }


# The fields of an overlapping VCF record that the genotype outputs actually use. Only these are kept, since
# holding every pysam VariantRecord for every locus until the output writers run would multiply the peak
# memory of a genome-wide run; write_genotypes_vcf() re-reads the full record from the input VCF anyway.
OverlappingVariant = collections.namedtuple("OverlappingVariant", ["chrom", "pos", "ref", "alts"])


def genotype_single_locus(tr_locus, vcf_file, fasta_obj, verbose=False, insertion_filter=None,
                          vcf_contig_lookup=None, fasta_contig_lookup=None, par_regions=None,
                          sex_chromosome_ploidy=None,
                          max_allele_length_for_edit_distance_purity=DEFAULT_MAX_ALLELE_LENGTH_FOR_EDIT_DISTANCE_PURITY):
    """Genotype a single tandem repeat locus using VCF variants.

    This function takes a TR locus from a catalog, fetches any overlapping VCF variants,
    extracts the diploid haplotype sequences, computes repeat counts, and returns a
    GenotypedTandemRepeat object with the complete genotyping results.

    Args:
        tr_locus (ReferenceTandemRepeat): The tandem repeat locus from the catalog
        vcf_file (pysam.VariantFile): An open pysam VariantFile object for fetching variants
        fasta_obj (pyfaidx.Fasta): Reference genome fasta object (with one_based_attributes=False)
        verbose (bool): If True, print detailed output for debugging
        insertion_filter (InsertionFilter): thresholds for deciding whether inserted bases belong to the repeat.
            An allele containing an insertion that isn't part of the repeat gets no call. When None, every
            inserted base between the locus start and end is counted. The genotype subcommand always passes a
            filter (see build_insertion_filter), so None is only for direct callers.
        vcf_contig_lookup (dict): the VCF's contig names indexed by build_contig_name_lookup, so a catalog
            and a VCF written in different naming conventions still line up. When None, the locus's contig
            name is used as given.
        fasta_contig_lookup (dict): the same thing for the reference FASTA.
        par_regions (dict): "X" and "Y" to lists of (start_0based, end) pseudoautosomal intervals (see
            get_PAR_region_coordinates). When None, no PARs are known.
        sex_chromosome_ploidy (dict): "X" and "Y" to the sample's ploidy of each outside the PARs (see
            detect_sex_chromosome_ploidy), or None to treat both as diploid. At a haploid locus a genotype that
            names only one allele is reported as HEMI rather than as a no call, and so is a locus whose two
            haplotypes come out identical.
        max_allele_length_for_edit_distance_purity (int): alleles longer than this get no edit-distance purity
            (see compute_repeat_counts_from_sequence)

    Returns:
        GenotypedTandemRepeat: A genotyped locus object containing:
            - Locus information (chrom, start, end, motif)
            - Overlapping variants
            - Allele sequences (or None if missing)
            - Repeat counts (or None if missing)
            - Zygosity classification (HOM, HET, HEMI, or None if missing)

    Notes:
        - No overlapping variants → both alleles equal reference (HOM), or a single reference allele (HEMI) at
          a haploid locus
        - A heterozygous variant at a haploid locus → the whole locus gets no call
        - Two or more heterozygous variants the VCF does not place on specific haplotypes → missing genotype
        - Insertion that isn't part of the tandem repeat on either allele → the whole locus gets no call
        - Variants that can't be applied to the reference on either allele → the whole locus gets no call
        - A genotype that leaves a haplotype uncalled → the whole locus gets no call, with the reason in
          no_call_reason, except at a haploid locus, where it's HEMI
        - A single-allele genotype such as "1" at a diploid locus → the whole locus gets no call
    """
    chrom = tr_locus.chrom
    start_0based = tr_locus.start_0based
    end = tr_locus.end_1based
    repeat_unit = tr_locus.repeat_unit

    if verbose:
        print(f"Genotyping locus: {tr_locus.locus_id}")

    no_call_reason = None
    no_call_detail = None

    # Fetch overlapping VCF variants under whatever name this VCF gives the contig. A contig the VCF simply
    # does not have yields no variants, which is the same answer as a contig it called and found nothing on.
    vcf_chrom = vcf_contig_lookup.get(normalize_chromosome_name(chrom), chrom) if vcf_contig_lookup else chrom
    variants = get_overlapping_vcf_variants(vcf_file, vcf_chrom, start_0based, end, repeat_unit)

    if verbose:
        print(f"  Found {len(variants)} overlapping variants")

    lightweight_variants = [OverlappingVariant(v.chrom, v.pos, v.ref, v.alts) for v in variants]

    # A chromosome the sample lacks (chrY in an XX sample) has nothing to genotype: with no records there, the
    # diploid path would report the reference on both haplotypes as a confident homozygous reference call.
    locus_ploidy = get_locus_ploidy(chrom, start_0based, end, par_regions, sex_chromosome_ploidy)
    if locus_ploidy == 0 and not variants:
        if verbose:
            print(f"  No call at {tr_locus.locus_id}: {NO_CALL_REASON_CHROMOSOME_ABSENT}")
        return GenotypedTandemRepeat(tr_locus=tr_locus, overlapping_variants=lightweight_variants,
                                     no_call_reason=NO_CALL_REASON_CHROMOSOME_ABSENT)

    # The chromosome was judged absent from the sample as a whole (see MIN_CHRY_RECORDS_FOR_HAPLOID_Y), yet
    # this locus has a called record on it, so the sample evidently carries the chromosome here. A VCF that
    # covers only part of the genome (an exome, a region subset) can fall under that whole-sample cutoff while
    # still calling real chrY variants, so the record wins and the locus is genotyped as the single copy it
    # would be on a haploid chromosome, or as the two copies a PAR locus has.
    if locus_ploidy == 0:
        locus_ploidy = 2 if overlaps_par(chrom, start_0based, end, par_regions) else 1

    # A haploid locus carries one copy, so a heterozygous call there contradicts the ploidy and there is no way
    # to tell which allele is real. DipCall produces these for a male whose assembly haplotypes aren't split by
    # parent, where one allele is Y-derived sequence or an assembly error (it flags them DIPX or DIPY).
    is_haploid_locus = locus_ploidy == 1
    if is_haploid_locus and any(is_heterozygous_genotype(v.samples[0].get("GT")) for v in variants):
        if verbose:
            print(f"  No call at {tr_locus.locus_id}: {NO_CALL_REASON_HET_AT_HAPLOID_LOCUS}")
        return GenotypedTandemRepeat(tr_locus=tr_locus, overlapping_variants=lightweight_variants,
                                     no_call_reason=NO_CALL_REASON_HET_AT_HAPLOID_LOCUS)

    # A diploid record with a "." allele leaves one haplotype uncalled. Reporting the surviving allele as a
    # hemizygous call is exactly what the build-error and non-repeat-insertion paths below refuse to do,
    # because that allele would then fill both the short and long allele columns and look like a real
    # hemizygous genotype. At a haploid locus the missing allele is expected rather than uncalled, so the call
    # is kept there and only the present haplotype is used.
    missing_genotype_positions = find_records_with_a_missing_genotype_allele(
        variants, require_every_allele_missing=is_haploid_locus)
    if missing_genotype_positions:
        detail = ", ".join(f"{pos:,d}" for pos in missing_genotype_positions)
        if verbose:
            print(f"  No call at {tr_locus.locus_id}: {NO_CALL_REASON_MISSING_GENOTYPE} ({detail})")
        return GenotypedTandemRepeat(tr_locus=tr_locus, overlapping_variants=lightweight_variants,
                                     no_call_reason=NO_CALL_REASON_MISSING_GENOTYPE,
                                     no_call_detail=f"at position(s) {detail}")

    # A single-allele GT such as "1" leaves the second haplotype just as undetermined as "1|." does, so at a
    # diploid locus it is a no call for the same reason. It gets its own reason because the allele isn't
    # uncalled: the caller treated the site as haploid, which on non-PAR chrX usually means
    # detect_sex_chromosome_ploidy disagreed with the caller about the sample's sex.
    if not is_haploid_locus:
        haploid_genotype_positions = find_records_with_a_haploid_genotype(variants)
        if haploid_genotype_positions:
            detail = ", ".join(f"{pos:,d}" for pos in haploid_genotype_positions)
            if verbose:
                print(f"  No call at {tr_locus.locus_id}: {NO_CALL_REASON_HAPLOID_GENOTYPE_AT_DIPLOID_LOCUS} "
                      f"({detail})")
            return GenotypedTandemRepeat(tr_locus=tr_locus, overlapping_variants=lightweight_variants,
                                         no_call_reason=NO_CALL_REASON_HAPLOID_GENOTYPE_AT_DIPLOID_LOCUS,
                                         no_call_detail=f"at position(s) {detail}")

    # Checked here as well as inside the extraction function so the reason reaches the output. Splitting these
    # variants onto haplotypes would mean guessing which alt sits with which.
    if not are_variants_unambiguously_phased(variants, verbose=verbose):
        return GenotypedTandemRepeat(tr_locus=tr_locus, overlapping_variants=lightweight_variants,
                                     no_call_reason=NO_CALL_REASON_AMBIGUOUS_PHASING)

    # Extract haplotype sequences
    # This handles phasing ambiguity (unresolved phase between heterozygous variants) by returning
    # missing genotypes
    haplotype0, haplotype1 = extract_haplotype_sequences_and_insertions_from_vcf(
        chrom, start_0based, end, fasta_obj, variants, verbose=verbose, insertion_filter=insertion_filter,
        fasta_contig_lookup=fasta_contig_lookup, repeat_unit=repeat_unit
    )

    # An insertion that isn't part of the repeat makes the whole locus a no call, not just the allele carrying it.
    # Dropping only that allele would leave a genotype indistinguishable from a real hemizygous call, and the
    # surviving allele would then be reported in both the short and long allele columns.
    if haplotype0.rejected_insertions or haplotype1.rejected_insertions:
        no_call_reason = NO_CALL_REASON_NON_REPEAT_INSERTION
        no_call_detail = (haplotype0.rejected_insertions or haplotype1.rejected_insertions)[0][1]
        if verbose:
            print(f"  No call at {tr_locus.locus_id}: "
                  f"{len(haplotype0.rejected_insertions) + len(haplotype1.rejected_insertions)} inserted "
                  f"sequence(s) are not sufficiently repetitive")
        haplotype0 = HaplotypeSequence(None, haplotype0.rejected_insertions, haplotype0.build_error,
                                        haplotype0.build_error_reason)
        haplotype1 = HaplotypeSequence(None, haplotype1.rejected_insertions, haplotype1.build_error,
                                        haplotype1.build_error_reason)

    # A haplotype whose variants couldn't be applied to the reference (a VCF ref allele that disagrees with the
    # fasta, a symbolic alt allele, two records phased onto the same haplotype with overlapping ref spans) makes
    # the whole locus a no call for the same reason: reporting only the surviving allele would look exactly like
    # a real hemizygous call, and that allele would fill both the short and long allele columns.
    if no_call_reason is None and (haplotype0.build_error or haplotype1.build_error):
        failed_haplotype = haplotype0 if haplotype0.build_error else haplotype1
        no_call_reason = failed_haplotype.build_error_reason or NO_CALL_REASON_HAPLOTYPE_BUILD_ERROR
        no_call_detail = failed_haplotype.build_error
        if verbose:
            print(f"  No call at {tr_locus.locus_id}: {no_call_reason} ({no_call_detail})")
        haplotype0 = HaplotypeSequence(None, haplotype0.rejected_insertions, haplotype0.build_error,
                                        haplotype0.build_error_reason)
        haplotype1 = HaplotypeSequence(None, haplotype1.rejected_insertions, haplotype1.build_error,
                                        haplotype1.build_error_reason)

    # A "1|." record and a ".|1" record at the same locus each leave a different haplotype uncalled, so neither
    # haplotype can be built. Without a reason here the locus would come out as a silent missing genotype.
    if (no_call_reason is None and is_haploid_locus
            and haplotype0.sequence is None and haplotype1.sequence is None):
        no_call_reason = NO_CALL_REASON_MISSING_GENOTYPE
        no_call_detail = "records leave different haplotypes uncalled"
        if verbose:
            print(f"  No call at {tr_locus.locus_id}: {no_call_reason} ({no_call_detail})")

    # At a haploid locus, two identical haplotypes are one chromosome copy reported twice: a locus with no
    # overlapping record gets the reference on both, and DipCall writes "1|1" for a male whose assembly
    # haplotypes aren't split by parent.
    # Keeping only one makes the locus HEMI like its neighbours with a ".|1" record. Heterozygous records were
    # already turned into no calls above, so both complete haplotypes always carry the same allele here.
    if (no_call_reason is None and is_haploid_locus
            and haplotype0.sequence is not None and haplotype0.sequence == haplotype1.sequence):
        haplotype1 = MISSING_HAPLOTYPE_SEQUENCES

    haplotype0_seq, haplotype1_seq = haplotype0.sequence, haplotype1.sequence

    if verbose:
        if haplotype0_seq is not None:
            print(f"  Haplotype 0: {len(haplotype0_seq)} bp")
        else:
            print(f"  Haplotype 0: missing")
        if haplotype1_seq is not None:
            print(f"  Haplotype 1: {len(haplotype1_seq)} bp")
        else:
            print(f"  Haplotype 1: missing")

    # Compute repeat counts for each haplotype that has a call
    allele1_counts = compute_repeat_counts_from_sequence(
        haplotype0_seq, repeat_unit, max_allele_length_for_edit_distance_purity)
    allele2_counts = compute_repeat_counts_from_sequence(
        haplotype1_seq, repeat_unit, max_allele_length_for_edit_distance_purity)

    if verbose:
        print(f"  Allele 1 counts: {allele1_counts}")
        print(f"  Allele 2 counts: {allele2_counts}")

    # Create and return the GenotypedTandemRepeat object
    genotyped = GenotypedTandemRepeat(
        tr_locus=tr_locus,
        overlapping_variants=lightweight_variants,
        allele1_sequence=haplotype0_seq,
        allele2_sequence=haplotype1_seq,
        num_repeats_allele1=allele1_counts["num_repeats"],
        num_repeats_allele2=allele2_counts["num_repeats"],
        allele1_purity=allele1_counts["purity"],
        allele2_purity=allele2_counts["purity"],
        allele1_purity_via_edit_distance=allele1_counts["purity_via_edit_distance"],
        allele2_purity_via_edit_distance=allele2_counts["purity_via_edit_distance"],
        num_alleles_with_non_repeat_insertions=(
            bool(haplotype0.rejected_insertions) + bool(haplotype1.rejected_insertions)),
        num_alleles_with_build_errors=(
            bool(haplotype0.build_error and haplotype0.build_error_reason not in REFERENCE_BUILD_ERROR_REASONS)
            + bool(haplotype1.build_error and haplotype1.build_error_reason not in REFERENCE_BUILD_ERROR_REASONS)),
        no_call_reason=no_call_reason,
        no_call_detail=no_call_detail,
    )

    if verbose:
        print(f"  Result: {genotyped}")

    return genotyped


def genotype_loci_using_open_files(loci, vcf_file, fasta_obj, args, vcf_contig_lookup, fasta_contig_lookup,
                                   par_regions, sex_chromosome_ploidy):
    """Genotype a list of tandem repeat loci using already-open VCF and reference fasta handles.

    This is the per-locus loop shared by the single-process and multi-process paths of genotype_all_loci().

    Args:
        loci (iterable): ReferenceTandemRepeat objects to genotype, in the order their results should be returned
        vcf_file (pysam.VariantFile): the open single-sample VCF
        fasta_obj (pyfaidx.Fasta): Reference genome fasta object (with one_based_attributes=False)
        args: Argument namespace (see genotype_all_loci)
        vcf_contig_lookup (dict): the VCF's contig names indexed by build_contig_name_lookup
        fasta_contig_lookup (dict): the reference's contig names indexed by build_contig_name_lookup
        par_regions (dict): pseudoautosomal regions from get_PAR_region_coordinates
        sex_chromosome_ploidy (dict): result of detect_sex_chromosome_ploidy, or None if no locus is on chrX/chrY

    Returns:
        tuple: (list, dict) containing:
            - List of GenotypedTandemRepeat objects, one per input locus, in input order
            - Dictionary of counters with genotyping statistics for these loci
    """
    verbose = getattr(args, 'verbose', False)
    counters = collections.defaultdict(int)
    genotyped_loci = []
    for tr_locus in loci:
        genotyped_locus = genotype_single_locus(
            tr_locus,
            vcf_file,
            fasta_obj,
            verbose=verbose,
            insertion_filter=build_insertion_filter(tr_locus.repeat_unit, args),
            vcf_contig_lookup=vcf_contig_lookup,
            fasta_contig_lookup=fasta_contig_lookup,
            par_regions=par_regions,
            sex_chromosome_ploidy=sex_chromosome_ploidy,
            max_allele_length_for_edit_distance_purity=getattr(
                args, "max_allele_length_for_edit_distance_purity", DEFAULT_MAX_ALLELE_LENGTH_FOR_EDIT_DISTANCE_PURITY),
        )
        genotyped_loci.append(genotyped_locus)

        # Update counters
        if genotyped_locus.num_overlapping_variants > 0:
            counters["loci_with_variants"] += 1
        else:
            counters["loci_without_variants"] += 1

        if genotyped_locus.num_alleles_with_non_repeat_insertions > 0:
            counters["loci_with_non_repeat_insertions"] += 1
            counters["alleles_with_non_repeat_insertions"] += (
                genotyped_locus.num_alleles_with_non_repeat_insertions)

        if genotyped_locus.num_alleles_with_build_errors > 0:
            counters["loci_with_haplotype_build_errors"] += 1
            counters["alleles_with_haplotype_build_errors"] += genotyped_locus.num_alleles_with_build_errors

        if genotyped_locus.no_call_reason is not None:
            counters[f"no_call: {genotyped_locus.no_call_reason}"] += 1

        if genotyped_locus.zygosity is None:
            counters["loci_with_missing_genotype"] += 1
        elif genotyped_locus.zygosity == "HOM":
            counters["loci_HOM"] += 1
        elif genotyped_locus.zygosity == "HET":
            counters["loci_HET"] += 1
        elif genotyped_locus.zygosity == "HEMI":
            counters["loci_HEMI"] += 1

    return genotyped_loci, counters


def genotype_loci_chunk_in_worker_process(loci_chunk, vcf_path, reference_fasta_path, args, par_regions,
                                          sex_chromosome_ploidy):
    """Genotype one chunk of catalog loci inside a multiprocessing.Pool worker.

    Open pysam and pyfaidx handles can't be sent to another process, so each chunk opens its own. The chunks
    are large enough (see MIN_LOCI_PER_GENOTYPING_CHUNK and NUM_GENOTYPING_CHUNKS_PER_WORKER_PROCESS) that the
    two opens cost a negligible fraction of the chunk's genotyping time.

    Args:
        loci_chunk (list): ReferenceTandemRepeat objects to genotype
        vcf_path (str): Path to the single-sample VCF file
        reference_fasta_path (str): Path to the reference genome fasta
        args: Argument namespace (see genotype_all_loci)
        par_regions (dict): pseudoautosomal regions from get_PAR_region_coordinates
        sex_chromosome_ploidy (dict): result of detect_sex_chromosome_ploidy, or None. It is computed once in
            the parent process since it requires scanning every chrX and chrY record in the VCF.

    Returns:
        tuple: (list, dict, dict or None) containing the genotyped loci and counters as returned by
            genotype_loci_using_open_files(), and the chunk's motif composition as returned by
            compute_motif_composition() if args.add_motif_composition is one of
            MOTIF_COMPOSITION_METHODS_COMPUTED_PER_GENOTYPING_CHUNK, or None otherwise.
    """
    vcf_file, _, vcf_contig_lookup = open_vcf_for_genotyping(vcf_path)
    fasta_obj = pyfaidx.Fasta(reference_fasta_path, one_based_attributes=False, as_raw=True)
    try:
        genotyped_loci, counters = genotype_loci_using_open_files(
            loci_chunk, vcf_file, fasta_obj, args, vcf_contig_lookup, build_contig_name_lookup(fasta_obj.keys()),
            par_regions, sex_chromosome_ploidy)
    finally:
        vcf_file.close()
        fasta_obj.close()

    if getattr(args, "add_motif_composition", None) in MOTIF_COMPOSITION_METHODS_COMPUTED_PER_GENOTYPING_CHUNK:
        return genotyped_loci, counters, compute_motif_composition(genotyped_loci, args)
    return genotyped_loci, counters, None


# How the catalog is split across worker processes when --threads > 1. Each chunk reopens the VCF and the
# reference fasta (roughly 50ms), so chunks must be large enough to make that negligible, while several chunks
# per worker let Pool.imap hand out work dynamically so one slow chunk doesn't leave the other workers idle.
NUM_GENOTYPING_CHUNKS_PER_WORKER_PROCESS = 4
MIN_LOCI_PER_GENOTYPING_CHUNK = 1000

# --add-motif-composition methods that run inside each genotyping worker process, on that worker's chunk of loci,
# so that they are parallelized by --threads too. "trf" is left out because it already runs --trf-threads TRF
# processes in parallel, and running it per worker would multiply that number by --threads.
MOTIF_COMPOSITION_METHODS_COMPUTED_PER_GENOTYPING_CHUNK = {"basic", MOTIF_DETECTION_METHOD_TRVIZ}


def genotype_all_loci(catalog_loci, vcf_path, fasta_obj, args):
    """Genotype all tandem repeat loci from a catalog using VCF variants.

    This function iterates through all TR loci in the catalog, genotypes each
    one using the VCF variants, and collects statistics about the genotyping
    results.

    Args:
        catalog_loci (list): List of ReferenceTandemRepeat objects from the catalog
        vcf_path (str): Path to the single-sample VCF file
        fasta_obj (pyfaidx.Fasta): Reference genome fasta object (with one_based_attributes=False)
        args: Argument namespace with optional attributes:
            - verbose (bool): If True, print detailed logging
            - show_progress_bar (bool): If True, display progress bar
            - threads (int): Number of worker processes. With more than 1, the loci are genotyped in chunks
              by a multiprocessing.Pool. The results are the same in either case, and in catalog order.
            - add_motif_composition (str or None): if one of MOTIF_COMPOSITION_METHODS_COMPUTED_PER_GENOTYPING_CHUNK,
              the motif composition is computed here too, by the same worker processes as the genotyping.

    Returns:
        tuple: (list, dict, dict or None) containing:
            - List of GenotypedTandemRepeat objects for all loci
            - Dictionary of counters with genotyping statistics
            - Motif composition of all loci as returned by compute_motif_composition(), or None if
              args.add_motif_composition is not one of MOTIF_COMPOSITION_METHODS_COMPUTED_PER_GENOTYPING_CHUNK
              (in which case the caller computes it, if requested)

    Notes:
        Memory usage scales with the number of loci. For very large catalogs
        (millions of loci), consider filtering to specific regions with -L.
    """
    # Get optional args with defaults
    verbose = getattr(args, 'verbose', False)
    show_progress_bar = getattr(args, 'show_progress_bar', False)
    num_threads = getattr(args, 'threads', 1)

    # Open VCF file once for all loci
    vcf_file, sample_name, vcf_contig_lookup = open_vcf_for_genotyping(vcf_path)

    # Index the reference's contig names once as well. Both lookups exist so that a "chr1" catalog can be
    # genotyped against a "1"-named VCF or reference, and a chrM locus against an MT-named one, without
    # re-deriving the spelling for every one of the catalog's millions of loci.
    fasta_contig_lookup = build_contig_name_lookup(fasta_obj.keys())
    par_regions = get_PAR_region_coordinates(fasta_obj, fasta_contig_lookup)

    # Detection scans every non-PAR chrX and chrY record in the VCF, so it is skipped when no locus being
    # genotyped (after -L) is on either chromosome, since its result would go unused.
    if any(normalize_chromosome_name(tr_locus.chrom) in ("X", "Y") for tr_locus in catalog_loci):
        sex_chromosome_ploidy = detect_sex_chromosome_ploidy(vcf_file, vcf_contig_lookup, par_regions)
    else:
        sex_chromosome_ploidy = None

    if verbose:
        print(f"Genotyping {len(catalog_loci):,d} TR loci using variants from sample: {sample_name}"
              + (f" with {num_threads} worker processes" if num_threads > 1 else ""))

    # Initialize counters for statistics
    counters = collections.defaultdict(int)
    counters["total_loci"] = len(catalog_loci)

    if num_threads == 1:
        loci_iterator = catalog_loci
        if show_progress_bar:
            loci_iterator = tqdm.tqdm(catalog_loci, unit=" loci", unit_scale=True)
        genotyped_loci, chunk_counters = genotype_loci_using_open_files(
            loci_iterator, vcf_file, fasta_obj, args, vcf_contig_lookup, fasta_contig_lookup, par_regions,
            sex_chromosome_ploidy)
        for key, value in chunk_counters.items():
            counters[key] += value
        motif_lists_by_locus = None
        if getattr(args, "add_motif_composition", None) in MOTIF_COMPOSITION_METHODS_COMPUTED_PER_GENOTYPING_CHUNK:
            motif_lists_by_locus = compute_motif_composition(genotyped_loci, args)
    else:
        chunk_size = max(MIN_LOCI_PER_GENOTYPING_CHUNK,
                         math.ceil(len(catalog_loci) / (num_threads * NUM_GENOTYPING_CHUNKS_PER_WORKER_PROCESS)))
        loci_chunks = [catalog_loci[i:i + chunk_size] for i in range(0, len(catalog_loci), chunk_size)]
        genotype_chunk = functools.partial(
            genotype_loci_chunk_in_worker_process,
            vcf_path=vcf_path, reference_fasta_path=fasta_obj.filename, args=args, par_regions=par_regions,
            sex_chromosome_ploidy=sex_chromosome_ploidy)
        progress_bar = tqdm.tqdm(total=len(catalog_loci), unit=" loci", unit_scale=True) if show_progress_bar else None
        genotyped_loci = []
        motif_lists_by_locus = (
            {} if getattr(args, "add_motif_composition", None) in MOTIF_COMPOSITION_METHODS_COMPUTED_PER_GENOTYPING_CHUNK
            else None)
        # Pool.imap returns the chunks in the order they were submitted, which keeps the output in catalog order
        # and lets each finished chunk be consumed as soon as it arrives rather than after the whole pool is done.
        with multiprocessing.Pool(min(num_threads, len(loci_chunks))) as pool:
            for chunk_genotyped_loci, chunk_counters, chunk_motif_lists_by_locus in pool.imap(
                    genotype_chunk, loci_chunks):
                genotyped_loci.extend(chunk_genotyped_loci)
                for key, value in chunk_counters.items():
                    counters[key] += value
                if chunk_motif_lists_by_locus is not None:
                    motif_lists_by_locus.update(chunk_motif_lists_by_locus)
                if progress_bar is not None:
                    progress_bar.update(len(chunk_genotyped_loci))
        if progress_bar is not None:
            progress_bar.close()

    # Close VCF file
    vcf_file.close()

    # Print summary if verbose
    if verbose:
        print(f"\nGenotyping complete:")
        print(f"  Total loci: {counters['total_loci']:,d}")
        print(f"  Loci with overlapping variants: {counters['loci_with_variants']:,d}")
        print(f"  Loci without overlapping variants: {counters['loci_without_variants']:,d}")
        print(f"  Loci with missing genotypes: {counters['loci_with_missing_genotype']:,d}")
        print(f"  Zygosity breakdown:")
        print(f"    HOM: {counters['loci_HOM']:,d}")
        print(f"    HET: {counters['loci_HET']:,d}")
        print(f"    HEMI: {counters['loci_HEMI']:,d}")
        # One line per reason, rather than a hardcoded line per reason: the per-reason counters are keyed on
        # the NoCallReason values themselves, so a new reason shows up here without touching this block.
        no_call_counter_keys = sorted(key for key in counters if key.startswith("no_call: "))
        if no_call_counter_keys:
            print(f"  Reasons loci were set to no call:")
            for counter_key in no_call_counter_keys:
                print(f"    {counters[counter_key]:10,d}  {counter_key[len('no_call: '):]}")
        # The counts above are loci. These are alleles, so they are reported separately rather than nested
        # under a locus count they can exceed: a locus has two alleles and both can be affected.
        if counters["alleles_with_non_repeat_insertions"] or counters["alleles_with_haplotype_build_errors"]:
            print(f"  Alleles at loci:")
            if counters["alleles_with_non_repeat_insertions"]:
                print(f"    {counters['alleles_with_non_repeat_insertions']:10,d}  carried an insertion "
                      f"that wasn't sufficiently repetitive")
            if counters["alleles_with_haplotype_build_errors"]:
                print(f"    {counters['alleles_with_haplotype_build_errors']:10,d}  had variants that "
                      f"couldn't be applied to the reference")

    return genotyped_loci, counters, motif_lists_by_locus


def do_catalog_subcommand(args):
    """Main function to parse arguments and run the tandem repeat detection pipeline."""

    fasta_obj = pyfaidx.Fasta(args.reference_fasta_path, one_based_attributes=False, as_raw=True)

    # parse input VCF
    counters = collections.defaultdict(int)
    filtered_alleles = {} if args.write_filtered_variants_to_vcf else None
    alleles_from_vcf = parse_input_vcf_file(args, counters, fasta_obj, filtered_alleles=filtered_alleles)

    # detect tandem repeats
    alleles_that_are_tandem_repeats, alleles_to_process_using_trf = detect_perfect_and_almost_perfect_tandem_repeats(
        alleles_from_vcf, counters, args, filtered_alleles=filtered_alleles)

    if not args.dont_run_trf:
        more_alleles_that_are_tandem_repeats = detect_tandem_repeats_using_trf(
            alleles_to_process_using_trf, counters, args, filtered_alleles=filtered_alleles,
            alleles_already_accepted={tr.allele for tr in alleles_that_are_tandem_repeats})
        alleles_that_are_tandem_repeats.extend(more_alleles_that_are_tandem_repeats)

    # write results to output file(s)
    if not args.output_prefix:
        args.output_prefix = args.input_vcf_prefix

    write_bed(alleles_that_are_tandem_repeats, args)

    if args.write_detailed_bed:
        write_bed(alleles_that_are_tandem_repeats, args, detailed=True)

    if args.write_vcf:
        write_vcf(alleles_that_are_tandem_repeats, args, only_write_filtered_out_alleles=False)

    if args.write_filtered_variants_to_vcf:
        write_vcf(alleles_that_are_tandem_repeats, args, only_write_filtered_out_alleles=True,
                  filtered_alleles=filtered_alleles)

    if args.write_tsv:
        write_tsv(alleles_that_are_tandem_repeats, args)

    if args.write_fasta:
        write_fasta(alleles_that_are_tandem_repeats, args)

    if args.verbose:
        print_stats(counters)
        print_tr_stats(alleles_that_are_tandem_repeats)


def detect_perfect_and_almost_perfect_tandem_repeats(alleles, counters, args, filtered_alleles=None):

    alleles_to_process_next = [(allele, DETECTION_MODE_PURE_REPEATS) for allele in alleles]
    alleles_to_process_next_using_trf = []
    tandem_repeat_alleles = []
    first_iteration = True
    while alleles_to_process_next:
        alleles_to_reprocess = []
        if args.verbose:
            print(f"Checking {len(alleles_to_process_next):,d} indel alleles for tandem repeats",
                   "after extending their flanking sequences" if not first_iteration else "")

        first_iteration = False

        if args.show_progress_bar:
            alleles_to_process_next = tqdm.tqdm(alleles_to_process_next, unit=" alleles", unit_scale=True)

        for allele, detection_mode in alleles_to_process_next:
            tandem_repeat_allele, filter_reason = check_if_allele_is_tandem_repeat(allele, args, detection_mode)

            if filter_reason:
                if allele.previously_increased_flanking_sequence_size:
                    raise ValueError(f"Logic error: allele {allele} was previously detected to have tandem repeats "
                                     f"that extended over the entire flanking sequence, but after extending their "
                                     f"flanking sequences, it was no longer found to be a tandem repeat using "
                                     f"detection mode {detection_mode} due to {filter_reason}")
                
                # this allele was not found to be a tandem repeat using the current detection mode
                if not args.dont_allow_interruptions and detection_mode == DETECTION_MODE_PURE_REPEATS:
                    # try the detection mode that allows for interruptions
                    alleles_to_reprocess.append((allele, DETECTION_MODE_ALLOW_INTERRUPTIONS))
                elif not args.dont_run_trf and len(allele.variant_bases) >= args.min_indel_size_to_run_trf:
                    # try using TRF
                    alleles_to_process_next_using_trf.append(allele)
                else:
                    counters[f"allele filter: {detection_mode}: {filter_reason}"] += 1
                    if filtered_alleles is not None:
                        filtered_alleles[(allele.chrom, allele.pos, allele.ref, allele.alt)] = filter_reason

                continue

            # reprocess the allele if the repeats were found to cover the entire left or right flanking sequence
            if need_to_reprocess_allele_with_extended_flanking_sequence(tandem_repeat_allele):
                counters[(f"allele op: increased flanking sequence size "
                          f"{tandem_repeat_allele.allele.number_of_times_flanking_sequence_size_was_increased}x for "
                          f"{detection_mode} repeats")] += 1
                alleles_to_reprocess.append((allele, detection_mode))  # reprocess the allele with the same detection mode
                #print(f"Detection mode [{detection_mode}]: Increasing flanking sequence size to {len(tandem_repeat_allele.allele.get_left_flanking_sequence()):,d}bp for {tandem_repeat_allele}")
                continue
            
            # this allele was found to be a tandem repeat using the current detection mode
            if tandem_repeat_allele.do_repeats_cover_entire_flanking_sequence() and not (
                    tandem_repeat_allele.allele.get_left_flank_stops_at_N() or
                    tandem_repeat_allele.allele.get_right_flank_stops_at_N()):
                print(f"WARNING: allele {allele} was found to be a tandem repeat using detection mode {detection_mode}, "
                      f"but the repeats cover the entire flanking sequence even though it is longer than "
                      f"{MAX_FLANKING_SEQUENCE_SIZE:,}bp. Skipping...")
            else:
                tandem_repeat_alleles.append(tandem_repeat_allele)

                if not args.dont_run_trf and tandem_repeat_allele.repeat_unit_length > 6 and len(allele.variant_bases) >= args.min_indel_size_to_run_trf:
                    # if this is a VNTR with a large motif, run TRF on it to see if it detects wider locus boundaries.
                    # The merge step can resolve redundant locus definitions.
                    alleles_to_process_next_using_trf.append(allele)

        alleles_to_process_next = alleles_to_reprocess

    if args.verbose:
        print(f"Found {sum(1 for tr in tandem_repeat_alleles if tr.is_pure_repeat):,d} perfect tandem repeat alleles and {sum(1 for tr in tandem_repeat_alleles if not tr.is_pure_repeat):,d} nearly-perfect tandem repeat alleles")

    return tandem_repeat_alleles, alleles_to_process_next_using_trf


def run_trf_batches_in_parallel(items, num_threads, worker_fn):
    """Distribute items round-robin across worker threads and run worker_fn on each batch in parallel.

    TRF is an external subprocess that releases the GIL, so running one TRF instance per thread provides
    real parallelism. This harness is shared by the 'catalog' and 'genotype' subcommands.

    Args:
        items (list): Items to distribute across threads (e.g. Allele objects or allele sequences).
        num_threads (int): Maximum number of worker threads. Capped at len(items).
        worker_fn (callable): Called as worker_fn(batch, thread_id), where batch is the sublist of items
            assigned to thread thread_id. Must return an iterable of result records.

    Yields:
        The result records returned by each worker_fn call, in thread order.
    """
    n_threads = max(1, min(num_threads, len(items)))
    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        futures = [
            executor.submit(
                worker_fn,
                [item for item_i, item in enumerate(items) if item_i % n_threads == thread_i],
                thread_i)
            for thread_i in range(n_threads)
        ]
        for future in futures:
            yield from future.result()


def detect_tandem_repeats_using_trf(alleles, counters, args, filtered_alleles=None, alleles_already_accepted=()):
    """Runs TandemRepeatFinder (TRF) on a list of indel alleles to detect tandem repeats.

    Args:
        alleles (list): Allele objects to run TRF on
        counters (dict): counter name to count, updated in place
        args (argparse.Namespace): command-line arguments parsed by parse_args()
        filtered_alleles (dict): (chrom, pos, ref, alt) to filter reason, updated in place for alleles TRF rejects,
            or None to skip recording them
        alleles_already_accepted (collection): Allele objects that pure or interrupted detection already
            accepted as tandem repeats and that are only run through TRF to look for wider locus boundaries.
            TRF finding nothing for one of these is not a filter, since the allele is in the output either way.

    Returns:
        list: TandemRepeatAllele objects for the alleles TRF found to be tandem repeats
    """

    tandem_repeat_alleles = []
    first_iteration = True
    alleles_to_process_next = alleles
    while alleles_to_process_next:
        before_counter = len(tandem_repeat_alleles)
        start_time = datetime.datetime.now()

        trf_working_dir = os.path.join(args.trf_working_dir, CURRENT_TIMESTAMP,
                                       f"{args.input_vcf_prefix}." + start_time.strftime("%Y%m%d_%H%M%S.%f"))
        if os.path.isdir(trf_working_dir):
            raise ValueError(f"ERROR: TRF working directory already exists: {trf_working_dir}. Each TRF run should be in a unique directory to avoid filename collisions.")

        if args.debug:
            print("-"*100)
            print(f"TRF working directory: {trf_working_dir}")
        os.makedirs(trf_working_dir)

        try:
            n_threads = min(args.trf_threads, len(alleles_to_process_next))
            if args.verbose:
                print(f"Launching {n_threads} TRF instance(s) to check {len(alleles_to_process_next):,d} indel alleles for tandem repeats",
                        "after extending their flanking sequences" if not first_iteration else "")
            first_iteration = False

            alleles_to_reprocess = []
            for tandem_repeat_allele, filter_reason, allele in run_trf_batches_in_parallel(
                    alleles_to_process_next, n_threads,
                    lambda batch, thread_i: run_trf(batch, args, thread_i, trf_working_dir)):
                if filter_reason:
                    if allele in alleles_already_accepted:
                        counters["allele op: TRF found no wider locus for an already accepted VNTR allele"] += 1
                    else:
                        counters[f"allele filter: TRF: {filter_reason}"] += 1
                        if filtered_alleles is not None:
                            filtered_alleles[(allele.chrom, allele.pos, allele.ref, allele.alt)] = filter_reason
                    continue

                # reprocess the allele if the repeats were found to cover the entire left or right flanking sequence
                if need_to_reprocess_allele_with_extended_flanking_sequence(tandem_repeat_allele):
                    counters[f"allele op: increased flanking sequence size {tandem_repeat_allele.allele.number_of_times_flanking_sequence_size_was_increased}x for TRF"] += 1
                    alleles_to_reprocess.append(allele)
                    continue

                # this allele was found to be a tandem repeat using TRF
                if tandem_repeat_allele.do_repeats_cover_entire_flanking_sequence() and not (
                        tandem_repeat_allele.allele.get_left_flank_stops_at_N() or
                        tandem_repeat_allele.allele.get_right_flank_stops_at_N()):
                    print(f"WARNING: allele {allele} was found to be a tandem repeat using TRF, but the repeats "
                          f"cover the entire flanking sequence even though it is longer than "
                          f"{MAX_FLANKING_SEQUENCE_SIZE:,}bp. Skipping...")
                else:
                    tandem_repeat_alleles.append(tandem_repeat_allele)

            alleles_to_process_next = alleles_to_reprocess

            elapsed = datetime.datetime.now() - start_time
            if args.verbose:
                print(f"Found {len(tandem_repeat_alleles) - before_counter:,d} additional tandem repeats after running TRF for {elapsed.seconds//60}m {elapsed.seconds%60}s"
                  + (f", and will recheck {len(alleles_to_process_next):,d} other alleles after extending their flanking sequences" if len(alleles_to_process_next) > 0 else ""))
        finally:
            if not args.debug:
                shutil.rmtree(trf_working_dir)

    return tandem_repeat_alleles


def parse_input_vcf_file(args, counters, fasta_obj, filtered_alleles=None):
    """Parse the input VCF file and return a list of Allele objects.

    Args:
        args (argparse.Namespace): command-line arguments parsed by parse_args()
        counters (dict): counter name to count, updated in place
        fasta_obj (pyfaidx.Fasta): the reference genome
        filtered_alleles (dict): (chrom, pos, ref, alt) to filter reason, updated in place for alleles dropped here,
            or None to skip recording them

    Returns:
        list: the Allele objects to check for tandem repeats
    """

    vcf_iterator = get_input_vcf_iterator(args, include_header=False)

    if args.show_progress_bar:
        vcf_iterator = tqdm.tqdm(vcf_iterator, unit=" variants", unit_scale=True)


    # iterate over all VCF rows
    alleles_from_vcf = []
    vcf_line_i = 0
    allele_order = 0
    for line in vcf_iterator:
        if line.startswith("#"):
            continue

        vcf_fields = line.strip().split("\t")
        if vcf_line_i < args.offset:
            vcf_line_i += 1
            continue

        if args.n is not None and vcf_line_i >= args.offset + args.n:
            break

        vcf_line_i += 1

        # parse the ALT allele(s)
        vcf_chrom = vcf_fields[0]
        vcf_pos = int(vcf_fields[1])
        vcf_ref = vcf_fields[3].upper()
        vcf_alt = vcf_fields[4].upper()
        alt_alleles = vcf_alt.split(",")

        info_field_dict = None
        if args.copy_info_field_keys_to_tsv and args.write_tsv and vcf_fields[7] and vcf_fields[7] != ".":
            info_field_dict = {}
            for info_field_value in vcf_fields[7].split(";"):
                info_field_key_value = info_field_value.split("=")
                key = info_field_key_value[0]
                if key in args.copy_info_field_keys_to_tsv:
                    value = info_field_key_value[1] if len(info_field_key_value) > 1 else True    
                    info_field_dict[key] = value
                    args.copy_info_field_keys_to_tsv[key] += 1

        if vcf_chrom not in fasta_obj:
            raise ValueError(f"Chromosome '{vcf_chrom}' not found in the reference fasta")

        # check for N's in the ref or alt sequences
        if "N" in vcf_ref or "N" in vcf_alt:
            counters[f"allele filter: {FILTER_ALLELE_WITH_N_BASES}"] += 1
            if filtered_alleles is not None:
                for alt_allele in alt_alleles:
                    filtered_alleles[(vcf_chrom, vcf_pos, vcf_ref, alt_allele)] = FILTER_ALLELE_WITH_N_BASES
            continue

        if not vcf_alt:
            raise ValueError(f"No ALT allele found in VCF row #{vcf_line_i + 1:,d}: {vcf_fields}")

        # Handle '*' alleles
        if "*" in alt_alleles:
            # if this variant has 1 regular allele and 1 "*" allele (which represents an overlapping deletion), discard the
            # "*" allele and recode the genotype as haploid
            alt_alleles = [a for a in alt_alleles if a != "*"]

        if len(alt_alleles) == 0:
            counters[f"WARNING: variant with no alt alleles"] += 1
            continue

        counters["variant counts: TOTAL variants"] += 1
        counters["allele counts: TOTAL alleles"] += len(alt_alleles)

        # check if the ALT alleles pass basic filters
        for alt_allele in alt_alleles:
            if len(vcf_ref) == len(alt_allele):
                counters[f"allele filter: {'SNV' if len(alt_allele) == 1 else 'MNV'}"] += 1
                continue
            elif len(vcf_ref) < len(alt_allele) and not alt_allele.startswith(vcf_ref):
                counters[f"allele filter: complex MNV deletion/insertion"] += 1
                continue
            elif len(alt_allele) < len(vcf_ref) and not vcf_ref.startswith(alt_allele):
                counters[f"allele filter: complex MNV insertion/deletion"] += 1
                continue

            allele = Allele(
                vcf_chrom, vcf_pos, vcf_ref, alt_allele, fasta_obj,
                order=allele_order,
                info_field_dict=info_field_dict)

            allele_order += 1

            if len(allele.variant_bases) > MAX_INDEL_SIZE:
                # this is a very large indel, so we don't want to process it
                counters[f"allele filter: {FILTER_ALLELE_TOO_BIG}"] += 1
                if filtered_alleles is not None:
                    filtered_alleles[(vcf_chrom, vcf_pos, vcf_ref, alt_allele)] = FILTER_ALLELE_TOO_BIG
                continue

            counters[f"allele counts: {allele.ins_or_del} alleles"] += 1

            alleles_from_vcf.append(allele)

    if args.verbose:
        print(f"Parsed {len(alleles_from_vcf):,d} indel alleles from {args.input_vcf_path}")
    
    return alleles_from_vcf


def check_if_allele_is_tandem_repeat(allele, args, detection_mode):
    """Determine if the given allele is a tandem repeat expansion or contraction or not.
    This is done by performing a brute-force scan for perfect (or nearly perfect) repeats in the allele sequence, and then extending the repeats
    into the flanking reference sequences.

    Args:
        allele (Allele): allele record
        args (argparse.Namespace): command-line arguments parsed by parse_args()
        detection_mode (str): Should be either DETECTION_MODE_PURE_REPEATS or DETECTION_MODE_ALLOW_INTERRUPTIONS

    Returns:
        2-tuple (TandemRepeatAllele, str):
            TandemRepeatAllele: if the allele represents a tandem repeat, this will be a TandemRepeatAllele object, otherwise it will be None.
            str: if the allele is not a tandem repeat, this will be a string describing the reason why the allele failed tandem repeat filters,
                or otherwise None if it passed all filters.
    """

    left_flanking_reference_sequence = allele.get_left_flanking_sequence()
    right_flanking_reference_sequence = allele.get_right_flanking_sequence()

    if detection_mode == DETECTION_MODE_PURE_REPEATS:
        # check whether this variant allele + flanking sequences represent a pure tandem repeat expansion or contraction
        (
            repeat_unit,
            num_total_repeats_in_variant_bases,
            _,
        ) = find_repeat_unit_without_allowing_interruptions(allele.variant_bases)

        num_total_repeats_left_flank = extend_repeat_into_sequence_without_allowing_interruptions(
            repeat_unit[::-1],
            left_flanking_reference_sequence[::-1])
        num_total_repeats_right_flank = extend_repeat_into_sequence_without_allowing_interruptions(
            repeat_unit,
            right_flanking_reference_sequence)

        num_repeat_bases_in_left_flank = num_total_repeats_left_flank * len(repeat_unit)
        num_repeat_bases_in_right_flank = num_total_repeats_right_flank * len(repeat_unit)

    elif detection_mode == DETECTION_MODE_ALLOW_INTERRUPTIONS:
        (
            repeat_unit,
            num_pure_repeats_in_variant_bases,
            num_total_repeats_in_variant_bases,
            repeat_unit_interruption_index,
            _
        ) = find_repeat_unit_allowing_interruptions(allele.variant_bases, allow_partial_repeats=False)

        reversed_repeat_unit_interruption_index = None
        if repeat_unit_interruption_index is not None:
            reversed_repeat_unit_interruption_index = (len(repeat_unit) - 1 - repeat_unit_interruption_index)

        num_pure_repeats_left_flank, num_total_repeats_left_flank, reversed_repeat_unit_interruption_index = extend_repeat_into_sequence_allowing_interruptions(
            repeat_unit[::-1],
            left_flanking_reference_sequence[::-1],
            repeat_unit_interruption_index=reversed_repeat_unit_interruption_index)

        if reversed_repeat_unit_interruption_index is not None:
            # reverse the repeat_unit_interruption_index
            repeat_unit_interruption_index = len(repeat_unit) - 1 - reversed_repeat_unit_interruption_index

        num_pure_repeats_right_flank, num_total_repeats_right_flank, repeat_unit_interruption_index = extend_repeat_into_sequence_allowing_interruptions(
            repeat_unit,
            right_flanking_reference_sequence,
            repeat_unit_interruption_index=repeat_unit_interruption_index)

        num_repeat_bases_in_left_flank = num_total_repeats_left_flank * len(repeat_unit)
        num_repeat_bases_in_right_flank = num_total_repeats_right_flank * len(repeat_unit)

        simplified_repeat_unit, _, _ = find_repeat_unit_without_allowing_interruptions(repeat_unit, allow_partial_repeats=False)
        repeat_unit = simplified_repeat_unit

    else:
        raise ValueError(f"Invalid detection_mode: '{detection_mode}'. It must be either '{DETECTION_MODE_PURE_REPEATS}' or '{DETECTION_MODE_ALLOW_INTERRUPTIONS}'.")

    tandem_repeat_allele = TandemRepeatAllele(
        allele,
        repeat_unit,
        adjust_repeat_unit=True,
        num_repeat_bases_in_left_flank=num_repeat_bases_in_left_flank,
        num_repeat_bases_in_variant=len(allele.variant_bases),
        num_repeat_bases_in_right_flank=num_repeat_bases_in_right_flank,
        detection_mode=detection_mode,
    )

    tandem_repeat_allele_failed_filters_reason = check_if_tandem_repeat_allele_failed_filters(args, tandem_repeat_allele)
    
    if args.debug: print(f"{detection_mode} repeats: {tandem_repeat_allele}, filter: {tandem_repeat_allele_failed_filters_reason}")
    if tandem_repeat_allele_failed_filters_reason is not None:
        return None, tandem_repeat_allele_failed_filters_reason
    else:
        return tandem_repeat_allele, None


def check_if_tandem_repeat_allele_failed_filters(args, tandem_repeat_allele, detected_by_trf=False):
    """Check if the given tandem repeat allele (represented by its repeat_unit and total_repeats) passes the filters
    specified in the command-line arguments.

    Args:
        args (argparse.Namespace): command-line arguments parsed by parse_args()
        tandem_repeat_allele (TandemRepeatAllele): The tandem repeat allele.
        detected_by_trf (bool): Whether this allele was detected by TRF (applies stricter filters if True).

    Returns:
        str: A string describing the reason why the allele failed filters, or None if the allele passed all filters.
    """
    total_repeats = tandem_repeat_allele.num_repeats_in_left_flank + tandem_repeat_allele.num_repeats_in_variant + tandem_repeat_allele.num_repeats_in_right_flank
    total_repeat_bases = tandem_repeat_allele.num_repeat_bases_in_left_flank + tandem_repeat_allele.num_repeat_bases_in_variant + tandem_repeat_allele.num_repeat_bases_in_right_flank
    repeat_unit = tandem_repeat_allele.repeat_unit

    if total_repeats == 1:
        # no repeat unit found in this allele
        return FILTER_ALLELE_INDEL_WITHOUT_REPEATS
    elif total_repeats < args.min_repeats:
        return FILTER_TR_ALLELE_NOT_ENOUGH_REPEATS.format(args.min_repeats)
    elif total_repeat_bases < args.min_tandem_repeat_length:
        return FILTER_TR_ALLELE_DOESNT_SPAN_ENOUGH_BASE_PAIRS.format(args.min_tandem_repeat_length)
    elif len(repeat_unit) < args.min_repeat_unit_length:
        return FILTER_TR_ALLELE_REPEAT_UNIT_TOO_SHORT.format(args.min_repeat_unit_length)
    elif len(repeat_unit) > args.max_repeat_unit_length:
        return FILTER_TR_ALLELE_REPEAT_UNIT_TOO_LONG.format(args.max_repeat_unit_length)

    if "N" in tandem_repeat_allele.variant_and_flanks_repeat_sequence:
        return FILTER_ALLELE_WITH_N_BASES

    if detected_by_trf:
        # apply extra criteria
        total_repeat_bases_in_reference = tandem_repeat_allele.end_1based - tandem_repeat_allele.start_0based
        if total_repeat_bases_in_reference < args.trf_min_repeats_in_reference * len(repeat_unit):
            return FILTER_TR_ALLELE_NOT_ENOUGH_REPEATS_IN_REFERENCE.format(args.trf_min_repeats_in_reference)
        if total_repeat_bases_in_reference > TRF_MAX_REPEATS_IN_REFERENCE_THRESHOLD * len(repeat_unit):
            return FILTER_TR_ALLELE_TOO_MANY_REPEATS.format(TRF_MAX_REPEATS_IN_REFERENCE_THRESHOLD)
        if total_repeat_bases_in_reference > TRF_MAX_SPAN_IN_REFERENCE_THRESHOLD:
            return FILTER_TR_ALLELE_SPANS_TOO_MANY_BASE_PAIRS.format(TRF_MAX_SPAN_IN_REFERENCE_THRESHOLD)
        if tandem_repeat_allele.repeat_purity < args.trf_min_purity:
            return FILTER_TR_ALLELE_PURITY_IS_TOO_LOW.format(args.trf_min_purity)

    return None  # did not fail filters


def run_trf(alleles, args, thread_id=0, trf_working_dir=None):
    """Run TRF on the given allele records.

    Args:
        alleles (list): List of Allele objects to run TRF on.
        args (argparse.Namespace): Command-line arguments parsed by parse_args().
        thread_id (int): ID of thread executing this function, starting from 0.
        trf_working_dir (str): Directory where TRF input/output files will be written.
            If None, uses the current working directory.

    Returns:
        list of 3-tuples: (tandem_repeat_allele, filter_reason, allele)
            tandem_repeat_allele (TandemRepeatAllele): tandem repeat allele object if the allele is a tandem repeat, None otherwise
            filter_reason (str): string describing the reason why the allele failed filters, or None if the allele passed all filters
            allele (Allele): Allele object that was processed
    """

    if trf_working_dir is None:
        trf_working_dir = os.getcwd()

    trf_fasta_filename = os.path.join(trf_working_dir, f"trf_input_sequences__thread{thread_id}.fa")
    with open(trf_fasta_filename, "wt") as f:
        # reverse the left flanking sequence + variant bases so that TRF starts detecting repeats from the
        # end of the variant sequence rather than some random point in the flanking region. This will
        # make it so the repeat unit is detected in the correct orientation (after being reversed back).
        # For example, if the variant was T > TCAGCAGCAGCAG , we want TRF to start from GAC ..

        for allele in alleles:

            left_flank_and_variant_bases = f"{allele.get_left_flanking_sequence()}{allele.variant_bases}"
            left_flank_and_variant_bases_reversed = left_flank_and_variant_bases[::-1]
            f.write(f">{allele.shortened_variant_id}$left\n")
            f.write(f"{left_flank_and_variant_bases_reversed}\n")

            variant_bases_and_right_flank = f"{allele.variant_bases}{allele.get_right_flanking_sequence()}"
            f.write(f">{allele.shortened_variant_id}$right\n")
            f.write(f"{variant_bases_and_right_flank}\n")


    # check that the results overlap with the variant bases
    trf_runner = TRFRunner(
        args.trf_executable_path,
        html_mode=True,
        min_motif_size=args.min_repeat_unit_length,
        max_motif_size=args.max_repeat_unit_length,
        match_score = 2,
        mismatch_penalty = args.trf_mismatch_penalty,
        indel_penalty = args.trf_indel_penalty,
        minscore = args.trf_min_score,
        debug=False,
        generate_motif_logo_plots=False,
    )

    trf_runner.run_trf_on_fasta_file(trf_fasta_filename)

    # parse the TRF output
    trf_allele_filter_counters = collections.Counter()
    results = []
    for allele_i, allele in enumerate(alleles):

        motif_size_to_matching_trf_results = collections.defaultdict(dict)
        for left_or_right in "left", "right":
            trf_results = trf_runner.parse_html_results(
                trf_fasta_filename,
                sequence_number=2 * allele_i + (0 if left_or_right == "left" else 1) + 1,
                total_sequences=2 * len(alleles))

            for trf_result in trf_results:
                if trf_result["sequence_name"] != f"{allele.shortened_variant_id}${left_or_right}":
                    raise ValueError(f"TRF result sequence name '{trf_result['sequence_name']}' does not match the expected variant ID '{allele.shortened_variant_id}${left_or_right}' in \n{pformat(trf_results)}")

                if trf_result["start_0based"] >= trf_result["repeat_unit_length"] or trf_result["end_1based"] <= len(allele.variant_bases) - trf_result["repeat_unit_length"]:
                    if args.debug:
                        print(f"TRF filtered out: {allele} because the repeat unit is {trf_result['repeat_unit_length']} bases long and starts at {trf_result['start_0based']} or because it ends at {trf_result['end_1based']} before the end of the variant ({len(allele.variant_bases)})")
                    # repeats must start close to the start of the variant bases and end near or after the end of the variant bases
                    continue

                if trf_result["repeat_unit_length"] > len(allele.variant_bases):
                    if args.debug:
                        print(f"TRF filtered out: {allele} because the repeat unit length ({trf_result['repeat_unit_length']}) is larger than the variant ({len(allele.variant_bases)})")
                    # to be a tandem repeat variant, it should represent an expansion or contraction by at least one repeat unit
                    continue

                # compute the number of tandem repeats in the variant bases
                variant_bases_sum = trf_result["start_0based"]
                repeat_count = 0
                while repeat_count < len(trf_result["repeats"]) and variant_bases_sum < len(allele.variant_bases):
                    current_repeat = trf_result["repeats"][repeat_count].replace("-", "")
                    variant_bases_sum += len(current_repeat)
                    repeat_count += 1

                trf_result["num_repeats_in_variant"] = repeat_count
                trf_result["tandem_repeat_bases_in_flank"]= max(0, trf_result["end_1based"] - len(allele.variant_bases))

                if left_or_right == "left":
                    trf_result["repeat_unit"] = trf_result["repeat_unit"][::-1]

                motif_size_to_matching_trf_results[trf_result["repeat_unit_length"]][left_or_right] = trf_result

        # filter out results where left and right have different repeat units
        motif_size_to_passing_trf_results = {}
        motif_size_to_tandem_repeat_allele = {}
        for motif_size, matching_trf_results in sorted(motif_size_to_matching_trf_results.items()):
            if matching_trf_results.get("left") and matching_trf_results.get("right") and not are_repeat_units_similar(
                compute_canonical_motif(matching_trf_results["left"]["repeat_unit"], include_reverse_complement=False), 
                compute_canonical_motif(matching_trf_results["right"]["repeat_unit"], include_reverse_complement=False)
            ):
                if args.debug: print(f"TRF filtered out: {allele}, filter reason: {matching_trf_results['left']['repeat_unit']} and {matching_trf_results['right']['repeat_unit']} are not similar")
                continue
                
            if matching_trf_results.get("right"):
                repeat_unit = matching_trf_results["right"]["repeat_unit"]

            elif matching_trf_results.get("left"):
                repeat_unit = matching_trf_results["left"]["repeat_unit"]

            else:
                raise ValueError(f"Logic error: No TRF results for allele {allele} with motif size {motif_size}")

            # check if the repeat unit itself consists of perfect repeats of a smaller repeat unit
            # (this happens in ~3% of TRs detected by TRF)
            simplified_repeat_unit, _, _ = find_repeat_unit_without_allowing_interruptions(repeat_unit, allow_partial_repeats=False)
            if len(simplified_repeat_unit) in motif_size_to_tandem_repeat_allele:
                # if a repeat unit has been simplified to a smaller repeat unit size that was already recorded, skip it.
                continue

            # update repeat unit to one that has the highest purity
            tandem_repeat_allele = TandemRepeatAllele(
                allele,
                repeat_unit=simplified_repeat_unit,
                adjust_repeat_unit=True,
                num_repeat_bases_in_left_flank=matching_trf_results.get("left", {}).get("tandem_repeat_bases_in_flank", 0),
                num_repeat_bases_in_variant=len(allele.variant_bases),
                num_repeat_bases_in_right_flank=matching_trf_results.get("right", {}).get("tandem_repeat_bases_in_flank", 0),
                detection_mode=DETECTION_MODE_TRF)

            if args.debug: print(f"TRF: checking if {tandem_repeat_allele} passes filters")
            filter_reason = check_if_tandem_repeat_allele_failed_filters(args, tandem_repeat_allele, detected_by_trf=True)

            if args.debug: print(f"TRF: {tandem_repeat_allele}, filter: {filter_reason}")
            if filter_reason is not None:
                trf_allele_filter_counters[f"TRF allele filter: {filter_reason}"] += 1
                continue

            motif_size_to_passing_trf_results[motif_size] = matching_trf_results
            motif_size_to_tandem_repeat_allele[motif_size] = tandem_repeat_allele

        if len(motif_size_to_passing_trf_results) == 0:
            if args.debug: print(f"TRF filtered out: {allele} because it has no repeat unit that passed filters")
            results.append((None, FILTER_ALLELE_INDEL_WITHOUT_REPEATS, allele))
            continue

        if args.allow_multiple_trf_results_per_locus:
            for motif_size, tr_allele in sorted(motif_size_to_tandem_repeat_allele.items()):
                results.append((tr_allele, None, allele))
        else:
            # get the entry with the smallest motif size
            best_motif_size = None
            for motif_size, matching_trf_results in sorted(motif_size_to_tandem_repeat_allele.items()):
                best_motif_size = motif_size
                break

            # if the smallest motif size is 2 or less, then pick the entry with the largest alignment score instead
            if best_motif_size <= 2 and len(motif_size_to_tandem_repeat_allele) > 1:
                # get the entry that has the max alignment score
                max_alignment_score = 0
                for motif_size, matching_trf_results in sorted(motif_size_to_passing_trf_results.items()):
                    alignment_score = 0
                    for _trf_result in matching_trf_results.values():
                        alignment_score += _trf_result["alignment_score"]
                    if alignment_score > max_alignment_score:
                        max_alignment_score = alignment_score
                        best_motif_size = motif_size

            tr_allele = motif_size_to_tandem_repeat_allele[best_motif_size]
            results.append((tr_allele, None, allele))

            """Other ideas for selecting the best TR allele:
            - drop definitions where the motif size is a larger multiple of another detected motif size ?
            - keep only definitions where the motif size is an exact multiple of the variant length or of the reference repeat length
            - keep definitions with the highest purity and/or quality score, or those above some threshold(s) 
            """

    if args.debug:
        print(f"thread {thread_id}: TRF allele filter reasons:")
        for key, count in sorted(trf_allele_filter_counters.items()):
            print(f"{count:10,d} {key}")

    return results


def compute_repeat_unit_id(canonical_repeat_unit):
    """Compute a unique identifier for a canonical repeat unit. Repeat units with the same id are treated as the same repeat unit.
    
    Args:
        canonical_repeat_unit (str): canonical repeat unit
    """
    
    if len(canonical_repeat_unit) <= 6:
        return canonical_repeat_unit
    else:
        # Intentional: for motifs longer than 6bp, only the motif length is used as the id, so any
        # two long motifs of the same length are treated as the same repeat unit when grouping loci
        # for merging. This is by design - large VNTR motifs vary between samples/assemblies, and
        # collapsing them by length keeps overlapping loci of the same period in a single merged group.
        return len(canonical_repeat_unit)


def are_repeat_units_similar(canonical_repeat_unit1, canonical_repeat_unit2):
    """Check if the two repeat units are similar enough to be considered the same repeat unit.
    
    Args:
        canonical_repeat_unit1 (str): canonical motif #1
        canonical_repeat_unit2 (str): canonical motif #2
    """
    
    return compute_repeat_unit_id(canonical_repeat_unit1) == compute_repeat_unit_id(canonical_repeat_unit2)
    

def merge_overlapping_tandem_repeat_loci(tandem_repeat_alleles, pyfaidx_fasta_obj, verbose=False):
    """Merge overlapping tandem repeats

    Args:
        tandem_repeat_alleles (list): TandemRepeatAllele or ReferenceTandemRepeat objects (the merge subcommand
            passes ReferenceTandemRepeat objects parsed from catalog BED files)
        pyfaidx_fasta_obj (pyfaidx_fasta.Fasta): reference fasta object
        verbose (bool): if True, print verbose output

    Returns:
        list: the input object unchanged for each locus that overlaps no other, plus a new ReferenceTandemRepeat
            for each group of overlapping loci that was merged
    """

    if verbose:
        print(f"Merging {len(tandem_repeat_alleles):,d} tandem repeat loci")

    before = len(tandem_repeat_alleles)

    tandem_repeat_alleles.sort(key=lambda x: (x.chrom, x.start_0based, x.end_1based, x.repeat_unit_length))
    
    # process alleles one chromosome at a time
    results = []
    for chrom, tr_alleles_for_chrom in itertools.groupby(tandem_repeat_alleles, key=lambda x: x.chrom):
        tr_allele_groups_to_merge = []
        
        motif_id_to_tr_allele_group = collections.defaultdict(list)
        motif_id_to_tr_allele_group_end_1based = collections.defaultdict(int)
        for tr_allele in tr_alleles_for_chrom:
            current_repeat_unit_id = compute_repeat_unit_id(
                compute_canonical_motif(tr_allele.repeat_unit, include_reverse_complement=True))
            if len(motif_id_to_tr_allele_group[current_repeat_unit_id]) == 0:
                motif_id_to_tr_allele_group_end_1based[current_repeat_unit_id] = tr_allele.end_1based
                motif_id_to_tr_allele_group[current_repeat_unit_id].append(tr_allele)
                continue

            # check if the current allele overlaps with the current group
            if tr_allele.start_0based <= motif_id_to_tr_allele_group_end_1based[current_repeat_unit_id] + 1:
                motif_id_to_tr_allele_group[current_repeat_unit_id].append(tr_allele)
                motif_id_to_tr_allele_group_end_1based[current_repeat_unit_id] = max(
                    motif_id_to_tr_allele_group_end_1based[current_repeat_unit_id], tr_allele.end_1based)
                continue

            # add the current group 
            tr_allele_groups_to_merge.append(motif_id_to_tr_allele_group[current_repeat_unit_id])

            # start a new group
            motif_id_to_tr_allele_group[current_repeat_unit_id] = [tr_allele]
            motif_id_to_tr_allele_group_end_1based[current_repeat_unit_id] = tr_allele.end_1based

        if len(motif_id_to_tr_allele_group) > 0:
            for current_repeat_unit_id, tr_allele_group in sorted(motif_id_to_tr_allele_group.items(), key=str):
                tr_allele_groups_to_merge.append(tr_allele_group)

        # merge tandem repeats in groups
        results_for_chrom = []
        for tr_alleles_group in tr_allele_groups_to_merge:
            #print(f"merging group_id: {group_id} which has {len(tr_alleles_with_similar_repeat_units):,d}
            #        tandem repeat alleles: {tr_alleles_with_similar_repeat_units}")
            tr_alleles_group = list(tr_alleles_group)
            if len(tr_alleles_group) == 1:
                results_for_chrom.append(tr_alleles_group[0])
                continue

            # merge the tandem repeats
            new_start_0based = min(tr_allele_i.start_0based for tr_allele_i in tr_alleles_group)
            new_end_1based = max(tr_allele_i.end_1based for tr_allele_i in tr_alleles_group)
            detection_modes = [tr_allele_i.detection_mode for tr_allele_i in tr_alleles_group if tr_allele_i.detection_mode is not None]
            if len(detection_modes) == 0:
                new_detection_mode = None
            else:
                new_detection_modes = set()
                for detection_mode in detection_modes:
                    for d in detection_mode.split(","):
                        new_detection_modes.add(d.replace("merged:", ""))
                new_detection_mode = "merged:" + ",".join(sorted(new_detection_modes))

            if new_end_1based - new_start_0based >= tr_alleles_group[0].repeat_unit_length:
                repeat_sequence = str(pyfaidx_fasta_obj[chrom][new_start_0based:new_end_1based]).upper()
                # An assembly N-gap can make the merged reference span all Ns; compute_most_common_motif rejects an
                # all-N motif, so drop the merged locus here instead of crashing (matches prior behavior where such
                # loci were absent from the combined catalog).
                if not repeat_sequence.strip("N"):
                    continue
                new_repeat_unit = compute_most_common_motif(repeat_sequence, tr_alleles_group[0].repeat_unit_length)

                simplified_motif, _, _ = find_repeat_unit_without_allowing_interruptions(new_repeat_unit, allow_partial_repeats=False)
                new_repeat_unit = simplified_motif

            else:
                new_repeat_unit = tr_alleles_group[0].repeat_unit

            merged_tr_allele = ReferenceTandemRepeat(
                chrom=chrom,
                start_0based=new_start_0based,
                end_1based=new_end_1based,
                repeat_unit=new_repeat_unit,
                detection_mode=new_detection_mode)
            results_for_chrom.append(merged_tr_allele)

        results_for_chrom.sort(key=lambda x: (x.start_0based, x.end_1based, x.repeat_unit_length))
        results.extend(results_for_chrom)

    if verbose:
        print(f"Dropped {before - len(results):,d} redundant tandem repeat loci, keeping {len(results):,d} tandem repeat loci")
    
    return results


def need_to_reprocess_allele_with_extended_flanking_sequence(tandem_repeat_allele):
    """Check if the tandem repeat allele needs to be reprocessed with an extended flanking sequence.

    Args:
        tandem_repeat_allele (TandemRepeatAllele): the tandem repeat allele to check
    """

    if tandem_repeat_allele.do_repeats_cover_entire_flanking_sequence():
        flanks_increased = False
        if tandem_repeat_allele.do_repeats_cover_entire_left_flanking_sequence() and not tandem_repeat_allele.allele.get_left_flank_stops_at_N():
            tandem_repeat_allele.allele.increase_left_flanking_sequence_size()
            flanks_increased = True
        if tandem_repeat_allele.do_repeats_cover_entire_right_flanking_sequence() and not tandem_repeat_allele.allele.get_right_flank_stops_at_N():
            tandem_repeat_allele.allele.increase_right_flanking_sequence_size()
            flanks_increased = True

        if flanks_increased \
            and tandem_repeat_allele.allele.get_expected_left_flanking_sequence_size() <= MAX_FLANKING_SEQUENCE_SIZE \
            and tandem_repeat_allele.allele.get_expected_right_flanking_sequence_size() <= MAX_FLANKING_SEQUENCE_SIZE:
            return True

    return False


def run_shell_command(command):
    """Run a shell command and raise if it fails.

    The writers below compress and index their output this way. A silent bgzip or tabix failure would leave a
    missing or unusable output file behind while the caller reports success, so the exit code is checked.

    Args:
        command (str): the shell command to run

    Raises:
        RuntimeError: if the command exits with a non-zero status
    """
    if os.system(command) != 0:
        raise RuntimeError(f"Command failed: {command}")


def write_bed(tandem_repeat_alleles, args, detailed=False):
    """Write the tandem repeat alleles to a BED file.
    
    Args:
        tandem_repeat_alleles (list): list of TandemRepeatAllele objects
        args (argparse.Namespace): command-line arguments parsed by parse_args()
        detailed (bool): if True, put extra info in the name field (ie. column 4) in addition to the repeat unit.
    """

    if detailed:
        bed_output_path = f"{args.output_prefix}.tandem_repeats.detailed.bed"
    else:
        bed_output_path = f"{args.output_prefix}.tandem_repeats.bed"


    tandem_repeat_alleles.sort(key=lambda x: (x.chrom, x.start_0based, x.end_1based, x.repeat_unit_length))

    with open(bed_output_path, "w") as f:
        for tandem_repeat_allele in tandem_repeat_alleles:
            if detailed:
                name_field = f"{tandem_repeat_allele.repeat_unit}:{tandem_repeat_allele.repeat_unit_length}bp"
                name_field += f":{(tandem_repeat_allele.end_1based - tandem_repeat_allele.start_0based)/tandem_repeat_allele.repeat_unit_length:0.1f}x"
                if tandem_repeat_allele.detection_mode is not None:
                    name_field += f":{tandem_repeat_allele.detection_mode}"
                # A reference span shorter than the motif (a zero-width locus, or a VNTR locus spanning a partial
                # copy) has no purity, and compute_repeat_purity returns nan for it
                if not math.isnan(tandem_repeat_allele.repeat_purity):
                    name_field += f":p{tandem_repeat_allele.repeat_purity:0.2f}"
            else:
                name_field = tandem_repeat_allele.repeat_unit

            f.write("\t".join(map(str, [
                tandem_repeat_allele.chrom,
                tandem_repeat_allele.start_0based,
                tandem_repeat_allele.end_1based,
                name_field,
                tandem_repeat_allele.repeat_unit_length,
            ])) + "\n")

    run_shell_command(f"bgzip -f {shlex.quote(bed_output_path)}")
    run_shell_command(f"tabix -f -p bed {shlex.quote(bed_output_path)}.gz")

    if args.verbose:
        print(f"Wrote {len(tandem_repeat_alleles):,d} tandem repeat alleles to {bed_output_path}.gz")


def write_tsv(tandem_repeat_alleles, args):
    """Write the tandem repeat alleles to a TSV file.
    
    Args:
        tandem_repeat_alleles (list): list of TandemRepeatAllele objects
        args (argparse.Namespace): command-line arguments parsed by parse_args()
    """

    tandem_repeat_alleles.sort(key=lambda x: (x.chrom, x.start_0based, x.end_1based, x.repeat_unit_length))

    tsv_output_path = f"{args.output_prefix}.tandem_repeats.tsv"

    header = [
        "Chrom",
        "Start0Based",
        "End1Based",
        "Locus",
        "LocusId",
        "INS_or_DEL",
        "Motif",
        "CanonicalMotif",
        "MotifSize",
        "NumRepeatsInReference",
        "VcfPos",
        "SummaryString",
        "IsFoundInReference",
        "IsPureRepeat",
        "DetectionMode",
    ]

    extra_header_columns = []
    if args.copy_info_field_keys_to_tsv:
        for key, count in sorted(args.copy_info_field_keys_to_tsv.items()):
            if count > 0:
                extra_header_columns.append(key)
            else:
                print(f"WARNING: INFO field key '{key}' was not found in any of rows of the input VCF. Skipping..")

    with open(tsv_output_path, "w") as f:
        f.write("\t".join(header + extra_header_columns) + "\n")
        
        for tandem_repeat_allele in tandem_repeat_alleles:
            # Handle ReferenceTandemRepeat objects which lack some TandemRepeatAllele properties
            if isinstance(tandem_repeat_allele, ReferenceTandemRepeat):
                ins_or_del = ""
                vcf_pos = ""
                is_pure_repeat = ""
                info_field_dict = {}
            else:
                ins_or_del = tandem_repeat_allele.ins_or_del
                vcf_pos = tandem_repeat_allele.allele.pos
                is_pure_repeat = tandem_repeat_allele.is_pure_repeat
                info_field_dict = tandem_repeat_allele.info_field_dict or {}

            output_row = [
                tandem_repeat_allele.chrom,
                tandem_repeat_allele.start_0based,
                tandem_repeat_allele.end_1based,
                f"{tandem_repeat_allele.chrom}:{tandem_repeat_allele.start_0based}-{tandem_repeat_allele.end_1based}",
                tandem_repeat_allele.locus_id,
                ins_or_del,
                tandem_repeat_allele.repeat_unit,
                tandem_repeat_allele.canonical_repeat_unit,
                tandem_repeat_allele.repeat_unit_length,
                tandem_repeat_allele.num_repeats_ref,
                vcf_pos,
                tandem_repeat_allele.summary_string,
                tandem_repeat_allele.end_1based > tandem_repeat_allele.start_0based,
                is_pure_repeat,
                tandem_repeat_allele.detection_mode,
            ]
            for key in extra_header_columns:
                output_row.append(info_field_dict.get(key, ""))

            f.write("\t".join(map(str, output_row)) + "\n")

    run_shell_command(f"bgzip -f {shlex.quote(tsv_output_path)}")
    if args.verbose:
        print(f"Wrote {len(tandem_repeat_alleles):,d} tandem repeat alleles to {tsv_output_path}.gz")


def write_fasta(tandem_repeat_alleles, args):
    """Write the tandem repeat alleles to a FASTA file. For insertion alleles, the alternate allele sequence is written while for deletion
    alleles, the reference allele sequence is written.
    
    Args:
        tandem_repeat_alleles (list): list of TandemRepeatAllele objects
        args (argparse.Namespace): command-line arguments parsed by parse_args()
    """

    tandem_repeat_alleles.sort(key=lambda x: (x.chrom, x.start_0based, x.end_1based, x.repeat_unit_length))

    fasta_output_path = f"{args.output_prefix}.tandem_repeats.fasta"
    with open(fasta_output_path, "w") as f:
        for tandem_repeat_allele in tandem_repeat_alleles:
            f.write(f">{tandem_repeat_allele.locus_id}__{tandem_repeat_allele.summary_string}\n")
            if tandem_repeat_allele.ins_or_del == "INS":
                f.write(f"{tandem_repeat_allele.alt_allele_repeat_sequence}\n")
            else:
                f.write(f"{tandem_repeat_allele.ref_allele_repeat_sequence}\n")

    run_shell_command(f"gzip -f {shlex.quote(fasta_output_path)}")
    if args.verbose:
        print(f"Wrote {len(tandem_repeat_alleles):,d} tandem repeat sequences to {fasta_output_path}.gz")


def get_input_vcf_iterator(args, include_header=False):
    if args.interval:
        if args.verbose:
            print(f"Parsing interval(s) {', '.join(args.interval)} from {args.input_vcf_path}")

        vcf_iterator = []
        tabix_file = pysam.TabixFile(args.input_vcf_path)
        # Normalize interval chromosomes to match the input VCF's contig naming convention so that
        # intervals work whether the VCF is indexed with 'chr'-prefixed contigs (e.g. chr1) or not (e.g. 1).
        normalize_chrom = create_normalize_chrom_function(
            has_chr_prefix=any(contig.startswith("chr") for contig in tabix_file.contigs))
        if include_header:
            vcf_iterator = (f"{line}\n" for line in tabix_file.header)

        intervals = []
        for interval_or_bed_file in args.interval:
            if ".bed" in interval_or_bed_file and file_exists(interval_or_bed_file):
                with open_file(interval_or_bed_file, is_text_file=True) as f:
                    for line in f:
                        chrom, start, end = line.strip().split("\t")[:3]
                        intervals.append(f"{chrom}:{int(start)}-{int(end)}")
            else:
                intervals.append(interval_or_bed_file)

        if intervals:
            interval_trees = collections.defaultdict(intervaltree.IntervalTree)
            for interval in intervals:
                chrom, start_0based, end = parse_interval(interval)
                chrom = normalize_chrom(chrom)
                interval_trees[chrom].addi(start_0based, end)
            normalized_intervals = []
            for chrom, interval_tree in sorted(interval_trees.items()):
                interval_tree.merge_overlaps()
                for interval in interval_tree:
                    normalized_intervals.append((chrom, interval.begin, interval.end))
            normalized_intervals.sort()

            intervals = [f"{chrom}:{start_0based}-{end}" for chrom, start_0based, end in normalized_intervals]

        def fetch_interval(tabix_file, interval):
            try:
                yield from tabix_file.fetch(interval)
            except Exception as e:
                print(f"WARNING: Unable to fetch interval {interval}: {e}. Skipping..")

        vcf_iterator = itertools.chain(
            vcf_iterator,
            (line for interval in intervals for line in fetch_interval(tabix_file, interval)),
        )
    else:
        if args.verbose:
            print(f"Parsing {args.input_vcf_path}")
        vcf_iterator = open_file(args.input_vcf_path, is_text_file=True)

    return vcf_iterator


def format_filter_id(filter_reason):
    """Turn a filter reason into a VCF FILTER ID.

    The FILTER column holds IDs declared in ##FILTER header lines, and an ID may not contain whitespace,
    semicolons or commas, while the filter reasons here read like "INDEL > 100,000bp" or "contains < 3 full
    repeats". The reason itself goes in the header line's Description.

    Args:
        filter_reason (str): one of the FILTER_* values, or "SNV", "MNV" or "not_TR"

    Returns:
        str: the reason with "<" and ">" spelled out and every other run of disallowed characters replaced by "_"
    """
    return re.sub("[^A-Za-z0-9_.]+", "_", filter_reason.replace("<", "lt").replace(">", "gt")).strip("_")


def write_vcf(tandem_repeat_alleles, args, only_write_filtered_out_alleles=False, filtered_alleles=None):
    """Write variants that either are or aren't tandem repeats to a VCF file.

    Args:
        tandem_repeat_alleles (list): list of TandemRepeatAllele objects
        args (argparse.Namespace): command-line arguments parsed by parse_args()
        only_write_filtered_out_alleles (bool): if True, only write the variants that are not in the tandem_repeats_alleles list
        filtered_alleles (dict): optional dict mapping (chrom, pos, ref, alt) to filter reason strings.
            Used when only_write_filtered_out_alleles=True to populate the FILTER column. A multi-allelic record
            whose ALT alleles were filtered for different reasons lists each reason once, ";"-separated.
    """

    vcf_iterator = get_input_vcf_iterator(args, include_header=True)

    # iterate over all VCF rows
    if only_write_filtered_out_alleles:
        output_vcf_path = f"{args.output_prefix}.not_tandem_repeats.vcf"
    else:
        output_vcf_path = f"{args.output_prefix}.tandem_repeats.vcf"

    tandem_repeat_alleles.sort(key=lambda x: x.order) # sort into their original order

    # Key by (chrom, pos, ref, alt) so that multiple ALT alleles at the same multi-allelic site are
    # tracked independently instead of collapsing to a single allele.
    vcf_alleles = {
        (tr.chrom, tr.allele.pos, tr.allele.ref, tr.allele.alt): tr for tr in tandem_repeat_alleles
    }

    # Declare the custom INFO fields appended below. Number=. (variable) because each is a
    # comma-separated list over only the tandem-repeat ALT alleles at a site, not all ALTs.
    # END is a VCF-reserved INFO key, so the tandem-repeat end coordinate is written as TR_END.
    tr_info_header_lines = [
        '##INFO=<ID=MOTIF,Number=.,Type=String,Description="Repeat unit of each tandem-repeat ALT allele at this site">',
        '##INFO=<ID=MOTIF_SIZE,Number=.,Type=Integer,Description="Repeat unit length of each tandem-repeat ALT allele">',
        '##INFO=<ID=START_0BASED,Number=.,Type=Integer,Description="0-based start coordinate of each tandem repeat">',
        '##INFO=<ID=TR_END,Number=.,Type=Integer,Description="1-based end coordinate of each tandem repeat">',
        '##INFO=<ID=DETECTED,Number=.,Type=String,Description="Motif detection method for each tandem-repeat ALT allele">',
    ]

    # The filtered-out VCF puts each record's filter reason in its FILTER column. Every ID written there must
    # be declared in the header, and IDs may not contain spaces, so the reasons are declared up front, each
    # with the readable reason as its Description. The three fallbacks cover records no filter reason was
    # recorded for.
    if only_write_filtered_out_alleles:
        filter_reasons = sorted(set(filtered_alleles.values()) if filtered_alleles else set())
        for filter_reason in filter_reasons + ["SNV", "MNV", "not_TR"]:
            tr_info_header_lines.append(
                f'##FILTER=<ID={format_filter_id(filter_reason)},Description="{filter_reason}">')

    with open(output_vcf_path, "w") as f:
        vcf_line_i = 0
        output_line_counter = 0
        for line in vcf_iterator:
            if line.startswith("#"):
                # Insert the custom INFO declarations just before the #CHROM column header line
                if line.startswith("#CHROM"):
                    for info_line in tr_info_header_lines:
                        f.write(info_line + "\n")
                f.write(line)
                continue

            vcf_fields = line.strip().split("\t")
            if vcf_line_i < args.offset:
                vcf_line_i += 1
                continue

            if args.n is not None and vcf_line_i >= args.offset + args.n:
                break

            vcf_line_i += 1

            # parse the ALT allele(s)
            vcf_chrom = vcf_fields[0]
            vcf_pos = int(vcf_fields[1])
            vcf_ref = vcf_fields[3].upper()

            # Evaluate tandem-repeat status per ALT allele so that multi-allelic sites with a mix of
            # tandem-repeat and non-tandem-repeat ALT alleles are handled correctly.
            alt_alleles = [a for a in vcf_fields[4].upper().split(",") if a != "*"]
            tr_alleles_at_site = [
                vcf_alleles[(vcf_chrom, vcf_pos, vcf_ref, alt)]
                for alt in alt_alleles
                if (vcf_chrom, vcf_pos, vcf_ref, alt) in vcf_alleles
            ]
            has_tandem_repeat_alt = len(tr_alleles_at_site) > 0
            has_non_tandem_repeat_alt = len(tr_alleles_at_site) < len(alt_alleles)

            if has_tandem_repeat_alt:
                # append TR info to the INFO field, aggregating across all tandem-repeat ALT alleles at this site
                if vcf_fields[7] == "." or not vcf_fields[7]:
                    vcf_fields[7] = ""
                else:
                    vcf_fields[7] += ";"
                vcf_fields[7] += "MOTIF=" + ",".join(tr.repeat_unit for tr in tr_alleles_at_site)
                vcf_fields[7] += ";MOTIF_SIZE=" + ",".join(str(tr.repeat_unit_length) for tr in tr_alleles_at_site)
                vcf_fields[7] += ";START_0BASED=" + ",".join(str(tr.start_0based) for tr in tr_alleles_at_site)
                vcf_fields[7] += ";TR_END=" + ",".join(str(tr.end_1based) for tr in tr_alleles_at_site)
                vcf_fields[7] += ";DETECTED=" + ",".join(tr.detection_mode for tr in tr_alleles_at_site)

            if only_write_filtered_out_alleles:
                # Write the record if it has at least one non-tandem-repeat ALT allele (its filter reason
                # would otherwise be lost). Records where every ALT is a tandem repeat are skipped here.
                # Intentional: for mixed multi-allelic sites (some TR and some non-TR ALTs) the original
                # record is written unchanged - we do not split it per-ALT or strip the tandem-repeat ALT
                # and its appended TR INFO. Keeping the whole record preserves the site's genotype fields
                # and the filtered-out VCF is only meant to record which sites had a non-TR allele.
                if not has_non_tandem_repeat_alt:
                    continue
                # set the FILTER column to the filter reason(s), one per filtered ALT allele
                filter_ids = list(dict.fromkeys(
                    format_filter_id(filtered_alleles[(vcf_chrom, vcf_pos, vcf_ref, alt)])
                    for alt in alt_alleles if filtered_alleles and (vcf_chrom, vcf_pos, vcf_ref, alt) in filtered_alleles))
                if filter_ids:
                    vcf_fields[6] = ";".join(filter_ids)
                else:
                    # determine filter reason for non-indel variants
                    if all(len(vcf_ref) == len(a) for a in alt_alleles):
                        vcf_fields[6] = "SNV" if all(len(a) == 1 for a in alt_alleles) else "MNV"
                    else:
                        vcf_fields[6] = "not_TR"
                f.write("\t".join(vcf_fields) + "\n")
                output_line_counter += 1
            elif has_tandem_repeat_alt:
                f.write("\t".join(vcf_fields) + "\n")
                output_line_counter += 1

    run_shell_command(f"bgzip -f {shlex.quote(output_vcf_path)}")
    run_shell_command(f"tabix -f -p vcf {shlex.quote(output_vcf_path)}.gz")

    if args.verbose:
        print(f"Wrote {output_line_counter:,d} variants to {output_vcf_path}.gz")
    

def print_stats(counters):
    """Print out all the counters"""

    key_prefixes = set()
    for key, _ in counters.items():
        tokens = key.split(":")
        key_prefixes.add(f"{tokens[0]}:")

    for print_totals_only in True, False:
        for key_prefix in sorted(key_prefixes):
            if print_totals_only ^ (key_prefix in ("variant counts:", "allele counts:")):
                continue

            current_counter = [(key, count) for key, count in counters.items() if key.startswith(key_prefix)]
            current_counter = sorted(current_counter, key=lambda x: (-x[1], x[0]))
            if current_counter:
                print("-"*15)
            for key, value in current_counter:
                if key_prefix.startswith("TR"):
                    total_key = "TR variant counts: TOTAL" if "variant" in key_prefix else "TR allele counts: TOTAL"
                else:
                    total_key = "variant counts: TOTAL variants" if "variant" in key_prefix else "allele counts: TOTAL alleles"

                total = counters[total_key]
                percent = f"{100*value / total:5.1f}%" if total > 0 else ""

                if print_totals_only:
                    print(f"{value:10,d}  {key}")
                else:
                    print(f"{value:10,d} out of {total:10,d} ({percent}) {key}")


def do_merge_subcommand(args):
    """Merge tandem repeat catalogs from two or more input BED files."""

    fasta_obj = pyfaidx.Fasta(args.reference_fasta_path, one_based_attributes=False, as_raw=True)

    all_trs = []

    if not args.output_prefix:
        if len(args.input_bed_paths) == 1:
            args.output_prefix = re.sub(".bed(.gz|.bgz)$", "", args.input_bed_paths[0]).replace(".tandem_repeats", "").replace(".detailed", "") + ".merged"
        else:
            args.output_prefix = f"combined.{len(args.input_bed_paths)}_catalogs"

    simplified_repeat_units_counter = 0
    input_bed_paths_iterator = args.input_bed_paths if not args.show_progress_bar else tqdm.tqdm(args.input_bed_paths, unit=" catalog")
    for path_i, input_bed_path in enumerate(input_bed_paths_iterator):
        if args.verbose:
            print("-"*100)

        input_files_to_close = []
        if args.interval:
            if args.verbose:
                print(f"Parsing {', '.join(args.interval)} from catalog #{path_i + 1}: {input_bed_path}")

            tabix_file = pysam.TabixFile(input_bed_path)
            bed_iterator = fetch_catalog_records_within_intervals(tabix_file, args.interval, input_bed_path)
            input_files_to_close.append(tabix_file)
        else:
            if args.verbose:
                print(f"Parsing catalog #{path_i + 1}: {input_bed_path}")
            bed_iterator = open_file(input_bed_path, is_text_file=True)
            input_files_to_close.append(bed_iterator)

        # parse the BED file into a list of ReferenceTandemRepeat objects
        current_catalog_trs = []
        for line_num, line in enumerate(bed_iterator, start=1):
            parsed_line = parse_catalog_bed_line(line, line_num, input_bed_path)
            if parsed_line is None:
                continue
            fields, repeat_unit = parsed_line
            name_field_tokens = fields[3].split(":")

            # A detailed BED name field reads "CAG:3bp:19.0x:pure:p0.95". The detection mode is everything between
            # the third token and the trailing purity token: write_bed leaves it out when a locus has none, and a
            # merged locus's mode ("merged:pure,trf") contains the ":" separator itself. Older detailed BEDs wrote
            # "pnan" for a locus shorter than its motif, which write_bed now leaves out.
            detection_mode = None
            if args.write_detailed_bed and len(name_field_tokens) >= 4:
                has_purity_token = re.match(r"^p(\d|nan)", name_field_tokens[-1]) is not None
                detection_mode = ":".join(name_field_tokens[3:-1] if has_purity_token else name_field_tokens[3:]) or None

            # check if the repeat unit itself consists of perfect repeats of a smaller repeat unit (this happens in ~3% of TRs detected by TRF)
            simplified_repeat_unit, _, _ = find_repeat_unit_without_allowing_interruptions(repeat_unit, allow_partial_repeats=False)
            if len(simplified_repeat_unit) != len(repeat_unit):
                simplified_repeat_units_counter += 1

            current_catalog_trs.append(ReferenceTandemRepeat(
                chrom=fields[0],
                start_0based=int(fields[1]),
                end_1based=int(fields[2]),
                repeat_unit=simplified_repeat_unit,
                detection_mode=detection_mode,
            ))
        
        if args.verbose:
            print_tr_stats(current_catalog_trs, title=f"Catalog #{path_i+1}: {input_bed_path}")

        all_trs.extend(current_catalog_trs)
        if len(all_trs) > args.batch_size or path_i == len(args.input_bed_paths) - 1:
            if args.verbose:
                print("="*100)

            all_trs = merge_overlapping_tandem_repeat_loci(all_trs, fasta_obj, verbose=args.verbose)

        for input_file in input_files_to_close:
            input_file.close()


    if args.verbose:
        if simplified_repeat_units_counter:
            print(f"Simplified {simplified_repeat_units_counter:,d} out of {len(all_trs):,d} ({100*simplified_repeat_units_counter/len(all_trs):5.1f}%) repeat units")
        print_tr_stats(all_trs, title=f"Merged catalog stats: ")

    write_bed(all_trs, args)    

    if args.write_detailed_bed:
        for tr in all_trs:
            repeat_sequence = str(fasta_obj[tr.chrom][tr.start_0based:tr.end_1based]).upper()
            tr.repeat_purity, _ = compute_repeat_purity(
                repeat_sequence, tr.repeat_unit, include_partial_repeats=True)

        write_bed(all_trs, args, detailed=True)


def print_tr_stats(tandem_repeat_alleles, title=None):
    """Print statistics about the tandem repeat alleles."""

    counters = collections.defaultdict(int)
    for tandem_repeat_allele in tandem_repeat_alleles:
        counters[f"total"] += 1
        ru_len = tandem_repeat_allele.repeat_unit_length
        if ru_len <= 6:
            counters[f"STRs"] += 1
            counters[f"STR{ru_len}"] += 1
        else:
            counters[f"VNTRs"] += 1

        if tandem_repeat_allele.detection_mode is not None:
            counters[f"detection_mode: {tandem_repeat_allele.detection_mode}"] += 1

    print("-"*15)
    if title:
        print(title)
    
    if counters['total'] > 0:
        for key, count in sorted(counters.items(), key=lambda x: (-x[1], x[0])):
            if key.startswith("detection_mode:"):
                print(f"{count:10,d} ({100*count/counters['total']:5.1f}%) {key}")

        print("-"*15)

        for ru_len in range(1, 7):
            print(f"{counters[f'STR{ru_len}']:10,d} ({100*counters[f'STR{ru_len}']/counters['total']:5.1f}%) {ru_len}bp motifs")

        str_stats = f"{counters['STRs']:10,d} ({100*counters['STRs']/counters['total']:5.1f}%)"
        vntr_stats = f"{counters['VNTRs']:10,d} ({100*counters['VNTRs']/counters['total']:5.1f}%)"
        print(f"{vntr_stats} 7+bp motifs")
        print("-"*15)
        print(f"{str_stats} STRs")
        print(f"{vntr_stats} VNTRs")

    print(f"{counters['total']:10,d} total TRs")


def compute_chrom_sort_key(chrom):
    """Return a sort key that orders chromosomes naturally and keeps every contig's records contiguous.

    Sorting on the raw chromosome name puts chr10 before chr2, so the genotype outputs would not be in genomic
    order. Unplaced, alt and decoy contigs have no natural rank, so they all sort after the standard
    chromosomes and are then ordered by name. Grouping them by name matters because tabix rejects a VCF whose
    records for one contig are not contiguous, which is what happens when every unplaced contig shares a rank
    and the records interleave by position.

    Args:
        chrom (str): chromosome name, with or without a "chr" prefix

    Returns:
        tuple: (rank, chrom) ordering key
    """
    chrom_without_prefix = chrom[3:] if chrom.startswith("chr") else chrom

    if chrom_without_prefix == "X":
        return 23, chrom
    if chrom_without_prefix == "Y":
        return 24, chrom
    if chrom_without_prefix in ("M", "MT"):
        return 25, chrom

    try:
        return int(chrom_without_prefix), chrom
    except ValueError:
        return 100, chrom  # unplaced/alt/decoy contigs sort last, grouped by name


def is_locus_genotyped_as_reference(genotyped_locus):
    """Whether every allele of a locus was genotyped as equal to the reference, so --skip-hom-ref-loci can drop it.

    A locus with no overlapping variants (after get_overlapping_vcf_variants drops hom-ref records and records
    that change nothing inside the locus) gets its alleles built from the reference: two at a diploid locus
    (HOM), one at a haploid locus such as chrX outside the PARs in a male (HEMI). A locus on a contig the VCF
    lacks also lands here, since the VCF's silence is read as "no variants".

    A locus set to no call before any variant was applied (its contig is missing from the reference FASTA, or
    it extends past the end of its contig) also has zero overlapping variants, but its alleles are unknown
    rather than reference, so it is not genotyped as reference and must survive into the outputs.

    Args:
        genotyped_locus (GenotypedTandemRepeat): the locus to test

    Returns:
        bool: True if the locus has no overlapping variants and was not set to no call
    """
    return genotyped_locus.num_overlapping_variants == 0 and genotyped_locus.no_call_reason is None


def write_genotypes_tsv(genotyped_loci, args, motif_lists_by_locus=None):
    """Write genotyped TR loci to a TSV file.

    This function writes the genotyped tandem repeat loci to a gzip-compressed
    TSV file with columns defined by GENOTYPE_TSV_OUTPUT_COLUMNS.

    Args:
        genotyped_loci (list): List of GenotypedTandemRepeat objects
        args (argparse.Namespace): Command-line arguments. Must have:
            - output_prefix (str): Prefix for output file path
            - verbose (bool): If True, print detailed output
            - skip_hom_ref_loci (bool): If True, skip loci with no overlapping variants.
        motif_lists_by_locus (dict): Optional dict mapping locus_id to {'allele1': [...], 'allele2': [...]}
            ordered motif lists (as returned by compute_motif_composition). If provided, the parsed motif
            sequence is written for each allele.

    Returns:
        str: Path to the output TSV file
    """
    # Get skip_hom_ref_loci with default of False
    skip_hom_ref_loci = getattr(args, 'skip_hom_ref_loci', False)

    # Filter loci based on skip_hom_ref_loci
    if skip_hom_ref_loci:
        filtered_loci = [locus for locus in genotyped_loci if not is_locus_genotyped_as_reference(locus)]
        if args.verbose:
            print(f"Skipped {len(genotyped_loci) - len(filtered_loci):,d} homozygous reference loci")
    else:
        filtered_loci = genotyped_loci

    # Sort by genomic position and motif, ordering chromosomes naturally so this matches the VCF output
    filtered_loci.sort(key=lambda x: (compute_chrom_sort_key(x.chrom), x.start_0based, x.end, x.motif))

    # Construct output filename
    tsv_output_path = f"{args.output_prefix}.tandem_repeat_genotypes.tsv.gz"

    # Write to gzip file
    with gzip.open(tsv_output_path, "wt") as f:
        # Write header
        f.write("\t".join(GENOTYPE_TSV_OUTPUT_COLUMNS) + "\n")

        # Write each genotyped locus
        for genotyped_locus in filtered_loci:
            tsv_dict = genotyped_locus.to_tsv_dict(
                motif_lists=motif_lists_by_locus.get(genotyped_locus.locus_id) if motif_lists_by_locus else None)
            row_values = [str(tsv_dict.get(col, "")) for col in GENOTYPE_TSV_OUTPUT_COLUMNS]
            f.write("\t".join(row_values) + "\n")

    if args.verbose:
        print(f"Wrote {len(filtered_loci):,d} genotyped TR loci to {tsv_output_path}")

    return tsv_output_path


def build_basic_split_motif_entry(sequence, motif_size):
    """Split a sequence into motif-sized chunks, keeping any trailing partial chunk as a suffix.

    Unlike a bare split_sequence_into_motifs() call, the trailing remainder is preserved (as the suffix)
    instead of being dropped.

    Args:
        sequence (str): The allele nucleotide sequence to split.
        motif_size (int): The motif size in base pairs.

    Returns:
        dict: A parsed-motif entry {"motifs": [...], "prefix": "", "suffix": "..."}, or None if the
            sequence is empty. For example, "CAGCAGCA" with motif_size 3 returns
            {"motifs": ["CAG", "CAG"], "prefix": "", "suffix": "CA"} which renders as "[CAG][CAG]CA".
    """
    if not sequence:
        return None
    motifs = split_sequence_into_motifs(sequence, motif_size)
    return {"motifs": motifs, "prefix": "", "suffix": sequence[len(motifs) * motif_size:]}


def build_trviz_motif_entry(sequence, motif, decomposer):
    """Split a sequence into motifs using the trviz decomposition algorithm.

    trviz aligns the sequence to repeated copies of the motif, so a copy may be longer or shorter than the motif
    where the allele has an insertion or deletion. A first or last piece shorter than the motif is kept as the
    prefix or suffix instead of being reported as a motif. The pieces are cut from the original sequence (not from
    trviz's uppercased copy) so the entry always reconstructs the allele exactly.

    Args:
        sequence (str): The allele nucleotide sequence to split.
        motif (str): The annotated locus motif.
        decomposer (trviz.decomposer.Decomposer): trviz decomposer instance.

    Returns:
        2-tuple (dict, str): A parsed-motif entry {"motifs": [...], "prefix": str, "suffix": str} and the method
            that produced it. Sequences or motifs containing bases other than A, C, G or T (eg. the degenerate
            GCN motif of a polyalanine repeat), which trviz rejects, fall back on the basic chunking method.
            Returns (None, None) if the sequence is empty. For example,
            "AGCAGCAGCA" with motif "CAG" returns
            ({"motifs": ["CAG", "CAG"], "prefix": "AG", "suffix": "CA"}, "trviz") which renders as
            "AG[CAG][CAG]CA".
    """
    if not sequence:
        return None, None
    if not set(sequence.upper()) <= set("ACGT") or not set(motif.upper()) <= set("ACGT"):
        return build_basic_split_motif_entry(sequence, len(motif)), MOTIF_DETECTION_METHOD_BASIC_SPLIT

    pieces = []
    offset = 0
    for piece in decomposer.decompose(sequence, [motif]):
        pieces.append(sequence[offset:offset + len(piece)])
        offset += len(piece)
    if offset != len(sequence):
        raise ValueError(f"trviz decomposition of {sequence} covers {offset} bases instead of {len(sequence)}")

    prefix = pieces.pop(0) if len(pieces) > 1 and len(pieces[0]) < len(motif) else ""
    suffix = pieces.pop() if pieces and len(pieces[-1]) < len(motif) else ""
    return {"motifs": pieces, "prefix": prefix, "suffix": suffix}, MOTIF_DETECTION_METHOD_TRVIZ


def compute_motif_lists_with_trviz(genotyped_loci, verbose=False):
    """Parse each allele's sequence into an ordered list of motifs using the trviz decomposition algorithm.

    trviz is imported here rather than at the top of the module, so that it is only required when
    --add-motif-composition trviz is used.

    Args:
        genotyped_loci (list): List of GenotypedTandemRepeat objects
        verbose (bool): If True, print progress information

    Returns:
        dict: Dictionary mapping locus_id to a dict with the same keys as compute_motif_composition() returns.
    """
    from trviz.decomposer import Decomposer

    if verbose:
        print("Parsing allele sequences into motifs using trviz...")
    decomposer = Decomposer()
    result = {}
    for locus in genotyped_loci:
        allele1_entry, allele1_method = build_trviz_motif_entry(locus.allele1_sequence, locus.motif, decomposer)
        allele2_entry, allele2_method = build_trviz_motif_entry(locus.allele2_sequence, locus.motif, decomposer)
        result[locus.locus_id] = {
            "allele1": allele1_entry,
            "allele2": allele2_entry,
            "allele1_method": allele1_method,
            "allele2_method": allele2_method,
        }
    return result


def format_motif_entry_as_sequence_string(entry):
    """Render a parsed-motif entry as a bracketed sequence string (eg. "CA[GCA][GCA][GCC]G").

    Args:
        entry (dict): Parsed-motif entry with keys "motifs" (list), "prefix" (str), and "suffix" (str),
            or None for a missing/empty allele.

    Returns:
        str: The rendered string, or None if entry is None or empty.
    """
    if not entry:
        return None
    return format_motifs_as_sequence_string(entry["motifs"], prefix=entry["prefix"], suffix=entry["suffix"])


def compute_motif_composition(genotyped_loci, args):
    """Compute the ordered list of motifs parsed from each allele's sequence, for each locus.

    Uses the basic chunking method (str_analysis.utils.find_motif_utils.split_sequence_into_motifs),
    TandemRepeatsFinder (TRF), or trviz, depending on args.add_motif_composition.

    Args:
        genotyped_loci (list): List of GenotypedTandemRepeat objects
        args: Argument namespace with attributes:
            - add_motif_composition (str or None): "basic", "trf", "trviz", or None (in which case nothing is
                computed)
            - trf_executable_path (str): Path to the TRF executable (required if add_motif_composition == "trf")
            - min_allele_length_for_trf_motif_splitting (int): Optional. Allele sequences shorter than
                max(this value, 2 * motif_size) are split using the basic method instead of TRF. Defaults to 12.
            - trf_threads (int): Optional. Number of TRF instances to run in parallel. Defaults to 1.
            - skip_hom_ref_loci (bool): If True, loci with no overlapping variants are skipped (matching the
                output writers) so that no motif composition is computed for loci that won't be written.
            - verbose (bool): If True, print progress information

    Returns:
        dict: Dictionary mapping locus_id to a dict with keys 'allele1' and 'allele2' (each a parsed-motif
            entry {"motifs": [...], "prefix": str, "suffix": str}, or None) and
            'allele1_method'/'allele2_method' (each naming the method that produced that allele's motifs:
            "trf", "trviz", "basic-split", or None). Returns an empty dict if args.add_motif_composition is not
            set. Example:
            {
                "chr1-100-150-CAG": {
                    "allele1": {"motifs": ["CAG", "CAG", "CCG", "CAG"], "prefix": "", "suffix": ""},
                    "allele2": {"motifs": ["CAG", "CAG"], "prefix": "", "suffix": "CA"},
                    "allele1_method": "basic-split",
                    "allele2_method": "basic-split",
                }
            }
    """
    if not args.add_motif_composition:
        return {}

    # Skip loci with no overlapping variants if requested, so motif composition (including TRF, which is
    # expensive) is not computed for loci that the output writers will omit.
    if getattr(args, "skip_hom_ref_loci", False):
        genotyped_loci = [locus for locus in genotyped_loci if not is_locus_genotyped_as_reference(locus)]

    if args.add_motif_composition == "basic":
        if args.verbose:
            print("Parsing allele sequences into motifs using the basic chunking method...")
        result = {}
        for locus in genotyped_loci:
            allele1_entry = build_basic_split_motif_entry(locus.allele1_sequence, locus.motif_size)
            allele2_entry = build_basic_split_motif_entry(locus.allele2_sequence, locus.motif_size)
            result[locus.locus_id] = {
                "allele1": allele1_entry,
                "allele2": allele2_entry,
                "allele1_method": MOTIF_DETECTION_METHOD_BASIC_SPLIT if allele1_entry else None,
                "allele2_method": MOTIF_DETECTION_METHOD_BASIC_SPLIT if allele2_entry else None,
            }
        return result

    if args.add_motif_composition == "trviz":
        return compute_motif_lists_with_trviz(genotyped_loci, verbose=args.verbose)

    return compute_motif_lists_with_trf(
        genotyped_loci, args.trf_executable_path,
        min_allele_length_for_trf=getattr(args, "min_allele_length_for_trf_motif_splitting", 12),
        num_threads=getattr(args, "trf_threads", 1),
        verbose=args.verbose)


def run_trf_motif_splitting(sequences_for_trf, trf_executable_path, trf_working_dir, thread_id=0):
    """Run TRF on a batch of allele sequences and split each into an ordered list of motifs.

    This is the per-thread worker used by compute_motif_lists_with_trf(). It writes the batch to its own
    FASTA file (named using thread_id so concurrent batches don't collide), runs TRF once on that file,
    and parses the result for each sequence. A TRF result is only accepted if it decomposes the allele
    using the annotated locus motif size and covers essentially the entire sequence (starts within one
    motif of the start and ends within one motif of the end); among accepted results, the one with the
    highest alignment score is used. Sequences with no qualifying TRF result fall back on the basic
    chunking method (build_basic_split_motif_entry).

    Args:
        sequences_for_trf (list): List of (seq_id, sequence, motif_size) tuples, where seq_id is
            "{locus_id}${allele_key}".
        trf_executable_path (str): Path to the TRF executable.
        trf_working_dir (str): Directory where this batch's TRF input/output files are written.
        thread_id (int): ID of the thread processing this batch, used to give the input FASTA a unique
            filename.

    Returns:
        list of 4-tuples: (locus_id, allele_key, entry, method), where entry is a parsed-motif dict
            ({"motifs": [...], "prefix": str, "suffix": str}) and method is "trf" or "basic-split".
    """
    if not sequences_for_trf:
        return []

    # Create TRFRunner with html_mode=True to get individual motif copies
    trf_runner = TRFRunner(trf_executable_path, html_mode=True)

    # Write this batch's sequences to a FASTA file unique to this thread
    trf_fasta_path = os.path.join(trf_working_dir, f"trf_motif_splitting_input__thread{thread_id}.fasta")
    with open(trf_fasta_path, "wt") as f:
        for seq_id, sequence, _ in sequences_for_trf:
            f.write(f">{seq_id}\n{sequence}\n")

    # Determine max period based on the longest sequence in this batch
    max_period = max(1, min(max(len(seq) for _, seq, _ in sequences_for_trf) // 2, 2000))

    # Run TRF once on this batch's FASTA
    trf_runner.run_trf_on_fasta_file(trf_fasta_path, max_period=max_period)

    total_sequences = len(sequences_for_trf)
    batch_results = []
    for seq_idx, (seq_id, sequence, motif_size) in enumerate(sequences_for_trf, start=1):
        # Parse the locus_id and allele from the sequence ID
        locus_id, allele_key = seq_id.rsplit("$", 1)

        # Parse TRF HTML results for this sequence
        trf_records = trf_runner.parse_html_results(
            trf_fasta_path,
            sequence_number=seq_idx,
            total_sequences=total_sequences,
            max_period=max_period,
        )

        # Only accept a TRF result that decomposes the allele using the annotated locus motif size and
        # covers essentially the entire sequence (starts within one motif of the start and ends within
        # one motif of the end). Among accepted results, prefer the one with the highest alignment score.
        accepted_records = [
            record for record in (trf_records or [])
            if record["repeat_unit_length"] == motif_size
            and record["start_0based"] < motif_size
            and record["end_1based"] > len(sequence) - motif_size
        ]

        entry = None
        if accepted_records:
            best_record = max(accepted_records, key=lambda x: x.get("alignment_score", 0))

            # Keep the ordered list of motifs from the 'repeats' list, dropping the dashes TRF writes where
            # the allele has a deletion relative to the consensus motif
            motifs = [m.replace("-", "") for m in best_record.get("repeats", [])]
            motifs = [m for m in motifs if m]
            if motifs:
                # Keep any bases before the first or after the last detected repeat as prefix/suffix
                candidate_entry = {
                    "motifs": motifs,
                    "prefix": sequence[:best_record["start_0based"]],
                    "suffix": sequence[best_record["end_1based"]:],
                }
                # Only accept a decomposition that puts back exactly the allele it came from. Anything else
                # means the parse dropped or duplicated bases, and reporting motifs that don't reconstruct the
                # allele is worse than falling back to the basic chunking method below.
                if candidate_entry["prefix"] + "".join(motifs) + candidate_entry["suffix"] == sequence:
                    entry = candidate_entry

        if entry:
            batch_results.append((locus_id, allele_key, entry, MOTIF_DETECTION_METHOD_TRF))
        else:
            # No qualifying TRF result for this allele - fall back to basic chunking rather than emitting None
            batch_results.append(
                (locus_id, allele_key, build_basic_split_motif_entry(sequence, motif_size),
                 MOTIF_DETECTION_METHOD_BASIC_SPLIT))

    return batch_results


def compute_motif_lists_with_trf(genotyped_loci, trf_executable_path, min_allele_length_for_trf=12,
                                 num_threads=1, verbose=False):
    """Parse each allele's sequence into an ordered list of motifs using TandemRepeatsFinder.

    This function runs TRF on allele sequences that are long enough to benefit from TRF analysis. A TRF
    result is only used if it decomposes the allele using the annotated locus motif size and covers
    essentially the entire sequence (it starts within one motif of the sequence start and ends within
    one motif of the sequence end). Any bases before the first or after the last detected repeat are kept
    as the entry's prefix/suffix. Sequences too short for TRF, and sequences where no TRF result qualifies,
    fall back on the basic chunking method (build_basic_split_motif_entry).

    The threshold for using TRF is: sequence length >= max(12bp, 2 * motif_size).
    This ensures TRF has enough sequence to work with (at least 2 repeats).

    Args:
        genotyped_loci (list): List of GenotypedTandemRepeat objects
        trf_executable_path (str): Path to the TRF executable
        min_allele_length_for_trf (int): Minimum allele sequence length to use TRF.
            Sequences shorter than max(min_allele_length_for_trf, 2 * motif_size)
            fall back to basic chunking. Default is 12.
        num_threads (int): Number of TRF instances to run in parallel (one per thread). Default is 1.
        verbose (bool): If True, print progress information

    Returns:
        dict: Dictionary mapping locus_id to a dict with keys 'allele1' and 'allele2' (each a parsed-motif
            entry {"motifs": [...], "prefix": str, "suffix": str}, or None for missing/empty alleles) and
            'allele1_method'/'allele2_method' (each naming the method that produced that allele's motifs:
            "trf", "basic-split", or None). Example:
            {
                "chr1-100-150-CAG": {
                    "allele1": {"motifs": ["GCA", "GCA", "GCC"], "prefix": "CA", "suffix": "G"},
                    "allele2": {"motifs": ["CAG", "CAG"], "prefix": "", "suffix": "CA"},
                    "allele1_method": "trf",
                    "allele2_method": "basic-split",
                }
            }
    """
    results = collections.defaultdict(dict)

    # Separate alleles into TRF vs basic based on sequence length
    # Use TRF if sequence length >= max(12bp, 2 * motif_size)
    sequences_for_trf = []  # (seq_id, sequence, motif_size)
    for locus in genotyped_loci:
        # Initialize both allele keys to None so loci with missing/empty allele sequences still appear in the
        # results (with None values), matching the basic method and ensuring the output writers emit the
        # motif fields (as null) rather than omitting them.
        results[locus.locus_id] = {
            "allele1": None, "allele2": None, "allele1_method": None, "allele2_method": None}
        for allele_key, sequence in [("allele1", locus.allele1_sequence),
                                     ("allele2", locus.allele2_sequence)]:
            if not sequence:
                continue

            if len(sequence) >= max(min_allele_length_for_trf, 2 * locus.motif_size):
                sequences_for_trf.append((f"{locus.locus_id}${allele_key}", sequence, locus.motif_size))
            else:
                results[locus.locus_id][allele_key] = build_basic_split_motif_entry(sequence, locus.motif_size)
                results[locus.locus_id][f"{allele_key}_method"] = MOTIF_DETECTION_METHOD_BASIC_SPLIT

    # If no sequences need TRF, return early
    if not sequences_for_trf:
        return results

    if verbose:
        print(f"Running TRF on {len(sequences_for_trf)} allele sequences "
              f"using {max(1, min(num_threads, len(sequences_for_trf)))} thread(s)...")

    # Distribute the sequences across threads, running one TRF instance per thread in its own FASTA file.
    trf_working_dir = tempfile.mkdtemp(prefix="motif_composition_trf__")
    try:
        for locus_id, allele_key, entry, method in run_trf_batches_in_parallel(
                sequences_for_trf, num_threads,
                lambda batch, thread_i: run_trf_motif_splitting(
                    batch, trf_executable_path, trf_working_dir, thread_i)):
            results[locus_id][allele_key] = entry
            results[locus_id][f"{allele_key}_method"] = method
    finally:
        shutil.rmtree(trf_working_dir, ignore_errors=True)

    return results


def write_genotypes_json(genotyped_loci, args, motif_lists_by_locus=None):
    """Write genotyped TR loci to a JSON file.

    This function writes the genotyped tandem repeat loci to a gzip-compressed
    JSON file. Unlike the TSV output, JSON preserves native types (int, float,
    bool, list) and can optionally include motif composition data.

    Args:
        genotyped_loci (list): List of GenotypedTandemRepeat objects
        args (argparse.Namespace): Command-line arguments. Must have:
            - output_prefix (str): Prefix for output file path
            - verbose (bool): If True, print detailed output
            - skip_hom_ref_loci (bool): If True, skip loci with no overlapping variants.
        motif_lists_by_locus (dict): Optional dict mapping locus_id to {'allele1': [...], 'allele2': [...]}
            ordered motif lists (as returned by compute_motif_composition). If provided, the parsed motif
            sequence and per-allele motif counts are written for each allele.

    Returns:
        str: Path to the output JSON file
    """
    # Get skip_hom_ref_loci with default of False
    skip_hom_ref_loci = getattr(args, 'skip_hom_ref_loci', False)

    # Filter loci based on skip_hom_ref_loci
    if skip_hom_ref_loci:
        filtered_loci = [locus for locus in genotyped_loci if not is_locus_genotyped_as_reference(locus)]
        if args.verbose:
            print(f"Skipped {len(genotyped_loci) - len(filtered_loci):,d} homozygous reference loci")
    else:
        filtered_loci = genotyped_loci

    # Sort by genomic position and motif, ordering chromosomes naturally so this matches the VCF output
    filtered_loci.sort(key=lambda x: (compute_chrom_sort_key(x.chrom), x.start_0based, x.end, x.motif))

    # Construct output filename
    json_output_path = f"{args.output_prefix}.tandem_repeat_genotypes.json.gz"

    # Write the records one at a time rather than building the whole list first, since for a catalog with millions
    # of loci that list is the largest object in memory. The text is identical to json.dump(list, f, indent=2):
    # each record is indented one level inside the enclosing list, and an empty list is written as "[]".
    with gzip.open(json_output_path, "wt") as f:
        f.write("[")
        for i, locus in enumerate(filtered_loci):
            json_dict = locus.to_json_dict(
                motif_lists=motif_lists_by_locus.get(locus.locus_id) if motif_lists_by_locus else None,
            )
            f.write(",\n  " if i > 0 else "\n  ")
            f.write(json.dumps(json_dict, indent=2).replace("\n", "\n  "))
        f.write("\n]" if filtered_loci else "]")

    if args.verbose:
        print(f"Wrote {len(filtered_loci):,d} genotyped TR loci to {json_output_path}")

    return json_output_path


def write_genotypes_vcf(genotyped_loci, input_vcf_path, args):
    """Write contributing variants to a VCF file with TR annotation.

    This function outputs a VCF file containing only the variants that overlapped
    TR loci during genotyping. Each variant is annotated with INFO fields
    indicating which TR locus/loci it contributed to.

    A variant that overlaps multiple TR loci is written once with comma-separated
    locus IDs and motifs in the INFO field.

    Args:
        genotyped_loci (list): List of GenotypedTandemRepeat objects
        input_vcf_path (str): Path to the input single-sample VCF file
        args (argparse.Namespace): Command-line arguments. Must have:
            - output_prefix (str): Prefix for output file path
            - verbose (bool): If True, print detailed output
            - skip_hom_ref_loci (bool): If True, skip loci with no overlapping variants.

    Returns:
        str: Path to the output VCF file (gzip-compressed)
    """
    # Get skip_hom_ref_loci with default of False
    skip_hom_ref_loci = getattr(args, 'skip_hom_ref_loci', False)

    # Filter loci based on skip_hom_ref_loci
    if skip_hom_ref_loci:
        filtered_loci = [locus for locus in genotyped_loci if not is_locus_genotyped_as_reference(locus)]
        if args.verbose:
            print(f"Skipped {len(genotyped_loci) - len(filtered_loci):,d} homozygous reference loci for VCF output")
    else:
        filtered_loci = genotyped_loci

    # Build a map of variant (chrom, pos, ref, alt_tuple) -> list of overlapping loci
    # We use (chrom, pos, ref, tuple(alts)) as key to uniquely identify variants
    variant_to_loci = collections.defaultdict(list)

    for locus in filtered_loci:
        for variant in locus.overlapping_variants:
            # Create a key that uniquely identifies this variant
            # variant.alts is a tuple of alt alleles
            variant_key = (variant.chrom, variant.pos, variant.ref, variant.alts)
            variant_to_loci[variant_key].append(locus)

    if not variant_to_loci:
        if args.verbose:
            print("No contributing variants found. Skipping VCF output.")
        return None

    # Open the input VCF to read the header
    input_vcf = pysam.VariantFile(input_vcf_path)

    # Create a new header with additional INFO fields
    new_header = input_vcf.header.copy()

    # Add INFO field definitions for TR annotations
    new_header.add_line(
        '##INFO=<ID=TR_LocusId,Number=.,Type=String,Description="Tandem repeat locus ID(s) that this variant overlaps (format: chr-start-end-motif, 0-based start). Comma-separated if multiple loci.">'
    )
    new_header.add_line(
        '##INFO=<ID=TR_Motif,Number=.,Type=String,Description="Motif(s) of the overlapping tandem repeat locus/loci. Comma-separated if multiple loci.">'
    )

    # Construct output filename (uncompressed initially, then bgzip)
    vcf_output_path = f"{args.output_prefix}.tandem_repeat_contributing_variants.vcf"

    # Open output VCF for writing
    output_vcf = pysam.VariantFile(vcf_output_path, "w", header=new_header)

    # Collect all variants we need to write, sorted by position
    variants_to_write = []
    for variant_key, loci in variant_to_loci.items():
        chrom, pos, ref, alts = variant_key
        variants_to_write.append((chrom, pos, ref, alts, loci))

    # Sort by chromosome and position, using a natural chromosome sort (chr1, chr2, ..., chr10, ...) that also
    # keeps each contig's records contiguous, which tabix requires in order to index the file.
    variants_to_write.sort(key=lambda x: (compute_chrom_sort_key(x[0]), x[1]))

    # Write each variant once
    output_variant_count = 0
    for chrom, pos, ref, alts, loci in variants_to_write:
        # Fetch the original variant record from input VCF to preserve all fields
        # We search for it by position
        try:
            original_records = list(input_vcf.fetch(chrom, pos - 1, pos))
        except ValueError:
            # Region not found (chromosome mismatch or no index)
            original_records = []

        # Find the exact matching record
        original_record = None
        for record in original_records:
            if record.pos == pos and record.ref == ref and record.alts == alts:
                original_record = record
                break

        if original_record is None:
            # If we can't find the original record, skip this variant
            # This shouldn't happen in normal operation
            if args.verbose:
                print(f"WARNING: Could not find original record for variant at {chrom}:{pos} {ref}->{alts}")
            continue

        # Create a new record based on the original
        new_record = output_vcf.new_record()
        new_record.chrom = original_record.chrom
        new_record.pos = original_record.pos
        new_record.id = original_record.id
        new_record.ref = original_record.ref
        new_record.alts = original_record.alts
        new_record.qual = original_record.qual
        # Copy all FILTER values from the original record. A record whose FILTER column is "." has no filter
        # keys, and the VCF spec distinguishes that ("filters not applied") from PASS ("all filters passed"),
        # so leave the new record's filter set empty rather than asserting PASS.
        for filter_key in original_record.filter.keys():
            new_record.filter.add(filter_key)

        # Copy all original INFO fields
        for key in original_record.info.keys():
            new_record.info[key] = original_record.info[key]

        # Add TR annotation INFO fields
        # Build comma-separated locus IDs and motifs
        locus_ids = [loc.locus_id for loc in loci]
        motifs = [loc.motif for loc in loci]

        new_record.info["TR_LocusId"] = ",".join(locus_ids)
        new_record.info["TR_Motif"] = ",".join(motifs)

        # Copy FORMAT and sample data
        for sample in original_record.samples:
            for key in original_record.samples[sample].keys():
                try:
                    new_record.samples[sample][key] = original_record.samples[sample][key]
                except (KeyError, TypeError):
                    pass  # Skip if format field not defined in header
            # Preserve genotype phasing (e.g. 0|1), which is stored separately from the
            # GT value and would otherwise default to unphased (0/1) on the new record
            new_record.samples[sample].phased = original_record.samples[sample].phased

        output_vcf.write(new_record)
        output_variant_count += 1

    # Close files
    output_vcf.close()
    input_vcf.close()

    run_shell_command(f"bgzip -f {shlex.quote(vcf_output_path)}")
    run_shell_command(f"tabix -f -p vcf {shlex.quote(vcf_output_path)}.gz")

    if args.verbose:
        print(f"Wrote {output_variant_count:,d} contributing variants to {vcf_output_path}.gz")
        if output_variant_count != len(variant_to_loci):
            print(f"  Note: {len(variant_to_loci) - output_variant_count} variants could not be written (not found in input VCF)")

    return f"{vcf_output_path}.gz"


def do_genotype_subcommand(args):
    """Genotype tandem repeat loci by looking at the genotypes of indels in the input single-sample VCF file.

    This function implements the main genotype subcommand which:
    1. Opens the reference fasta
    2. Parses the catalog BED file to get TR loci
    3. Genotypes each locus using variants from the input VCF
    4. Writes output files (TSV and optionally VCF)

    Args:
        args: Argument namespace with:
            - reference_fasta_path (str): Path to reference genome fasta
            - catalog_bed (str): Path to catalog BED file with TR loci
            - input_vcf_path (str): Path to single-sample VCF with variants
            - output_prefix (str): Output file prefix (optional)
            - interval (list): Optional list of genomic intervals to filter to
            - verbose (bool): If True, print detailed logs
            - show_progress_bar (bool): If True, display progress bar
            - write_vcf (bool): If True, output VCF with contributing variants
    """
    # Set up output prefix
    if args.output_prefix is None:
        args.output_prefix = args.input_vcf_prefix

    # Open reference fasta
    if not file_exists(args.reference_fasta_path):
        raise ValueError(f"Reference fasta not found: {args.reference_fasta_path}")
    fasta_obj = pyfaidx.Fasta(args.reference_fasta_path, one_based_attributes=False, as_raw=True)

    if args.verbose:
        print(f"Reference fasta: {args.reference_fasta_path}")

    # Parse catalog BED file
    if not file_exists(args.catalog_bed):
        raise ValueError(f"Catalog BED file not found: {args.catalog_bed}")

    intervals = args.interval if hasattr(args, 'interval') and args.interval else None
    catalog_loci = parse_catalog_bed_file(
        args.catalog_bed,
        intervals=intervals,
        verbose=args.verbose,
    )

    if len(catalog_loci) == 0:
        print("WARNING: No TR loci found in catalog. Nothing to genotype.")
        return

    # Genotype all loci
    genotyped_loci, counters, motif_lists_by_locus = genotype_all_loci(
        catalog_loci,
        args.input_vcf_path,
        fasta_obj,
        args,
    )

    # Parse each allele sequence into an ordered list of motifs if requested, unless genotype_all_loci already did
    # so in its worker processes. This is computed once and shared by both the TSV and JSON outputs.
    if motif_lists_by_locus is None:
        motif_lists_by_locus = compute_motif_composition(genotyped_loci, args)

    # Write TSV output
    write_genotypes_tsv(genotyped_loci, args, motif_lists_by_locus)

    # Write JSON output if requested
    if args.write_json:
        write_genotypes_json(genotyped_loci, args, motif_lists_by_locus)

    # Write VCF output if requested
    if args.write_vcf:
        write_genotypes_vcf(genotyped_loci, args.input_vcf_path, args)

    # Print summary
    if args.verbose:
        print(f"\nGenotype subcommand complete.")
        print(f"  Output prefix: {args.output_prefix}")


def main():
    """Main function to parse arguments and run the tandem repeat detection pipeline."""

    args = parse_args()

    if args.subcommand == "catalog":
        do_catalog_subcommand(args)

    elif args.subcommand == "merge":
        do_merge_subcommand(args)

    elif args.subcommand == "genotype":
        do_genotype_subcommand(args)

if __name__ == "__main__":
    main()


