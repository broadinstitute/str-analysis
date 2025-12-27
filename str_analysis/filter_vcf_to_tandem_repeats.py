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
the results to output files. Unlike the original filter_vcf_to_STR_variants.py script, it now separates tandem repeat
locus discovery from genotyping (with genotyping now an optional downstream step than can be performed using the
filter_vcf_to_genotype_tandem_repeats.py script).

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

  3. merge tandem repeat alleles that overlap and have similar motifs
  4. write results to output files
"""

import argparse
import collections
import datetime
import gzip
import itertools
import multiprocessing
import os
import pyfaidx
import pysam
import re
import shutil
import tqdm

from concurrent.futures import ThreadPoolExecutor
from pprint import pformat

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.find_repeat_unit import find_repeat_unit_allowing_interruptions
from str_analysis.utils.find_repeat_unit import find_repeat_unit_without_allowing_interruptions
from str_analysis.utils.find_repeat_unit import extend_repeat_into_sequence_allowing_interruptions
from str_analysis.utils.find_repeat_unit import extend_repeat_into_sequence_without_allowing_interruptions
from str_analysis.utils.find_repeat_unit import extend_repeat_into_sequence_base_by_base
from str_analysis.utils.file_utils import open_file, file_exists
from str_analysis.utils.trf_runner import TRFRunner
from str_analysis.utils.find_motif_utils import compute_repeat_purity, compute_most_common_motif

DETECTION_MODE_PURE_REPEATS = "pure"
DETECTION_MODE_ALLOW_INTERRUPTIONS = "interrupted"
DETECTION_MODE_TRF = "trf"

DETECTION_MODE_ORDER = [
    DETECTION_MODE_PURE_REPEATS,
    DETECTION_MODE_ALLOW_INTERRUPTIONS,
    DETECTION_MODE_TRF,
]

CURRENT_TIMESTAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S.%f")
TRF_WORKING_DIR = f"trf_working_dir"


MAX_INDEL_SIZE = 100_000  # bp
MAX_FLANKING_SEQUENCE_SIZE = 1_000_000  # bp


FILTER_ALLELE_WITH_N_BASES = "contains Ns in the variant sequence"
FILTER_ALLELE_SNV_OR_MNV = "SNV/MNV"
FILTER_ALLELE_MNV_INDEL = "complex multinucleotide indel"
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

#FILTER_TR_ALLELE_PARTIAL_REPEAT = "ends in partial repeat"

FILTER_TR_ALLELE_REPEAT_UNIT_TOO_SHORT = "repeat unit < %d bp"
FILTER_TR_ALLELE_REPEAT_UNIT_TOO_LONG = "repeat unit > %d bp"


def parse_args():
    """Parse command-line arguments."""

    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparsers = p.add_subparsers(dest="subcommand", required=True, title="Subcommand")
    catalog_p = subparsers.add_parser("catalog", help="step1: Discover tandem repeat loci in a VCF file.")

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

    catalog_p.add_argument("--trf-min-repeats-in-reference", default=2, help="For TRF results, require a locus to span "
        "at least this many repeats in the reference genome. This helps filter out noisy, low-quality TRF results.")
    catalog_p.add_argument("--trf-min-purity", default=0.2, help="For TRF results, filter out locus definitions where the "
        "repeat purity is below this threshold (defined as the fraction of bases that correspond to perfect repeats "
        "of the locus motif).")
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
    catalog_p.add_argument("--write-filtered-out-variants-to-vcf", help="Output a VCF file with variants where one allele was found to be an TR, "
                   "but that were still filtered out for reasons such as being multiallelic and having alleles with different motifs, "
                   "or because one allele was an TR while the other was an SNV. These types of variants are filtered out to reduce complexity "
                   "in downstream analyses.",
                   action="store_true")
    catalog_p.add_argument("--write-fasta", help="Output a FASTA file containing all TR alleles", action="store_true")
    catalog_p.add_argument("--write-tsv", help="Output a TSV file containing all TR alleles", action="store_true")
    catalog_p.add_argument("-ik", "--copy-info-field-keys-to-tsv", help="Copy the values of these INFO field keys from the input "
                                                                 "VCF to the output TSV files.", action="append")
    catalog_p.add_argument("-L", "--interval", help="Only process variants in this genomic interval (format: chrom:start-end) or BED file",
                   action="append")
    
    catalog_p.add_argument("input_vcf_path", help="Input single-sample VCF file path. This script was designed and tested on VCFs produced by DipCall, "
                   "but should work with any single-sample VCF.")

    # merge subcommand
    merge_p = subparsers.add_parser("merge", help="step2: Optionally merge catalogs produced in step1 from two or more VCF files.")
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
    genotype_p = subparsers.add_parser("genotype", help="step3: Genotype tandem repeat loci by looking at the genotypes of indels in the input single-sample VCF file.")
    genotype_p.add_argument("--catalog-bed", help="Input BED file containing a tandem repeat catalog where the name field contains "
                            "(or starts with) the locus repeat unit. This catalog can be generated by the 'catalog' subcommand, or can be from other sources.", required=True)
    genotype_p.add_argument("-o", "--output-prefix", help="Output file prefix. If not specified, it will be computed based on "
                            "the input vcf filename")
    genotype_p.add_argument("--write-vcf", help="Output a VCF file with the subset of variants that contributed to TR genotyping.", action="store_true")
    genotype_p.add_argument("--write-motif-composition", help="Output a JSON file containing the distribution of observed motifs at each TR locus.", action="store_true")
    genotype_p.add_argument("input_vcf_path", help="Input VCF single-sample VCF file containing variant genotypes from which to compute the TR genotypes.")

    args = p.parse_args()

    if args.subcommand == "catalog":
        if args.dont_allow_interruptions:
            args.dont_run_trf = True

        if not args.dont_run_trf and not args.trf_executable_path:
            p.error(f"Must specify --trf-executable-path or --dont-run-trf")

        if args.copy_info_field_keys_to_tsv:
            args.copy_info_field_keys_to_tsv = {key: 0 for key in args.copy_info_field_keys_to_tsv}
        
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
        elif left_or_right == "right":
            chrom_size = len(self._fasta_obj[self._chrom])
            self._right_flank_end = min(self._right_flank_start_0based + num_flanking_bases, chrom_size)
            self._right_flanking_reference_sequence = str(self._fasta_obj[self._chrom][self._right_flank_start_0based : self._right_flank_end]).upper()
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
            adjust_repeat_unit_and_boundaries_to_maximize_purity,
            num_repeat_bases_in_left_flank,
            num_repeat_bases_in_variant, 
            num_repeat_bases_in_right_flank, 
            detection_mode,
    ):
        """Initialize a TandemRepeatAllele object.

        Args:
            allele (Allele): the allele record that this TandemRepeatAllele object is based on
            repeat_unit (str): the repeat unit of the tandem repeat allele
            adjust_repeat_unit_and_boundaries_to_maximize_purity (bool): whether to set the repeat unit to the most common motif in the variant sequence of the same length as the given repeat unit
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

        if adjust_repeat_unit_and_boundaries_to_maximize_purity:
            self._adjust_repeat_unit_and_boundaries_to_maximize_purity()

        self._repeat_unit_length = len(self._repeat_unit)
        self._detection_mode = detection_mode
        self._summary_string = None

        self._canonical_repeat_unit = None
        self._repeat_purity = None

        self._start_0based = self._allele.get_left_flank_end() - self._num_repeat_bases_in_left_flank
        self._end_1based = self._allele.get_right_flank_start_0based() + self._num_repeat_bases_in_right_flank

        if self._start_0based > self._end_1based:
            raise ValueError(f"Logic error: start_0based ({self._start_0based}) > end_1based ({self._end_1based})")

    def _adjust_repeat_unit_and_boundaries_to_maximize_purity(self):
        if self._num_repeat_bases_in_left_flank + self._num_repeat_bases_in_variant + self._num_repeat_bases_in_right_flank < len(self._repeat_unit):
            return

        most_common_motif = compute_most_common_motif(self.variant_and_flanks_repeat_sequence, len(self._repeat_unit))
        self.repeat_unit_adjusted = most_common_motif != self._repeat_unit  # for debugging
        self._repeat_unit = most_common_motif

        left_flanking_sequence = self._allele.get_left_flanking_sequence()
        extra_bases_in_left_flank = extend_repeat_into_sequence_base_by_base(
            self._repeat_unit[::-1], left_flanking_sequence[::-1][self._num_repeat_bases_in_left_flank:])
        self._num_repeat_bases_in_left_flank += extra_bases_in_left_flank

        right_flanking_sequence = self._allele.get_right_flanking_sequence()
        extra_bases_in_right_flank = extend_repeat_into_sequence_base_by_base(
            self._repeat_unit, right_flanking_sequence[self._num_repeat_bases_in_right_flank:])
        self._num_repeat_bases_in_right_flank += extra_bases_in_right_flank

        self.added_extra_bases_to_left_flank = bool(extra_bases_in_left_flank)
        self.added_extra_bases_to_right_flank = bool(extra_bases_in_right_flank)
        if self.added_extra_bases_to_left_flank or self.added_extra_bases_to_right_flank:
            most_common_motif = compute_most_common_motif(self.variant_and_flanks_repeat_sequence, len(self._repeat_unit))
            self.repeat_unit_adjusted = self.repeat_unit_adjusted or most_common_motif != self._repeat_unit
            self._repeat_unit = most_common_motif

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
            self._summary_string += f"{self.num_repeats_in_variant_and_flanks/self.repeat_unit_length:0.1f}x:"
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
    """Represents a tandem repeat locus in the reference genome"""

    def __init__(
            self,
            chrom,
            start_0based,
            end_1based,
            repeat_unit,
            detection_mode=None,
        ):
        """Initialize a TandemRepeatAllele object.

        Args:
            chrom (str): repeat locus chromosome name
            start_0based (int): repeat locus start coordinate
            end_1based (int): repeat locus end coordinate
            repeat_unit (str): the repeat unit
            detection_mode (str): the detection mode used to find the tandem repeat allele
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
        self._repeat_sequence = None
        self._summary_string = None

    @property
    def chrom(self):
        return self._chrom
    
    @property
    def start_0based(self):
        return self._start_0based
    
    @property
    def end_1based(self):
        return self._end_1based

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
    def ref_interval_size(self):
        return self._end_1based - self._start_0based

    @property
    def locus_id(self):
        return f"{self._chrom}-{self._start_0based}-{self._end_1based}-{self.repeat_unit}"

    @property
    def num_repeats_ref(self):
        return self.ref_interval_size // self.repeat_unit_length

    @property
    def summary_string(self):
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


def do_catalog_subcommand(args):
    """Main function to parse arguments and run the tandem repeat detection pipeline."""

    fasta_obj = pyfaidx.Fasta(args.reference_fasta_path, one_based_attributes=False, as_raw=True)

    # parse input VCF
    counters = collections.defaultdict(int)
    alleles_from_vcf = parse_input_vcf_file(args, counters, fasta_obj)

    # detect tandem repeats
    alleles_that_are_tandem_repeats, alleles_to_process_using_trf = detect_perfect_and_almost_perfect_tandem_repeats(
        alleles_from_vcf, counters, args)
    
    if not args.dont_run_trf:
        more_alleles_that_are_tandem_repeats = detect_tandem_repeats_using_trf(
            alleles_to_process_using_trf, counters, args)
        alleles_that_are_tandem_repeats.extend(more_alleles_that_are_tandem_repeats)

    # write results to output file(s)
    if not args.output_prefix:
        args.output_prefix = args.input_vcf_prefix

    write_bed(alleles_that_are_tandem_repeats, args)

    if args.write_detailed_bed:
        write_bed(alleles_that_are_tandem_repeats, args, detailed=True)

    if args.write_vcf:
        write_vcf(alleles_that_are_tandem_repeats, args, only_write_filtered_out_alleles=False)
    
    if args.write_filtered_out_variants_to_vcf:
        write_vcf(alleles_that_are_tandem_repeats, args, only_write_filtered_out_alleles=True)

    if args.write_tsv:
        write_tsv(alleles_that_are_tandem_repeats, args)

    if args.write_fasta:
        write_fasta(alleles_that_are_tandem_repeats, args)

    if args.verbose:
        print_stats(counters)
        print_tr_stats(alleles_that_are_tandem_repeats)


def detect_perfect_and_almost_perfect_tandem_repeats(alleles, counters, args):

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
            if tandem_repeat_allele.do_repeats_cover_entire_flanking_sequence():
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


def detect_tandem_repeats_using_trf(alleles, counters, args):
    """Runs TandemRepeatFinder (TRF) on a list of indel alleles to detect tandem repeats."""

    # set up TRF working dir
    original_working_dir = os.getcwd()

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
        os.chdir(trf_working_dir)
        
        n_threads = min(args.trf_threads, len(alleles_to_process_next))
        if args.verbose:
            print(f"Launching {n_threads} TRF instance(s) to check {len(alleles_to_process_next):,d} indel alleles for tandem repeats", 
                    "after extending their flanking sequences" if not first_iteration else "")
        first_iteration = False

        alleles_to_reprocess = []
        with ThreadPoolExecutor(max_workers=n_threads) as ex:
            
            futures = []
            for thread_i in range(0, n_threads):
                thread_input_alleles = [allele for allele_i, allele in enumerate(alleles_to_process_next) if allele_i % n_threads == thread_i]
                futures.append(ex.submit(run_trf, thread_input_alleles, args, thread_i))
                
            # collect and process results from all threads
            for thread_i in range(0, n_threads):
                for tandem_repeat_allele, filter_reason, allele in futures[thread_i].result():
                    if filter_reason:
                        counters[f"allele filter: TRF: {filter_reason}"] += 1
                        #if args.debug:
                        #    print(f"TRF filtered out: {allele}, filter reason: {filter_reason}")
                        continue
                    
                    # reprocess the allele if the repeats were found to cover the entire left or right flanking sequence
                    if need_to_reprocess_allele_with_extended_flanking_sequence(tandem_repeat_allele):
                        counters[f"allele op: increased flanking sequence size {tandem_repeat_allele.allele.number_of_times_flanking_sequence_size_was_increased}x for TRF"] += 1
                        alleles_to_reprocess.append(allele)
                        continue

                    # this allele was found to be a tandem repeat using TRF
                    if tandem_repeat_allele.do_repeats_cover_entire_flanking_sequence():
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
        os.chdir(original_working_dir)
        if not args.debug:
            shutil.rmtree(trf_working_dir)

    return tandem_repeat_alleles


def parse_input_vcf_file(args, counters, fasta_obj):
    """Parse the input VCF file and return a list of Allele objects."""

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
            continue

        if not vcf_alt:
            raise ValueError(f"No ALT allele found in VCF row #{vcf_line_i + 1:,d}: {vcf_fields}")

        if len(alt_alleles) > 2:
            if args.verbose:
                vcf_chrom_without_chr_prefix = vcf_chrom.replace("chr", "")
                print(f"WARNING: VCF row #{vcf_line_i:,d}: {vcf_chrom_without_chr_prefix}-{vcf_pos}-{vcf_ref}-{vcf_alt}: multi-allelic variant has "
                    f"{len(alt_alleles)} alt alleles. This script doesn't support more than 2 alt alleles. Skipping...")
            counters[f"WARNING: multi-allelic variant with {len(alt_alleles)} alt alleles"] += 1
            continue

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
                #allele_filter_reason = FILTER_ALLELE_SNV_OR_MNV
                continue
            elif len(vcf_ref) < len(alt_allele) and not alt_allele.startswith(vcf_ref):
                counters[f"allele filter: complex MNV deletion/insertion"] += 1
                #allele_filter_reason = FILTER_ALLELE_MNV_INDEL
                continue
            elif len(alt_allele) < len(vcf_ref) and not vcf_ref.startswith(alt_allele):
                counters[f"allele filter: complex MNV insertion/deletion"] += 1
                #allele_filter_reason = FILTER_ALLELE_MNV_INDEL
                continue

            allele = Allele(
                vcf_chrom, vcf_pos, vcf_ref, alt_allele, fasta_obj,
                order=allele_order,
                info_field_dict=info_field_dict)

            allele_order += 1

            if len(allele.variant_bases) > MAX_INDEL_SIZE:
                # this is a very large indel, so we don't want to process it
                counters[f"allele filter: {FILTER_ALLELE_TOO_BIG}"] += 1
                #allele_filter_reason = FILTER_ALLELE_TOO_BIG
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
        adjust_repeat_unit_and_boundaries_to_maximize_purity=True,
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


def run_trf(alleles, args, thread_id=0):
    """Run TRF on the given allele records.

    Args:
        alleles (list): List of Allele objects to run TRF on.
        args (argparse.Namespace): Command-line arguments parsed by parse_args().
        thread_id (int): ID of thread executing this function, starting from 0.

    Returns:
        list of 3-tuples: (tandem_repeat_allele, filter_reason, allele)
            tandem_repeat_allele (TandemRepeatAllele): tandem repeat allele object if the allele is a tandem repeat, None otherwise
            filter_reason (str): string describing the reason why the allele failed filters, or None if the allele passed all filters
            allele (Allele): Allele object that was processed
    """


    trf_fasta_filename = f"trf_input_sequences__thread{thread_id}.fa"
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
        mismatch_penalty = 7,
        indel_penalty = 7,
        minscore = 20,
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
                adjust_repeat_unit_and_boundaries_to_maximize_purity=True,
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
                for motif_size, matching_trf_results in motif_size_to_passing_trf_results.items():
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

    if args.verbose:
        print(f"thread {thread_id}: TRF allele filter reasons:")
        for key, count in trf_allele_filter_counters.items():
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
        tandem_repeat_alleles (list): list of TandemRepeatAllele objects
        pyfaidx_fasta_obj (pyfaidx_fasta.Fasta): reference fasta object
        verbose (bool): if True, print verbose output

    Returns:
        list: list of TandemRepeatAllele objects
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
            for current_repeat_unit_id, tr_allele_group in motif_id_to_tr_allele_group.items():
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
                new_repeat_unit = compute_most_common_motif(repeat_sequence, tr_alleles_group[0].repeat_unit_length)
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
        if tandem_repeat_allele.do_repeats_cover_entire_left_flanking_sequence():
            tandem_repeat_allele.allele.increase_left_flanking_sequence_size()
        if tandem_repeat_allele.do_repeats_cover_entire_right_flanking_sequence():
            tandem_repeat_allele.allele.increase_right_flanking_sequence_size()
        
        if tandem_repeat_allele.allele.get_expected_left_flanking_sequence_size() <= MAX_FLANKING_SEQUENCE_SIZE \
            and tandem_repeat_allele.allele.get_expected_right_flanking_sequence_size() <= MAX_FLANKING_SEQUENCE_SIZE:
            return True
    
    return False


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
                name_field += f":p{tandem_repeat_allele.repeat_purity:0.2}"
            else:
                name_field = tandem_repeat_allele.repeat_unit

            f.write("\t".join(map(str, [
                tandem_repeat_allele.chrom,
                tandem_repeat_allele.start_0based,
                tandem_repeat_allele.end_1based,
                name_field,
                tandem_repeat_allele.repeat_unit_length,
            ])) + "\n")

    os.system(f"bgzip -f {bed_output_path}")
    os.system(f"tabix -p bed {bed_output_path}.gz")

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
        for key, count in args.copy_info_field_keys_to_tsv.items():
            if count > 0:
                extra_header_columns.append(key)
            else:
                print(f"WARNING: INFO field key '{key}' was not found in any of rows of the input VCF. Skipping..")

    with open(tsv_output_path, "w") as f:
        f.write("\t".join(header + extra_header_columns) + "\n")
        
        for tandem_repeat_allele in tandem_repeat_alleles:
            output_row = [
                tandem_repeat_allele.chrom,
                tandem_repeat_allele.start_0based,
                tandem_repeat_allele.end_1based,
                f"{tandem_repeat_allele.chrom}:{tandem_repeat_allele.start_0based}-{tandem_repeat_allele.end_1based}",
                tandem_repeat_allele.locus_id,
                tandem_repeat_allele.ins_or_del,
                tandem_repeat_allele.repeat_unit,
                tandem_repeat_allele.canonical_repeat_unit,
                tandem_repeat_allele.repeat_unit_length,
                tandem_repeat_allele.num_repeats_ref,
                tandem_repeat_allele.allele.pos,
                tandem_repeat_allele.summary_string,
                tandem_repeat_allele.end_1based > tandem_repeat_allele.start_0based,
                tandem_repeat_allele.is_pure_repeat,
                tandem_repeat_allele.detection_mode,
            ]
            for key in extra_header_columns:
                output_row.append(tandem_repeat_allele.info_field_dict.get(key, ""))

            f.write("\t".join(map(str, output_row)) + "\n")

    os.system(f"bgzip -f {tsv_output_path}")
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

    os.system(f"gzip -f {fasta_output_path}")
    if args.verbose:
        print(f"Wrote {len(tandem_repeat_alleles):,d} tandem repeat sequences to {fasta_output_path}.gz")


def get_input_vcf_iterator(args, include_header=False):
    if args.interval:
        if args.verbose:
            print(f"Parsing interval(s) {', '.join(args.interval)} from {args.input_vcf_path}")

        vcf_iterator = []
        tabix_file = pysam.TabixFile(args.input_vcf_path)
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

        vcf_iterator = itertools.chain(
            vcf_iterator,
            (line for interval in intervals for line in tabix_file.fetch(interval)),
        )
    else:
        if args.verbose:
            print(f"Parsing {args.input_vcf_path}")
        vcf_iterator = open_file(args.input_vcf_path, is_text_file=True)

    return vcf_iterator


def write_vcf(tandem_repeat_alleles, args, only_write_filtered_out_alleles=False):
    """Write variants that either are or aren't tandem repeats to a VCF file.
    
    Args:
        tandem_repeat_alleles (list): list of TandemRepeatAllele objects
        args (argparse.Namespace): command-line arguments parsed by parse_args()
        only_write_filtered_out_alleles (bool): if True, only write the variants that are not in the tandem_repeats_alleles list
    """

    vcf_iterator = get_input_vcf_iterator(args, include_header=True)

    # iterate over all VCF rows
    if only_write_filtered_out_alleles:
        output_vcf_path = f"{args.output_prefix}.not_tandem_repeats.vcf"
    else:
        output_vcf_path = f"{args.output_prefix}.tandem_repeats.vcf"

    tandem_repeat_alleles.sort(key=lambda x: x.order) # sort into their original order

    vcf_alleles = {
        (tr.chrom, tr.allele.pos, tr.allele.ref): tr for tr in tandem_repeat_alleles
    }

    with open(output_vcf_path, "w") as f:
        vcf_line_i = 0
        output_line_counter = 0
        for line in vcf_iterator:
            if line.startswith("#"):
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


            key = (vcf_chrom, vcf_pos, vcf_ref)
            is_tandem_repeat = key in vcf_alleles
            if is_tandem_repeat:
                # append TR info to the INFO field
                tr = vcf_alleles[key]
                if vcf_fields[7] == ".":
                    vcf_fields[7] = ""
                else:
                    vcf_fields[7] += ";"
                vcf_fields[7] += f"MOTIF={tr.repeat_unit}"
                vcf_fields[7] += f";MOTIF_SIZE={tr.repeat_unit_length}"
                vcf_fields[7] += f";START_0BASED={tr.start_0based}"
                vcf_fields[7] += f";END={tr.end_1based}"
                vcf_fields[7] += f";DETECTED={tr.detection_mode}"

                line = "\t".join(vcf_fields) + "\n"

            if (not only_write_filtered_out_alleles and is_tandem_repeat) or (only_write_filtered_out_alleles and not is_tandem_repeat):
                f.write(line)
                output_line_counter += 1

    os.system(f"bgzip -f {output_vcf_path}")
    os.system(f"tabix -p vcf {output_vcf_path}.gz")

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
    input_catalogs_have_details = False

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
            bed_iterator = (line for interval in args.interval for line in tabix_file.fetch(interval))
            input_files_to_close.append(tabix_file)
        else:
            if args.verbose:
                print(f"Parsing catalog #{path_i + 1}: {input_bed_path}")
            bed_iterator = open_file(input_bed_path, is_text_file=True)
            input_files_to_close.append(bed_iterator)

        # parse the BED file into a list of ReferenceTandemRepeat objects
        current_catalog_trs = []
        for line in bed_iterator:
            fields = line.strip().split("\t")
            if len(fields) < 4:
                raise ValueError(f"Invalid BED file format in {input_bed_path}: {line}")
            
            name_field_tokens = fields[3].split(":")
            repeat_unit = name_field_tokens[0]
            detection_mode = None
            if args.write_detailed_bed and len(name_field_tokens) >= 3:
                #motif_size = int(name_field_tokens[1])
                detection_mode = name_field_tokens[3]
                input_catalogs_have_details = True

            # check if the repeat unit itself consists of perfect repeats of a smaller repeat unit (this happends in ~3% of TRs detected by TRF)
            simplified_repeat_unit, _, _ = find_repeat_unit_without_allowing_interruptions(repeat_unit, allow_partial_repeats=False)
            if len(simplified_repeat_unit) != len(repeat_unit):
                simplified_repeat_units_counter += 1

            repeat_unit = simplified_repeat_unit

            current_catalog_trs.append(ReferenceTandemRepeat(
                chrom=fields[0],
                start_0based=int(fields[1]),
                end_1based=int(fields[2]),
                repeat_unit=repeat_unit,
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

    if args.write_detailed_bed and input_catalogs_have_details:
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


def do_genotype_subcommand(args):
    """Genotype tandem repeat loci by looking at the genotypes of indels in the input single-sample VCF file."""

    args.catalog_bed
    args.input_vcf_path

    if args.output_prefix is None:
        args.output_prefix = args.input_vcf_prefix

    
    raise NotImplementedError("Not implemented yet")


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


