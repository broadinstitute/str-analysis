#!/usr/bin/env python3

"""
This script takes a VCF (either single-sample or multi-sample) and filters it to the subset of insertions and deletions
that represent tandem repeat (TR) expansions or contractions. It does this by checking each indel to see if the inserted or deleted sequence
consists entirely of repeats of some motif, and if yes, whether these repeats can be extended into the flanking reference sequences immediately 
to the left or right of the variant. The output is a set of tandem repeat loci (including their motifs and reference start and end coordinates) 
that can then be used for downstream analyses, such as genotyping.

This script is the next iteration of the filter_vcf_to_STR_variants.py script. It implements multiple approaches to detecting repeat sequences
within each variant - first doing a fast, brute-force scan for perfect (or nearly perfect) repeats. If not repeats are detected in this first 
step, it runs TandemRepeatFinder to discover more imperfect repeats (particularly VNTRs). The script then merges overlapping tandem repeat alleles 
that have very similar motifs and writes the results to output files. Unlike the original filter_vcf_to_STR_variants.py script, it now 
separates tandem repeat locus discovery from genotyping (with genotyping now an optional downstream step than can be performed using the 
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

  3. merge tandem repeat alleles that overlap and have very similar motifs
  4. write results to output files
"""

import argparse
import collections
import datetime
import itertools
import math
import multiprocessing
import os
import pyfaidx
import pysam
import re
import shutil
import tqdm

from concurrent.futures import ThreadPoolExecutor
from pprint import pformat, pprint

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.find_repeat_unit import find_repeat_unit_allowing_interruptions
from str_analysis.utils.find_repeat_unit import find_repeat_unit_without_allowing_interruptions
from str_analysis.utils.find_repeat_unit import extend_repeat_into_sequence_allowing_interruptions
from str_analysis.utils.find_repeat_unit import extend_repeat_into_sequence_without_allowing_interruptions
from str_analysis.utils.file_utils import open_file
from str_analysis.utils.trf_runner import TRFRunner


DETECTION_MODE_PURE_REPEATS = "pure"
DETECTION_MODE_ALLOW_INTERRUPTIONS = "interrupted"
DETECTION_MODE_TRF = "trf"

DETECTION_MODE_ORDER = [
    DETECTION_MODE_PURE_REPEATS,
    DETECTION_MODE_ALLOW_INTERRUPTIONS,
    DETECTION_MODE_TRF,
]

CURRENT_TIMESTAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S.%f")
TRF_WORKING_DIR = f"trf_working_dir/{CURRENT_TIMESTAMP}"


MAX_INDEL_SIZE = 100_000  # bp
MAX_FLANKING_SEQUENCE_SIZE = 1_000_000  # bp


FILTER_ALLELE_WITH_N_BASES = "contains Ns in the variant sequence"
FILTER_ALLELE_SNV_OR_MNV = "SNV/MNV"
FILTER_ALLELE_MNV_INDEL = "complex multinucleotide indel"
FILTER_ALLELE_INDEL_WITHOUT_REPEATS = "INDEL without repeats"
FILTER_ALLELE_TOO_BIG = f"INDEL > {MAX_INDEL_SIZE}bp"
FILTER_TR_ALLELE_NOT_ENOUGH_REPEATS = "contains only %d full repeats"
FILTER_TR_ALLELE_DOESNT_SPAN_ENOUGH_BASE_PAIRS = "spans < %d bp"

#FILTER_TR_ALLELE_PARTIAL_REPEAT = "ends in partial repeat"

FILTER_TR_ALLELE_REPEAT_UNIT_TOO_SHORT = "repeat unit < %d bp"
FILTER_TR_ALLELE_REPEAT_UNIT_TOO_LONG = "repeat unit > %d bp"


def parse_args():
    """Parse command-line arguments."""

    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-R", "--reference-fasta-path", help="Reference genome fasta path.", required=True)
    p.add_argument("--dont-allow-interruptions", action="store_true", help="Only detect perfect repeats. This implicitly "
                   "also enables --dont-run-trf since detection of pure repeats does not require running TandemRepeatFinder (TRF).")
    p.add_argument("--dont-run-trf", action="store_true", help="Don't use TandemRepeatFinder (TRF) to help detect imperfect TRs. Instead, only "
                   "use the simpler algorithm that allows one position within the repeat unit to vary across repeats.")
    p.add_argument("--trf-executable-path", help="Path to the TandemRepeatFinder (TRF) executable. This is "
                   "required unless --dont-run-trf is specified.")
    p.add_argument("-t", "--trf-threads", default=max(1, multiprocessing.cpu_count() - 2), type=int, help="Number of TandemRepeatFinder (TRF) "
                   "instances to run in parallel.")
    p.add_argument("--min-indel-size-to-run-trf", default=7, type=int, help="Only run TandemRepeatFinder (TRF) "
                    "on insertions and deletions that are at least this many base pairs.")

    p.add_argument("--min-tandem-repeat-length", type=int, default=9, help="Only detect tandem repeat variants that are at least this long (in base pairs). "
                   "This threshold will be applied to the total repeat sequence including any repeats in the flanking sequence to the left "
                   "and right of the variant in addition to the inserted or deleted bases")
    p.add_argument("--min-repeats", type=int, default=3, help="Only detect tandem repeat loci that consist of at least this many repeats. "
                   "This threshold will be applied to the total repeat sequence including any repeats in the flanking sequence to "
                   "the left and right of the variant in addition to the inserted or deleted bases")
    p.add_argument("--min-repeat-unit-length", type=int, default=1, help="Minimum repeat unit length in base pairs.")
    p.add_argument("--max-repeat-unit-length", type=int, default=10**9, help="Max repeat unit length in base pairs.")
    p.add_argument("--show-progress-bar", help="Show a progress bar in the terminal when processing variants.",
                   action="store_true")
    p.add_argument("-v", "--verbose", help="Print detailed logs.", action="store_true")
    p.add_argument( "--debug", help="Print any debugging info and don't delete intermediate files.", action="store_true")

    p.add_argument("--offset", default=0, type=int, help="Skip the first N variants in the VCF file. This is useful for testing ")
    p.add_argument("-n", type=int, help="Only process N rows from the VCF (after applying --offset). Useful for testing.")

    p.add_argument("-o", "--output-prefix", help="Output file prefix. If not specified, it will be computed based on "
                   "the input vcf filename")
    
    p.add_argument("--write-detailed-bed", help="Output a BED file with all TR alleles where the name field (ie. column 4) contains additional info besides the repeat unit.", action="store_true")
    p.add_argument("--write-vcf", help="Output a VCF file with all variants that were found to be TRs.", action="store_true")
    p.add_argument("--write-filtered-out-variants-to-vcf", help="Output a VCF file with variants where one allele was found to be an TR, "
                   "but that were still filtered out for reasons such as being multiallelic and having alleles with different motifs, "
                   "or because one allele was an TR while the other was an SNV. These types of variants are filtered out to reduce complexity "
                   "in downstream analyses.",
                   action="store_true")
    p.add_argument("--write-fasta", help="Output a FASTA file containing all TR alleles", action="store_true")
    p.add_argument("--write-tsv", help="Output a TSV file containing all TR alleles", action="store_true")
    p.add_argument("-L", "--interval", help="Only process variants in this genomic interval (format: chrom:start-end)",
                   action="append")
    
    p.add_argument("input_vcf_path", help="Input single-sample VCF file path. This script was designed and tested on VCFs produced by DipCall, "
                   "but should work with any single-sample VCF.")

    args = p.parse_args()

    if not args.dont_run_trf and not args.trf_executable_path:
        p.error(f"Must specify --trf-executable-path or --dont-run-trf")

    args.input_vcf_prefix = re.sub(".vcf(.gz|.bgz)?$", "", os.path.basename(args.input_vcf_path))

    return args


class Allele:
    """Represents a single VCF allele."""

    def __init__(self, chrom, pos, ref, alt, fasta_obj):

        self._chrom = chrom
        self._pos = pos
        self._ref = ref
        self._alt = alt
        self._fasta_obj = fasta_obj

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
        self._previously_increased_flanking_sequence_size = False  # tracks whether either the left or right flanking sequence size was increased from it's starting size
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
        return f"{self._chrom}:{self._pos:,d} {self._ref}>{self._alt} ({self._ins_or_del})"
    
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


class TandemRepeatAllele:
    """Stores additional information about a VCF insertion or deletion allele that represents a tandem repeat expansion or contraction."""

    def __init__(
            self, 
            allele, 
            repeat_unit, 
            num_repeat_bases_in_left_flank, 
            num_repeat_bases_in_variant, 
            num_repeat_bases_in_right_flank, 
            detection_mode,
            is_pure_repeat=True,
            motif_interruption_indices=None,
    ):
        """Initialize a TandemRepeatAllele object.

        Args:
            allele (Allele): the allele record that this TandemRepeatAllele object is based on
            repeat_unit (str): the repeat unit of the tandem repeat allele
            num_repeat_bases_in_left_flank (int): the number of repeat bases in the left flanking sequence
            num_repeat_bases_in_variant (int): the number of repeat bases in the variant
            num_repeat_bases_in_right_flank (int): the number of repeat bases in the right flanking sequence
            detection_mode (str): the detection mode used to find the tandem repeat allele
            is_pure_repeat (bool): True if the allele is a pure repeat (no interruptions), False otherwise
            motif_interruption_indices (list): indices of interruptions in the repeat unit, or None if no interruptions were found
        """

        self._allele = allele
        self._repeat_unit = repeat_unit
        self._canonical_repeat_unit = None
        self._num_repeat_bases_in_left_flank = num_repeat_bases_in_left_flank
        self._num_repeat_bases_in_variant = num_repeat_bases_in_variant
        self._num_repeat_bases_in_right_flank = num_repeat_bases_in_right_flank
        self._detection_mode = detection_mode
        self._is_pure_repeat = is_pure_repeat

        if motif_interruption_indices is not None:
            if not isinstance(motif_interruption_indices, list) or not all(isinstance(i, int) for i in motif_interruption_indices):
                raise ValueError(f"Logic error: motif_interruption_indices must be a list or integers, or None rather than '{motif_interruption_indices}'.")
            self._motif_interruption_indices = motif_interruption_indices 
        else:
            self._motif_interruption_indices = None

        self._summary_string = None

        self._start_0based = self._allele.get_left_flank_end() - num_repeat_bases_in_left_flank
        self._end_1based = self._allele.get_right_flank_start_0based() + num_repeat_bases_in_right_flank

        self._num_repeats_ref = (num_repeat_bases_in_left_flank + num_repeat_bases_in_right_flank) // len(repeat_unit)
        self._num_repeats_alt = self._num_repeats_ref
        if self._allele.ins_or_del == "INS":
            self._num_repeats_alt += num_repeat_bases_in_variant // len(repeat_unit)
        elif self._allele.ins_or_del == "DEL":
            self._num_repeats_ref += num_repeat_bases_in_variant // len(repeat_unit)
        else:
            raise ValueError(f"Logic error: variant {allele} is a complex MNV insertion/deletion")

        self._ref_allele_repeat_sequence = ""
        self._alt_allele_repeat_sequence = ""
        if num_repeat_bases_in_left_flank:
            left_flanking_sequence = self._allele.get_left_flanking_sequence()[-num_repeat_bases_in_left_flank:]
            self._ref_allele_repeat_sequence += left_flanking_sequence
            self._alt_allele_repeat_sequence += left_flanking_sequence

        if num_repeat_bases_in_variant == 0:
            raise ValueError(f"{self._allele} repeat_bases_in_variant specified as 0, implying that this is not a tandem repeat")
    
        if self._allele.ins_or_del == "INS":
            self._alt_allele_repeat_sequence += self._allele.variant_bases
        elif self._allele.ins_or_del == "DEL":
            self._ref_allele_repeat_sequence += self._allele.variant_bases

        if num_repeat_bases_in_right_flank:
            right_flanking_sequence = self._allele.get_right_flanking_sequence()[:num_repeat_bases_in_right_flank]
            self._ref_allele_repeat_sequence += right_flanking_sequence
            self._alt_allele_repeat_sequence += right_flanking_sequence

        self._ref_allele_repeat_sequence = ""
        self._alt_allele_repeat_sequence = ""
        
        if num_repeat_bases_in_left_flank:
            left_flanking_sequence = self._allele.get_left_flanking_sequence()[-num_repeat_bases_in_left_flank:]
            self._ref_allele_repeat_sequence += left_flanking_sequence
            self._alt_allele_repeat_sequence += left_flanking_sequence

        if num_repeat_bases_in_variant == 0:
            raise ValueError(f"{self._allele} repeat_bases_in_variant specified as 0, implying that this is not a tandem repeat")
    
        if self._allele.ins_or_del == "INS":
            self._alt_allele_repeat_sequence += self._allele.variant_bases
        elif self._allele.ins_or_del == "DEL":
            self._ref_allele_repeat_sequence += self._allele.variant_bases

        if num_repeat_bases_in_right_flank:
            right_flanking_sequence = self._allele.get_right_flanking_sequence()[:num_repeat_bases_in_right_flank]
            self._ref_allele_repeat_sequence += right_flanking_sequence
            self._alt_allele_repeat_sequence += right_flanking_sequence

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
        return self._num_repeats_ref

    @property
    def num_repeats_alt(self):
        return self._num_repeats_alt

    @property
    def ref_allele_repeat_sequence(self):
        return self._ref_allele_repeat_sequence

    @property
    def alt_allele_repeat_sequence(self):
        return self._alt_allele_repeat_sequence

    @property
    def variant_and_flanks_repeat_sequence(self):
        return self._alt_allele_repeat_sequence if self.ins_or_del == "INS" else self._ref_allele_repeat_sequence
    
    @property
    def allele(self):
        return self._allele

    @property
    def repeat_unit(self):
        return self._repeat_unit

    @property
    def canonical_repeat_unit(self):
        if self._canonical_repeat_unit is None:
            self._canonical_repeat_unit = compute_canonical_motif(self._repeat_unit, include_reverse_complement=True)
        return self._canonical_repeat_unit

    @property
    def detection_mode(self):
        return self._detection_mode

    @property
    def locus_id(self):
        return f"{self._allele.chrom}-{self._start_0based}-{self._end_1based}-{self._repeat_unit}"

    @property
    def summary_string(self):
        if self._summary_string is None:                
            self._summary_string = f"{len(self.repeat_unit)}bp:"
            if len(self.repeat_unit) > 30:
                self._summary_string += f"{self.repeat_unit[:30]}...:"
            else:
                self._summary_string += f"{self.repeat_unit}:"
            self._summary_string += f"{self.ins_or_del}:{self.detection_mode}"

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
        return self._num_repeat_bases_in_left_flank // len(self._repeat_unit)

    @property
    def num_repeats_in_variant(self):
        return self._num_repeat_bases_in_variant // len(self._repeat_unit)

    @property
    def num_repeats_in_right_flank(self):
        return self._num_repeat_bases_in_right_flank // len(self._repeat_unit)
    
    @property
    def ins_or_del(self):
        return self._allele.ins_or_del

    @property
    def is_pure_repeat(self):
        return self._is_pure_repeat
    
    @property
    def motif_interruption_indices(self):
        return self._motif_interruption_indices

    @property
    def motif_interruption_indices_string(self):
        return ",".join(map(str, self._motif_interruption_indices)) if self._motif_interruption_indices is not None else ""

    @property
    def variant_id(self):
        return f"{self._allele.chrom}-{self._allele.pos}-{self._allele.ref}-{self._allele.alt}"

    def do_repeats_cover_entire_left_flanking_sequence(self):
        return self._num_repeat_bases_in_left_flank > len(self._allele.get_left_flanking_sequence()) - len(self._repeat_unit)

    def do_repeats_cover_entire_right_flanking_sequence(self):
        return self._num_repeat_bases_in_right_flank > len(self._allele.get_right_flanking_sequence()) - len(self._repeat_unit)

    def do_repeats_cover_entire_flanking_sequence(self):
        return self.do_repeats_cover_entire_left_flanking_sequence() or self.do_repeats_cover_entire_right_flanking_sequence()
    

    def __str__(self):
        return f"{self._allele.chrom}:{self._start_0based:,d}-{self._end_1based:,d} {self._repeat_unit}x{self._num_repeats_ref} [{self._detection_mode}]"
    
    def __repr__(self):
        return self.__str__()


def main():
    """Main function to parse arguments and run the tandem repeat detection pipeline."""

    args = parse_args()

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


    # merge overlapping tandem repeats if they have the same (or similar) repeat units
    before = len(alleles_that_are_tandem_repeats)
    alleles_that_are_tandem_repeats = merge_overlapping_tandem_repeat_loci(alleles_that_are_tandem_repeats, counters)

    alleles_that_are_tandem_repeats.sort(key=lambda x: (x.chrom, x.start_0based, x.end_1based))
    print(f"Merged {before - len(alleles_that_are_tandem_repeats):,d} overlapping tandem repeat alleles with very similar repeat units")

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

    

def detect_perfect_and_almost_perfect_tandem_repeats(alleles, counters, args):

    alleles_to_process_next = [(allele, DETECTION_MODE_PURE_REPEATS) for allele in alleles]
    alleles_to_process_next_using_trf = []
    tandem_repeat_alleles = []
    first_iteration = True
    while alleles_to_process_next:
        alleles_to_reprocess = []
        print(f"Checking {len(alleles_to_process_next):,d} indel alleles for tandem repeats", 
              "after extending their flanking sequences" if not first_iteration else "")
        first_iteration = False

        if args.show_progress_bar:
            alleles_to_process_next = tqdm.tqdm(alleles_to_process_next, unit=" alleles", unit_scale=True)

        for allele, detection_mode in alleles_to_process_next:
            tandem_repeat_allele, filter_reason = check_if_allele_is_tandem_repeat(allele, args, counters, detection_mode)

            if filter_reason:
                if allele.previously_increased_flanking_sequence_size:
                    raise ValueError(f"Logic error: allele {allele} was previously detected to have tandem repeats that extended over the entire "
                                     f"flanking sequence, but after extending their flanking sequences, it was no longer found to be a tandem "
                                     f"repeat using detection mode {detection_mode} due to {filter_reason}")
                
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

            # reprocess the allele if the repeats were found tocover the entire left or right flanking sequence
            if need_to_reprocess_allele_with_extended_flanking_sequence(tandem_repeat_allele):
                counters[f"allele op: increased flanking sequence size {tandem_repeat_allele.allele.number_of_times_flanking_sequence_size_was_increased}x for {detection_mode} repeats"] += 1
                alleles_to_reprocess.append((allele, detection_mode))  # reprocess the allele with the same detection mode
                #print(f"Detection mode [{detection_mode}]: Increasing flanking sequence size to {len(tandem_repeat_allele.allele.get_left_flanking_sequence()):,d}bp for {tandem_repeat_allele}")
                continue
            
            # this allele was found to be a tandem repeat using the current detection mode
            if tandem_repeat_allele.do_repeats_cover_entire_flanking_sequence():
                print(f"WARNING: allele {allele} was found to be a tandem repeat using detection mode {detection_mode}, but the repeats cover the entire flanking sequence even though it is longer than {MAX_FLANKING_SEQUENCE_SIZE:,}bp. Skipping...")
            else:
                tandem_repeat_alleles.append(tandem_repeat_allele)

                if not args.dont_run_trf and len(tandem_repeat_allele.repeat_unit) > 6 and len(allele.variant_bases) >= args.min_indel_size_to_run_trf:
                    # if this is a VNTR with a large motif, run TRF on it to see if it detects wider locus boundaries. The merge step can resolve redudant locus definitions.
                    alleles_to_process_next_using_trf.append(allele)

        alleles_to_process_next = alleles_to_reprocess

    print(f"Found {len(tandem_repeat_alleles):,d} indel alleles that represent perfect (or nearly perfect) tandem repeat expansions or contractions")

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

        trf_working_dir = os.path.join(TRF_WORKING_DIR, f"{args.input_vcf_prefix}." + start_time.strftime("%Y%m%d_%H%M%S.%f"))
        if os.path.isdir(trf_working_dir):
            raise ValueError(f"ERROR: TRF working directory already exists: {trf_working_dir}. Each TRF run should be in a unique directory to avoid filename collisions.")

        if args.debug: 
            print("-"*100)
            print(f"TRF working directory: {trf_working_dir}")
        os.makedirs(trf_working_dir)
        os.chdir(trf_working_dir)
        
        n_threads = min(args.trf_threads, len(alleles_to_process_next))
        print(f"Launching {n_threads} TRF instance(s) to check {len(alleles_to_process_next):,d} indel alleles for tandem repeats", 
                "after extending their flanking sequences" if not first_iteration else "")
        first_iteration = False

        alleles_to_reprocess = []
        with ThreadPoolExecutor(max_workers=n_threads) as ex:
            
            futures = []
            for thread_i in range(0, n_threads):
                thread_input_alleles = [allele for allele_i, allele in enumerate(alleles_to_process_next) if allele_i % n_threads == thread_i]
                futures.append(ex.submit(run_trf, thread_input_alleles, args, counters, thread_i))
                
            # collect and process results from all threads
            for thread_i in range(0, n_threads):
                for tandem_repeat_allele, filter_reason, allele in futures[thread_i].result():
                    if filter_reason:
                        counters[f"allele filter: TRF: {filter_reason}"] += 1
                        if args.debug:
                            print(f"TRF filtered out: {allele}, filter reason: {filter_reason}")
                        continue
                    
                    # reprocess the allele if the repeats were found tocover the entire left or right flanking sequence
                    if need_to_reprocess_allele_with_extended_flanking_sequence(tandem_repeat_allele):
                        counters[f"allele op: increased flanking sequence size {tandem_repeat_allele.allele.number_of_times_flanking_sequence_size_was_increased}x for TRF"] += 1
                        alleles_to_reprocess.append(allele)
                        continue

                    # this allele was found to be a tandem repeat using TRF
                    if tandem_repeat_allele.do_repeats_cover_entire_flanking_sequence():
                        print(f"WARNING: allele {allele} was found to be a tandem repeat using TRF, but the repeats cover the entire flanking sequence even though it is longer than {MAX_FLANKING_SEQUENCE_SIZE:,}bp. Skipping...")
                    else:
                        tandem_repeat_alleles.append(tandem_repeat_allele)

        alleles_to_process_next = alleles_to_reprocess

        elapsed = datetime.datetime.now() - start_time
        print(f"Found {len(tandem_repeat_alleles) - before_counter:,d} additional tandem repeats after running TRF for {elapsed.seconds//60}m {elapsed.seconds%60}s"
              + (", but will recheck " + str(len(alleles_to_process_next)) + " alleles after extending their flanking sequences" if len(alleles_to_process_next) > 0 else ""))
        os.chdir(original_working_dir)
        if not args.debug:
            shutil.rmtree(trf_working_dir)

    return tandem_repeat_alleles

    
def parse_input_vcf_file(args, counters, fasta_obj):
    """Parse the input VCF file and return a list of Allele objects."""

    input_files_to_close = []
    if args.interval:
        print(f"Parsing interval(s) {', '.join(args.interval)} from {args.input_vcf_path}")

        tabix_file = pysam.TabixFile(args.input_vcf_path)
        vcf_iterator = (line for interval in args.interval for line in tabix_file.fetch(interval))
        input_files_to_close.append(tabix_file)
    else:
        print(f"Parsing {args.input_vcf_path}")
        vcf_iterator = open_file(args.input_vcf_path, is_text_file=True)
        input_files_to_close.append(vcf_iterator)

    if args.show_progress_bar:
        vcf_iterator = tqdm.tqdm(vcf_iterator, unit=" variants", unit_scale=True)


    # iterate over all VCF rows
    alleles_from_vcf = []
    vcf_line_i = 0
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

            allele = Allele(vcf_chrom, vcf_pos, vcf_ref, alt_allele, fasta_obj)
            if len(allele.variant_bases) > MAX_INDEL_SIZE:
                # this is a very large indel, so we don't want to process it
                counters[f"allele filter: {FILTER_ALLELE_TOO_BIG}"] += 1
                #allele_filter_reason = FILTER_ALLELE_TOO_BIG
                continue

            counters[f"allele counts: {allele.ins_or_del} alleles"] += 1

            alleles_from_vcf.append(allele)

    print(f"Parsed {len(alleles_from_vcf):,d} indel alleles from {args.input_vcf_path}")
    
    return alleles_from_vcf



def check_if_allele_is_tandem_repeat(allele, args, counters, detection_mode):
    """Determine if the given allele is a tandem repeat expansion or contraction or not. 
    This is done by performing a brute-force scan for perfect (or nearly perfect) repeats in the allele sequence, and then extending the repeats
    into the flanking reference sequences.

    Args:
        allele (Allele): allele record
        args (argparse.Namespace): command-line arguments parsed by parse_args()
        counters (dict): Dictionary of counters to collect summary stats about the number of TR variants found, etc.
        detection_mode (str): Should be either DETECTION_MODE_PURE_REPEATS or DETECTION_MODE_ALLOW_INTERRUPTIONS
    Return:
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

        is_pure_repeat = True
        motif_interruption_indices = None

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

        motif_interruption_indices = [repeat_unit_interruption_index] if repeat_unit_interruption_index is not None else None
        is_pure_repeat = False

        num_repeat_bases_in_left_flank = num_total_repeats_left_flank * len(repeat_unit)
        num_repeat_bases_in_right_flank = num_total_repeats_right_flank * len(repeat_unit)

    else:
        raise ValueError(f"Invalid detection_mode: '{detection_mode}'. It must be either '{DETECTION_MODE_PURE_REPEATS}' or '{DETECTION_MODE_ALLOW_INTERRUPTIONS}'.")


    tandem_repeat_allele = TandemRepeatAllele(
        allele,
        repeat_unit,
        num_repeat_bases_in_left_flank,
        len(allele.variant_bases),
        num_repeat_bases_in_right_flank,
        detection_mode,
        is_pure_repeat,
        motif_interruption_indices)
    
    tandem_repeat_allele_failed_filters_reason = check_if_tandem_repeat_allele_failed_filters(args, tandem_repeat_allele)
    
    if args.debug: print(f"{detection_mode} repeats: {tandem_repeat_allele}, filter: {tandem_repeat_allele_failed_filters_reason}")
    if tandem_repeat_allele_failed_filters_reason is not None:
        return None, tandem_repeat_allele_failed_filters_reason
    else:
        return tandem_repeat_allele, None



def check_if_tandem_repeat_allele_failed_filters(args, tandem_repeat_allele):
    """Check if the given tandem repeat allele (represented by its repeat_unit and total_repeats) passes the filters
    specified in the command-line arguments.

    Args:
        args (argparse.Namespace): command-line arguments parsed by parse_args()
        tandem_repeat_allele (TandemRepeatAllele): The tandem repeat allele.

    Return:
        str: A string describing the reason why the allele failed filters, or None if the allele passed all filters.
    """
    total_repeats = tandem_repeat_allele.num_repeats_in_left_flank + tandem_repeat_allele.num_repeats_in_variant + tandem_repeat_allele.num_repeats_in_right_flank
    total_repeat_bases = tandem_repeat_allele.num_repeat_bases_in_left_flank + tandem_repeat_allele.num_repeat_bases_in_variant + tandem_repeat_allele.num_repeat_bases_in_right_flank
    repeat_unit = tandem_repeat_allele.repeat_unit



    if total_repeats == 1:
        # no repeat unit found in this allele
        return FILTER_ALLELE_INDEL_WITHOUT_REPEATS
    elif total_repeats < args.min_repeats:
        return FILTER_TR_ALLELE_NOT_ENOUGH_REPEATS % int(total_repeats)
    elif total_repeat_bases < args.min_tandem_repeat_length:
        return FILTER_TR_ALLELE_DOESNT_SPAN_ENOUGH_BASE_PAIRS % args.min_tandem_repeat_length
    elif len(repeat_unit) < args.min_repeat_unit_length:
        return FILTER_TR_ALLELE_REPEAT_UNIT_TOO_SHORT % args.min_repeat_unit_length
    elif len(repeat_unit) > args.max_repeat_unit_length:
        return FILTER_TR_ALLELE_REPEAT_UNIT_TOO_LONG % args.max_repeat_unit_length
    
    if "N" in tandem_repeat_allele.variant_and_flanks_repeat_sequence:
        return FILTER_ALLELE_WITH_N_BASES
    
    return None  # did not fail filters


def run_trf(alleles, args, counters, thread_id=1):
    """Run TRF on the given allele records.
    
    Args:
        alleles (list): List of Allele objects to run TRF on.
        args (argparse.Namespace): Command-line arguments parsed by parse_args().
        counters (dict): Dictionary of counters to collect summary stats about the number of TR variants found, etc.
        thread_id (int): ID of thread executing this function, starting from 0.

    Return:
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
    trf_runner = TRFRunner(args.trf_executable_path, html_mode=True,
                    min_motif_size=args.min_repeat_unit_length,
                    max_motif_size=args.max_repeat_unit_length)

    trf_runner.run_trf_on_fasta_file(trf_fasta_filename)

    # parse the TRF output
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

                # check that repeat boundary matches the variant breakpoint
                #if variant_bases_sum + fuzz < len(allele.variant_bases) or variant_bases_sum - fuzz > len(allele.variant_bases):
                #    continue

                trf_result["num_repeats_in_variant"] = repeat_count
                trf_result["tandem_repeat_bases_in_flank"]= max(0, trf_result["end_1based"] - len(allele.variant_bases))

                if left_or_right == "left":
                    trf_result["repeat_unit"] = trf_result["repeat_unit"][::-1]

                motif_size_to_matching_trf_results[trf_result["repeat_unit_length"]][left_or_right] = trf_result

        # filter out results where left and right have different repeat units
        motif_size_to_passing_trf_results = {}
        motif_size_to_tandem_repeat_allele = {}
        for motif_size, matching_trf_results in motif_size_to_matching_trf_results.items():
            if matching_trf_results.get("left") and matching_trf_results.get("right") and not are_repeat_units_similar(
                matching_trf_results["left"]["repeat_unit"], matching_trf_results["right"]["repeat_unit"]
            ):
                if args.debug: print(f"TRF filtered out: {allele}, filter reason: {matching_trf_results['left']['repeat_unit']} and {matching_trf_results['right']['repeat_unit']} are not similar")
                continue
                
            if matching_trf_results.get("right"):
                repeat_unit = matching_trf_results["right"]["repeat_unit"]
                motif_interruption_indices = list(sorted([position_in_motif for position_in_motif, _ in matching_trf_results["right"]["motif_positions_with_interruptions"].items()]))

            elif matching_trf_results.get("left"):
                repeat_unit = matching_trf_results["left"]["repeat_unit"]
                motif_interruption_indices = list(sorted([len(repeat_unit) - position_in_motif + 1 for position_in_motif, _ in matching_trf_results["left"]["motif_positions_with_interruptions"].items()]))

            else:
                raise ValueError(f"Logic error: No TRF results for allele {allele} with motif size {motif_size}")

            tandem_repeat_allele = TandemRepeatAllele(
                allele,
                repeat_unit,
                matching_trf_results.get("left", {}).get("tandem_repeat_bases_in_flank", 0),
                len(allele.variant_bases),
                matching_trf_results.get("right", {}).get("tandem_repeat_bases_in_flank", 0),
                detection_mode=DETECTION_MODE_TRF,
                is_pure_repeat=False,
                motif_interruption_indices=motif_interruption_indices or None,
            )
            
            if args.debug: print(f"TRF: checking if {tandem_repeat_allele} passes filters")
            filter_reason = check_if_tandem_repeat_allele_failed_filters(args, tandem_repeat_allele)
            if args.debug: print(f"TRF: {tandem_repeat_allele}, filter: {filter_reason}")
            if filter_reason is not None:
                continue

            motif_size_to_passing_trf_results[motif_size] = matching_trf_results
            motif_size_to_tandem_repeat_allele[motif_size] = tandem_repeat_allele

        if len(motif_size_to_passing_trf_results) == 0:
            if args.debug: print(f"TRF filtered out: {allele} because it has no repeat unit that passed filters")
            results.append((None, FILTER_ALLELE_INDEL_WITHOUT_REPEATS, allele))
            continue

        # get the entry with the smallest motif size
        best_motif_size = None
        for motif_size, matching_trf_results in sorted(motif_size_to_passing_trf_results.items()):
            best_motif_size = motif_size
            break

        # if the smallest motif size is 2 or less, then pick the entry with the largest alignment score instead
        if best_motif_size <= 2 and len(motif_size_to_passing_trf_results) > 1:
            # get the entry that has the max alignment score
            max_alignment_score = 0
            for motif_size, matching_trf_results in motif_size_to_passing_trf_results.items():
                alignment_score = 0
                for _trf_result in matching_trf_results.values():
                    alignment_score += _trf_result["alignment_score"]
                if alignment_score > max_alignment_score:
                    max_alignment_score = alignment_score
                    best_motif_size = motif_size

        tandem_repeat_allele = motif_size_to_tandem_repeat_allele[best_motif_size]
        results.append((tandem_repeat_allele, None, allele))
        
    return results


def are_repeat_units_similar(canonical_repeat_unit1, canonical_repeat_unit2):
    """Check if the two repeat units are similar enough to be considered the same repeat unit.
    
    Args:
        canonical_repeat_unit1 (str): canonical motif #1
        canonical_repeat_unit2 (str): canonical motif #2
    """
    
    if canonical_repeat_unit1 == canonical_repeat_unit2:
        return True

    if len(canonical_repeat_unit1) > 6 and len(canonical_repeat_unit2) > 6:
        fuzz = int(math.log10(2 * len(canonical_repeat_unit1))) + 1
        return abs(len(canonical_repeat_unit1) - len(canonical_repeat_unit2)) <= fuzz

    return False


def merge_overlapping_tandem_repeat_loci(tandem_repeat_alleles, counters):
    """Merge overlapping tandem repeats
    
    Args:
        tandem_repeat_alleles (list): list of TandemRepeatAllele objects
    """

    # For each chromosome, create an interval tree of the tandem repeat alleles
    tandem_repeat_alleles.sort(key=lambda x: (x.chrom, x.start_0based, x.end_1based))
    results = []
    for chrom, alleles in itertools.groupby(tandem_repeat_alleles, key=lambda x: x.chrom):
        alleles = list(alleles)

        current_i = 0
        current_end_1based = alleles[current_i].end_1based
        next_i = 1
        while next_i < len(alleles):
            # check if they have similar motifs
            merge_loci = False
            if alleles[next_i].start_0based < current_end_1based - len(alleles[next_i].repeat_unit) + 1 and (
                are_repeat_units_similar(alleles[current_i].canonical_repeat_unit, alleles[next_i].canonical_repeat_unit)):
                    merge_loci = True

            if merge_loci:
                #print(f"Merging {current_i}: {alleles[current_i]} and {next_i}: {alleles[next_i]}")
                if alleles[current_i].ref_interval_size < alleles[next_i].ref_interval_size:
                    current_i = next_i  # keep the locus definiton that has the larger interval
                current_end_1based = max(current_end_1based, alleles[next_i].end_1based)
            else:
                results.append(alleles[current_i])
                current_i = next_i
                current_end_1based = alleles[current_i].end_1based

            next_i += 1

        # append the last one
        results.append(alleles[current_i])
            

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

    with open(bed_output_path, "w") as f:
        for tandem_repeat_allele in tandem_repeat_alleles:
            f.write("\t".join(map(str, [
                tandem_repeat_allele.chrom,
                tandem_repeat_allele.start_0based,
                tandem_repeat_allele.end_1based,
                f"{tandem_repeat_allele.repeat_unit}:{len(tandem_repeat_allele.repeat_unit)}bp:{tandem_repeat_allele.detection_mode}" if detailed else tandem_repeat_allele.repeat_unit,
                len(tandem_repeat_allele.repeat_unit),
            ])) + "\n")

    os.system(f"bgzip -f {bed_output_path}")
    os.system(f"tabix -p bed {bed_output_path}.gz")

    print(f"Wrote {len(tandem_repeat_alleles):,d} tandem repeat alleles to {bed_output_path}.gz")


def write_tsv(tandem_repeat_alleles, args):
    """Write the tandem repeat alleles to a TSV file.
    
    Args:
        tandem_repeat_alleles (list): list of TandemRepeatAllele objects
        args (argparse.Namespace): command-line arguments parsed by parse_args()
    """

    tsv_output_path = f"{args.output_prefix}.tandem_repeats.tsv"

    with open(tsv_output_path, "w") as f:
        f.write("\t".join([
            "Chrom",
            "Start0Based",
            "End1Based",
            "Locus",
            "LocusId",
            "INS_or_DEL",
            "Motif",
            "MotifInterruptionIndices",
            "CanonicalMotif",
            "MotifSize",
            "NumRepeatsInReference",
            "VcfPos",
            "SummaryString",
            "IsFoundInReference",
            "IsPureRepeat",
            "DetectionMode",
        ]) + "\n")
        
        for tandem_repeat_allele in tandem_repeat_alleles:
            f.write("\t".join(map(str, [
                tandem_repeat_allele.chrom,
                tandem_repeat_allele.start_0based,
                tandem_repeat_allele.end_1based,
                f"{tandem_repeat_allele.chrom}:{tandem_repeat_allele.start_0based}-{tandem_repeat_allele.end_1based}",
                tandem_repeat_allele.locus_id,
                tandem_repeat_allele.ins_or_del,
                tandem_repeat_allele.repeat_unit,
                tandem_repeat_allele.motif_interruption_indices_string,
                tandem_repeat_allele.canonical_repeat_unit,
                len(tandem_repeat_allele.repeat_unit),
                tandem_repeat_allele.num_repeats_ref,
                tandem_repeat_allele.allele.pos,
                tandem_repeat_allele.summary_string,
                tandem_repeat_allele.end_1based > tandem_repeat_allele.start_0based,
                tandem_repeat_allele.is_pure_repeat,
                tandem_repeat_allele.detection_mode,
            ])) + "\n")

    os.system(f"bgzip -f {tsv_output_path}")
    print(f"Wrote {len(tandem_repeat_alleles):,d} tandem repeat alleles to {tsv_output_path}.gz")


def write_fasta(tandem_repeat_alleles, args):
    """Write the tandem repeat alleles to a FASTA file. For insertion alleles, the alternate allele sequence is written while for deletion
    alleles, the reference allele sequence is written.
    
    Args:
        tandem_repeat_alleles (list): list of TandemRepeatAllele objects
        args (argparse.Namespace): command-line arguments parsed by parse_args()
    """

    fasta_output_path = f"{args.output_prefix}.tandem_repeats.fasta"
    with open(fasta_output_path, "w") as f:
        for tandem_repeat_allele in tandem_repeat_alleles:
            f.write(f">{tandem_repeat_allele.locus_id}__{tandem_repeat_allele.summary_string}\n")
            if tandem_repeat_allele.ins_or_del == "INS":
                f.write(f"{tandem_repeat_allele.alt_allele_repeat_sequence}\n")
            else:
                f.write(f"{tandem_repeat_allele.ref_allele_repeat_sequence}\n")

    os.system(f"gzip -f {fasta_output_path}")
    print(f"Wrote {len(tandem_repeat_alleles):,d} tandem repeat sequences to {fasta_output_path}.gz")


def write_vcf(tandem_repeat_alleles, args, only_write_filtered_out_alleles=False):
    """Write variants that either are or aren't tandem repeats to a VCF file.
    
    Args:
        tandem_repeat_alleles (list): list of TandemRepeatAllele objects
        args (argparse.Namespace): command-line arguments parsed by parse_args()
        only_write_filtered_out_alleles (bool): if True, only write the variants that are not in the tandem_repeats_alleles list
    """

    input_files_to_close = []
    if args.interval:
        tabix_file = pysam.TabixFile(args.input_vcf_path)
        vcf_iterator = itertools.chain(
            (f"{line}\n" for line in tabix_file.header), 
            (line for interval in args.interval for line in tabix_file.fetch(interval)),
        )
        input_files_to_close.append(tabix_file)
    else:
        vcf_iterator = open_file(args.input_vcf_path, is_text_file=True)
        input_files_to_close.append(vcf_iterator)

    # iterate over all VCF rows
    if only_write_filtered_out_alleles:
        output_vcf_path = f"{args.output_prefix}.not_tandem_repeats.vcf"
    else:
        output_vcf_path = f"{args.output_prefix}.tandem_repeats.vcf"

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
                vcf_fields[7] += f";MOTIF_SIZE={len(tr.repeat_unit)}"
                vcf_fields[7] += f";START_0BASED={tr.start_0based}"
                vcf_fields[7] += f";END={tr.end_1based}"
                if tr.motif_interruption_indices:
                    vcf_fields[7] += f";INTERUPTIONS={tr.motif_interruption_indices_string}"
                vcf_fields[7] += f";DETECTED={tr.detection_mode}"

                line = "\t".join(vcf_fields) + "\n"

            if (not only_write_filtered_out_alleles and is_tandem_repeat) or (only_write_filtered_out_alleles and not is_tandem_repeat):
                f.write(line)
                output_line_counter += 1

    os.system(f"bgzip -f {output_vcf_path}")
    os.system(f"tabix -p vcf {output_vcf_path}.gz")

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
                print("--------------")
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


if __name__ == "__main__":
    main()


