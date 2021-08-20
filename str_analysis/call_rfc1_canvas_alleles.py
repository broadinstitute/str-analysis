#!/usr/bin/env python3

"""
This script takes a bam/cram file and outputs a .json file with informatino
"""

import argparse
import collections
import json
import os
import peewee as pw
import pysam
from pprint import pprint
import re
import subprocess

from str_analysis.utils.canonical_repeat_unit import compute_canonical_repeat_unit
from str_analysis.utils.ehdn_info_for_locus import parse_ehdn_info_for_locus
from str_analysis.utils.most_frequent_repeat_unit import compute_most_frequent_repeat_unit

CANVAS_REPEAT_UNIT_SIZE = 5
MARGIN = 7   # base pairs
FLANK_SIZE = 2000   # base pairs
MIN_MAPQ = 3
MIN_FRACTION_OF_BASES_COVERED = 0.7

RFC1_LOCUS_COORDS_0BASED = {
    "37": ("4", 39350044, 39350099),
    "38": ("chr4", 39348424, 39348479),   # chr4:39348425-39348479
}

RFC1_LOCUS_KNOWN_ALLELES_BY_CATEGORY = {
        "BENIGN": {"AAAAG", "AAAGG"},
    "PATHOGENIC": {"AAGGG", "ACAGG"}, # ACAGG is from "A Māori specific RFC1 pathogenic repeat..." [Beecroft 2021]
    #"UNCERTAIN": {"AAGAG", "AGAGG",}, # from [Akcimen 2019]
}

ALL_RFC1_LOCUS_KNOWN_ALLELES = {a for allele_set in RFC1_LOCUS_KNOWN_ALLELES_BY_CATEGORY.values() for a in allele_set}

GENOME_VERSION_ALIASES = {
    "GRCh37": "37", "hg19": "37", "hg37": "37", "37": "37",
    "GRCh38": "38", "hg38": "38", "38": "38",
}


OFFTARGET_REGIONS = {
    "38": {
        "AAAAG": [
            "chr1:157959795-157960160",
            "chr10:38112579-38112976",
            "chr11:114394023-114394461",
            "chr12:125823177-125823548",
            "chr18:6909097-6909460",
            "chr18:46255244-46255614",
            "chr2:162752749-162753183",
            "chr4:104453000-104453524",
            "chr5:70672067-70672497",
            "chr5:84539198-84539654",
            "chr6:51553892-51554235",
            "chr8:69433867-69434260",
            "chr9:19473778-19474188",
            "chr9:99458785-99459128",
            "chrUn_JTFH01001554v1_decoy:488-852",
        ],
        "AAGGG":  [
            "chr10:25874193-25874766",
            "chr21:39583664-39584274",
            "chr4:43136889-43137319",
        ],
        "ACAGG":  [
            "chr12:12047082-12047450",
            "chr12:15936726-15937132",
            "chr16:81823861-81824253",
            "chr17:73224663-73225192",
            "chr2:109248715-109249087",
            "chr7:18430857-18431231",
            "chr7:137196425-137196784",
            "chr7:149793911-149794272",
            "chr8:105878913-105879266",
            "chr9:136061962-136062334",
            "chrX:68262640-68263034",
        ],

    }
}

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-g", "--genome-version", choices=GENOME_VERSION_ALIASES.keys(), required=True)
    p.add_argument("-R", "--reference", help="Reference fasta path. The reference fasta is sometimes necessary for "
                                             "decoding cram files.")
    p.add_argument("-o", "--output-prefix", help="Output filename prefix")
    p.add_argument("-s", "--sample-id", help="The sample id to put in the output json file. If not specified, it "
        "will be retrieved from the bam/cram header or filename prefix.")
    p.add_argument("-e", "--ehdn-profile", help="If specified, information relevant to the RFC1 locus will be "
        "transferred from this ExpansionHunterDenovo profile to the output json")
    p.add_argument("-r", "--run-expansion-hunter", action="store_true", help="If this option is specified, this "
         "script will run ExpansionHunter once for each of the motif(s) it detects at the RFC1 locus. "
         "ExpansionHunter doesn't currently support genotyping multiallelic repeats such as RFC1 where "
         "an individual may have 2 alleles with motifs that differ from eachother (and from the reference motif). "
         "Running ExpansionHunter separately for each motif provides a work-around.")
    p.add_argument("--expansion-hunter-path", help="The path of the ExpansionHunter executable to use if -r is "
        "specified. This must be ExpansionHunter v3 or v4.", default="ExpansionHunter")
    p.add_argument("--reference-fasta-path", help="The path of the reference genome fasta to use if -r is "
        "specified")
    p.add_argument("-t", "--temp-dir", help="Directory for intermediate files such as those generated when "
                                            "running ExpansionHunter", default=".")
    p.add_argument("-v", "--verbose", action="store_true", help="Print detailed log messages")
    p.add_argument("bam_or_cram_path", help="bam or cram path")

    args = p.parse_args()
    args.genome_version = GENOME_VERSION_ALIASES[args.genome_version]

    if not os.path.isfile(args.bam_or_cram_path):
        p.error(f"{args.bam_or_cram_path} not found")

    if args.reference and not os.path.isfile(args.reference):
        p.error(f"{args.reference} not found")

    if args.run_expansion_hunter and not args.reference_fasta_path:
        p.error("--reference-fasta-path is required when --run-expansion-hunter is used")

    return args


def generate_variant_catalog(locus_id, repeat_unit, chrom, start_1based, end_1based, offtarget_regions=None):
    return {
        "LocusId": locus_id,
        "LocusStructure": f"({repeat_unit})*",
        "ReferenceRegion": f"{chrom}:{start_1based - 1}-{end_1based}",
        "VariantType": "RareRepeat",
        "OfftargetRegions": [] if not offtarget_regions else offtarget_regions,
    }


def generate_ehv2_repeat_spec(locus_id, repeat_unit, chrom, start_1based, end_1based, offtarget_regions=None):
    return {
        "RepeatId": locus_id,
        "RepeatUnit": repeat_unit,
        "TargetRegion": f"{chrom}:{start_1based - 1}-{end_1based}",
        "CommonUnit": "false",
        "OffTargetRegions": [] if not offtarget_regions else offtarget_regions,
    }


def run_expansion_hunter(
        sample_id,
        expansion_hunter_path,
        genome_version,
        reference_fasta_path,
        bam_or_cram_path,
        well_supported_repeat_units,
        result,
        directory=".",
        verbose=False,
):
    """
    :param well_supported_repeat_units:
    :param expansion_hunter_path:
    :param result:
    :param verbose:
    :return:
    """

    if genome_version not in ("37", "38"):
        raise ValueError(f"Unexpected genome version: {genome_version}. Must be '37' or '38'")

    chrom, start_1based, end_1based = RFC1_LOCUS_COORDS_0BASED[genome_version]

    for repeat_unit in well_supported_repeat_units:
        # get off-target regions
        """
        db = pw.SqliteDatabase(f"../data/{genome_version}")
        cursor = exons_db.execute_sql(
            f"SELECT annotation_source, gene_id, transcript_id, gene_name, gene_type, transcript_type FROM exons "
            f"WHERE chrom='{repeat.chrom}' AND start_1based <= {repeat.end_1based} AND end_1based >= {repeat.start_1based}")
        #  f"AND transcript_type NOT IN ('transcribed_unprocessed_pseudogene', 'unprocessed_pseudogene', 'processed_pseudogene', 'lncRNA')"

        result = cursor.fetchone()
        cursor.close()
    
        if result:
            repeat.is_exonic = True
            repeat.gene_annotation_source = result[0]
            repeat.gene_id = result[1]
            repeat.transcript_id = result[2]
            repeat.gene_name = result[3]
            repeat.gene_type = result[4]
            repeat.transcript_type = result[5]
    
            yield repeat
    
        """

        if False:
            # generate variant catalog
            variant_catalog = generate_variant_catalog(f"RFC1_{repeat_unit}",
                repeat_unit, chrom, start_1based, end_1based, offtarget_regions=OFFTARGET_REGIONS[genome_version][repeat_unit])

            variant_catalog_path = os.path.join(directory, f"{repeat_unit}.variant_catalog.json")
            with open(variant_catalog_path, "wt") as f:
                json.dump([variant_catalog], f)

            # run expansion hunter
            print("--"*10)
            print(f"Running ExpansionHunter on {sample_id} for repeat unit {repeat_unit}")
            if verbose:
                print("Using variant catalog: ")
                pprint(variant_catalog)

            filename_prefix = f"{sample_id}.{repeat_unit}"
            output_prefix = f"{os.path.join(directory, filename_prefix)}.expansion_hunter4"
            subprocess.check_call(f"""{expansion_hunter_path} \
                --sex male \
                --aligner path-aligner \
                --reference {reference_fasta_path} \
                --reads {bam_or_cram_path} \
                --variant-catalog {variant_catalog_path} \
                --output-prefix {output_prefix}
            """, shell=True)

            # parse result
            with open(f"{output_prefix}.json", "rt") as f:
                expansion_hunter_output_json = json.load(f)

            if verbose:
                pprint(expansion_hunter_output_json)

        # run EHv2
        if True:
            os.system(f"mkdir -p repeat_spec_{repeat_unit}")
            repeat_spec = generate_ehv2_repeat_spec(f"RFC1_{repeat_unit}",
                repeat_unit, chrom, start_1based, end_1based, offtarget_regions=OFFTARGET_REGIONS[genome_version][repeat_unit])
            with open(os.path.join(directory, f"repeat_spec_{repeat_unit}/{repeat_unit}.json"), "wt") as f:
                json.dump(repeat_spec, f)

            #bam_dir = os.path.abspath(os.path.dirname(bam_or_cram_path))
            #reference_fasta_dir = os.path.abspath(os.path.dirname(reference_fasta_path))
            #current_dir = os.path.abspath(os.getcwd())
            output_prefix = f"{sample_id}.{repeat_unit}.expansion_hunter2"
            subprocess.check_call(f"""
                /Users/weisburd/p1/bin/ExpansionHunter \
                    --ref-fasta {reference_fasta_path} \
                    --bam {bam_or_cram_path} \
                    --repeat-specs repeat_spec_{repeat_unit} \
                    --read-depth {result['left_flank_coverage']} \
                    --sex male \
                    --vcf {output_prefix}.vcf \
                    --json {output_prefix}.json \
                    --log {output_prefix}.log
            """, shell=True)
            # parse result
            with open(f"{output_prefix}.json", "rt") as f:
                expansion_hunter_output_json = json.load(f)

            if verbose:
                pprint(expansion_hunter_output_json)


    # reorder well_supported_repeat_units based on EH results
    return well_supported_repeat_units


def main():
    args = parse_args()

    bam_cram_prefix = re.sub(".bam$|.cram$", "", os.path.basename(args.bam_or_cram_path))
    args.output_prefix = args.output_prefix or bam_cram_prefix

    chrom, locus_start_0based, locus_end = RFC1_LOCUS_COORDS_0BASED[args.genome_version]

    result = {}

    # process bam/cram
    print(f"Processing {args.bam_or_cram_path}")
    with pysam.Samfile(args.bam_or_cram_path, reference_filename=args.reference) as f:

        # try to get sample id from bam/cram header
        if not args.sample_id:
            read_groups = f.header.as_dict().get("RG", [])
            if read_groups:
                args.sample_id = read_groups[0].get("SM")
                if args.sample_id:
                    print(f"Using sample id '{args.sample_id}' from the bam/cram header")

            if not args.sample_id:
                args.sample_id = bam_cram_prefix
                print(f"Using sample id '{args.sample_id}' based on the input filename prefix")

        read_length = None
        while read_length is None:
            try:
                r = next(f)
                read_length = r.infer_read_length()
            except StopIteration:
                raise ValueError(f"{args.bam_or_cram_path} contains too few reads")

        # count reads in the left & right flanks to estimate read depth
        # NOTE: f.fetch retrieves all reads that *overlap* the given interval
        left_flank_n_well_aligned_bases = sum((r.query_alignment_length for r in f.fetch(
                chrom,
                locus_start_0based - MARGIN - FLANK_SIZE,
                locus_start_0based - MARGIN,
            ) if r.mapq >= MIN_MAPQ))
        right_flank_n_well_aligned_bases = sum((r.query_alignment_length for r in f.fetch(
                chrom,
                locus_end + MARGIN,
                locus_end + MARGIN + FLANK_SIZE,
            ) if r.mapq >= MIN_MAPQ))

        # get all sequences that overlap the RFC1 locus (regardless of whether they're soft-clipped)
        overlapping_sequences = []
        for r in f.fetch(chrom, locus_start_0based - MARGIN, locus_end + MARGIN):
            # see https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.query_alignment_sequence
            if r.mapq < MIN_MAPQ:
                continue

            read_sequence = r.seq
            #has_soft_clipped_bases_on_left = r.query_alignment_start > 0
            #has_soft_clipped_bases_on_right = r.query_alignment_end < len(read_sequence)
            read_start_pos_including_soft_clips = r.reference_start - r.query_alignment_start
            read_end_pos_including_soft_clips = read_start_pos_including_soft_clips + len(read_sequence)
            start_offset = 0
            if read_start_pos_including_soft_clips < locus_start_0based:
                start_offset = locus_start_0based - read_start_pos_including_soft_clips

            end_offset = len(read_sequence)
            if read_end_pos_including_soft_clips > locus_end:
                end_offset = locus_end - read_end_pos_including_soft_clips

            relevant_bases = read_sequence[start_offset:end_offset]
            if len(relevant_bases) >= CANVAS_REPEAT_UNIT_SIZE:
                overlapping_sequences.append(relevant_bases)

    left_flank_coverage = left_flank_n_well_aligned_bases / FLANK_SIZE
    right_flank_coverage = right_flank_n_well_aligned_bases / FLANK_SIZE
    result.update({
        "sample_id": args.sample_id,
        "left_flank_coverage": left_flank_coverage,
        "right_flank_coverage": right_flank_coverage,
    })

    # compute the repeat unit(s) found in the soft-clipped reads, and how many times each one occurs
    repeat_unit_to_n_occurrences = collections.defaultdict(int)
    repeat_unit_to_read_count = collections.defaultdict(int)
    for overlapping_sequence in overlapping_sequences:
        # in gnomAD, EHdn sometimes finds 6bp repeat units (eg. AAAGGG), so check for those as well
        for repeat_unit_size in (CANVAS_REPEAT_UNIT_SIZE, CANVAS_REPEAT_UNIT_SIZE + 1):
            if len(overlapping_sequence) < repeat_unit_size:
                continue

            repeat_unit, count = compute_most_frequent_repeat_unit(
                overlapping_sequence,
                repeat_unit_size=repeat_unit_size,
                min_occurrences=3,
                min_fraction_bases_covered=MIN_FRACTION_OF_BASES_COVERED)

            if args.verbose:
                if repeat_unit:
                    print(f"Found {compute_canonical_repeat_unit(repeat_unit)} occurs {count}x in read bases that "
                          f"overlap the RFC1 locus: {overlapping_sequence}")
                else:
                    if repeat_unit_size == CANVAS_REPEAT_UNIT_SIZE:
                        print(f"Didn't find a consistent {repeat_unit_size}bp repeat unit in read bases "
                            f"that overlap the RFC1 locus: {overlapping_sequence}")

            if repeat_unit is not None:
                canonical_repeat_unit = compute_canonical_repeat_unit(repeat_unit)
                repeat_unit_to_read_count[canonical_repeat_unit] += 1
                repeat_unit_to_n_occurrences[canonical_repeat_unit] += count

    result.update({
        "found_n_reads_overlap_rfc1_locus": len(overlapping_sequences),
        "found_repeats_in_n_reads": sum(repeat_unit_to_read_count.values()),
        "found_repeats_in_fraction_of_reads": sum(repeat_unit_to_read_count.values())/len(overlapping_sequences),
    })

    # evaluate the repeat units
    well_supported_repeat_units = []
    for repeat_unit, read_count in repeat_unit_to_read_count.items():
        if "N" in repeat_unit:
            continue
        if read_count < 3:
            continue
        well_supported_repeat_units.append(repeat_unit)

    well_supported_repeat_units.sort(key=lambda repeat_unit: repeat_unit_to_n_occurrences[repeat_unit], reverse=True)

    if args.run_expansion_hunter:
        well_supported_repeat_units = run_expansion_hunter(
            args.sample_id,
            args.expansion_hunter_path,
            args.genome_version,
            args.reference_fasta_path,
            args.bam_or_cram_path,
            well_supported_repeat_units,
            result,
            directory=args.temp_dir,
            verbose=args.verbose)

    n_pathogenic_alleles = 0
    n_benign_alleles = 0
    n_total_well_supported_alleles = 0
    for i in 0, 1:
        allele_number = i + 1
        result.update({
            f"allele{allele_number}_repeat_unit": None,
            f"allele{allele_number}_read_count": None,
            f"allele{allele_number}_n_occurrences": None,
        })

        if len(well_supported_repeat_units) <= i:
            continue
        n_total_well_supported_alleles += 1
        repeat_unit = well_supported_repeat_units[i]
        read_count = repeat_unit_to_read_count.get(repeat_unit)
        n_occurrences = repeat_unit_to_n_occurrences.get(repeat_unit)
        if repeat_unit in RFC1_LOCUS_KNOWN_ALLELES_BY_CATEGORY["PATHOGENIC"]:
            n_pathogenic_alleles += 1
        elif repeat_unit in RFC1_LOCUS_KNOWN_ALLELES_BY_CATEGORY["BENIGN"]:
            n_benign_alleles += 1

        result.update({
            f"allele{allele_number}_repeat_unit": repeat_unit,
            f"allele{allele_number}_read_count": read_count,
            f"allele{allele_number}_n_occurrences": n_occurrences,
        })

    # decide which combination of alleles is supported by the data
    # NOTE: there's no attempt to determine the size of the expansion and whether it's in the pathogenic range
    final_call = None
    if n_total_well_supported_alleles > 0 and left_flank_coverage >= 5 and right_flank_coverage >= 5:
        if n_pathogenic_alleles == n_total_well_supported_alleles:
            # reads support only known pathogenic allele(s)
            final_call = "PATHOGENIC MOTIF / PATHOGENIC MOTIF"
        elif n_benign_alleles == n_total_well_supported_alleles:
            # reads support only known benign allele(s)
            final_call = "BENIGN MOTIF / BENIGN MOTIF"
        elif n_benign_alleles == 0 and n_pathogenic_alleles == 0:
            # reads support one or more non-reference alleles of unknown significance
            final_call = "MOTIF OF UNCERTAIN SIGNIFICANCE / MOTIF OF UNCERTAIN SIGNIFICANCE"
        elif n_benign_alleles > 0 and n_pathogenic_alleles > 0:
            # reads support one known benign allele and one pathogenic allele
            final_call = "BENIGN MOTIF / PATHOGENIC MOTIF"
        elif n_pathogenic_alleles > 0:
            # reads support one pathogenic allele and at least one other allele of unknown significance
            final_call = "PATHOGENIC MOTIF / MOTIF OF UNCERTAIN SIGNIFICANCE"
        elif n_benign_alleles > 0:
            # reads support one known benign allele and at least one other allele of unknown significance
            final_call = "BENIGN MOTIF / MOTIF OF UNCERTAIN SIGNIFICANCE"

    result.update({
        "n_total_well_supported_alleles": n_total_well_supported_alleles,
        "n_benign_alleles": n_benign_alleles,
        "n_pathogenic_alleles": n_pathogenic_alleles,
    })

    result["call"] = final_call

    print(f"Final call: {final_call}")

    # process EHdn profile if one was provided
    if args.ehdn_profile:
        print(f"Parsing {args.ehdn_profile}")
        with open(args.ehdn_profile, "rt") as f:
            data = json.load(f)

        records, sample_read_depth, _ = parse_ehdn_info_for_locus(data, chrom, locus_start_0based, locus_end)
        result["ehdn_sample_read_depth"] = sample_read_depth

        # get the 2 alleles with the most read support
        records.sort(key=lambda r: (
            -r["anchored_irr_count_for_this_repeat_unit_and_region"],
            -r["total_irr_count_for_this_repeat_unit_and_region"],
            r["repeat_unit"],
        ))

        for i in 0, 1:
            record = records[i] if len(records) > i else {}
            allele_number = i + 1
            result.update({
                f"ehdn_allele{allele_number}_repeat_unit": record.get("repeat_unit"),
                f"ehdn_allele{allele_number}_anchored_irr_count": record.get("anchored_irr_count_for_this_repeat_unit_and_region"),
                f"ehdn_allele{allele_number}_n_anchored_regions": record.get("n_anchored_regions_for_this_repeat_unit"),
                f"ehdn_allele{allele_number}_paired_irr_count": record.get("paired_irr_count_for_this_repeat_unit"),
                f"ehdn_allele{allele_number}_total_irr_count": record.get("total_irr_count_for_this_repeat_unit_and_region"),
            })

    # generate output
    output_filename = f"{args.output_prefix}.rfc1_canvas_alleles.json"
    with open(output_filename, "wt") as f:
        json.dump(result, f, indent=2)
    print(f"Wrote results to {output_filename}")
    if args.verbose:
        pprint(result)


if __name__ == "__main__":
    main()