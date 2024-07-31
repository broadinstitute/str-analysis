"""This script converts a TRGT output VCF to the ExpansionHunter output .json format
which I use as a common input format in downstream scripts.
"""

"""
TRGT output vcf example 1: hom-ref call

chr1	
674824	
.	
AGGAGCAGG	
.	
0	
.	T
RID=chr1_674823_674832;END=674832;MOTIFS=AGG;STRUC=(AGG)n	
GT:AL:ALLR:SD:MC:MS:AP:AM	
0/0:9,9:9-9,9-9:3,2:3,3:0(0-9),0(0-9):0.888889,0.888889:.,.

--------------------
TRGT output vcf example 2: multi-allelic call

chr1	
3668516	
.	
GTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT	
GTGTGTGTGTGTGTGTGTGTGT,GTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT	
0	
.	
TRID=chr1_3668515_3668547;END=3668547;MOTIFS=GT;STRUC=(GT)n	
GT:AL:ALLR:SD:MC:MS:AP:AM	
1/2:22,36:22-22,34-38:7,14:11,18:0(0-22),0(0-36):1,1:.,.
"""

"""
ExpansionHunter output format:

  "LocusResults": {
        "chr12-57610122-57610131-GCA": {
          "AlleleCount": 2,
          "Coverage": 50.469442942130875,
          "FragmentLength": 433,
          "LocusId": "chr12-57610122-57610131-GCA",
          "ReadLength": 151,
          "Variants": {
            "chr12-57610122-57610131-GCA": {
              "CountsOfFlankingReads": "(1, 1), (2, 4)",
              "CountsOfInrepeatReads": "()",
              "CountsOfSpanningReads": "(2, 1), (3, 48), (6, 1)",
              "Genotype": "3/3",
              "GenotypeConfidenceInterval": "3-3/3-3",
              "ReferenceRegion": "chr12:57610122-57610131",
              "RepeatUnit": "GCA",
              "VariantId": "chr12-57610122-57610131-GCA",
              "VariantType": "Repeat"
            }
          }
        },

  "SampleParameters": {
        "SampleId": "NA19239",
        "Sex": "Female"
  }
"""


import argparse
import gzip
import simplejson as json
import re
from tqdm import tqdm

def main():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--discard-hom-ref", action="store_true", help="Discard hom-ref calls")
    p.add_argument("--dont-output-REF-ALT-fields", action="store_true", help="Exclude the VCF REF and ALT fields from "
                   "the output as they can take up a lot of space.")
    p.add_argument("--parse-genotype-from-AL-field", action="store_true", help="By default, the genotype is taken from "
                   "the MC field. If this option is specified, it will be taken from the AL field for loci with only "
                   "1 motif, and from the MC field only for loci with multiple motifs.")
    grp = p.add_mutually_exclusive_group()
    grp.add_argument("--parse-reference-region-from-locus-id", action="store_true", help="Parse the reference region from the " 
                   "locus ID (which is expected to have the format '{chrom}-{start_0based}-{end}-{motif}') instead of "
                   "from the VCF row position and END fields.")
    grp.add_argument("--set-locus-id", action="store_true", help="If specified, the locus id will be set to "
                     "'{chrom}-{start_0based}-{end}-{motif} based on the VCF row position and END fields.")
    p.add_argument("--verbose", action="store_true", help="Print verbose output")
    p.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    p.add_argument("--sample-id",
                   help="If not specified, the sample id will be parsed from the last column of the vcf header.")
    p.add_argument("vcf_path", help="TRGT vcf path")
    args = p.parse_args()

    print(f"Processing {args.vcf_path}")
    locus_results = process_trgt_vcf(
        args.vcf_path,
        sample_id=args.sample_id,
        discard_hom_ref=args.discard_hom_ref,
        use_trgt_locus_id=not args.set_locus_id,
        parse_genotype_from_AL_field=args.parse_genotype_from_AL_field,
        parse_reference_region_from_locus_id=args.parse_reference_region_from_locus_id,
        dont_output_REF_ALT_fields=args.dont_output_REF_ALT_fields,
        verbose=args.verbose,
        show_progress_bar=args.show_progress_bar,
    )

    output_json_path = re.sub(".vcf(.gz)?$", "", args.vcf_path) + ".json"
    print(f"Writing {len(locus_results['LocusResults']):,d} loci to {output_json_path}")
    with open(output_json_path, "wt") as f:
        json.dump(locus_results, f, indent=3, ignore_nan=True)


def process_trgt_vcf(vcf_path, sample_id=None, discard_hom_ref=True, use_trgt_locus_id=False,
                     parse_genotype_from_AL_field=False,
                     parse_reference_region_from_locus_id=False,
                     dont_output_REF_ALT_fields=False,
                     verbose=False, show_progress_bar=False):
    locus_results = {
        "LocusResults": {},
        "SampleParameters": {
            "SampleId": sample_id,
            "Sex": None,
        },
    }

    fopen = gzip.open if vcf_path.endswith("gz") else open

    with fopen(vcf_path, "rt") as vcf:
        for line in vcf:
            if line.startswith("#CHROM"):
                header_fields = line.strip().split("\t")
                if sample_id is None and len(header_fields) == 10:
                    print(f"Got sample id '{header_fields[9]}' from the VCF header")
                    locus_results["SampleParameters"]["SampleId"] = header_fields[9]
            if not line.startswith("#"):
                break

    with fopen(vcf_path, "rt") as vcf:
        line_counter = 0
        if show_progress_bar:
            vcf = tqdm(vcf, unit=" vcf records", unit_scale=True, unit_divisor=1000)

        for line in vcf:
            if line.startswith("#"):
                continue

            line_counter += 1
            fields = line.strip().split("\t")
            chrom = fields[0]
            start_1based = int(fields[1])
            info = fields[7]
            if not fields[9] or fields[9] == ".":  # no genotype
                continue

            info_dict = dict([key_value.split("=") for key_value in info.split(";")])
            genotype_fields = fields[8].split(":")
            genotype_values = fields[9].split(":")
            genotype_dict = dict(zip(genotype_fields, genotype_values))

            if discard_hom_ref and genotype_dict["GT"] == "0/0":
                continue

            # GT:AL:ALLR:SD:MC:MS:AP:AM
            if genotype_dict["AL"] == ".":   # no genotype
                continue

            try:
                end_1based = int(info_dict["END"])

                motifs = info_dict["MOTIFS"].split(",")
                allele_sizes_bp = [int(allele_size_bp) for allele_size_bp in genotype_dict["AL"].split(",")]
                flip_alleles = len(allele_sizes_bp) == 2 and allele_sizes_bp[0] > allele_sizes_bp[1]
                if flip_alleles:
                    for key in "AL", "ALLR", "SD", "MC", "MS", "AP", "AM":
                        genotype_dict[key] = ",".join(genotype_dict[key].split(",")[::-1])

                parse_genotype_from_AL_field = parse_genotype_from_AL_field and len(motifs) == 1
                if parse_genotype_from_AL_field:
                    locus_id = info_dict["TRID"] if use_trgt_locus_id else f"{chrom}-{start_1based - 1}-{end_1based}-{repeat_unit}"
                else:
                    locus_id = info_dict["TRID"]

                if parse_reference_region_from_locus_id:
                    locus_id_fields = re.split("[_-]", locus_id)
                    if len(locus_id_fields) < 3:
                        raise ValueError(f"Unable to parse chrom, start, end from locus ID field: {locus_id}")
                    chrom = locus_id_fields[0]
                    start_1based = int(locus_id_fields[1]) + 1
                    end_1based = int(locus_id_fields[2])

                if parse_genotype_from_AL_field:
                    repeat_unit = motifs[0]

                    genotype = []
                    for allele_size_bp in genotype_dict["AL"].split(","):
                        genotype.append(str(int(allele_size_bp)//len(repeat_unit)))
                    genotype = "/".join(genotype)

                    genotype_CI = []
                    for ci in genotype_dict["ALLR"].split(","):
                        ci = ci.split("-")
                        ci = [str(int(ci_value)//len(repeat_unit)) for ci_value in ci]
                        genotype_CI.append("-".join(ci))
                    genotype_CI = "/".join(genotype_CI)

                    locus_results["LocusResults"][locus_id] = {
                        "AlleleCount": genotype_dict["AL"].count(",") + 1,
                        "LocusId": locus_id,
                        #"Coverage": None,
                        #"ReadLength": None,
                        #"FragmentLength": None,
                        "Variants": {
                            locus_id: info_dict | genotype_dict | {
                                "Genotype": genotype,   #"17/17",
                                "GenotypeConfidenceInterval": genotype_CI, #"17-17/17-17",
                                "ReferenceRegion": f"{chrom}:{start_1based - 1}-{end_1based}",
                                "RepeatUnit": repeat_unit,
                                "VariantId": locus_id,
                                "VariantType": "Repeat",
                                "Ref": None if dont_output_REF_ALT_fields else fields[3],
                                "Alt": None if dont_output_REF_ALT_fields else fields[4],
                            }
                        }
                    }
                else:
                    motif_counts = [[int(c) for c in mc.split("_")] for mc in genotype_dict["MC"].split(",")]

                    motif_count_genotypes = []
                    motif_count_genotypes_CIs = []
                    for motif_i, motif in enumerate(motifs):
                        genotype = []
                        for motif_counts_for_haplotype in motif_counts:
                            genotype.append(str(motif_counts_for_haplotype[motif_i]))

                        motif_count_genotypes.append("/".join(genotype))
                        motif_count_genotypes_CIs.append("/".join([f"{g}-{g}" for g in genotype]))
                    # parse allele sizes from string like '0(0-9)_1(9-24),0(0-9)_1(9-24)'
                    #interval_sizes_bp = [
                    #    re.findall(r"(\d+)[(](\d+)-(\d+)[)]", interval_size_bp) for interval_size_bp in genotype_dict["MS"].split(",")
                    #]
                    #interval_sizes_bp = [
                    #    {
                    #        int(motif_i): (int(end) - int(start)) for motif_i, start, end in interval_sizes_bp_for_haplotype
                    #    }
                    #    for interval_sizes_bp_for_haplotype in interval_sizes_bp
                    #]

                    #interval_size_genotypes = []
                    #for motif_i, motif in enumerate(motifs):
                    #    genotype = []
                    #    for interval_sizes_bp_for_haplotype in interval_sizes_bp:
                    #        span_of_current_motif = interval_sizes_bp_for_haplotype.get(motif_i, 0)
                    #        genotype.append(span_of_current_motif//len(motif))
                    #    interval_size_genotypes.append("/".join(map(str, genotype)))

                    #for motif_count_genotype, interval_size_genotype in zip(motif_count_genotypes, interval_size_genotypes):
                    #    if motif_count_genotype != interval_size_genotype:
                    #        raise ValueError(f"Motif count genotype {motif_count_genotype} != interval size genotype {interval_size_genotype} in locus {locus_id}: {line}")

                    variant_records = {}
                    for motif_i, motif in enumerate(motifs):
                        variant_id = f"{locus_id}-m{motif_i}-{motif}"
                        variant_records[variant_id] = info_dict | genotype_dict | {
                            "Genotype": motif_count_genotypes[motif_i],   #"17/17",
                            "GenotypeConfidenceInterval": motif_count_genotypes_CIs[motif_i], #"17-17/17-17",
                            "ReferenceRegion": f"{chrom}:{start_1based - 1}-{end_1based}",
                            "RepeatUnit": motif,
                            "VariantId": variant_id,
                            "VariantType": "Repeat",
                            "Ref": None if dont_output_REF_ALT_fields else fields[3],
                            "Alt": None if dont_output_REF_ALT_fields else fields[4],
                        }


                    locus_results["LocusResults"][locus_id] = {
                        "AlleleCount": genotype_dict["AL"].count(",") + 1,
                        "LocusId": locus_id,
                        #"Coverage": None,
                        #"ReadLength": None,
                        #"FragmentLength": None,
                        "Variants": variant_records,
                    }

            except Exception as e:
                print(f"Error on vcf record #{line_counter}: {e}")
                print(line)
                print(genotype_dict)
                # print stack trace
                import traceback
                traceback.print_exc()


    return locus_results


if __name__ == "__main__":
    main()
