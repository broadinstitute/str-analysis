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
import os
from pprint import pprint
import re


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sample-id",
                   help="If not specified, the sample id will be parsed from the last column of the vcf header.")
    p.add_argument("vcf_path", help="TRGT vcf path")
    args = p.parse_args()

    print(f"Processing {args.vcf_path}")
    locus_results = process_trgt_vcf(args.vcf_path, sample_id=args.sample_id)

    output_json_path = re.sub(".vcf(.gz)?$", "", args.vcf_path) + ".json"
    print(f"Writing results for {len(locus_results['LocusResults']):,d} loci to {output_json_path}")
    with open(output_json_path, "wt") as f:
        json.dump(locus_results, f, indent=3)


def process_trgt_vcf(vcf_path, sample_id=None):
    locus_results = {
        "LocusResults": {},
        "SampleParameters": {
            "SampleId": sample_id,
            "Sex": None,
        },
    }

    fopen = gzip.open if vcf_path.endswith("gz") else open
    with fopen(vcf_path, "rt") as vcf:
        line_counter = 0
        for line in vcf:
            if line.startswith("#CHROM"):
                header_fields = line.strip().split("\t")
                if sample_id is None and len(header_fields) == 10:
                    print(f"Got sample id '{header_fields[9]}' from the VCF header")
                    locus_results["SampleParameters"]["SampleId"] = header_fields[9]

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

            if genotype_dict["AL"] == ".":   # no genotype
                continue

            try:

                if "," in info_dict["MOTIFS"]:
                    raise ValueError(f"Support for multiple motifs not yet implemented: {info_dict['MOTIFS']}")

                repeat_unit = info_dict["MOTIFS"]
                end_1based = int(info_dict["END"])

                locus_id = f"{chrom}-{start_1based - 1}-{end_1based}-{repeat_unit}"

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
                    "Coverage": None,
                    "ReadLength": None,
                    "FragmentLength": None,
                    "Variants": {
                        locus_id: info_dict | genotype_dict | {
                            "Genotype": genotype,   #"17/17",
                            "GenotypeConfidenceInterval": genotype_CI, #"17-17/17-17",
                            "ReferenceRegion": f"{chrom}:{start_1based - 1}-{end_1based}",
                            "RepeatUnit": repeat_unit,
                            "VariantId": locus_id,
                            "VariantType": "Repeat",
                            "Ref": fields[3],
                            "Alt": fields[4],
                        }
                    }
                }
            except Exception as e:
                print(f"Error on vcf record #{line_counter}: {e}")
                print(line)
                print(genotype_dict)

    return locus_results


if __name__ == "__main__":
    main()
