"""This script converts a HipSTR output VCF to the .json format ExpansionHunter uses to output results.
This makes it easier to pass HipSTR results to downstream scripts.
"""

"""
HipSTR output vcf format:

chrX    
155292771       
chrX-155292770-155292798-AAAC   
AAACAAACAAACAAACAAACAAACAAAC    
AAACAAACAAACAAACAAACAAACAAACAAACAAAC    
.       
.       
INFRAME_PGEOM=0.95;INFRAME_UP=0.05;INFRAME_DOWN=0.05;OUTFRAME_PGEOM=0.95;OUTFRAME_UP=0.01;OUTFRAME_DOWN=0.01;START=155292771;END=155292798;PERIOD=4;NSKIP=0;NFILT=0;BPDIFFS=8;DP=57;DSNP=0;DSTUTTER=0;DFLANKINDEL=3;AN=2;REFAC=1;AC=1 
GT:GB:Q:PQ:DP:DSNP:DSTUTTER:DFLANKINDEL:PDP:PSNP:GLDIFF:AB:DAB:FS:ALLREADS:MALLREADS    
0|1:0|8:1.00:0.50:57:0:0:3:28.96|28.04:0|0:33.39:-0.00:37:-0.13:0|16;8|12:0|14;8|13
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
import re


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sample-id",
                   help="If not specified, the sample id will be parsed from the last column of the vcf header.")
    p.add_argument("vcf_path", help="HipSTR vcf path(s)")
    args = p.parse_args()

    print(f"Processing {args.vcf_path}")
    locus_results = process_hipstr_vcf(args.vcf_path, sample_id=args.sample_id)

    output_json_path = re.sub(".vcf(.gz)?$", "", args.vcf_path) + ".json"
    print(f"Writing results for {len(locus_results['LocusResults']):,d} loci to {output_json_path}")
    with open(output_json_path, "wt") as f:
        json.dump(locus_results, f, indent=4)


def process_hipstr_vcf(vcf_path, sample_id=None):
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
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    header_fields = line.strip().split("\t")
                    if sample_id is None and len(header_fields) > 9:
                        print(f"Got sample id '{header_fields[9]}' from the VCF header")
                        locus_results["SampleParameters"]["SampleId"] = header_fields[9]

                continue

            line_counter += 1
            fields = line.strip().split("\t")
            if not fields[9] or fields[9] == ".":
                continue  # no genotype

            genotype_fields = fields[8].split(":")
            genotype_values = fields[9].split(":")
            genotype_dict = dict(zip(genotype_fields, genotype_values))

            genotype = genotype_dict["GT"].replace("|", "/").replace("\\", "/")
            if genotype == ".":
                continue   # no genotype

            chrom = fields[0]
            start_1based = int(fields[1])
            locus_id = fields[2]
            _, _, _, repeat_unit = locus_id.split("-")
            ref = fields[3]
            alts = fields[4].split(",")
            alleles = [ref] + alts

            info = fields[7]
            info_dict = dict([key_value.split("=") for key_value in info.split(";")])
            for key, value in info_dict.items():
                try:
                    info_dict[key] = float(value.strip('"'))
                except ValueError:
                    continue

            for key, value in genotype_dict.items():
                try:
                    genotype_dict[key] = float(value.strip('"'))
                except ValueError:
                    continue
            end_1based = int(info_dict["END"])

            try:
                left_allele, right_allele = map(int, genotype.split("/"))
                allele1 = alleles[left_allele]
                allele2 = alleles[right_allele]
                if len(allele1) < len(allele2):
                    short_allele, long_allele = allele1, allele2
                else:
                    short_allele, long_allele = allele2, allele1

                num_repeats1 = int(len(short_allele)/len(repeat_unit))
                num_repeats2 = int(len(long_allele)/len(repeat_unit))

                locus_results["LocusResults"][locus_id] = {
                    "AlleleCount": 2,
                    "LocusId": locus_id,
                    "Coverage": float(genotype_dict["DP"]),
                    "ReadLength": None,
                    "FragmentLength": None,
                    "Variants": {
                        locus_id: info_dict | genotype_dict | {
                            "Genotype": f"{num_repeats1}/{num_repeats2}",
                            "GenotypeConfidenceInterval": f"{num_repeats1}-{num_repeats1}/{num_repeats2}-{num_repeats2}",
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
