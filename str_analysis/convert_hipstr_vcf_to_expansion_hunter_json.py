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
    p.add_argument("--skip-loci-with-invalid-locus-ids", action="store_true",
                   help="Locus Ids are expected to have the format 'chrom-start-end-motif'. Skip loci with ids"
                        "that don't match this format instead of throwing an error")
    p.add_argument("--skip-hom-ref-loci", action="store_true", help="Filter ou loci that were called as "
                                                                         "homozygous reference")

    p.add_argument("vcf_path", help="HipSTR vcf path(s)")
    args = p.parse_args()

    print(f"Processing {args.vcf_path}")
    locus_results = process_hipstr_vcf(
        args.vcf_path,
        sample_id=args.sample_id,
        skip_hom_ref_loci=args.skip_hom_ref_loci,
        skip_loci_with_invalid_locus_ids=args.skip_loci_with_invalid_locus_ids,
    )

    output_json_path = re.sub(".vcf(.gz)?$", "", args.vcf_path) + ".json"
    print(f"Writing {len(locus_results['LocusResults']):,d} loci to {output_json_path}")
    with open(output_json_path, "wt") as f:
        json.dump(locus_results, f, indent=4)


def process_hipstr_vcf(vcf_path, sample_id=None, skip_hom_ref_loci=False, skip_loci_with_invalid_locus_ids=False):
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
        locus_pos_adjust_counter = 0

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

            if skip_hom_ref_loci and genotype == "0/0":
                continue

            chrom = fields[0]
            pos = int(fields[1])
            locus_id = fields[2]
            locus_id_fields = locus_id.split("-")
            if len(locus_id_fields) >= 4 and len(set(locus_id_fields[3]) - set("ACGTRN")) == 0:
                repeat_unit = locus_id_fields[3]
            else:
                if skip_loci_with_invalid_locus_ids:
                    continue
                else:
                    raise ValueError(
                        f"VCF row #{line_counter}: Unable to parse repeat unit from locus id '{locus_id}'. "
                        f"Expected locus id to have the format 'chrom-start-end-motif'. Instead, found '{locus_id}'")

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

            start_1based = int(info_dict["START"])
            end_1based = int(info_dict["END"])
            period = int(info_dict["PERIOD"])

            try:
                left_allele, right_allele = map(int, genotype.split("/"))
                allele1 = alleles[left_allele]
                allele2 = alleles[right_allele]
                if len(allele1) < len(allele2):
                    short_allele, long_allele = allele1, allele2
                else:
                    short_allele, long_allele = allele2, allele1

                if pos < start_1based:
                    if start_1based - pos >= len(short_allele):
                        raise ValueError(f"start_1based ({start_1based:,d}) - pos ({pos:,d}) >= len(short_allele) ({len(short_allele):,d}). The short allele is '{short_allele}'")
                    if start_1based - pos >= len(long_allele):
                        raise ValueError(f"start_1based ({start_1based:,d}) - pos ({pos:,d}) >= len(long_allele) ({len(long_allele):,d}). The long allele is '{long_allele}'")
                    locus_pos_adjust_counter += 1
                    print(f"WARNING: {locus_id} VCF row has pos < start_1based: {pos:,d} < {start_1based:,d}. Diff: {start_1based - pos}bp. "
                          f"Trimming alleles: {short_allele} => {short_allele[start_1based - pos:]} and "
                          f"{long_allele} => {long_allele[start_1based - pos:]}")
                    short_allele = short_allele[start_1based - pos:]
                    long_allele = long_allele[start_1based - pos:]

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

        if locus_pos_adjust_counter:
            print(f"Adjusted the start position of {locus_pos_adjust_counter:,d} out of {line_counter:,d} ({locus_pos_adjust_counter/line_counter:.2%}) loci.")

    return locus_results


if __name__ == "__main__":
    main()
