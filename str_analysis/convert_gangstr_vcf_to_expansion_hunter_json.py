"""This script converts a GangSTR output VCF to the .json format ExpansionHunter uses to output results.
This makes it easier to pass GangSTR results to downstream scripts.
"""

"""
GangSTR output vcf format:

chr14   
39042968        
.       
tttatttatttatttattta    
.       
.       
.       
END=39042990;RU=ttta;PERIOD=4;REF=5;GRID=2,8;STUTTERUP=0.05;STUTTERDOWN=0.05;STUTTERP=0.9;EXPTHRESH=-1  
GT:DP:Q:REPCN:REPCI:RC:ENCLREADS:FLNKREADS:ML:INS:STDERR:QEXP   
0/0:18:0.893606:5,5:5-5,5-5:7,11,0,0:5,7:NULL:124.576:570.176,153.445:0,0:-1,-1,-1
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
    p.add_argument("--variant-catalog", help="ExpansionHunter variant catalog. If specified, fields from this"
                                             " variant catalog will be added to the output .json.")
    p.add_argument("vcf_path", help="GangSTR vcf path")
    args = p.parse_args()

    variant_catalog = None
    if args.variant_catalog:
        if not os.path.isfile(args.variant_catalog):
            p.error(f"{args.variant_catalog} not found")

        with open(args.variant_catalog) as f:
            variant_catalog = json.load(f)

    print(f"Processing {args.vcf_path}")
    locus_results = process_gangstr_vcf(args.vcf_path, variant_catalog=variant_catalog, sample_id=args.sample_id)

    output_json_path = re.sub(".vcf(.gz)?$", "", args.vcf_path) + ".json"
    print(f"Writing results for {len(locus_results['LocusResults']):,d} loci to {output_json_path}")
    with open(output_json_path, "wt") as f:
        json.dump(locus_results, f, indent=4)


def create_variant_catalog_lookup(variant_catalog):
    variant_catalog_lookup = {}
    for record in variant_catalog:
        reference_regions = record["ReferenceRegion"]
        variant_types = record["VariantType"]
        variant_ids = record.get("VariantId", record["LocusId"])
        if not isinstance(reference_regions, list):
            reference_regions = [reference_regions]
            variant_types = [variant_types]
            variant_ids = [variant_ids]
        for variant_id, variant_type, reference_region in zip(variant_ids, variant_types, reference_regions):
            variant_catalog_lookup[reference_region] = dict(record)
            variant_catalog_lookup[reference_region]["ReferenceRegion"] = reference_region
            variant_catalog_lookup[reference_region]["VariantType"] = variant_type
            variant_catalog_lookup[reference_region]["VariantId"] = variant_id

    return variant_catalog_lookup


def process_gangstr_vcf(vcf_path, variant_catalog=None, sample_id=None):
    locus_results = {
        "LocusResults": {},
        "SampleParameters": {
            "SampleId": sample_id,
            "Sex": None,
        },
    }

    variant_catalog_lookup = create_variant_catalog_lookup(variant_catalog or [])

    fopen = gzip.open if vcf_path.endswith("gz") else open
    with fopen(vcf_path, "rt") as vcf:
        line_counter = 0
        for line in vcf:
            if line.startswith("#"):
                if line.startswith("##command="):
                    # Parse the "--bam-samps RGP_1126_1 --samp-sex F" args from the command
                    sample_id_match = re.search("bam-samps ([^-]+)", line)
                    if sample_id_match:
                        locus_results["SampleParameters"]["SampleId"] = sample_id_match.group(1).strip()
                    sample_sex_match = re.search("samp-sex ([^- ]+)", line)
                    if sample_sex_match:
                        sample_sex = sample_sex_match.group(1).strip()
                        locus_results["SampleParameters"]["Sex"] = "Female" if sample_sex.upper().startswith("F") else ("Male" if sample_sex.upper().startswith("M") else None)
                elif line.startswith("#CHROM"):
                    header_fields = line.strip().split("\t")
                    if sample_id is None and len(header_fields) == 10:
                        print(f"Got sample id '{header_fields[9]}' from the VCF header")
                        locus_results["SampleParameters"]["SampleId"] = header_fields[9]

                continue

            line_counter += 1
            fields = line.strip().split("\t")
            chrom = fields[0]
            start_1based = int(fields[1])
            info = fields[7]
            if not fields[9] or fields[9] == ".":  # no genotype
                continue

            # example1: chr8	141029449	.	gtggtggtg	.	.	.
            #    END=141029457;RU=gtg;PERIOD=3;REF=3;GRID=1,6;STUTTERUP=0.05;STUTTERDOWN=0.05;STUTTERP=0.9;EXPTHRESH=-1
            #    GT:DP:Q:REPCN:REPCI:RC:ENCLREADS:FLNKREADS:ML:INS:STDERR:QEXP
            #    0/0:47:1:3,3:3-3,3-3:38,9,0,0:3,38:NULL:232.623:342,114:0,0:-1,-1,-1
            #
            # example2: chr8	140697417	.	tgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtg	tgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtg	.	.
            #   END=140697454;RU=tg;PERIOD=2;REF=19;GRID=13,23;STUTTERUP=0.05;STUTTERDOWN=0.05;STUTTERP=0.9;EXPTHRESH=-1
            #   GT:DP:Q:REPCN:REPCI:RC:ENCLREADS:FLNKREADS:ML:INS:STDERR:QEXP
            #   0/1:41:0.999999:19,20:19-19,19-20:24,7,0,9:19,14|20,10:6,1|8,1|9,2|11,1|13,1|14,1|15,1|17,1:247.512:342,114:0.29703,0.564617:-1,-1,-1
            #
            # example3: chr8	140594826	.	tgtgtgtgtgtgtgtgtgtgtgtgtgtg	tgtgtgtgtgtgtgtgtgtgtgtgtgtgtg,tgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtg	.	.
            #   END=140594853;RU=tg;PERIOD=2;REF=14;GRID=12,26;STUTTERUP=0.05;STUTTERDOWN=0.05;STUTTERP=0.9;EXPTHRESH=-1
            #   GT:DP:Q:REPCN:REPCI:RC:ENCLREADS:FLNKREADS:ML:INS:STDERR:QEXP
            #   1/2:36:0.994974:15,23:15-15,20-23:21,7,0,8:15,19|23,2:5,1|8,1|9,1|10,1|11,1|14,1|20,2:227.086:342,114:0,1.56141:-1,-1,-1

            info_dict = dict([key_value.split("=") for key_value in info.split(";")])
            genotype_fields = fields[8].split(":")
            genotype_values = fields[9].split(":")
            genotype_dict = dict(zip(genotype_fields, genotype_values))

            try:
                repeat_unit = info_dict["RU"].upper()
                end_1based = int(info_dict["END"])
                ref_repeat_count = int(info_dict["REF"])

                variant_catalog_record = variant_catalog_lookup.get(f"{chrom}:{start_1based - 1}-{end_1based}", {})
                locus_id = variant_catalog_record.get("LocusId", f"{chrom}-{start_1based - 1}-{end_1based}-{repeat_unit}")
                variant_id = variant_catalog_record.get("VariantId", locus_id)
                if locus_id != variant_id:
                    print(f"Skipping adjacent locus {locus_id}, variant_id: {variant_id}")
                    continue

                locus_results["LocusResults"][locus_id] = {
                    "AlleleCount": genotype_dict["REPCN"].count(",") + 1,
                    "LocusId": locus_id,
                    "Coverage": float(genotype_dict["DP"]),  #10.757737459978655,
                    "ReadLength": None,
                    "FragmentLength": None,
                    "Variants": {
                        locus_id: info_dict | genotype_dict | {
                            "Genotype": genotype_dict["REPCN"].replace(",", "/"), #"17/17",
                            "GenotypeConfidenceInterval": genotype_dict["REPCI"].replace(",", "/"), #"17-17/17-17",
                            "ReferenceRegion": f"{chrom}:{start_1based - 1}-{end_1based}",
                            "RepeatUnit": repeat_unit,
                            "VariantId": variant_id,
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
