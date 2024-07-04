"""This script converts a Straglr BED output file to the ExpansionHunter output .json format
which I use as a common input format in downstream scripts.
"""


import argparse
import gzip
import simplejson as json
import re
from tqdm import tqdm

def main():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--verbose", action="store_true", help="Print verbose output")
    p.add_argument("--sample-id", help="If not specified, the sample id will be based on the filename.")
    p.add_argument("bed_path", help="TRGT vcf path")
    args = p.parse_args()

    print(f"Processing {args.bed_path}")
    locus_results = process_straglr_bed(
        args.bed_path,
        sample_id=args.sample_id,
        verbose=args.verbose,
    )

    output_json_path = re.sub(".bed(.gz)?$", "", args.bed_path) + ".json"
    print(f"Writing {len(locus_results['LocusResults']):,d} loci to {output_json_path}")
    with open(output_json_path, "wt") as f:
        json.dump(locus_results, f, indent=3, ignore_nan=True)


def process_straglr_bed(bed_path, sample_id=None, verbose=False):
    locus_results = {
        "LocusResults": {},
        "SampleParameters": {
            "SampleId": sample_id,
            "Sex": None,
        },
    }

    fopen = gzip.open if bed_path.endswith("gz") else open

    with fopen(bed_path, "rt") as bed_file:
        # #chrom, start, end, repeat_unit, allele1:size,  allele1:copy_number, allele1:support allele2:size, allele2:copy_number, allele2:support
        header_fields = bed_file.readline().strip().strip("#").split("\t")

        if verbose:
            bed_file = tqdm(bed_file, unit=" BED lines", unit_scale=True, unit_divisor=1000)

        line_counter = 0
        for line in bed_file:
            line_counter += 1
            fields = line.strip().split("\t")

            record = dict(zip(header_fields, fields))
            chrom = record["chrom"]
            start_0based = int(record["start"])
            end_1based = int(record["end"])
            repeat_unit = record["repeat_unit"]
            reference_allele_num_repeats = (end_1based - start_0based)//len(repeat_unit)

            locus_id = f"{chrom}-{start_0based}-{end_1based}-{repeat_unit}"
            if record["allele1:copy_number"] == "-":
                # skip missing genotypes
                continue

            allele1_num_repeats = int(float(record["allele1:copy_number"]))
            allele2_num_repeats = int(float(record["allele2:copy_number"])) if record["allele2:copy_number"] != "-" else None
            if int(allele1_num_repeats) == reference_allele_num_repeats and allele2_num_repeats is None:
                # skip hom-ref
                continue

            allele1_read_support = int(record["allele1:support"])
            allele2_read_support = int(record["allele2:support"]) if record["allele2:support"] != "-" else 0

            short_allele = allele1_num_repeats if allele2_num_repeats is None else min(allele1_num_repeats, allele2_num_repeats)
            long_allele = allele1_num_repeats if allele2_num_repeats is None else max(allele1_num_repeats, allele2_num_repeats)
            genotype = f"{short_allele}/{long_allele}"
            genotype_CI = f"{short_allele}-{short_allele}/{long_allele}-{long_allele}"

            locus_results["LocusResults"][locus_id] = {
                "AlleleCount": 1 if allele2_num_repeats is None else 2,
                "LocusId": locus_id,
                "Coverage": allele1_read_support + allele2_read_support,

                #"ReadLength": None,
                #"FragmentLength": None,
                "Variants": {
                    locus_id:{
                        "Genotype": genotype,   #"17/17",
                        "GenotypeConfidenceInterval": genotype_CI, #"17-17/17-17",
                        "ReferenceRegion": f"{chrom}:{start_0based}-{end_1based}",
                        "RepeatUnit": repeat_unit,
                        "VariantId": locus_id,
                        "VariantType": "Repeat",
                    }
                }
            }

    return locus_results


if __name__ == "__main__":
    main()
