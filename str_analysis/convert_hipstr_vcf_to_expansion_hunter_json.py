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
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--sample-id",
                   help="If not specified, the sample id will be parsed from the last column of the vcf header.")
    p.add_argument("--skip-hom-ref-loci", action="store_true", help="Filter ou loci that were called as "
                                                                         "homozygous reference")
    p.add_argument("--output-path", help="Output path for the json file. By default, it is based on the input filename.")
    p.add_argument("vcf_path", help="HipSTR vcf path(s)")
    args = p.parse_args()

    print(f"Processing {args.vcf_path}")
    locus_results = process_hipstr_vcf(
        args.vcf_path,
        sample_id=args.sample_id,
        skip_hom_ref_loci=args.skip_hom_ref_loci,
    )

    output_json_path = args.output_path or (re.sub(".vcf(.gz)?$", "", os.path.basename(args.vcf_path)) + ".json")
    print(f"Writing {len(locus_results['LocusResults']):,d} loci to {output_json_path}")
    with open(output_json_path, "wt") as f:
        json.dump(locus_results, f, indent=4, ignore_nan=True)


def trim_allele(allele, left_trim, right_trim):
    """Trim left_trim bases off the start of an allele and right_trim bases off its end.

    The same number of flanking bases is removed from every allele of a record, which leaves each allele's own
    length intact. An empty allele (a symbolic ALT, which spans zero bases of the locus) is already in trimmed form.

    Args:
        allele (str): the allele sequence, or "" for a symbolic allele.
        left_trim (int): bases to remove from the start.
        right_trim (int): bases to remove from the end.

    Returns:
        str: the trimmed allele.
    """
    if not allele:
        return ""
    return allele[left_trim: len(allele) - right_trim]


def process_hipstr_vcf(vcf_path, sample_id=None, skip_hom_ref_loci=False):
    locus_results = {
        "LocusResults": {},
        "SampleParameters": {
            "SampleId": sample_id,
            "Sex": None,
        },
    }

    fopen = gzip.open if vcf_path.endswith((".gz", ".bgz")) else open
    with fopen(vcf_path, "rt") as vcf:
        line_counter = 0
        locus_coordinates_adjusted_counter = 0
        locus_coordinates_adjusted_became_hom_ref_counter = 0
        record_does_not_span_locus_counter = 0
        allele_shorter_than_flanks_counter = 0
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

            ref = fields[3]
            # A symbolic ALT (LongTR writes <DEL>) is an allele that spans zero bases of the locus. Keeping the
            # literal text would measure it as len("<DEL>") // period repeats, so represent it as what it is.
            alts = ["" if alt.startswith("<") else alt for alt in fields[4].split(",")]
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

                # HipSTR and LongTR both extend the genotyped region past the requested bed interval, so the record
                # can span more than the locus. The flanking bases are shared by REF and every ALT, so the same
                # number of bases comes off each end of each allele, which leaves the allele's own length intact.
                # Slicing to a fixed width instead would truncate every expanded allele back to the reference
                # length and report it as homozygous reference.
                left_trim = start_1based - pos
                right_trim = (pos + len(ref) - 1) - end_1based
                if left_trim < 0 or right_trim < 0:
                    # the record doesn't span the whole locus interval, so its alleles can't be trimmed to it
                    record_does_not_span_locus_counter += 1
                    continue

                locus_coordinates_changed = left_trim > 0 or right_trim > 0
                if locus_coordinates_changed:
                    # Check both alleles, not just the shorter one: a symbolic allele is already "" and so always
                    # sorts first, which would otherwise leave a real allele shorter than the flanks unchecked and
                    # silently trimmed away to "" (0 repeats).
                    if any(allele and left_trim + right_trim > len(allele)
                           for allele in (short_allele, long_allele)):
                        # an allele shorter than the flanks being removed can't be trimmed to the locus interval
                        allele_shorter_than_flanks_counter += 1
                        continue

                    locus_coordinates_adjusted_counter += 1
                    print(f"WARNING: {locus_id} VCF row #{line_counter:,d} has " +
                          (f"pos < start_1based: {pos:,d} < {start_1based:,d}" if left_trim > 0 else "") +
                          (" and " if left_trim > 0 and right_trim > 0 else "") +
                          (f"pos + len(ref) - 1 > end_1based: {pos + len(ref) - 1:,d} > {end_1based:,d}" if right_trim > 0 else "") + " "
                          f"Diff: {left_trim}bp on the left, {right_trim}bp on the right. Trimming alleles: "
                          f"{short_allele} => {trim_allele(short_allele, left_trim, right_trim)} and "
                          f"{long_allele} => {trim_allele(long_allele, left_trim, right_trim)}")

                    short_allele = trim_allele(short_allele, left_trim, right_trim)
                    long_allele = trim_allele(long_allele, left_trim, right_trim)

                num_repeats_in_reference = int((end_1based - start_1based + 1) / period)
                num_repeats1 = int(len(short_allele)/period)
                num_repeats2 = int(len(long_allele)/period)
                if locus_coordinates_changed and num_repeats1 == num_repeats_in_reference and num_repeats2 == num_repeats_in_reference:
                    locus_coordinates_adjusted_became_hom_ref_counter += 1
                    if skip_hom_ref_loci:
                        continue

                repeat_unit_candidates = [long_allele[i:i+period] for i in range(0, len(long_allele), period)]
                repeat_unit_candidates += [short_allele[i:i+period] for i in range(0, len(short_allele), period)]
                if not repeat_unit_candidates:
                    # both alleles are symbolic, so the reference allele is the only sequence left to read it from
                    trimmed_ref = trim_allele(ref, left_trim, right_trim)
                    repeat_unit_candidates = [trimmed_ref[i:i+period] for i in range(0, len(trimmed_ref), period)]
                repeat_unit = max(repeat_unit_candidates, key=repeat_unit_candidates.count)

                locus_results["LocusResults"][locus_id] = {
                    "AlleleCount": 2,
                    "LocusId": locus_id,
                    "Coverage": float(genotype_dict["DP"]),
                    #"ReadLength": None,
                    #"FragmentLength": None,
                    "Variants": {
                        locus_id: {**info_dict, **genotype_dict, **{
                            "Genotype": f"{num_repeats1}/{num_repeats2}",
                            "GenotypeConfidenceInterval": f"{num_repeats1}-{num_repeats1}/{num_repeats2}-{num_repeats2}",
                            "ReferenceRegion": f"{chrom}:{start_1based - 1}-{end_1based}",
                            "RepeatUnit": repeat_unit,
                            "VariantId": locus_id,
                            "VariantType": "Repeat",
                            "Ref": fields[3],
                            "Alt": fields[4],
                        }}
                    }
                }
            except Exception as e:
                print(f"Error on vcf record #{line_counter}: {e}")
                print(line)
                print(genotype_dict)

        if locus_coordinates_adjusted_counter:
            print(f"Adjusted the start position of {locus_coordinates_adjusted_counter:,d} out of {line_counter:,d} ({locus_coordinates_adjusted_counter/line_counter:.2%}) loci.")
        if locus_coordinates_adjusted_became_hom_ref_counter:
            print(f"Adjusted loci that became homozygous reference after coordinate adjustment: {locus_coordinates_adjusted_became_hom_ref_counter:,d} out of {locus_coordinates_adjusted_counter:,d} ({locus_coordinates_adjusted_became_hom_ref_counter/locus_coordinates_adjusted_counter:.2%}) loci.")
        if record_does_not_span_locus_counter:
            print(f"Skipped {record_does_not_span_locus_counter:,d} out of {line_counter:,d} ({record_does_not_span_locus_counter/line_counter:.2%}) loci where the vcf record didn't span the whole locus interval.")
        if allele_shorter_than_flanks_counter:
            print(f"Skipped {allele_shorter_than_flanks_counter:,d} out of {line_counter:,d} ({allele_shorter_than_flanks_counter/line_counter:.2%}) loci with an allele shorter than the flanking bases being trimmed off.")

    return locus_results


if __name__ == "__main__":
    main()
