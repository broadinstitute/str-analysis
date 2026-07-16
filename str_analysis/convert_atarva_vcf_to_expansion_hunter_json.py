"""This script converts an ATaRVa (https://github.com/SowpatiLab/ATaRVa) `genotype` output VCF to the
ExpansionHunter output .json format which is used as a common input format in downstream scripts.

ATaRVa writes a single-sample VCF whose FORMAT column is always

    GT:AL:CN:LPM:AR:SD:DP:SN:SQ:MA:MR:DS:MV

and whose INFO column carries the repeat coordinates and motif:

    MOTIF   the repeat motif
    START   start position of the repeat region (0-based)
    END     end position of the repeat region (1-based)
    REFCN   reference-allele copy number ((END - START) // len(MOTIF))

Per-allele genotype is read from the FORMAT fields:

    AL   allele length in base pairs, one value per allele
    CN   motif copy number for each allele (== AL // len(MOTIF))

This converter uses CN as the per-allele repeat count, matching the integer repeat-count Genotype that the other
tool converters (vamos, TRGT, inquiSTR, ...) emit into the ExpansionHunter json. The truth-set locus id
"{chrom}-{start_0based}-{end_1based}-{motif}" is reconstructed from CHROM, INFO START, INFO END and INFO MOTIF so
the result merges against the truth set the same way the other tools do (ATaRVa is fed the single-motif truth-set
loci bed, so each locus has exactly one motif).

Records that ATaRVa could not genotype (FILTER=LESS_READS, whose SAMPLE column is all ".") are skipped, as are
"multizygous" loci where ATaRVa reported more than two alleles (which can't be expressed as a diploid genotype and
would not match a diploid truth-set locus anyway).
"""

"""
ATaRVa genotype output vcf example (FORMAT is GT:AL:CN:LPM:AR:SD:DP:SN:SQ:MA:MR:DS:MV; 1 homozygous, 1 heterozygous,
and 1 failed (LESS_READS) data line shown):

#CHROM  POS  ID  REF   ALT   QUAL  FILTER      INFO                                                      FORMAT             HG002
chr1    11227 .  NNN.. AGG.. 0     PASS        AC=2;AN=2;MOTIF=AGG;START=11226;END=11462;ID=.;REFCN=78   GT:AL:CN:...:MV    1/1:234,234:78,78:...
chr1    15797 .  CTT.. CTT.. 0     PASS        AC=1;AN=2;MOTIF=CTT;START=15796;END=15849;ID=.;REFCN=17   GT:AL:CN:...:MV    0|1:53,51:17,17:...
chr2    20000 .  AT..  .     0     LESS_READS  AC=0;AN=0;MOTIF=AT;START=19998;END=20040;ID=.;REFCN=21    GT:AL:CN:...:MV    .:.:.:...:.
"""

"""
ExpansionHunter output format:

  "LocusResults": {
        "chr12-57610122-57610131-GCA": {
          "AlleleCount": 2,
          "LocusId": "chr12-57610122-57610131-GCA",
          "Variants": {
            "chr12-57610122-57610131-GCA": {
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
        "SampleId": "HG002",
        "Sex": null
  }
"""


import argparse
import gzip
import re
import simplejson as json
from tqdm import tqdm


def parse_info(info):
    """Parse a VCF INFO column string into a dict. Flag fields (no '=') map to True."""
    info_dict = {}
    for entry in info.split(";"):
        if not entry:
            continue
        if "=" in entry:
            key, value = entry.split("=", 1)
            info_dict[key] = value
        else:
            info_dict[entry] = True
    return info_dict


def parse_sample(format_field, sample_field):
    """Zip a VCF FORMAT column string with its SAMPLE column into a dict of {field: value}."""
    return dict(zip(format_field.split(":"), sample_field.split(":")))


def main():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--discard-hom-ref", action="store_true", help="Discard loci where both alleles equal the number "
                   "of repeats in the reference.")
    p.add_argument("--sample-id", help="If not specified, the sample id is parsed from the last column of the VCF "
                   "header.")
    p.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    p.add_argument("--verbose", action="store_true", help="Print extra info about discarded or skipped loci")
    p.add_argument("atarva_vcf_path", help="ATaRVa genotype output VCF path (.vcf or .vcf.gz)")
    args = p.parse_args()

    print(f"Processing {args.atarva_vcf_path}")
    locus_results = process_atarva_vcf(
        args.atarva_vcf_path,
        sample_id=args.sample_id,
        discard_hom_ref=args.discard_hom_ref,
        show_progress_bar=args.show_progress_bar,
        verbose=args.verbose,
    )

    output_json_path = re.sub(".vcf(.gz)?$", "", args.atarva_vcf_path) + ".json"
    print(f"Writing {len(locus_results['LocusResults']):,d} loci to {output_json_path}")
    with open(output_json_path, "wt") as f:
        json.dump(locus_results, f, indent=3, ignore_nan=True)


def process_atarva_vcf(atarva_vcf_path, sample_id=None, discard_hom_ref=True, show_progress_bar=False, verbose=False):
    locus_results = {
        "LocusResults": {},
        "SampleParameters": {
            "SampleId": sample_id,
            "Sex": None,
        },
    }

    fopen = gzip.open if atarva_vcf_path.endswith("gz") else open

    no_genotype_count = skipped_multizygous_count = discarded_hom_ref_count = 0
    with fopen(atarva_vcf_path, "rt") as vcf:
        if show_progress_bar:
            vcf = tqdm(vcf, unit=" vcf records", unit_scale=True, unit_divisor=1000)

        line_counter = 0
        for line in vcf:
            if line.startswith("#CHROM"):
                header_fields = line.rstrip("\n").split("\t")
                if sample_id is None and len(header_fields) >= 10:
                    locus_results["SampleParameters"]["SampleId"] = header_fields[9]
                    print(f"Got sample id '{header_fields[9]}' from the VCF header")
                continue
            if line.startswith("#"):
                continue

            line_counter += 1
            try:
                fields = line.rstrip("\n").split("\t")
                chrom = fields[0]
                info_dict = parse_info(fields[7])
                sample_dict = parse_sample(fields[8], fields[9])

                # ATaRVa echoes the repeat coordinates in INFO (START is 0-based, END is 1-based). Prefer these over
                # POS so the locus id matches the truth-set catalog exactly.
                start_0based = int(info_dict["START"])
                end_1based = int(info_dict["END"])
                motif = info_dict["MOTIF"]
                motif_size = len(motif)

                # CN gives the per-allele motif copy number. Failed loci (FILTER=LESS_READS) have CN "." -> no-call.
                cn_field = sample_dict.get("CN", ".")
                allele_sizes = [int(cn) for cn in cn_field.split(",") if cn not in (".", "")]
                if len(allele_sizes) == 0:
                    no_genotype_count += 1
                    continue
                if len(allele_sizes) > 2:
                    # ATaRVa reported >2 alleles (multizygous); can't be expressed as a diploid genotype.
                    skipped_multizygous_count += 1
                    continue

                allele_sizes.sort()
                num_repeats_in_reference = (end_1based - start_0based) // motif_size
                if discard_hom_ref and all(a == num_repeats_in_reference for a in allele_sizes):
                    discarded_hom_ref_count += 1
                    continue

                locus_id = f"{chrom}-{start_0based}-{end_1based}-{motif}"
                locus_results["LocusResults"][locus_id] = {
                    "AlleleCount": len(allele_sizes),
                    "LocusId": locus_id,
                    "Variants": {
                        locus_id: {
                            "Genotype": "/".join(str(a) for a in allele_sizes),
                            "GenotypeConfidenceInterval": "/".join(f"{a}-{a}" for a in allele_sizes),
                            "ReferenceRegion": f"{chrom}:{start_0based}-{end_1based}",
                            "RepeatUnit": motif,
                            "VariantId": locus_id,
                            "VariantType": "Repeat",
                        }
                    }
                }
            except Exception as e:
                print(f"Error while parsing ATaRVa vcf record #{line_counter}: {e}")
                print(f"    {line}")

    if verbose or no_genotype_count or skipped_multizygous_count or discarded_hom_ref_count:
        print(f"Parsed {len(locus_results['LocusResults']):,d} loci, skipped {no_genotype_count:,d} loci without a "
              f"genotype (no-call), {skipped_multizygous_count:,d} multizygous loci, and discarded "
              f"{discarded_hom_ref_count:,d} hom-ref loci")

    return locus_results


if __name__ == "__main__":
    main()
