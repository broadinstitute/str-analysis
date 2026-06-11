"""This script converts a vamos (https://github.com/ChaissonLab/vamos) `--read` mode output VCF to the
ExpansionHunter output .json format which is used as a common input format in downstream scripts.

vamos `--read` writes a single-sample diploid VCF with a symbolic <VNTR> ALT allele. Each record carries these
INFO fields:

    END         end position of the locus (1-based, inclusive) echoed from the vamos catalog
    RU          comma-separated list of motifs ("repeat units") used to annotate the locus
    SVTYPE      always "VNTR"
    ALTANNO_H1  comma-separated indices into RU giving the motif composition of haplotype 1
    LEN_H1      total count of motifs in the haplotype 1 allele (== len(ALTANNO_H1))
    ALTANNO_H2  same as ALTANNO_H1 for haplotype 2 (present only when the genotype is heterozygous, GT=1/2)
    LEN_H2      total count of motifs in the haplotype 2 allele

The vamos catalog this pipeline feeds it has exactly one motif per locus (it is generated from the truth set's
single-motif ExpansionHunter catalog), so RU has a single motif M and the per-haplotype repeat count is simply the
motif count LEN_Hx. This converter therefore sets RepeatUnit = M and the genotype to the LEN_H1/LEN_H2 motif counts.
Loci whose RU contains more than one motif are true VNTRs that can't be expressed as an integer repeat count of a
single motif, so they are skipped (they would not match a truth set STR locus anyway).

vamos is run in its default 1-based mode (no -Z flag), so VCF POS equals the catalog start (1-based) and the 0-based
start is POS - 1 (see vamos src/vcf.cpp: POS = ref_start + 1 - oneOffset, with oneOffset defaulting to 1). The truth
set locus id "{chrom}-{start_0based}-{end}-{motif}" is reconstructed from chrom, POS, END and M so the result merges
against the truth set the same way the other tools do.
"""

"""
vamos --read output vcf example (1 header line shown, then 1 homozygous and 1 heterozygous data line):

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG002
chr1	11226	.	N	<VNTR>	.	PASS	END=11462;RU=AGG;SVTYPE=VNTR;ALTANNO_H1=0,0,4,0;LEN_H1=4;	GT	1/1
chr1	15796	.	N	<VNTR>	.	PASS	END=15849;RU=CTT;SVTYPE=VNTR;ALTANNO_H1=0,0,0;LEN_H1=3;ALTANNO_H2=0,0,0,0;LEN_H2=4;	GT	1/2
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


def haplotype_num_repeats(info_dict, suffix):
    """Return the number of motif copies for one haplotype of a single-motif locus, or None if it is absent.

    Args:
        info_dict (dict): parsed VCF INFO fields for the locus.
        suffix (str): the haplotype suffix, "H1" or "H2".

    Returns:
        int or None: the motif count (LEN_<suffix>), or None when this haplotype has no annotation (which for vamos
            `--read` means the genotype is homozygous and only H1 is reported).
    """
    if f"LEN_{suffix}" in info_dict:
        return int(info_dict[f"LEN_{suffix}"])
    # Fall back to counting the ALTANNO indices if an older vamos build didn't emit LEN. vamos separates the indices
    # with "," in --read mode and "-" in --contig mode.
    altanno = info_dict.get(f"ALTANNO_{suffix}")
    if altanno is None:
        return None
    return len([x for x in re.split("[,-]", altanno.strip()) if x != ""])


def main():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--discard-hom-ref", action="store_true", help="Discard loci where both alleles equal the number "
                   "of repeats in the reference.")
    p.add_argument("--sample-id", help="If not specified, the sample id is parsed from the last column of the VCF "
                   "header.")
    p.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    p.add_argument("--verbose", action="store_true", help="Print extra info about discarded or skipped loci")
    p.add_argument("vamos_vcf_path", help="vamos --read output VCF path (.vcf or .vcf.gz)")
    args = p.parse_args()

    print(f"Processing {args.vamos_vcf_path}")
    locus_results = process_vamos_vcf(
        args.vamos_vcf_path,
        sample_id=args.sample_id,
        discard_hom_ref=args.discard_hom_ref,
        show_progress_bar=args.show_progress_bar,
        verbose=args.verbose,
    )

    output_json_path = re.sub(".vcf(.gz)?$", "", args.vamos_vcf_path) + ".json"
    print(f"Writing {len(locus_results['LocusResults']):,d} loci to {output_json_path}")
    with open(output_json_path, "wt") as f:
        json.dump(locus_results, f, indent=3, ignore_nan=True)


def process_vamos_vcf(vamos_vcf_path, sample_id=None, discard_hom_ref=True, show_progress_bar=False, verbose=False):
    locus_results = {
        "LocusResults": {},
        "SampleParameters": {
            "SampleId": sample_id,
            "Sex": None,
        },
    }

    fopen = gzip.open if vamos_vcf_path.endswith("gz") else open

    skipped_multi_motif_count = no_genotype_count = discarded_hom_ref_count = 0
    with fopen(vamos_vcf_path, "rt") as vcf:
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
                start_0based = int(fields[1]) - 1  # vamos runs in default 1-based mode (no -Z), so POS == catalog start
                info_dict = parse_info(fields[7])
                end_1based = int(info_dict["END"])

                motifs = info_dict["RU"].split(",")
                if len(motifs) != 1:
                    skipped_multi_motif_count += 1
                    continue
                motif = motifs[0]
                motif_size = len(motif)

                # ALTANNO_H2/LEN_H2 are absent for a homozygous call (GT=1/1), where vamos reports only H1.
                allele_sizes = [n for n in (haplotype_num_repeats(info_dict, "H1"),
                                            haplotype_num_repeats(info_dict, "H2")) if n is not None]
                if len(allele_sizes) == 0:
                    no_genotype_count += 1
                    continue
                if len(allele_sizes) == 1:
                    allele_sizes = [allele_sizes[0], allele_sizes[0]]

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
                print(f"Error while parsing vamos vcf record #{line_counter}: {e}")
                print(f"    {line}")

    if verbose or skipped_multi_motif_count or no_genotype_count or discarded_hom_ref_count:
        print(f"Parsed {len(locus_results['LocusResults']):,d} loci, skipped {skipped_multi_motif_count:,d} "
              f"multi-motif loci, {no_genotype_count:,d} loci without a genotype, and {discarded_hom_ref_count:,d} "
              f"hom-ref loci")

    return locus_results


if __name__ == "__main__":
    main()
