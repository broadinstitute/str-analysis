"""This script converts an EnsembleTR consensus output VCF to the ExpansionHunter output .json format which is used as
a common input format by downstream str-truth-set-v2 comparison scripts (combine_str_json_to_tsv.py, etc.).

EnsembleTR (https://github.com/gymrek-lab/EnsembleTR) merges per-caller TR genotypes into a single consensus call. Its
output VCF:
  - drops the per-locus ID (ID column is "."),
  - rewrites INFO RU to a *canonical* rotation of the repeat unit (which generally differs from the truth-set motif),
  - reports the consensus diploid genotype as FORMAT NCOPY (comma-separated repeat copy numbers per allele),
  - reports a consensus FORMAT SCORE.

Because the downstream comparison merges tool results to the truth set on (LocusId, Motif, MotifSize), this script
reconciles each EnsembleTR record back to the truth-set locus by genomic coordinate using the ExpansionHunter variant
catalog (--variant-catalog) that the upstream callers were run against. The reconciled truth-set LocusId, motif, and
reference region are written to the json so that combine_str_json_to_tsv.py reproduces the exact merge keys. The repeat
COUNT is rotation-invariant (count = allele_bp / period and the motif length is unchanged), so the NCOPY copy numbers
remain valid under the truth-set motif.

The EnsembleTR SCORE is stored as the json "Q" field; run combine_str_json_to_tsv.py with
--include-extra-ensembletr-fields to surface it as the tool's Q column.
"""

import argparse
import bisect
import collections
import gzip
import simplejson as json
import re


def load_catalog_loci(variant_catalog_path):
    """Load an ExpansionHunter variant catalog (a json list of {LocusId, ReferenceRegion, ...}) and return:
      - exact_lookup: dict (chrom_no_chr, start_0based, end_1based) -> (locus_id, motif)
      - per_chrom_intervals: dict chrom_no_chr -> sorted list of (start_0based, end_1based, locus_id, motif)
    The motif is parsed from the trailing field of the LocusId ('{chrom}-{start0}-{end1}-{motif}').
    """
    fopen = gzip.open if variant_catalog_path.endswith("gz") else open
    with fopen(variant_catalog_path, "rt") as f:
        catalog = json.load(f)

    exact_lookup = {}
    per_chrom_intervals = collections.defaultdict(list)
    for record in catalog:
        locus_id = record["LocusId"]
        reference_region = record["ReferenceRegion"]
        # multi-locus structures store a list of regions; those aren't part of the single-locus truth set, so skip them
        if isinstance(reference_region, list):
            continue
        chrom, coords = reference_region.split(":")
        chrom = chrom[3:] if chrom.startswith("chr") else chrom
        start_0based, end_1based = (int(x) for x in coords.split("-"))
        motif = locus_id.rsplit("-", 1)[-1]
        exact_lookup[(chrom, start_0based, end_1based)] = (locus_id, motif)
        per_chrom_intervals[chrom].append((start_0based, end_1based, locus_id, motif))

    for chrom in per_chrom_intervals:
        per_chrom_intervals[chrom].sort()
    return exact_lookup, per_chrom_intervals


def reconcile_locus(chrom, start_0based, end_1based, period, exact_lookup, per_chrom_intervals):
    """Return (locus_id, motif, start_0based, end_1based) of the truth-set locus matching the EnsembleTR record, or None
    if no catalog locus overlaps. Tries an exact coordinate match first, then the catalog locus with the largest overlap
    whose motif length equals the EnsembleTR PERIOD."""
    exact = exact_lookup.get((chrom, start_0based, end_1based))
    if exact is not None:
        return exact[0], exact[1], start_0based, end_1based

    intervals = per_chrom_intervals.get(chrom)
    if not intervals:
        return None
    # candidate window: catalog loci can start anywhere before this record's end; scan from the first locus whose start
    # is within a motif-period of this record's start, walking until catalog starts move past this record's end
    lo = bisect.bisect_left(intervals, (start_0based - period - 1,))
    best = None
    best_overlap = 0
    for i in range(max(0, lo), len(intervals)):
        cat_start, cat_end, cat_locus_id, cat_motif = intervals[i]
        if cat_start >= end_1based:
            break
        if len(cat_motif) != period:
            continue
        overlap = min(end_1based, cat_end) - max(start_0based, cat_start)
        if overlap > best_overlap:
            best_overlap = overlap
            best = (cat_locus_id, cat_motif, cat_start, cat_end)
    return best


def process_ensembletr_vcf(vcf_path, exact_lookup, per_chrom_intervals, sample_id=None):
    locus_results = {
        "LocusResults": {},
        "SampleParameters": {"SampleId": sample_id, "Sex": None},
    }
    counters = collections.Counter()

    fopen = gzip.open if vcf_path.endswith("gz") else open
    with fopen(vcf_path, "rt") as vcf:
        sample_column_index = 9
        for line in vcf:
            if line.startswith("#CHROM"):
                header_fields = line.rstrip("\n").split("\t")
                if sample_id is None and len(header_fields) > 9:
                    locus_results["SampleParameters"]["SampleId"] = header_fields[9]
                continue
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            chrom_no_chr = chrom[3:] if chrom.startswith("chr") else chrom
            info_dict = dict(kv.split("=", 1) for kv in fields[7].split(";") if "=" in kv)
            format_keys = fields[8].split(":")
            sample_values = fields[sample_column_index].split(":")
            sample_dict = dict(zip(format_keys, sample_values))

            counters["total records"] += 1
            ncopy = sample_dict.get("NCOPY", ".")
            genotype_field = sample_dict.get("GT", ".")
            if ncopy in (".", "") or genotype_field in (".", "./.", ".|."):
                counters["no-call records skipped"] += 1
                continue

            try:
                period = int(info_dict["PERIOD"])
                start_1based = int(info_dict["START"])
                end_1based = int(info_dict["END"])
            except (KeyError, ValueError):
                counters["records missing START/END/PERIOD"] += 1
                continue
            start_0based = start_1based - 1
            end_0based_exclusive = end_1based - 1  # EnsembleTR END is 1-based-inclusive-of-the-last-ref-base + 1

            match = reconcile_locus(chrom_no_chr, start_0based, end_0based_exclusive, period,
                                    exact_lookup, per_chrom_intervals)
            if match is None:
                counters["records with no catalog locus (skipped)"] += 1
                continue
            locus_id, motif, ref_start_0based, ref_end_1based = match
            counters["records reconciled to a truth-set locus"] += 1

            allele_counts = sorted(int(round(float(c))) for c in ncopy.split(",") if c not in (".", ""))
            if not allele_counts:
                counters["no-call records skipped"] += 1
                continue
            genotype = "/".join(str(c) for c in allele_counts)
            genotype_ci = "/".join(f"{c}-{c}" for c in allele_counts)  # EnsembleTR has no per-allele CI; use zero width

            variant_record = {
                "Genotype": genotype,
                "GenotypeConfidenceInterval": genotype_ci,
                "ReferenceRegion": f"chr{chrom_no_chr}:{ref_start_0based}-{ref_end_1based}",
                "RepeatUnit": motif,
                "VariantId": locus_id,
                "Q": float(sample_dict.get("SCORE", "nan")),
            }
            locus_results["LocusResults"][locus_id] = {
                "AlleleCount": len(allele_counts),
                "LocusId": locus_id,
                "Variants": {locus_id: variant_record},
            }

    return locus_results, counters


def main():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-c", "--variant-catalog", required=True, help="ExpansionHunter variant catalog json (the catalog "
                   "the upstream callers were run against) used to reconcile EnsembleTR records to truth-set LocusIds.")
    p.add_argument("--sample-id", help="Sample id. If not specified, it is taken from the EnsembleTR VCF #CHROM line.")
    p.add_argument("-o", "--output-json", help="Output json path. Defaults to the input VCF path with the .vcf[.gz] "
                   "suffix replaced by .json.")
    p.add_argument("vcf_path", help="EnsembleTR consensus output VCF (or VCF.gz).")
    args = p.parse_args()

    print(f"Loading variant catalog {args.variant_catalog}")
    exact_lookup, per_chrom_intervals = load_catalog_loci(args.variant_catalog)
    print(f"Loaded {len(exact_lookup):,d} catalog loci")

    print(f"Processing {args.vcf_path}")
    locus_results, counters = process_ensembletr_vcf(
        args.vcf_path, exact_lookup, per_chrom_intervals, sample_id=args.sample_id)

    for key, count in counters.most_common():
        print(f"  {count:,d}  {key}")

    output_json_path = args.output_json or (re.sub(r".vcf(.gz)?$", "", args.vcf_path) + ".json")
    print(f"Writing {len(locus_results['LocusResults']):,d} loci to {output_json_path}")
    with open(output_json_path, "wt") as f:
        json.dump(locus_results, f, indent=3, ignore_nan=True)


if __name__ == "__main__":
    main()
