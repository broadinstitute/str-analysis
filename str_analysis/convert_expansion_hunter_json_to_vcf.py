"""This script converts an ExpansionHunter output .json file (one or more shards) into a single-sample
ExpansionHunter-style VCF. It exists because the str-truth-set-v2 ExpansionHunter pipeline only persists the EH .json
output (not the .vcf), but EnsembleTR consumes per-caller VCFs and auto-detects ExpansionHunter input by the
'##ALT=<ID=STR\\d+...' header lines plus the INFO VARID/RU/RL and FORMAT REPCN fields. The reconstructed VCF carries
exactly those fields so that TRTools' TRRecordHarmonizer (used internally by EnsembleTR) parses it as ExpansionHunter
output.

The REF allele is written as a single placeholder base ("N") because EnsembleTR obtains reference allele sequences from
its --ref FASTA and TRTools' ExpansionHunter harmonizer derives allele sizes from INFO RL / RU rather than the REF/ALT
sequences.
"""

import argparse
import collections
import gzip
import simplejson as json
import os


def parse_reference_region(reference_region):
    """Parse an ExpansionHunter ReferenceRegion string 'chrom:start_0based-end_1based' into (chrom, start0, end1)."""
    chrom, coords = reference_region.split(":")
    start_0based, end_1based = coords.split("-")
    return chrom, int(start_0based), int(end_1based)


def parse_expansion_hunter_json(json_path):
    """Yield (chrom, pos_1based, end_1based, locus_id, motif, ref_count, allele_counts, ci_strings) tuples from an
    ExpansionHunter output json file.

    allele_counts is the list of integer repeat counts per allele (length 1 for hemizygous calls), and ci_strings is the
    matching list of '{start}-{end}' confidence-interval strings (one per allele).
    """
    fopen = gzip.open if json_path.endswith("gz") else open
    with fopen(json_path, "rt") as f:
        contents = json.load(f)

    sample_id = contents.get("SampleParameters", {}).get("SampleId")
    for locus_record in contents.get("LocusResults", {}).values():
        for variant_json in locus_record.get("Variants", {}).values():
            if "Genotype" not in variant_json or not variant_json.get("RepeatUnit"):
                continue
            chrom, start_0based, end_1based = parse_reference_region(variant_json["ReferenceRegion"])
            motif = variant_json["RepeatUnit"]
            ref_count = (end_1based - start_0based) // len(motif)
            allele_counts = [int(c) for c in variant_json["Genotype"].split("/")]
            confidence_intervals = (variant_json.get("GenotypeConfidenceInterval") or "").split("/")
            ci_strings = []
            for i, count in enumerate(allele_counts):
                ci = confidence_intervals[i] if i < len(confidence_intervals) and confidence_intervals[i] else f"{count}-{count}"
                ci_strings.append(ci)
            # Represent hemizygous calls (1 allele, e.g. male chrX/chrY) as a duplicated diploid genotype. This matches
            # the truth set's hemi representation and, critically, avoids a crash in EnsembleTR's ExpansionHunter score
            # function (utils.GetEHScore indexes REPCN[1]/REPCI[1], assuming two alleles).
            if len(allele_counts) == 1:
                allele_counts = allele_counts * 2
                ci_strings = ci_strings * 2
            yield (sample_id, chrom, start_0based + 1, end_1based, variant_json["VariantId"], motif, ref_count,
                   allele_counts, ci_strings)


def main():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-o", "--output-vcf", required=True, help="Output VCF path (will be written uncompressed; bgzip + "
                   "tabix separately if a .gz path is needed).")
    p.add_argument("--sample-id", help="Sample id to write in the VCF #CHROM line. If not specified, it is taken from "
                   "the SampleParameters.SampleId field of the first input json.")
    p.add_argument("json_paths", nargs="+", help="One or more ExpansionHunter output json (or json.gz) files. Records "
                   "from all files are combined into a single sorted VCF.")
    args = p.parse_args()

    records = []
    distinct_counts = set()
    sample_id = args.sample_id
    for json_path in args.json_paths:
        print(f"Parsing {json_path}")
        for (json_sample_id, chrom, pos, end, locus_id, motif, ref_count, allele_counts, ci_strings) in \
                parse_expansion_hunter_json(json_path):
            if sample_id is None:
                sample_id = json_sample_id
            # build the symbolic ALT alleles for the distinct non-reference repeat counts, in first-seen order
            alt_counts = []
            gt_indices = []
            for count in allele_counts:
                if count == ref_count:
                    gt_indices.append(0)
                else:
                    if count not in alt_counts:
                        alt_counts.append(count)
                    gt_indices.append(1 + alt_counts.index(count))
            distinct_counts.add(ref_count)
            distinct_counts.update(alt_counts)
            alt = ",".join(f"<STR{count}>" for count in alt_counts) if alt_counts else "."
            info = f"END={end};REF={ref_count};REPID={locus_id};RL={end - (pos - 1)};RU={motif};SVTYPE=STR;VARID={locus_id}"
            genotype = "/".join(str(i) for i in gt_indices)
            repcn = "/".join(str(c) for c in allele_counts)
            repci = "/".join(ci_strings)
            records.append((chrom, pos, locus_id, alt, info, genotype, repcn, repci))

    def chrom_sort_key(chrom):
        c = chrom[3:] if chrom.startswith("chr") else chrom
        return (0, int(c)) if c.isdigit() else (1, c)

    records.sort(key=lambda r: (chrom_sort_key(r[0]), r[1]))

    # gather the set of contigs actually used, so the header stays minimal but valid
    contigs = collections.OrderedDict()
    for chrom, *_ in records:
        contigs.setdefault(chrom, True)

    header_lines = ["##fileformat=VCFv4.1"]
    header_lines += [f'##ALT=<ID=STR{count},Description="Allele comprised of {count} repeat units">'
                     for count in sorted(distinct_counts)]
    header_lines += [
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">',
        '##INFO=<ID=REF,Number=1,Type=Integer,Description="Reference copy number">',
        '##INFO=<ID=REPID,Number=1,Type=String,Description="Repeat identifier as specified in the variant catalog">',
        '##INFO=<ID=RL,Number=1,Type=Integer,Description="Reference length in bp">',
        '##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit in the reference orientation">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=VARID,Number=1,Type=String,Description="Variant identifier as specified in the variant catalog">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=REPCN,Number=1,Type=String,Description="Number of repeat units spanned by the allele">',
        '##FORMAT=<ID=REPCI,Number=1,Type=String,Description="Confidence interval for REPCN">',
    ]
    header_lines += [f"##contig=<ID={chrom}>" for chrom in contigs]
    header_lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + (sample_id or "SAMPLE"))

    with open(args.output_vcf, "wt") as out:
        out.write("\n".join(header_lines) + "\n")
        for chrom, pos, locus_id, alt, info, genotype, repcn, repci in records:
            out.write(f"{chrom}\t{pos}\t{locus_id}\tN\t{alt}\t.\t.\t{info}\tGT:REPCN:REPCI\t{genotype}:{repcn}:{repci}\n")

    print(f"Wrote {len(records):,d} records ({len(contigs)} contigs) for sample '{sample_id}' to {args.output_vcf}")


if __name__ == "__main__":
    main()
