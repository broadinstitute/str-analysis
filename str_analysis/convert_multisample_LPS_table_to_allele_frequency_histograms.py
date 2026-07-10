"""
This script takes a wide-format TSV file with columns:

trid     (example: "10-100000859-100000887-A" or "10-100000859-100000887-A,10-100001413-100001429-T")
motif    (example: "A")
<sample1> <sample2> ...  (example: "3,3" meaning two alleles of size 3; "." for no-call)

and outputs a per-locus summary table with allele frequency histograms and statistics.

When ``--vcf-interval-tsv`` is provided (the small TSV produced by
extract_vcf_interval_metadata.py with columns ``trid, locus_id, motif,
interval, vc``), each output row also carries:

    LocusId  = the resolved single locus_id (chrom-start-end-motif) for this row
    Interval = "{chrom}:{vcf_start_0based}-{vcf_end_1based}"  (always set)
    VC       = the inner span from INFO/STRUC if the row was genotyped as
               part of a variation cluster (<VC:...>), or "" for an isolated TR
               (<TR:...>).

These columns disambiguate the rows that share a LocusId because the same
LocusId was genotyped under multiple TRGT catalog intervals (e.g. once as a
standalone TR and once inside a VC).

For each LPS row's (trid, motif) key, the script pops one chunk from the
pre-built (trid, motif) -> deque[chunk] map, where each chunk groups all
LocusIds genotyped by the same VCF record (so a compound TRID with three
LocusIds ending in ``-AGAA`` for motif ``AGAA`` produces a single chunk
containing three locus_ids, all sharing the same interval/vc, and the
script emits three output rows from the LPS row's allele data).

The script enforces a process-wide uniqueness invariant: no two output rows
share the same (LocusId, Interval, VC) tuple.
"""

import argparse
import collections
import gzip
import json
import os
import numpy as np
import tqdm


HEADER_FIELDS = [
    "LocusId",
    "Motif",
    "Interval",
    "VC",
    "AlleleSizeHistogram",
    "BiallelicHistogram",
    "Min",
    "Mode",
    "Mean",
    "Stdev",
    "Median",
    "99thPercentile",
    "Max",
    "ShortAllele99thPercentile",
    "ShortAlleleMax",
    "HemiAllele99thPercentile",
    "HemiAlleleMax",
    "UniqueAlleleLengths",
    "NumCalledAlleles",
]


def load_vcf_interval_metadata(tsv_path):
    """Loads the small interval-metadata TSV into ``(trid, motif) -> deque[chunk]``.

    The TSV is produced by ``data-prep/hprc-lps/extract_vcf_interval_metadata.py``
    with columns ``trid, locus_id, motif, interval, vc``. Consecutive rows
    with the same ``(trid, motif, interval, vc)`` come from the same VCF
    record and are grouped into a single chunk
    ``(interval, vc, [locus_id, ...])``. The deque preserves VCF order across
    distinct VCF records that share a ``(trid, motif)`` key (e.g. one TR plus
    one VC genotyping of the same LocusId).

    Args:
        tsv_path: Path to the gzipped or plain TSV.

    Returns:
        ``dict[(trid, motif), collections.deque[(interval, vc, [locus_id, ...])]]``.
    """
    opener = gzip.open if str(tsv_path).endswith(".gz") else open
    metadata = {}
    current_key = None
    current_chunk_key = None
    current_locus_ids = None

    def flush():
        if current_key is None or current_chunk_key is None:
            return
        interval, vc = current_chunk_key
        metadata.setdefault(current_key, collections.deque()).append(
            (interval, vc, current_locus_ids)
        )

    with opener(tsv_path, "rt") as f:
        header = next(f).rstrip("\n").split("\t")
        expected = ["trid", "locus_id", "motif", "interval", "vc"]
        if header != expected:
            raise ValueError(
                f"--vcf-interval-tsv header must be {expected!r}; got {header!r}"
            )
        for line in f:
            trid, locus_id, motif, interval, vc = line.rstrip("\n").split("\t")
            key = (trid, motif)
            chunk_key = (interval, vc)
            if key != current_key or chunk_key != current_chunk_key:
                flush()
                current_key = key
                current_chunk_key = chunk_key
                current_locus_ids = [locus_id]
            else:
                current_locus_ids.append(locus_id)
        flush()
    return metadata


def _format_decimal(value):
    """Format a numeric value as an integer if whole, otherwise round to 2 decimal places."""
    if value == int(value):
        return int(value)
    return round(float(value), 2)


def compute_histograms(allele_sizes, alleles_by_sample_id):
    """Compute allele size and biallelic histogram strings.

    Args:
        allele_sizes: list of allele sizes
        alleles_by_sample_id: dict mapping sample_id to list of allele sizes

    Returns:
        tuple of (allele_size_histogram_str, biallelic_histogram_str), or ("", "") if empty
    """
    if not allele_sizes:
        return "", ""

    allele_counts = collections.Counter(allele_sizes)
    genotype_counts = collections.defaultdict(int)
    for allele_list in alleles_by_sample_id.values():
        if len(allele_list) == 1:
            allele_list = allele_list * 2
        genotype_counts[tuple(sorted(allele_list))] += 1

    return (
        ",".join(f"{size}x:{count}" for size, count in sorted(allele_counts.items())),
        ",".join(f"{g[0]}/{g[1]}:{count}" for g, count in sorted(genotype_counts.items(), key=lambda x: (x[0][0], x[0][1]))),
    )


def compute_row(locus_id, motif, allele_sizes, alleles_by_sample_id, interval="", vc="", sample_id_to_sex=None):
    """Compute statistics for a group of allele sizes.

    Args:
        locus_id (str): the single LocusId (chrom-start-end-motif) for this row
        motif (str): the motif
        allele_sizes (list): the allele sizes for all samples
        alleles_by_sample_id (dict): the alleles for the current key by sample id
        interval (str): TRGT interval ``"{chrom}:{vcf_start_0based}-{vcf_end_1based}"`` or ``""``
        vc (str): inner ``<VC:...>`` span or ``""`` for an isolated TR
        sample_id_to_sex (dict): sample_id -> "male"/"female", or None/empty if unavailable

    Returns:
        dict: a dictionary mapping HEADER_FIELDS keys to values, or None if allele_sizes is empty
    """

    if not allele_sizes:
        return None

    allele_sizes = sorted(allele_sizes)
    allele_counts = collections.Counter(allele_sizes)

    short_alleles = []
    for sample_id, allele_list in alleles_by_sample_id.items():
        if len(allele_list) == 1:
            short_alleles.append(allele_list[0])
        elif len(allele_list) == 2:
            short_alleles.append(min(allele_list))
        else:
            raise ValueError(f"Found {len(allele_list)} alleles for {sample_id} in {locus_id} {motif}")

    allele_histogram, biallelic_histogram = compute_histograms(allele_sizes, alleles_by_sample_id)

    # Hemi* columns cover male-only allele sizes at chrX/chrY loci (hemizygous calls),
    # so they aren't diluted by the female diploid calls that ShortAllele*/AlleleSize* mix in for chrX.
    hemi_allele_99th_percentile = ""
    hemi_allele_max = ""
    chrom = locus_id.split("-", 1)[0]
    if sample_id_to_sex and chrom in ("X", "Y"):
        hemi_alleles = [
            allele_size
            for sample_id, allele_list in alleles_by_sample_id.items()
            if sample_id_to_sex.get(sample_id) == "male"
            for allele_size in allele_list
        ]
        if hemi_alleles:
            hemi_allele_99th_percentile = _format_decimal(np.percentile(hemi_alleles, 99))
            hemi_allele_max = int(max(hemi_alleles))

    return {
        "LocusId": locus_id,
        "Motif": motif,
        "Interval": interval,
        "VC": vc,
        "AlleleSizeHistogram": allele_histogram,
        "BiallelicHistogram": biallelic_histogram,
        "Min": int(min(allele_sizes)),
        "Mode": min(size for size, count in allele_counts.items() if count == allele_counts.most_common(1)[0][1]),
        "Mean": f"{np.mean(allele_sizes):.2f}",
        "Stdev": f"{np.std(allele_sizes):.2f}",
        "Median": _format_decimal(np.median(allele_sizes)),
        "99thPercentile": _format_decimal(np.percentile(allele_sizes, 99)),
        "Max": int(max(allele_sizes)),
        "ShortAllele99thPercentile": _format_decimal(np.percentile(short_alleles, 99)),
        "ShortAlleleMax": int(max(short_alleles)),
        "HemiAllele99thPercentile": hemi_allele_99th_percentile,
        "HemiAlleleMax": hemi_allele_max,
        "UniqueAlleleLengths": len(set(allele_sizes)),
        "NumCalledAlleles": len(allele_sizes),
    }


"""
Example row in 1kgp metadata TSV file:

SampleId                                                            HG00459
Sex                                                                    male
BiosampleId                                                      SAME125269
Subpopulation                                                           CHS
SubpopulationName                                      Southern Han Chinese
Population                                                              EAS
PopulationName                                          East Asian Ancestry
SubpopulationElasticId                                                  CHS
DataSource                1000 Genomes 30x on GRCh38,1000 Genomes phase ...

Population counts:
   1427 AFR
    535 AMR
    623 EAS
    672 EUR
      1 EUR,AFR
    661 SAS


 2343 female
   2653 male

Input table example row:

trid       1-19175-19184-TCC
motif                    TCC
HG00096                  3,3
HG00097                  3,3
HG00099                  3,3
...


"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-metadata-tsv", default="https://storage.googleapis.com/tandem-repeat-catalog/1kGP_metadata.tsv", help="Sample ancestry metadata TSV file (local path or URL)")
    parser.add_argument("--input-table", default="hprc-lps_2025-12-06/hprc-lps.txt.gz", help="Combine HPRC LPS dataset")
    parser.add_argument("--vcf-interval-tsv", help="Path to a small TSV.gz produced by data-prep/hprc-lps/extract_vcf_interval_metadata.py (columns: trid, locus_id, motif, interval, vc). Required to populate the LocusId/Interval/VC output columns and resolve compound TRIDs.")
    parser.add_argument("--no-header", action="store_true", help="If set, assume the first row is data (not a header) and generate synthetic sample names (_s1, _s2, ...)")
    parser.add_argument("--population", choices=["AFR", "AMR", "EAS", "EUR", "SAS"], help="If specified, only process samples from this population")
    parser.add_argument("--sex", choices=["male", "female"], help="If specified, only process samples from this sex")
    parser.add_argument("--stratify-by-population", action="store_true", help="If set, add per-population histogram columns to the output")
    parser.add_argument("--stratify-by-sex", action="store_true", help="If set, add per-sex histogram columns to the output")
    parser.add_argument("--output-format", choices=["TSV", "JSON"], default="TSV", help="Output format (default: TSV)")
    parser.add_argument("-n", "--num-samples", type=int, default=None, help="Number of samples to process")
    parser.add_argument("-l", "--num-loci", type=int, default=None, help="Number of loci to process")
    args = parser.parse_args()

    if args.stratify_by_population and args.population:
        parser.error("--stratify-by-population and --population are mutually exclusive")
    if args.stratify_by_sex and args.sex:
        parser.error("--stratify-by-sex and --sex are mutually exclusive")
    if (args.stratify_by_population or args.stratify_by_sex) and args.no_header:
        parser.error("--stratify-by-population and --stratify-by-sex require a header row with real sample ids")

    if not os.path.isfile(args.input_table):
        parser.error(f"Input file {args.input_table} does not exist")

    if args.vcf_interval_tsv and not os.path.isfile(args.vcf_interval_tsv):
        parser.error(f"--vcf-interval-tsv {args.vcf_interval_tsv} does not exist")

    fopen = gzip.open if args.input_table.endswith("gz") else open

    if args.no_header:
        # Read first line to determine number of columns, then generate synthetic header
        with fopen(args.input_table, "rt") as infile:
            first_line = next(infile).strip().split("\t")
        num_sample_columns = len(first_line) - 2
        header_fields = ["trid", "motif"] + [f"_s{i+1}" for i in range(num_sample_columns)]
        sample_ids_in_input_table = header_fields[2:]
    else:
        with fopen(args.input_table, "rt") as infile:
            header_fields = next(infile).strip().split("\t")
            if header_fields[0:2] != ["trid", "motif"]:
                parser.error(f"Expected header to start with 'trid', 'motif' columns in {args.input_table}. Found: {header_fields}")
        sample_ids_in_input_table = header_fields[2:]
        print(f"Parsed {len(sample_ids_in_input_table):,d} sample ids from {args.input_table}")

    if args.no_header:
        if args.population or args.sex:
            parser.error("--population and --sex filters require a header row with real sample ids and a --sample-metadata-tsv")
        sample_ids_to_include_list = sample_ids_in_input_table
        if args.num_samples is not None and len(sample_ids_to_include_list) > args.num_samples:
            sample_ids_to_include_list = sample_ids_to_include_list[:args.num_samples]
        sample_id_to_strata = {}
        strata_labels = []
        sample_id_to_sex = {}
    else:
        import pandas as pd
        df_metadata = pd.read_table(args.sample_metadata_tsv, dtype={"SampleId": str})
        all_sample_ids = set(df_metadata.SampleId)
        unexpected_sample_ids = set(sample_ids_in_input_table) - all_sample_ids
        if unexpected_sample_ids:
            parser.error(f"{len(unexpected_sample_ids):,d} sample id(s) not found in {args.sample_metadata_tsv}: {', '.join(sorted(unexpected_sample_ids))}")

        df_metadata = df_metadata[df_metadata.SampleId.isin(set(sample_ids_in_input_table))]
        if args.population:
            df_metadata = df_metadata[df_metadata.Population == args.population]
            print(f"Kept {len(df_metadata):,d} samples from population {args.population}")
        if args.sex:
            df_metadata = df_metadata[df_metadata.Sex == args.sex]
            print(f"Kept {len(df_metadata):,d} samples from sex {args.sex}")

        valid_ids = set(df_metadata.SampleId)
        sample_ids_to_include_list = [s for s in sample_ids_in_input_table if s in valid_ids]
        if args.num_samples is not None and len(sample_ids_to_include_list) > args.num_samples:
            sample_ids_to_include_list = sample_ids_to_include_list[:args.num_samples]
            df_metadata = df_metadata[df_metadata.SampleId.isin(set(sample_ids_to_include_list))]

        sample_id_to_sex = dict(zip(df_metadata.SampleId, df_metadata.Sex))

        # Build sample_id -> list of stratum labels. When both stratify flags are
        # set, each sample contributes to its Pop_Sex cell plus the Pop row-marginal
        # (across both sexes) and the Sex column-marginal (across all populations).
        sample_id_to_strata = {}
        strata_labels = []
        if args.stratify_by_population or args.stratify_by_sex:
            df_included = df_metadata[df_metadata.SampleId.isin(set(sample_ids_to_include_list))]
            if args.stratify_by_population and args.stratify_by_sex:
                sample_id_to_strata = {
                    sid: [f"{pop}_{sex}", pop, sex]
                    for sid, pop, sex in zip(df_included.SampleId, df_included.Population, df_included.Sex)
                }
            elif args.stratify_by_population:
                sample_id_to_strata = {sid: [pop] for sid, pop in zip(df_included.SampleId, df_included.Population)}
            else:
                sample_id_to_strata = {sid: [sex] for sid, sex in zip(df_included.SampleId, df_included.Sex)}
            strata_labels = sorted({label for labels in sample_id_to_strata.values() for label in labels})

    sample_ids_to_include = set(sample_ids_to_include_list)

    # Build output header with stratified histogram columns
    output_header = list(HEADER_FIELDS)
    for label in strata_labels:
        output_header.append(f"AlleleSizeHistogram__{label}")
    for label in strata_labels:
        output_header.append(f"BiallelicHistogram__{label}")

    # Build output path
    output_dir = os.path.dirname(args.input_table)
    output_path = os.path.join(output_dir, os.path.basename(args.input_table).replace(".tsv", "").replace(".txt", "").replace(".gz", ""))
    output_path += ".per_locus_and_motif"
    if args.population:
        output_path += f".only_{args.population}"
    if args.sex:
        output_path += f".only_{args.sex}"
    if args.stratify_by_population:
        output_path += ".by_population"
    if args.stratify_by_sex:
        output_path += ".by_sex"
    output_path += f".{len(sample_ids_to_include_list)}_samples"
    if args.output_format == "JSON":
        output_path += ".json.gz"
    else:
        output_path += ".tsv.gz"

    # Build the full (trid, motif) -> deque[chunk] map from the small interval
    # TSV when --vcf-interval-tsv is given. Each chunk groups all LocusIds
    # genotyped by one VCF record (same interval/vc, differing locus_id).
    vcf_metadata = {}
    if args.vcf_interval_tsv:
        print(f"Loading interval metadata from {args.vcf_interval_tsv}")
        vcf_metadata = load_vcf_interval_metadata(args.vcf_interval_tsv)
        chunks_total = sum(len(v) for v in vcf_metadata.values())
        print(f"Loaded {chunks_total:,d} VCF-record chunks across {len(vcf_metadata):,d} unique (trid, motif) keys")

    # Process-wide uniqueness check for emitted (LocusId, Interval, VC) tuples.
    seen_output_keys = set()

    # Atomic write: stream into a .tmp path next to the destination so a
    # mid-run exception (malformed cell, end-of-input drain mismatch, etc.)
    # doesn't destroy a previous good output file. os.replace promotes it
    # at the end of the function, after the end-of-processing assertion has
    # passed. A try/finally guarantees the .tmp is removed on any failure
    # path, matching extract_vcf_interval_metadata.py and the purity/methylation
    # script's pattern.
    tmp_output_path = output_path + ".tmp"
    print(f"Writing data from {len(sample_ids_to_include_list):,d} samples to {output_path}")
    promoted = False
    try:
     with fopen(args.input_table, "rt") as infile, gzip.open(tmp_output_path, "wt") as outfile:
        if not args.no_header:
            next(infile)  # skip header

        if args.output_format == "TSV":
            outfile.write("\t".join(output_header) + "\n")
        else:
            outfile.write("[\n")
            json_first_row = True

        line_count = 0
        rows_written = 0
        for line_number, line in tqdm.tqdm(enumerate(infile), unit=" lines", unit_scale=True):
            fields = line.strip().split("\t")
            if len(fields) < 2:
                continue  # skip blank or malformed lines
            if fields[0] == "trid" and fields[1] == "motif":
                continue  # skip internal header rows in concatenated files

            if args.num_loci is not None and line_count >= args.num_loci:
                break
            line_count += 1

            if len(fields) != len(header_fields):
                raise ValueError(f"Line {line_number + 1} has {len(fields)} fields, expected {len(header_fields)}: {line.strip()[:200]}")

            trid = fields[0]
            motif = fields[1]

            # Resolve the (trid, motif) key against the pre-loaded VCF metadata
            # to determine which LocusIds and which (interval, vc) this LPS row
            # represents. If --vcf-interval-tsv was not provided, fall back to
            # the legacy single-LocusId / empty Interval&VC behavior.
            if vcf_metadata:
                chunks = vcf_metadata.get((trid, motif))
                if chunks is None:
                    raise ValueError(
                        f"Line #{line_number + 1}: no VCF record found for "
                        f"(trid={trid!r}, motif={motif!r}) in --vcf-interval-tsv"
                    )
                if not chunks:
                    raise ValueError(
                        f"Line #{line_number + 1}: --vcf-interval-tsv has fewer "
                        f"chunks than the LPS table has rows for "
                        f"(trid={trid!r}, motif={motif!r}); all chunks for this "
                        f"key were already consumed by earlier LPS rows."
                    )
                interval, vc, chunk_locus_ids = chunks.popleft()
            else:
                if "," in trid:
                    raise ValueError(
                        f"Line #{line_number + 1}: compound TRID {trid!r} cannot "
                        f"be resolved without --vcf-interval-tsv"
                    )
                interval = ""
                vc = ""
                chunk_locus_ids = [trid]

            alleles = []
            alleles_by_sample_id = collections.defaultdict(list)
            stratum_alleles = collections.defaultdict(list)
            stratum_alleles_by_sample_id = collections.defaultdict(lambda: collections.defaultdict(list))
            for sample_id, allele_sizes in zip(header_fields[2:], fields[2:]):
                if sample_id not in sample_ids_to_include:
                    continue
                if allele_sizes == ".":
                    continue
                for allele_size in allele_sizes.split(","):
                    try:
                        allele_size = int(allele_size)
                    except ValueError:
                        raise ValueError(f"Expected integer allele size, got {allele_size} in line #{line_number + 1}: {line}")
                    alleles.append(allele_size)
                    alleles_by_sample_id[sample_id].append(allele_size)
                    for stratum in sample_id_to_strata.get(sample_id, ()):
                        stratum_alleles[stratum].append(allele_size)
                        stratum_alleles_by_sample_id[stratum][sample_id].append(allele_size)

            # Emit one output row per LocusId in the chunk. All rows share the
            # same per-sample allele data (computed once from the LPS row) but
            # differ in LocusId. The (LocusId, Interval, VC) tuple must be
            # process-globally unique.
            stratified_columns = {}
            for label in strata_labels:
                allele_histogram, biallelic_histogram = compute_histograms(
                    stratum_alleles.get(label, []),
                    stratum_alleles_by_sample_id.get(label, {}),
                )
                stratified_columns[f"AlleleSizeHistogram__{label}"] = allele_histogram
                stratified_columns[f"BiallelicHistogram__{label}"] = biallelic_histogram

            for locus_id in chunk_locus_ids:
                row = compute_row(locus_id, motif, alleles, alleles_by_sample_id, interval=interval, vc=vc, sample_id_to_sex=sample_id_to_sex)
                if row is None:
                    continue
                key = (row["LocusId"], row["Interval"], row["VC"])
                if key in seen_output_keys:
                    raise ValueError(
                        f"Line #{line_number + 1}: duplicate output tuple "
                        f"(LocusId={key[0]!r}, Interval={key[1]!r}, VC={key[2]!r})"
                    )
                seen_output_keys.add(key)
                row.update(stratified_columns)
                if args.output_format == "JSON":
                    if not json_first_row:
                        outfile.write(",\n")
                    json_first_row = False
                    outfile.write("  " + json.dumps(row, indent=2).replace("\n", "\n  "))
                else:
                    outfile.write("\t".join(str(row[field]) for field in output_header) + "\n")
                rows_written += 1

        if args.output_format == "JSON":
            outfile.write("\n]\n")

     # End-of-processing assertion: every chunk popped from the deque must
     # correspond to an LPS row. Leftover chunks indicate a mismatch between
     # the LPS table and --vcf-interval-tsv (e.g. extract emitted records that
     # the LPS table doesn't have, which would mean LPS rows got paired with
     # the wrong (Interval, VC) chunks via FIFO).
     #
     # Skip the check when --num-loci is in effect: the main loop intentionally
     # breaks early, so leftover chunks are expected and not a sign of mismatch.
     if vcf_metadata and args.num_loci is None:
         unconsumed = sum(len(deq) for deq in vcf_metadata.values())
         if unconsumed:
             examples = [key for key, deq in vcf_metadata.items() if deq][:3]
             raise ValueError(
                 f"{unconsumed:,d} VCF-record chunks in --vcf-interval-tsv were "
                 f"never consumed by an LPS row (e.g. {examples!r}). The two "
                 f"inputs are out of sync; re-extract --vcf-interval-tsv from "
                 f"the same VCF the LPS table was generated from."
             )

     # Atomic promote: only if the run succeeded (no exception, drain check passed).
     os.replace(tmp_output_path, output_path)
     promoted = True
    finally:
        # Clean up the .tmp file on any failure path (mid-loop exception,
        # drain-check failure, anything else); the os.replace above only
        # fires when the full run succeeds.
        if not promoted and os.path.exists(tmp_output_path):
            try:
                os.remove(tmp_output_path)
            except OSError:
                pass
    print(f"Wrote {rows_written:9,d} rows to {output_path}")

if __name__ == "__main__":
    main()
