"""Extract per-allele repeat sequences from bw2 ExpansionHunter's output JSON.

Sibling of extract_allele_sequences_from_vcf.py, for the one tool whose sequence has to come from a different
source format: ExpansionHunter's VCF ALT is a symbolic allele ("<STR12>"), never an actual sequence, so it can't be
parsed the way TRGT/ATaRVa/HipSTR's VCFs are. What IS available is the "ConsensusSequences" field the bw2 fork can
write to its JSON output (docs @ https://github.com/Illumina/ExpansionHunter/blob/master/docs/05_OutputJsonFiles.md,
plus the bw2-only reviewer/ConsensusSequence.cpp) -- a real base-level consensus built from realigned reads, one
string per allele, ordered short-then-long (same order as "Genotype"). It's off by default (kept off in this
project's genotyping runs to keep the regenerated truth-set JSON lean -- see expansion_hunter_pipeline.py's
enable_consensus_sequences); a run for the sequence-accuracy benchmark must pass enable_consensus_sequences=True.

ExpansionHunter marks ambiguous bases it couldn't resolve from the reads as 'N'. Those are real "we don't know"
positions, not evidence of a match, so add_sequence_accuracy_columns.py's edit-distance calculation never lets an
'N' count as equal to anything -- including another 'N' in the truth sequence.

Output schema matches extract_allele_sequences_from_vcf.py's {prefix}.allele_sequences.tsv.gz exactly (see
OUTPUT_COLUMNS there), so add_sequence_accuracy_columns.py consumes either one unmodified:

    LocusId                    locus id with any leading "chr" stripped, so it joins against the tool comparison tables
    AlleleSequence: Allele 1   the shorter allele's consensus sequence
    AlleleSequence: Allele 2   the longer allele's consensus sequence
    AlleleStatus: Allele 1     "called" or "no_call"
    AlleleStatus: Allele 2     "called" or "no_call"
    ExtractionStatus           "ok" or "multiallelic" (per locus)

Both of optimized-streaming mode's genotyping paths emit consensus sequences: the fast spanning-read heuristic
(HtsLowMemStreamingHelpers.cpp, which handles most loci) as well as the full graph-alignment genotyper the harder
loci fall through to. The two derive the consensus differently -- the fast path uses only whole-repeat-spanning
reads and does no realignment -- so a fast-path sequence will not always match what full genotyping would have
produced for the same locus; see "Consensus sequences in fast genotyping mode" in ExpansionHunter's
docs/05_OutputJsonFiles.md for the full list of differences. Loci genotyped by the fast path are flagged with
"QuickGenotype": true in the json.

Allele shapes: a heterozygous call carries one consensus per allele, while a homozygous call carries a SINGLE
consensus covering both of its identical alleles (both paths do this) -- that one is duplicated into both output
slots, exactly as a haploid call is. A locus contributes no_call/no_call when ExpansionHunter didn't genotype it at
all (no "Genotype" key), when consensus sequences were disabled or couldn't be built, or when the shape is anything
other than those two. A locus with more than one repeat variant is also always no_call/no_call, since
locus/LocusAnalyzer.cpp explicitly skips consensus-building for those ("multiple repeat variants not supported").
ref_mismatch (from the VCF extractor) doesn't apply here: there's no REF/ALT trimmed against a reference interval,
just the tool's own consensus.
"""

import argparse
import collections
import gzip
import json
import os
import re

from str_analysis.extract_allele_sequences_from_vcf import (
    CALLED, NO_CALL, EXTRACTION_OK, EXTRACTION_MULTIALLELIC, OUTPUT_COLUMNS, make_row)


def extract_allele_sequences_from_locus_json(locus_json):
    """Convert one LocusResults entry into an allele-sequence row.

    Args:
        locus_json (dict): one value from json_contents["LocusResults"].

    Returns:
        dict: a row with the OUTPUT_COLUMNS keys, or None if the locus has no usable LocusId or no repeat variant
        with a "Genotype" (ExpansionHunter never emitted a call for it).
    """
    locus_id = locus_json.get("LocusId")
    if not locus_id:
        return None

    variant_json = None
    for candidate in locus_json.get("Variants", {}).values():
        if "Genotype" in candidate:
            if variant_json is not None:
                # more than one genotyped repeat variant at this locus -- LocusAnalyzer.cpp doesn't build a
                # consensus for these either, so there's nothing to extract
                return make_row(locus_id, extraction_status=EXTRACTION_MULTIALLELIC)
            variant_json = candidate

    if variant_json is None:
        return None

    genotype_fields = variant_json["Genotype"].split("/")
    if len(genotype_fields) > 2:
        return make_row(locus_id, extraction_status=EXTRACTION_MULTIALLELIC)

    # A homozygous call carries ONE consensus for its two identical alleles (both the full genotyper and the fast
    # path do this -- see "Homozygous calls get a single consensus" in HtsLowMemStreamingHelpers.cpp), so
    # "one sequence per allele" and "one sequence for an N-ploid homozygous call" are both valid shapes. Anything
    # else is unexpected and is treated as no data rather than guessed at.
    consensus_sequences = variant_json.get("ConsensusSequences")
    is_homozygous = len(set(genotype_fields)) == 1
    if not consensus_sequences or not (
            len(consensus_sequences) == len(genotype_fields) or (len(consensus_sequences) == 1 and is_homozygous)):
        # no consensus available for this locus (e.g. --dont-output-consensus-sequences was used, or the read
        # support was too poor to build one), or a malformed/unexpected ConsensusSequences shape
        return make_row(locus_id)

    sequences = list(consensus_sequences)
    statuses = [CALLED] * len(sequences)
    if len(sequences) == 1:
        # Either a homozygous diploid call (one consensus covering both identical alleles) or a haploid call (male
        # chrX/chrY outside the PAR). The truth set stores HEMI loci as a duplicated diploid genotype, and
        # add_tool_results_columns.py mirrors that for the repeat-count columns, so duplicate here too rather than
        # leaving Allele 2 as a no-call -- which would drop every homozygous locus from the benchmark.
        sequences.append(sequences[0])
        statuses.append(CALLED)

    # defensively re-sort short-first (ExpansionHunter already writes Genotype/ConsensusSequences in that order --
    # see JsonWriter.cpp's alleleRank<=0 -> shortAlleleSizeInUnitsCi() -- but matching how the VCF extractor orders
    # its output costs nothing and doesn't depend on that ordering holding in every ExpansionHunter version)
    sequences, statuses = map(list, zip(*sorted(zip(sequences, statuses), key=lambda pair: len(pair[0]))))

    return make_row(locus_id, sequences[0], sequences[1], statuses[0], statuses[1])


def process_expansion_hunter_json_files(json_paths, output_path, show_progress_bar=False):
    """Extract allele sequences from one or more ExpansionHunter JSON output files and write them to a single tsv.

    Args:
        json_paths (list): paths of the tool's output json(.gz) file(s). This pipeline always runs ExpansionHunter
            unsharded, so there's normally just one, but multiple are accepted for parity with the VCF extractor.
        output_path (str): path of the output .tsv.gz.
        show_progress_bar (bool): whether to show a tqdm progress bar.

    Returns:
        collections.Counter: counts of each ExtractionStatus, each AlleleStatus, and "total" loci.
    """
    counters = collections.Counter()
    with gzip.open(output_path, "wt") as output_file:
        output_file.write("\t".join(OUTPUT_COLUMNS) + "\n")
        for json_path in json_paths:
            print(f"Processing {json_path}")
            with (gzip.open if json_path.endswith("gz") else open)(json_path, "rt") as f:
                json_contents = json.load(f)

            locus_results = list(json_contents.get("LocusResults", {}).values())
            if show_progress_bar:
                from tqdm import tqdm
                locus_results = tqdm(locus_results, unit=" loci", unit_scale=True, unit_divisor=1000)

            for locus_json in locus_results:
                counters["total"] += 1
                row = extract_allele_sequences_from_locus_json(locus_json)
                if row is None:
                    counters["no_genotype"] += 1
                    continue
                counters[row["ExtractionStatus"]] += 1
                counters[row["AlleleStatus: Allele 1"]] += 1
                counters[row["AlleleStatus: Allele 2"]] += 1
                output_file.write("\t".join(str(row[column]) for column in OUTPUT_COLUMNS) + "\n")

    print(f"Parsed {counters['total']:,d} loci from {len(json_paths):,d} file(s):")
    print(f"{counters[EXTRACTION_OK]:12,d} loci extracted ok")
    print(f"{counters[EXTRACTION_MULTIALLELIC]:12,d} loci skipped: multiple repeat variants (no consensus built)")
    print(f"{counters['no_genotype']:12,d} loci skipped: ExpansionHunter didn't genotype them")
    print(f"{counters[CALLED]:12,d} alleles called, {counters[NO_CALL]:,d} alleles not called "
          f"(no ConsensusSequences, e.g. --dont-output-consensus-sequences)")
    print(f"Wrote {counters['total'] - counters['no_genotype']:,d} rows to {output_path}")

    return counters


def main():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    p.add_argument("-o", "--output-path",
                   help="Output .tsv.gz path. Defaults to {input json name}.allele_sequences.tsv.gz")
    p.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    p.add_argument("json_paths", nargs="+", help="ExpansionHunter output json(.gz) path(s).")
    args = p.parse_args()

    for json_path in args.json_paths:
        if not os.path.isfile(json_path):
            p.error(f"{json_path} not found")

    output_path = args.output_path or (
        re.sub(r"\.json(\.gz)?$", "", os.path.basename(args.json_paths[0])) + ".allele_sequences.tsv.gz")

    process_expansion_hunter_json_files(args.json_paths, output_path, show_progress_bar=args.show_progress_bar)


if __name__ == "__main__":
    main()
