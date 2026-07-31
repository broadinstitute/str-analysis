"""Extract per-allele repeat sequences from a tandem repeat genotyping tool's output VCF.

This is the input to the sequence-accuracy benchmark, which scores a tool by the edit distance between the allele
sequence it reports and the assembly-derived truth allele sequence, rather than by repeat count (which is what
plot_tool_accuracy_by_allele_size.py measures).

Supported tools: TRGT (v3/v5), ATaRVa, and HipSTR - the three tools in the comparison whose VCF carries an actual
allele sequence in REF/ALT.

The output is a {prefix}.allele_sequences.tsv.gz table with one row per VCF record:

    LocusId                    locus id with any leading "chr" stripped, so it joins against the tool comparison tables
    AlleleSequence: Allele 1   the shorter allele's sequence, spanning exactly the locus interval
    AlleleSequence: Allele 2   the longer allele's sequence
    AlleleStatus: Allele 1     "called" or "no_call"
    AlleleStatus: Allele 2     "called" or "no_call"
    ExtractionStatus           "ok", "ref_mismatch", or "multiallelic" (per locus)

An empty AlleleSequence is a genuine zero-length allele ONLY when its AlleleStatus is "called" (ATaRVa writes <DEL>
for these). Every other empty value is explained by the row's statuses, so downstream code must never treat an empty
string as a zero-length allele without first checking the status.

Coordinates: each supported tool's locus id has the truth set's "{chrom}-{start_0based}-{end}-{motif}" form, so the
locus interval is parsed from the locus id itself. That interval - not the VCF record's own POS/REF span - is what the
truth allele sequences span, so REF and every ALT are trimmed to it. Both TRGT and HipSTR emit records that extend
past the repeat: TRGT writes a leading padding base (POS is one before the repeat start), and HipSTR can extend the
genotyped region beyond the requested bed interval. With --reference-fasta, the trimmed REF is checked against the
reference genome, and any locus that doesn't match is reported as ExtractionStatus=ref_mismatch rather than being
silently mis-scored.
"""

import argparse
import collections
import gzip
import os
import re

TOOL_CHOICES = ("TRGTv3", "TRGTv5", "ATaRVa", "HipSTR")

CALLED = "called"
NO_CALL = "no_call"

EXTRACTION_OK = "ok"
EXTRACTION_REF_MISMATCH = "ref_mismatch"
EXTRACTION_MULTIALLELIC = "multiallelic"

OUTPUT_COLUMNS = [
    "LocusId",
    "AlleleSequence: Allele 1",
    "AlleleSequence: Allele 2",
    "AlleleStatus: Allele 1",
    "AlleleStatus: Allele 2",
    "ExtractionStatus",
]


def parse_info(info):
    """Parse a VCF INFO column string into a dict. Flag fields (no '=') map to True.

    Args:
        info (str): the INFO column of one VCF record.

    Returns:
        dict: {key: value} for "key=value" entries, {key: True} for flags.
    """
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


def get_locus_id(tool, chrom, id_field, info_dict):
    """Return one VCF record's locus id, in the truth set's "{chrom}-{start_0based}-{end}-{motif}" form.

    Args:
        tool (str): one of TOOL_CHOICES.
        chrom (str): the record's CHROM field.
        id_field (str): the record's ID field.
        info_dict (dict): the record's parsed INFO field.

    Returns:
        str: the locus id.
    """
    if tool in ("TRGTv3", "TRGTv5"):
        return info_dict["TRID"]
    if tool == "HipSTR":
        return id_field
    if tool == "ATaRVa":
        # ATaRVa doesn't echo the catalog's locus id, so rebuild it from the repeat coordinates it does echo
        # (START is 0-based, END is 1-based), the same way convert_atarva_vcf_to_expansion_hunter_json.py does.
        return f"{chrom}-{info_dict['START']}-{info_dict['END']}-{info_dict['MOTIF']}"
    raise ValueError(f"Unexpected tool: {tool}")


def parse_locus_interval(locus_id):
    """Parse a "{chrom}-{start_0based}-{end}-{motif}" locus id into its repeat interval.

    Args:
        locus_id (str): the locus id.

    Returns:
        tuple: (start_0based, end), or None if the locus id isn't in that form.
    """
    fields = re.split("[_-]", locus_id)
    if len(fields) < 3 or not fields[1].isdigit() or not fields[2].isdigit():
        return None
    return int(fields[1]), int(fields[2])


def make_row(locus_id, sequence1="", sequence2="", status1=NO_CALL, status2=NO_CALL, extraction_status=EXTRACTION_OK):
    """Assemble one output row. See OUTPUT_COLUMNS for what each field means."""
    return {
        "LocusId": re.sub("^chr", "", locus_id),
        "AlleleSequence: Allele 1": sequence1,
        "AlleleSequence: Allele 2": sequence2,
        "AlleleStatus: Allele 1": status1,
        "AlleleStatus: Allele 2": status2,
        "ExtractionStatus": extraction_status,
    }


def extract_allele_sequences_from_vcf_record(line, tool, get_reference_sequence=None):
    """Convert one VCF data line into an allele-sequence row.

    Resolves the GT indices to REF/ALT sequences, trims them to the locus interval encoded in the locus id, and
    size-sorts the two alleles so Allele 1 is the shorter one (matching how add_tool_results_columns.py orders the
    repeat-count columns).

    Args:
        line (str): one non-header VCF line.
        tool (str): one of TOOL_CHOICES; selects how the locus id is obtained (see get_locus_id).
        get_reference_sequence (callable): optional f(chrom, start_0based, end) -> uppercase reference sequence. When
            given, a trimmed REF that doesn't match the reference genome makes the locus ExtractionStatus=ref_mismatch.

    Returns:
        dict: a row with the OUTPUT_COLUMNS keys, or None if the record is malformed (no genotype column, an
        unparseable locus id, or a GT index with no corresponding allele). Callers count the None rows as loci the
        tool didn't produce a usable call for.
    """
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 10:
        return None

    chrom = fields[0]
    pos_1based = int(fields[1])
    ref = fields[3].upper()
    info_dict = parse_info(fields[7])

    locus_id = get_locus_id(tool, chrom, fields[2], info_dict)
    locus_interval = parse_locus_interval(locus_id)
    if locus_interval is None:
        return None
    start_0based, end = locus_interval

    genotype_indices = dict(zip(fields[8].split(":"), fields[9].split(":"))).get("GT", ".").replace("|", "/").split("/")
    if len(genotype_indices) > 2:
        # e.g. ATaRVa's "multizygous" loci, where it reported more than 2 alleles
        return make_row(locus_id, extraction_status=EXTRACTION_MULTIALLELIC)

    alt_alleles = [] if fields[4] == "." else fields[4].split(",")
    if len(alt_alleles) > 2:
        return make_row(locus_id, extraction_status=EXTRACTION_MULTIALLELIC)

    # Trim the flanking bases that make the record span more than the locus interval (TRGT's leading padding base,
    # HipSTR's extended genotyping region). The padding is shared by REF and every ALT, so the same offsets apply.
    left_trim = start_0based - (pos_1based - 1)
    right_trim = (pos_1based - 1 + len(ref)) - end
    if left_trim < 0 or right_trim < 0 or left_trim + right_trim > len(ref):
        return make_row(locus_id, extraction_status=EXTRACTION_REF_MISMATCH)

    alleles = [ref[left_trim: len(ref) - right_trim]]
    if get_reference_sequence is not None and alleles[0] != get_reference_sequence(chrom, start_0based, end):
        return make_row(locus_id, extraction_status=EXTRACTION_REF_MISMATCH)

    for alt in alt_alleles:
        if alt.startswith("<"):
            # a symbolic allele: ATaRVa writes <DEL> for an allele that spans zero bases of the locus interval
            alleles.append("")
        elif len(alt) < left_trim + right_trim:
            return make_row(locus_id, extraction_status=EXTRACTION_REF_MISMATCH)
        else:
            alt = alt.upper()
            alleles.append(alt[left_trim: len(alt) - right_trim])

    sequences = []
    statuses = []
    for genotype_index in genotype_indices:
        if not genotype_index.isdigit():
            sequences.append("")
            statuses.append(NO_CALL)
        elif int(genotype_index) >= len(alleles):
            return None
        else:
            sequences.append(alleles[int(genotype_index)])
            statuses.append(CALLED)

    if len(sequences) == 1:
        # a haploid call (male chrX/chrY outside the PAR). The truth set stores these HEMI loci as a duplicated
        # diploid genotype, and add_tool_results_columns.py mirrors that for the repeat-count columns, so do the
        # same here rather than leaving Allele 2 as a no-call.
        sequences.append(sequences[0])
        statuses.append(statuses[0])

    if statuses[0] == CALLED and statuses[1] == CALLED:
        sequences.sort(key=len)
    elif statuses[1] == CALLED:
        # only one of the two alleles was called; keep it in the Allele 1 slot
        sequences.reverse()
        statuses.reverse()

    return make_row(locus_id, sequences[0], sequences[1], statuses[0], statuses[1])


def create_reference_sequence_func(reference_fasta_path):
    """Return f(chrom, start_0based, end) -> the uppercase reference sequence of that interval.

    The sequence is upper-cased because hg38 soft-masks tandem repeats in lower case, which would otherwise make
    every REF comparison in the repeat regions of interest fail.

    Args:
        reference_fasta_path (str): path of the reference genome fasta (its .fai must be next to it).

    Returns:
        callable: f(chrom, start_0based, end) -> str
    """
    import pysam

    fasta = pysam.FastaFile(reference_fasta_path)
    reference_has_chr_prefix = fasta.references[0].startswith("chr")

    def get_reference_sequence(chrom, start_0based, end):
        chrom = chrom.replace("chr", "")
        if reference_has_chr_prefix:
            chrom = f"chr{chrom}"
        return fasta.fetch(chrom, start_0based, end).upper()

    return get_reference_sequence


def process_vcfs(vcf_paths, tool, output_path, get_reference_sequence=None, show_progress_bar=False):
    """Extract allele sequences from one or more VCFs and write them to a single tsv.

    Args:
        vcf_paths (list): paths of the tool's output VCF(s). HipSTR is sharded, so it has more than one.
        tool (str): one of TOOL_CHOICES.
        output_path (str): path of the output .tsv.gz.
        get_reference_sequence (callable): optional, see extract_allele_sequences_from_vcf_record.
        show_progress_bar (bool): whether to show a tqdm progress bar.

    Returns:
        collections.Counter: counts of each ExtractionStatus, each AlleleStatus, and "malformed" / "total" records.
    """
    counters = collections.Counter()
    with gzip.open(output_path, "wt") as output_file:
        output_file.write("\t".join(OUTPUT_COLUMNS) + "\n")
        for vcf_path in vcf_paths:
            print(f"Processing {vcf_path}")
            with (gzip.open if vcf_path.endswith("gz") else open)(vcf_path, "rt") as vcf:
                if show_progress_bar:
                    from tqdm import tqdm
                    vcf = tqdm(vcf, unit=" vcf records", unit_scale=True, unit_divisor=1000)
                for line in vcf:
                    if line.startswith("#"):
                        continue
                    counters["total"] += 1
                    row = extract_allele_sequences_from_vcf_record(line, tool, get_reference_sequence)
                    if row is None:
                        counters["malformed"] += 1
                        continue
                    counters[row["ExtractionStatus"]] += 1
                    counters[row["AlleleStatus: Allele 1"]] += 1
                    counters[row["AlleleStatus: Allele 2"]] += 1
                    output_file.write("\t".join(str(row[column]) for column in OUTPUT_COLUMNS) + "\n")

    print(f"Parsed {counters['total']:,d} vcf records from {len(vcf_paths):,d} file(s):")
    print(f"{counters[EXTRACTION_OK]:12,d} loci extracted ok")
    print(f"{counters[EXTRACTION_REF_MISMATCH]:12,d} loci skipped: REF does not match the reference interval")
    print(f"{counters[EXTRACTION_MULTIALLELIC]:12,d} loci skipped: more than 2 alleles")
    print(f"{counters['malformed']:12,d} loci skipped: malformed record (no usable locus id or genotype)")
    print(f"{counters[CALLED]:12,d} alleles called, {counters[NO_CALL]:,d} alleles not called")
    print(f"Wrote {counters['total'] - counters['malformed']:,d} rows to {output_path}")

    return counters


def main():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    p.add_argument("--tool", choices=TOOL_CHOICES, required=True,
                   help="Which tool produced the input VCF(s). Selects how the locus id is read from each record.")
    p.add_argument("-R", "--reference-fasta", help="Reference genome fasta. When specified, each record's trimmed "
                   "REF is checked against it, and loci that don't match are reported as ref_mismatch instead of "
                   "being scored against a mis-aligned truth sequence.")
    p.add_argument("-o", "--output-path", help="Output .tsv.gz path. Defaults to {input vcf name}.allele_sequences.tsv.gz")
    p.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    p.add_argument("vcf_paths", nargs="+", help="Tool output VCF path(s). HipSTR is sharded, so it has more than one.")
    args = p.parse_args()

    for vcf_path in args.vcf_paths:
        if not os.path.isfile(vcf_path):
            p.error(f"{vcf_path} not found")

    output_path = args.output_path or (
        re.sub(".vcf(.gz)?$", "", os.path.basename(args.vcf_paths[0])) + ".allele_sequences.tsv.gz")

    process_vcfs(
        args.vcf_paths,
        args.tool,
        output_path,
        get_reference_sequence=create_reference_sequence_func(args.reference_fasta) if args.reference_fasta else None,
        show_progress_bar=args.show_progress_bar)


if __name__ == "__main__":
    main()
