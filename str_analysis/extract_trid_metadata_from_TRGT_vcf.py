"""Extract per-record interval metadata from a TRGT multisample VCF.

For each VCF record × motif, emits one row per matching LocusId:

    trid  locus_id  motif  interval  vc

where ``trid`` is the full ``INFO/TRID`` (a single LocusId, a comma-separated
list of LocusIds, or ``VC:chrom:start-end`` for a cluster under the newer
convention), ``locus_id`` is one of the LocusIds the record covers whose motif
suffix matches the row's ``motif`` -- read from ``INFO/TRID``, or from
``INFO/STRUC`` when the TRID names the cluster itself (so a compound TRID listing three LocusIds that end in ``-AGAA``
with motif ``AGAA`` emits three rows, one per LocusId, all sharing the same
``interval``/``vc``), ``interval`` is ``{chrom_no_chr}:{POS}-{END}``, and
``vc`` is the variation cluster's own span in that same form (e.g.
``13:102161564-102161724``) or ``""`` for an isolated TR row. A cluster is
recognized under any of the catalog conventions in circulation: ``STRUC``
starting with ``<VC``, which covers both the bare counter ``<VC1>`` and
``<VC:chrom:start-end>``, or ``TRID`` starting with ``VC:``.

Rows from the same VCF record share the same ``(trid, motif, interval, vc)``
and appear consecutively in the output, so downstream consumers can group
them back into per-VCF-record chunks.

Output is a small gzipped TSV (91 MB for the 2026-05 HPRC256 VCF)
intended to be loaded into memory by downstream LPS processing in lieu of
re-parsing the 16 GB VCF.

Parallelism: one worker per chromosome (using ``tabix`` to fetch records).
By default every contig listed in the VCF's tabix index is extracted, including
alt contigs, and the output follows the order of the index (``tabix -l``), whether
or not the contig names carry a ``chr`` prefix. ``--contig`` restricts the run to the named contigs,
which is what lets a per-contig pipeline shard run this against the same
subset VCF its other per-contig steps read; the output then follows the order
the flags were given.

Gene-named repeat ids: a TRGT catalog built from the known disease loci catalog can name repeats
after genes rather than coordinates, e.g. ``TMEM185A_CGCCGT`` or ``HTT``. Those can't be matched to
a motif by their ending, so they're first replaced with the chrom-start-end-motif LocusId that
``--known-loci-catalog`` defines for them (see ``load_known_locus_ids``).

Caching: if the output exists and its mtime is newer than the input VCF's mtime,
this script's own mtime and the known loci catalog's mtime, the script skips
extraction. Pass ``--force`` to override.
"""

import argparse
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

DEFAULT_KNOWN_LOCI_CATALOG = Path(__file__).parent / "variant_catalogs" / "variant_catalog_without_offtargets.GRCh38.json"


def load_known_locus_ids(catalog_path):
    """Returns a dict mapping each repeat id in a known disease loci catalog to its chrom-start-end-motif LocusId.

    Catalogs like variant_catalog_without_offtargets.GRCh38.json name a single-repeat locus by its
    LocusId (e.g. ``EP400``) and each repeat of a multi-repeat locus by its VariantId (e.g.
    ``TMEM185A_CGCCGT`` and ``TMEM185A``). A TRGT catalog built from one carries those names into
    INFO/TRID. Each is mapped to the LocusId built from that repeat's own ReferenceRegion and the
    motif its LocusStructure gives it, so ``TMEM185A_CGCCGT`` becomes ``X-149631723-149631735-CGCCGT``.

    Args:
        catalog_path (str): path of an ExpansionHunter-style variant catalog JSON (optionally gzipped).

    Returns:
        dict: repeat id -> chrom-start-end-motif LocusId, with no "chr" prefix.

    Raises:
        ValueError: if a record's VariantIds, ReferenceRegions and LocusStructure motifs don't line up.
    """
    fopen = gzip.open if str(catalog_path).endswith((".gz", ".bgz")) else open
    with fopen(catalog_path, "rt") as f:
        records = json.load(f)

    known_locus_ids = {}
    for record in records:
        motifs = re.findall(r"[(]([A-Z]+)[)]", record["LocusStructure"])
        reference_regions = record["ReferenceRegion"]
        if not isinstance(reference_regions, list):
            reference_regions = [reference_regions]
        repeat_ids = record.get("VariantId") or [record["LocusId"]]
        if not isinstance(repeat_ids, list):
            repeat_ids = [repeat_ids]
        if not len(motifs) == len(reference_regions) == len(repeat_ids):
            raise ValueError(f"{catalog_path} record {record['LocusId']} has {len(repeat_ids)} VariantIds, "
                             f"{len(reference_regions)} ReferenceRegions and {len(motifs)} LocusStructure motifs")
        for repeat_id, reference_region, motif in zip(repeat_ids, reference_regions, motifs):
            chrom, start_0based, end_1based = re.fullmatch(r"(.+):(\d+)-(\d+)", reference_region).groups()
            known_locus_ids[repeat_id] = f"{_strip_chr(chrom)}-{start_0based}-{end_1based}-{motif}"

    return known_locus_ids


def list_vcf_contigs(vcf_path):
    """Returns the contigs that have records in a VCF, in the order its tabix index lists them.

    Args:
        vcf_path (str): path of a bgzipped, tabix-indexed VCF.

    Returns:
        list: contig names, exactly as the VCF writes them.
    """
    return subprocess.run(["tabix", "-l", vcf_path], capture_output=True, text=True, check=True).stdout.split()


def parse_struc_vc_span(struc, trid, chrom, vcf_start_0based, vcf_end_1based):
    """Returns the VC span for a variation cluster record, or ``""`` for an isolated TR row.

    Three catalog conventions are in circulation. The v1.0.1 variation_clusters_and_isolated_TRs
    catalog, the one the gnomAD-LR pipeline ran TRGT against, puts the comma-joined repeat ids in
    the TRID and writes ``STRUC`` as a bare counter, ``<VC1>``, ``<VC2>``, which carries no span at
    all. A later convention keeps the ids in the TRID but writes the span into ``STRUC``,
    ``<VC:chrom:start-end>``. The newest one, which fixes the duplicate-TRID problem in
    https://github.com/PacificBiosciences/trgt-lps/issues/5, gives the cluster its own TRID
    ``VC:chrom:start-end`` and moves the repeat ids into ``STRUC``.

    All three put the span on the record itself, so derive it from POS and END rather than from
    whichever field happens to carry it. That leaves only recognition, and testing ``STRUC`` for
    the ``<VC`` prefix covers the first two forms at once: nothing else in these catalogs starts
    with it, since isolated repeats write ``<TR1>`` or ``<TR:...>`` and vamos writes ``(CTT)n``.
    """
    if not struc.startswith("<VC") and not trid.startswith("VC:"):
        return ""
    return f"{_strip_chr(chrom)}:{vcf_start_0based}-{vcf_end_1based}"


def parse_constituent_locus_ids(trid, struc):
    """Returns the ids of the repeats a record covers, under any of the catalog conventions.

    A cluster written in the newer convention names itself in TRID and lists its repeats in
    ``STRUC``; every other row lists them in TRID.
    """
    if trid.startswith("VC:") and struc.startswith("<VC:"):
        return struc[4:-1].split(",") if struc.endswith(">") else struc[4:].split(",")
    return trid.split(",")


def _strip_chr(chrom):
    return chrom[3:] if chrom.startswith("chr") else chrom


def no_sample_has_first_allele_called(fmt_fields, sample_fields):
    """Returns True when no sample's GT has a called first allele.

    Matches the records that trgt-lps drops from its LPS output. Tested against trgt-lps 0.11.0, it
    keeps a record only when at least one sample's first allele is called: a record where every
    sample is ``.``, ``./.``, ``.|.`` or ``./0`` gets no LPS row, while one where a sample is
    ``0/.`` does. The extract step must skip the same records to keep ``--vcf-trid-metadata-tsv``
    aligned with the LPS table.

    Args:
        fmt_fields (list): the FORMAT column split on ":".
        sample_fields (list): the sample columns.

    Returns:
        bool: True if the record should be skipped.
    """
    if "GT" not in fmt_fields:
        # No GT field in this record; can't decide, so keep the record.
        return False
    gt_index = fmt_fields.index("GT")
    for s in sample_fields:
        sample_values = s.split(":")
        gt = sample_values[gt_index] if gt_index < len(sample_values) else ""
        if gt.replace("|", "/").split("/")[0] not in ("", "."):
            return False
    return True


def parse_vcf_line(line, known_locus_ids=None):
    """Parses a non-header VCF line and yields ``(trid, locus_id, motif, interval, vc)`` rows.

    For each unique motif in ``INFO/MOTIFS``, emits one row per LocusId the record
    covers whose suffix matches ``-{motif}``. Those LocusIds come from ``INFO/TRID``,
    or from ``INFO/STRUC`` when the TRID names the cluster itself; see
    ``parse_constituent_locus_ids``. A gene-named id found in ``known_locus_ids`` is replaced with
    its chrom-start-end-motif LocusId first, and a LocusId that no motif matches gets a warning,
    since it gets no row. Standalone TRs emit one
    row (the TRID is a single LocusId); compound TRGT rows can emit multiple
    rows per motif when several of the listed LocusIds share that motif.

    Records where no sample has a called first allele are skipped to match trgt-lps's LPS output
    cardinality (see ``no_sample_has_first_allele_called``).

    Returns an empty tuple if the line is malformed, has no TRID or MOTIFS, is skipped by the rule
    above, has no END, or no LocusId in the TRID matches any motif. A record skipped for having no
    END gets a warning, since trgt-lps still writes LPS rows for it, which will then have no
    metadata.

    Args:
        line (str): a VCF record line.
        known_locus_ids (dict): gene-named repeat id -> chrom-start-end-motif LocusId, as returned by
            ``load_known_locus_ids``, or None to use the ids as they are.

    Returns:
        list or tuple: the rows, or an empty tuple when the record gets none.
    """
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 8:
        return ()
    chrom = fields[0]
    pos = fields[1]
    info = fields[7]

    trid = None
    motifs = None
    end = None
    struc = ""
    for kv in info.split(";"):
        if kv.startswith("TRID="):
            trid = kv[5:]
        elif kv.startswith("MOTIFS="):
            motifs = kv[7:]
        elif kv.startswith("END="):
            end = kv[4:]
        elif kv.startswith("STRUC="):
            struc = kv[6:]

    if trid is None or motifs is None:
        return ()

    if len(fields) >= 10 and no_sample_has_first_allele_called(fields[8].split(":"), fields[9:]):
        return ()

    if end is None:
        print(f"WARNING: skipping the VCF record at {chrom}:{pos} with TRID={trid} because it has no END "
              f"in its INFO field")
        return ()
    try:
        vcf_start_0based = int(pos)
        vcf_end_1based = int(end)
    except ValueError:
        return ()

    interval = f"{_strip_chr(chrom)}:{vcf_start_0based}-{vcf_end_1based}"
    vc = parse_struc_vc_span(struc, trid, chrom, vcf_start_0based, vcf_end_1based)

    known_locus_ids = known_locus_ids or {}
    locus_ids = [known_locus_ids.get(locus_id, locus_id) for locus_id in parse_constituent_locus_ids(trid, struc)]

    rows = []
    for motif in dict.fromkeys(motifs.split(",")):
        suffix = f"-{motif}"
        for locus_id in locus_ids:
            if locus_id.endswith(suffix):
                rows.append((trid, locus_id, motif, interval, vc))

    matched_locus_ids = {row[1] for row in rows}
    unmatched_locus_ids = [locus_id for locus_id in locus_ids if locus_id not in matched_locus_ids]
    if unmatched_locus_ids:
        print(f"WARNING: the VCF record at {chrom}:{pos} with TRID={trid} has repeat id(s) "
              f"{', '.join(unmatched_locus_ids)} that no motif in MOTIFS={motifs} matches, so they get no "
              f"metadata row")
    return rows


def extract_chrom_to_tempfile(vcf_path, chrom, tmpdir, known_locus_ids):
    """Streams VCF records for one chromosome via ``tabix`` and writes a per-chrom TSV.

    Also verifies that each ``(TRID, LocusId, Interval)`` tuple appears at most once
    on this chromosome (a global uniqueness invariant: LocusId/Interval both
    encode chrom, so per-chrom uniqueness implies global uniqueness). Two records that cover the
    same repeat over the same span under different TRIDs both get their rows, since their TRIDs
    tell them apart; the converter decides what to do with them when it writes output rows.

    Returns ``(chrom, tmp_path, row_count)``.
    """
    tmp_path = os.path.join(tmpdir, f"{chrom}.tsv")
    rows_written = 0
    seen_trid_locus_interval = set()
    proc = subprocess.Popen(
        ["tabix", vcf_path, chrom],
        stdout=subprocess.PIPE,
        text=True,
        bufsize=1 << 20,
    )
    try:
        with open(tmp_path, "wt") as out:
            for line in proc.stdout:
                for row in parse_vcf_line(line, known_locus_ids):
                    trid, locus_id, _, interval, _ = row
                    key = (trid, locus_id, interval)
                    if key in seen_trid_locus_interval:
                        raise RuntimeError(
                            f"Duplicate (TRID, LocusId, Interval) on {chrom}: TRID={trid!r} "
                            f"LocusId={locus_id!r} Interval={interval!r} — the VCF contains the same "
                            f"record twice."
                        )
                    seen_trid_locus_interval.add(key)
                    out.write("\t".join(row) + "\n")
                    rows_written += 1
    finally:
        # Ensure the tabix subprocess is reaped even when an exception fires
        # mid-stream (e.g. the duplicate-key RuntimeError above).
        if proc.stdout is not None:
            proc.stdout.close()
        rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"tabix for {chrom} returned non-zero exit code {rc}")

    return chrom, tmp_path, rows_written


def output_is_fresh(output_path, vcf_path, known_loci_catalog):
    """Returns True if output_path exists and is newer than vcf_path, known_loci_catalog and this script.

    The VCF is a fixed published artifact, so its mtime alone would never invalidate the
    cache; an edit to the parsing below (which convention a cluster is recognized by, which
    LocusIds a record expands to, which records are skipped) would leave a TSV produced by
    the previous version of this code looking current, and the pipeline would join against
    it. Include this file's own mtime so changing the parser re-extracts, and the known loci
    catalog's so changing which LocusIds gene-named ids map to does too.

    The check does not look at --contig, so re-running with a different contig set but the same
    --output-tsv would be skipped. Give each contig its own output path (which is what a sharded
    pipeline does anyway), or pass --force.
    """
    if not os.path.isfile(output_path):
        return False
    newest_input = max(os.path.getmtime(p) for p in (vcf_path, __file__, known_loci_catalog))
    return os.path.getmtime(output_path) >= newest_input


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input-vcf", required=True, help="TRGT multisample VCF (gzipped + tabix-indexed).")
    parser.add_argument("--output-tsv", required=True, help="Output TSV.gz path.")
    parser.add_argument("--contig", action="append", metavar="CONTIG",
                        help="Restrict extraction to this contig. May be repeated. Defaults to "
                             "every contig in the VCF's tabix index. Duplicates are ignored, and the "
                             "output follows the order given.")
    parser.add_argument("--workers", type=int, default=8, help="Number of parallel chromosome workers.")
    parser.add_argument("--known-loci-catalog", default=str(DEFAULT_KNOWN_LOCI_CATALOG),
                        help="Variant catalog JSON that defines the gene-named repeat ids a TRGT catalog may use "
                             "(e.g. TMEM185A_CGCCGT or HTT), used to replace them with chrom-start-end-motif "
                             "LocusIds. Point this at the version the VCF's TRGT catalog was built from if it "
                             "differs from the one bundled with str_analysis.")
    parser.add_argument("--force", action="store_true", help="Re-extract even if the output is newer than the VCF.")
    args = parser.parse_args()

    if not os.path.isfile(args.input_vcf):
        parser.error(f"--input-vcf does not exist: {args.input_vcf}")
    tbi = args.input_vcf + ".tbi"
    if not os.path.isfile(tbi):
        parser.error(f"Tabix index not found: {tbi}; run `tabix -p vcf {args.input_vcf}` first.")
    if not os.path.isfile(args.known_loci_catalog):
        parser.error(f"--known-loci-catalog does not exist: {args.known_loci_catalog}")

    if not args.force and output_is_fresh(args.output_tsv, args.input_vcf, args.known_loci_catalog):
        print(f"Output {args.output_tsv} is newer than {args.input_vcf}; skipping extraction (--force to override).")
        return

    # Preserve the caller's order rather than sorting into canonical order: a caller that named
    # its contigs asked for those, and every consumer either filters by contig or relies only on
    # rows from one VCF record staying consecutive, which holds within a chromosome either way.
    vcf_contigs = list_vcf_contigs(args.input_vcf)
    if args.contig:
        chroms = list(dict.fromkeys(args.contig))
        # Only a warning, since a shard's contig can legitimately have no records (e.g. chrM), but
        # it also catches a naming mismatch such as passing "chr1" for a VCF that writes "1".
        contigs_without_records = [c for c in chroms if c not in vcf_contigs]
        if contigs_without_records:
            print(f"WARNING: {args.input_vcf} has no records on --contig {', '.join(contigs_without_records)}")
    else:
        chroms = vcf_contigs

    known_locus_ids = load_known_locus_ids(args.known_loci_catalog)
    print(f"Loaded {len(known_locus_ids):,d} gene-named repeat ids from {args.known_loci_catalog}")

    print(f"Extracting interval metadata from {args.input_vcf} -> {args.output_tsv} "
          f"({len(chroms)} contig(s), workers={args.workers})")
    tmpdir = tempfile.mkdtemp(prefix="vcf_interval_metadata_")
    chrom_results = {}
    # Atomic write: stream into a tmp path next to the destination, then
    # os.replace at the end. Prevents a Ctrl-C / OOM / network glitch mid-
    # concat from leaving a partial gzip whose mtime is fresher than the VCF
    # (which would silently short-circuit future runs via output_is_fresh).
    tmp_output = args.output_tsv + ".tmp"
    # Not a `with` block: its __exit__ shuts down without cancel_futures, so one worker's
    # failure would still cost the full wall time of every chromosome queued behind it.
    ex = ThreadPoolExecutor(max_workers=args.workers)
    futures = {}
    try:
        try:
            futures = {
                ex.submit(extract_chrom_to_tempfile, args.input_vcf, chrom, tmpdir, known_locus_ids): chrom
                for chrom in chroms
            }
            for fut in as_completed(futures):
                chrom, path, count = fut.result()
                chrom_results[chrom] = (path, count)
                print(f"  {chrom}: {count:,d} rows -> {path}")
        finally:
            # Drops the chromosomes still queued; the ones already running finish. Cancelling
            # each future by hand rather than passing shutdown(cancel_futures=True), which needs
            # Python 3.9 while setup.py still declares support back to 3.7.
            for fut in futures:
                fut.cancel()
            ex.shutdown(wait=True)

        # Concatenate the per-chrom temp files into the tmp output in the same order they were
        # fanned out in, so the output order is the one the module docstring promises.
        total = 0
        with gzip.open(tmp_output, "wt") as out:
            out.write("trid\tlocus_id\tmotif\tinterval\tvc\n")
            for chrom in chroms:
                if chrom not in chrom_results:
                    continue
                path, count = chrom_results[chrom]
                with open(path) as f:
                    for line in f:
                        out.write(line)
                total += count
        # Promote tmp -> final atomically; do this AFTER the gzip writer
        # closed cleanly so the final path is never a partial gzip.
        os.replace(tmp_output, args.output_tsv)
        print(f"Wrote {total:,d} rows to {args.output_tsv}")
    finally:
        # Clean up the tmp output if it's still around (e.g. exception fired
        # before os.replace).
        if os.path.isfile(tmp_output):
            try:
                os.remove(tmp_output)
            except OSError:
                pass
        # rmtree rather than removing the paths in chrom_results: on the error path that
        # dict stopped being filled at the raise, so workers still running at that moment
        # wrote temp files it never recorded.
        shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    main()
