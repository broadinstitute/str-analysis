"""This script computes statistics about the overlap between the STR truth set vcf and other STR catalogs and
genomic regions of interest.
"""

import collections
import sys
from intervaltree import Interval, IntervalTree
import tqdm

from str_analysis.utils.file_utils import open_file

GENE_MODELS = {
    "gencode": "gs://str-truth-set/hg38/ref/other/gencode.v50.annotation.gtf.gz",
    "mane": "gs://str-truth-set/hg38/ref/other/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz",
    "refseq": "gs://str-truth-set/hg38/ref/other/hg38.ncbiRefSeq.gtf.gz", # downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/
}


# default approx. max size of a gene promoter region upstream of a gene
PROMOTER_SIZE = 1000  # bp

# maps (gtf_path, promoter_size) to a dict of per-chromosome IntervalTrees. The promoter size is part of the key
# because the promoter intervals stored in the trees depend on it.
GENCODE_INTERVAL_TREES = {}

TRANSCRIPT_ID_TO_CDS_COORDS_MAP = {}

# The fields of a GTF record that compute_genomic_region_of_interval and compute_UTR_type use. The IntervalTrees store
# these instead of whole records, since a whole-genome GTF has millions of records.
CompactGtfRecord = collections.namedtuple("CompactGtfRecord", [
    "feature_type", "chrom", "start_1based", "end_1based", "strand", "gene_id", "transcript_id", "gene_name",
    "UTR_type_from_gtf"])


def convert_gtf_record_to_compact_gtf_record(record):
    """Returns the CompactGtfRecord of a record from generate_gtf_records. Its strings are interned, so that the many
    records of a gene or transcript share one copy of its ids and name."""
    return CompactGtfRecord(
        feature_type=sys.intern(record["feature_type"]),
        chrom=sys.intern(record["chrom"]),
        start_1based=record["start_1based"],
        end_1based=record["end_1based"],
        strand=sys.intern(record["strand"]),
        gene_id=sys.intern(record["gene_id"]),
        transcript_id=sys.intern(record["transcript_id"]),
        gene_name=sys.intern(record["gene_name"]),
        UTR_type_from_gtf=record["UTR_type_from_gtf"],
    )


def parse_gtf_to_interval_trees(gtf_path=GENE_MODELS["gencode"], verbose=False, show_progress_bar=False, n=None,
                                promoter_size=PROMOTER_SIZE):
    """Converts GTF file to interval trees for fast overlap queries.

    compute_genomic_region_of_interval uses, of the records of each feature type that overlap an interval, only the one
    with the lowest transcript id. So records that have the same chromosome, start, end and feature type are stored
    once, as the one with the lowest transcript id. GENCODE, for example, lists an exon once for every transcript that
    has it, so this stores about a quarter of its records. Each interval's data is a CompactGtfRecord.

    Args:
        gtf_path (str): path to GTF file. This can be a local file or a gs:// path and can be compressed or not.
        verbose (bool): if True, print verbose output
        show_progress_bar (bool): if True, show a progress bar
        n (int): number of records to process. If None, process all records.
        promoter_size (int): number of base pairs upstream of each transcript to annotate as its promoter
    """

    # parse the GTF file one record at a time
    iterator = generate_gtf_records(gtf_path, promoter_size=promoter_size)
    if verbose:
        print(f"Parsing {gtf_path}...")
    if show_progress_bar:
        iterator = tqdm.tqdm(iterator, unit=" GTF records", unit_scale=True)

    counters = collections.defaultdict(int)

    # for check_that_transcript_ids_of_unlabeled_UTR_records_have_only_one_transcript_record, which checks all records
    # even when n limits the records that are stored
    transcript_record_counts = collections.Counter()
    transcript_ids_of_unlabeled_UTR_records = set()

    # (chrom without "chr", 0-based start, end, feature type) => the CompactGtfRecord with the lowest transcript id
    record_with_lowest_transcript_id_by_interval_and_feature_type = {}

    # cache all CDS coords to use later for converting "UTR" records into "5' UTR" and "3' UTR"
    transcript_id_to_cds_coords_map = {}

    for record in iterator:
        if record["feature_type"] == "transcript":
            transcript_record_counts[record["transcript_id"]] += 1
        elif record["feature_type"] == "UTR" and record["UTR_type_from_gtf"] is None:
            transcript_ids_of_unlabeled_UTR_records.add(record["transcript_id"])

        if n is not None and counters["total"] >= n:
            continue

        key = (record["chrom"].replace("chr", ""), record["start_1based"] - 1, record["end_1based"],
               record["feature_type"])
        stored_record = record_with_lowest_transcript_id_by_interval_and_feature_type.get(key)
        if stored_record is None or record["transcript_id"] < stored_record.transcript_id:
            record_with_lowest_transcript_id_by_interval_and_feature_type[key] = convert_gtf_record_to_compact_gtf_record(
                record)

        if record["feature_type"] == "CDS" and record["transcript_id"] not in transcript_id_to_cds_coords_map:
            transcript_id_to_cds_coords_map[record["transcript_id"]] = (
                record["chrom"],
                record["start_1based"],
                record["end_1based"],
                record["strand"],
            )

        counters[record["feature_type"]] += 1
        counters[f"total"] += 1

    # check before anything is cached, so that a later call can't reuse the results of a GTF that failed this check
    check_that_transcript_ids_of_unlabeled_UTR_records_have_only_one_transcript_record(
        transcript_record_counts, transcript_ids_of_unlabeled_UTR_records, gtf_path)

    if verbose:
        print(f"Adding {len(record_with_lowest_transcript_id_by_interval_and_feature_type):,d} gtf records with distinct "
              f"intervals and feature types to IntervalTree for overlap detection...")

    # create IntervalTree for each chromosome
    intervals_by_chrom = collections.defaultdict(list)
    for (chrom, start_0based, end_1based, _), record in record_with_lowest_transcript_id_by_interval_and_feature_type.items():
        intervals_by_chrom[chrom].append(Interval(start_0based, end_1based, data=record))
    record_with_lowest_transcript_id_by_interval_and_feature_type.clear()

    interval_trees = collections.defaultdict(IntervalTree)
    for chrom in list(intervals_by_chrom):
        interval_trees[chrom] = IntervalTree(intervals_by_chrom.pop(chrom))

    GENCODE_INTERVAL_TREES[(gtf_path, promoter_size)] = interval_trees
    TRANSCRIPT_ID_TO_CDS_COORDS_MAP[gtf_path] = transcript_id_to_cds_coords_map

    if verbose:
        for key, count in sorted(counters.items(), key=lambda s: s[1], reverse=True):
            print(f"{count:10,d} {key}")


def check_that_transcript_ids_of_unlabeled_UTR_records_have_only_one_transcript_record(
        transcript_record_counts, transcript_ids_of_unlabeled_UTR_records, gtf_path):
    """Raise an error if a UTR record that isn't labeled as 5UTR or 3UTR in the GTF (as in GENCODE GTFs) has a
    transcript id that appears in more than one transcript record. compute_UTR_type labels such UTRs as 5' or 3' by
    comparing them with the first CDS record of their transcript id, so if several copies of a transcript share an id,
    it would compare the UTRs of the other copies with the wrong CDS.

    Args:
        transcript_record_counts (collections.Counter): the number of transcript records of each transcript id
        transcript_ids_of_unlabeled_UTR_records (set): the transcript ids of the UTR records that aren't labeled as
            5UTR or 3UTR in the GTF
        gtf_path (str): path of the GTF file the records came from, used in the error message

    Raises:
        ValueError: if any such transcript id is found
    """
    transcript_ids_with_more_than_one_transcript_record = sorted(
        transcript_id for transcript_id in transcript_ids_of_unlabeled_UTR_records
        if transcript_record_counts[transcript_id] > 1)
    if transcript_ids_with_more_than_one_transcript_record:
        raise ValueError(
            f"{len(transcript_ids_with_more_than_one_transcript_record):,d} transcript id(s) in {gtf_path} have "
            f"UTR records that aren't labeled as 5UTR or 3UTR, and also appear in more than one transcript record, so "
            f"these UTRs can't be labeled as 5' or 3' by comparing them with the transcript's CDS. Examples: "
            f"{', '.join(transcript_ids_with_more_than_one_transcript_record[:10])}")


def generate_gtf_records(gtf_path, promoter_size=PROMOTER_SIZE):
    """Parse the GTF file. For each transcript record, also yield a "promoter" record that covers the promoter_size
    base pairs upstream of the transcript.

    Args:
        gtf_path (str): path to GTF file. This can be a local file or a gs:// path and can be compressed or not.
        promoter_size (int): number of base pairs upstream of each transcript to annotate as its promoter. If 0, no
            promoter records are generated.

    Yields:
        dict: one record per GTF line with a relevant feature type, plus one promoter record per transcript
    """

    if promoter_size < 0:
        raise ValueError(f"promoter_size must be >= 0: {promoter_size}")

    gtf_file = open_file(gtf_path, is_text_file=True, download_local_copy_before_opening=True)
    for i, line in enumerate(gtf_file):
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        feature_type = fields[2]
        if feature_type not in {
            "transcript",
            "exon",
            "CDS",
            "UTR",
            "5UTR",
            "3UTR",
            "promoter",
        }:
            if feature_type.lower() not in {"gene", "start_codon", "stop_codon", "selenocysteine"}:
                print(f"WARNING: Unexpected feature type '{feature_type}' in {gtf_path} line #{i}: {line.strip()}")
            continue

        # GTFs such as UCSC's label each UTR as a 5UTR or 3UTR. Keep that label, but treat these records as UTR records
        # so that they're picked the same way as the UTR records of GENCODE GTFs.
        UTR_type_from_gtf = {"5UTR": "5' UTR", "3UTR": "3' UTR"}.get(feature_type)
        if UTR_type_from_gtf is not None:
            feature_type = "UTR"

        annotation_source = fields[1]
        chrom = fields[0]
        start_1based = int(fields[3])
        end_1based = int(fields[4])
        if end_1based <= start_1based:
            continue

        meta_fields = {}
        for meta_field in fields[8].strip("; ").split(";"):
            key, value = meta_field.strip().replace('"', '').split()
            meta_fields[key] = value

        strand = fields[6]

        gene_id = meta_fields["gene_id"].split(".")[0]
        transcript_id = meta_fields["transcript_id"].split(".")[0]

        record = {
            "feature_type": feature_type,
            "chrom": chrom,
            "start_1based": start_1based,
            "end_1based": end_1based,
            "annotation_source": annotation_source,
            "strand": strand,
            "gene_id": gene_id,
            "transcript_id": transcript_id,
            "gene_name": meta_fields["gene_name"],
            "gene_type": meta_fields.get("transcript_type"),
            "transcript_type":  meta_fields.get("transcript_type"),
            "UTR_type_from_gtf": UTR_type_from_gtf,
        }

        yield record

        # create a separate record for the transcript's promoter
        if feature_type == "transcript":
            if strand == "+":
                promoter_start_1based = max(1, start_1based - promoter_size)
                promoter_end_1based = start_1based - 1
            elif strand == "-":
                promoter_start_1based = end_1based + 1
                promoter_end_1based = end_1based + promoter_size
            else:
                raise ValueError(f"Unexpected strand value in {gtf_path} line #{i}: {line}")

            if promoter_end_1based <= promoter_start_1based:
                continue

            promoter_record = dict(record)
            promoter_record.update({
                "feature_type": "promoter",
                "start_1based": promoter_start_1based,
                "end_1based": promoter_end_1based,
            })

            yield promoter_record

    gtf_file.close()


def compute_UTR_type(utr_record, transcript_id_to_cds_coords):

    if utr_record["feature_type"] != "UTR":
        raise ValueError(f"Expected a UTR record, but got a {utr_record['feature_type']} record: {utr_record}")

    # Use the GTF's own 5UTR or 3UTR label when it has one. The lookup below takes the first CDS record of the
    # transcript id, so it compares the UTR with the wrong CDS when several copies of a gene share a transcript id.
    if utr_record.get("UTR_type_from_gtf"):
        return utr_record["UTR_type_from_gtf"]

    cds_coords = transcript_id_to_cds_coords.get(utr_record["transcript_id"])
    if cds_coords is None:
        print(f"WARNING: CDS not found for {utr_record['transcript_id']}. Unable to determine if UTR is 5' or 3'")
        return "UTR"

    cds_chrom, cds_start_1based, cds_end_1based, cds_strand = cds_coords
    if utr_record["chrom"] != cds_chrom:
        print(f"WARNING: {utr_record['transcript_id']} chromosome in CDS record != chromosome in UTR record: {cds_chrom} vs {utr_record['chrom']}. Unable to determine if UTR is 5' or 3'")
        return "UTR"

    if utr_record["strand"] != cds_strand:
        print(f"WARNING: {utr_record['transcript_id']} strand in CDS record != strand in UTR record: {cds_strand} vs {utr_record['strand']}")
        return "UTR"

    if (utr_record["strand"] == "+" and utr_record["end_1based"] < cds_start_1based) \
            or (utr_record["strand"] == "-" and utr_record["start_1based"] > cds_end_1based):
        return "5' UTR"
    elif (utr_record["strand"] == "-" and utr_record["end_1based"] < cds_start_1based) \
            or (utr_record["strand"] == "+" and utr_record["start_1based"] > cds_end_1based):
        return "3' UTR"
    else:
        print(f"WARNING: Something wrong with UTR info: {utr_record['start_1based']}-{utr_record['end_1based']} {utr_record['strand']}, "
              f"or CDS info: {cds_coords}. Unable to determine if UTR is 5' or 3'")
        return "UTR"


def compute_genomic_region_of_interval(chrom, start_1based, end_1based, genes_gtf_path=None, verbose=False, show_progress_bar=False, n=None, promoter_size=PROMOTER_SIZE):
    """This method computes the genomic region of a given interval. The first time it is called for a given GTF file
    and promoter size, it will parse the GTF file and create an interval tree for fast overlap queries. Subsequent
    calls with the same GTF file and promoter size will reuse the interval tree.

    Args:
        chrom (str): chromosome
        start_1based (int): start position (1-based)
        end_1based (int): end position (1-based)
        genes_gtf_path (str): path to GTF file. This can be a local file or a gs:// path and can be compressed or not.
        verbose (bool): if True, print verbose output
        show_progress_bar (bool): if True, show a progress bar
        n (int): number of records to process. If None, process all records.
        promoter_size (int): number of base pairs upstream of each transcript to annotate as its promoter

    Returns: 4-tuple containing
        region ("CDS", "5' UTR", "3' UTR", "UTR", "exon", "intron", or "promoter")
        gene name         (or None if region is not in a transcript)
        gene id           (or None if region is not in a transcript)
        transcript id     (or None if region is not in a transcript)
    """

    if genes_gtf_path is None:
        genes_gtf_path = GENE_MODELS["gencode"]

    if (genes_gtf_path, promoter_size) not in GENCODE_INTERVAL_TREES:
        parse_gtf_to_interval_trees(genes_gtf_path, verbose=verbose, show_progress_bar=show_progress_bar, n=n,
                                    promoter_size=promoter_size)

    interval_trees = GENCODE_INTERVAL_TREES[(genes_gtf_path, promoter_size)]
    transcript_id_to_cds_coords_map = TRANSCRIPT_ID_TO_CDS_COORDS_MAP[genes_gtf_path]

    input_chrom = str(chrom).replace("chr", "")
    input_locus_interval = Interval(start_1based - 1, end_1based)

    # The given interval may overlap more than one region in the GTF. Return the most biologically important region.
    overlapping_intervals_from_gtf = interval_trees[input_chrom].overlap(input_locus_interval)
    overlapping_intervals_feature_types = {}
    for i in sorted(overlapping_intervals_from_gtf, key=lambda i: i.data.transcript_id):
        feature_type = i.data.feature_type
        if feature_type not in overlapping_intervals_feature_types:
            overlapping_intervals_feature_types[feature_type] = i.data

    if "transcript" in overlapping_intervals_feature_types and "exon" not in overlapping_intervals_feature_types:
        overlapping_intervals_feature_types["intron"] = overlapping_intervals_feature_types["transcript"]

    for feature_type in "CDS", "5' UTR", "3' UTR", "UTR", "exon", "intron", "promoter":  # handle both coding and non-coding transcripts
        if feature_type not in overlapping_intervals_feature_types:
            continue
        record = overlapping_intervals_feature_types[feature_type]
        if feature_type == "UTR":
            feature_type = compute_UTR_type(record._asdict(), transcript_id_to_cds_coords_map)
        return feature_type, record.gene_name, record.gene_id, record.transcript_id

    return "intergenic", None, None, None




