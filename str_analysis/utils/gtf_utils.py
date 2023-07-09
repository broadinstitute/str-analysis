"""This script computes statistics about the overlap between the STR truth set vcf and other STR catalogs and
genomic regions of interest.
"""

import collections
from intervaltree import Interval, IntervalTree
import tqdm

from str_analysis.utils.file_utils import open_file

GENE_MODELS = {
    "gencode": "gs://str-truth-set/hg38/ref/other/gencode.v43.annotation.gtf.gz",
    "mane": "gs://str-truth-set/hg38/ref/other/MANE.GRCh38.v1.0.ensembl_genomic.gtf.gz",
}


# approx. max size of a gene promoter region upstream of a gene
PROMOTER_SIZE = 1000  # bp

GENCODE_INTERVAL_TREES = {}

TRANSCRIPT_ID_TO_CDS_COORDS_MAP = {}

def parse_gtf_to_interval_trees(gtf_path=GENE_MODELS["gencode"], verbose=False):
    """Converts GTF file to interval trees for fast overlap queries"""

    # parse the GTF file into IntervalTrees
    iterator = generate_gtf_records(gtf_path)
    if verbose:
        print(f"Parsing {gtf_path}...")
        iterator = tqdm.tqdm(iterator, unit=" GTF records", unit_scale=True)

    counters = collections.defaultdict(int)
    gtf_records = list(iterator)

    # create IntervalTree for each chromosome
    interval_trees = GENCODE_INTERVAL_TREES[gtf_path] = collections.defaultdict(IntervalTree)

    # cache all CDS coords to use later for converting "UTR" records into "5' UTR" and "3' UTR"
    transcript_id_to_cds_coords_map = TRANSCRIPT_ID_TO_CDS_COORDS_MAP[gtf_path] = {}

    if verbose:
        print(f"Adding gtf records to IntervalTree for overlap detection...")
        iterator = tqdm.tqdm(gtf_records, unit=" GTF records", unit_scale=True)
    else:
        iterator = gtf_records

    for record in iterator:
        chrom = record["chrom"].replace("chr", "")
        interval_trees[chrom].add(Interval(
            record["start_1based"] - 1,
            record["end_1based"],
            data=record,
        ))

        if record["feature_type"] == "CDS" and record["transcript_id"] not in transcript_id_to_cds_coords_map:
            transcript_id_to_cds_coords_map[record["transcript_id"]] = (
                record["chrom"],
                record["start_1based"],
                record["end_1based"],
                record["strand"],
            )

        counters[record["feature_type"]] += 1
        counters[f"total"] += 1

    if verbose:
        for key, count in sorted(counters.items(), key=lambda s: s[1], reverse=True):
            print(f"{count:10,d} {key}")


def generate_gtf_records(gtf_path):
    """Parse the GTF file"""

    gtf_file = open_file(gtf_path)
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
            "promoter",
        }:
            continue

        annotation_source = fields[1]
        chrom = fields[0]
        start_1based = int(fields[3])
        end_1based = int(fields[4])

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
            "gene_type": meta_fields["transcript_type"],
            "transcript_type":  meta_fields["transcript_type"],
        }

        yield record

        # create a separate record for the transcript's promoter
        if feature_type == "transcript":
            if strand == "+":
                promoter_start_1based = max(1, start_1based - PROMOTER_SIZE)
                promoter_end_1based = start_1based - 1
            elif strand == "-":
                promoter_start_1based = end_1based + 1
                promoter_end_1based = end_1based + PROMOTER_SIZE
            else:
                raise ValueError(f"Unexpected strand value in {gtf_path} line #{i}: {line}")

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
        raise utr_record["feature_type"]

    cds_coords = transcript_id_to_cds_coords.get(utr_record["transcript_id"])
    if cds_coords is None:
        print("WARNING: CDS not found for", utr_record["transcript_id"])
        print(utr_record)
        return "UTR"

    cds_chrom, cds_start_1based, cds_end_1based, cds_strand = cds_coords
    if utr_record["chrom"] != cds_chrom:
        print("ERROR:", utr_record["transcript_id"], "chrom in CDS utr_record != chrom in UTR utr_record:", cds_chrom, "vs", utr_record["chrom"])
        return "UTR"

    if utr_record["strand"] != cds_strand:
        print("ERROR:", utr_record["transcript_id"], "strand in CDS utr_record != strand in UTR utr_record:", cds_strand, "vs", utr_record["strand"])
        return "UTR"

    if (utr_record["strand"] == "+" and utr_record["end_1based"] < cds_start_1based) \
            or (utr_record["strand"] == "-" and utr_record["start_1based"] > cds_end_1based):
        return "5' UTR"
    elif (utr_record["strand"] == "-" and utr_record["end_1based"] < cds_start_1based) \
            or (utr_record["strand"] == "+" and utr_record["start_1based"] > cds_end_1based):
        return "3' UTR"
    else:
        print("ERROR: Something wrong with UTR info:", utr_record["strand"], utr_record["start_1based"],  utr_record["end_1based"],
              "or CDS info:", cds_coords)
        return "UTR"


def compute_genomic_region_of_interval(chrom, start_1based, end_1based, genes_gtf_path=None, verbose=False):
    """This method computes the genomic region of a given interval. The first time it is called, it will parse the
    GTF file and create an interval tree for fast overlap queries. Subsequent calls will reuse the interval tree.

    Args:
        chrom (str): chromosome
        start_1based (int): start position (1-based)
        end_1based (int): end position (1-based)
        genes_gtf_path (str): path to GTF file. This can be a local file or a gs:// path and can be compressed or not.
        verbose (bool): if True, print verbose output

    Returns: 4-tuple containing
        region ("CDS", "5' UTR", "3' UTR", "UTR", "exon", "intron", or "promoter")
        gene name         (or None if region is not in a transcript)
        gene id           (or None if region is not in a transcript)
        transcript id     (or None if region is not in a transcript)
    """

    if genes_gtf_path is None:
        genes_gtf_path = GENE_MODELS["gencode"]

    if genes_gtf_path not in GENCODE_INTERVAL_TREES:
        parse_gtf_to_interval_trees(genes_gtf_path, verbose=verbose)

    interval_trees = GENCODE_INTERVAL_TREES[genes_gtf_path]
    transcript_id_to_cds_coords_map = TRANSCRIPT_ID_TO_CDS_COORDS_MAP[genes_gtf_path]

    input_chrom = str(chrom).replace("chr", "")
    input_locus_interval = Interval(start_1based - 1, end_1based)

    # The truth set STR may overlap more than one region in the GTF. Return the most biologically important region.
    overlapping_intervals_from_gtf = interval_trees[input_chrom].overlap(input_locus_interval)
    overlapping_intervals_feature_types = {}
    for i in sorted(overlapping_intervals_from_gtf, key=lambda i: i.data["transcript_id"]):
        feature_type = i.data["feature_type"]
        if feature_type not in overlapping_intervals_feature_types:
            overlapping_intervals_feature_types[feature_type] = i.data

    if "transcript" in overlapping_intervals_feature_types and "exon" not in overlapping_intervals_feature_types:
        overlapping_intervals_feature_types["intron"] = overlapping_intervals_feature_types["transcript"]

    for feature_type in "CDS", "5' UTR", "3' UTR", "UTR", "exon", "intron", "promoter":  # handle both coding and non-coding transcripts
        if feature_type not in overlapping_intervals_feature_types:
            continue
        record = overlapping_intervals_feature_types[feature_type]
        if feature_type == "UTR":
            feature_type = compute_UTR_type(overlapping_intervals_feature_types["UTR"], transcript_id_to_cds_coords_map)
        return feature_type, record["gene_name"], record["gene_id"], record["transcript_id"]

    return "intergenic", None, None, None




