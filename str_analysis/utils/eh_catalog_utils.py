import ijson
import re
from str_analysis.utils.file_utils import open_file
from str_analysis.utils.misc_utils import parse_interval

ACGT_REGEX = re.compile("^[ACGT]+$", re.IGNORECASE)

def parse_motifs_from_locus_structure(locus_structure):
    """Takes an ExpansionHunter LocusStructure like "(CAG)*AGAC(GCC)+" or a TRGT STRUC like "TCCAG(CAG)nAGAC(GCC)nGCACGG"
    and returns a list of motifs ["CAG", "GCC"].
    """
    return [
        locus_structure_part.split("(")[-1] for locus_structure_part in locus_structure.split(")")[:-1]
    ]


def convert_json_records_to_bed_format_tuples(json_records):
    """Takes an iterator over json records in ExpansionHunter format, and coverts them to BED format tuples

    Yield:
        5-tuples: (chrom, start_0based, end_1based, motif, repeat_count)
    """
    for record in json_records:
        locus_structure = record["LocusStructure"]
        reference_regions = record["ReferenceRegion"]
        if not isinstance(reference_regions, list):
            reference_regions = [reference_regions]

        motifs = parse_motifs_from_locus_structure(locus_structure)
        if len(motifs) != len(reference_regions):
            raise ValueError(f"Number of motifs ({len(motifs)}) != number of reference regions "
                             f"({len(reference_regions)}) in record: {record}")

        for reference_region, motif in zip(reference_regions, motifs):
            chrom, start_0based, end_1based = parse_interval(reference_region)
            yield chrom, start_0based, end_1based, motif, f"{(end_1based - start_0based) / len(motif):0.1f}"


def get_variant_catalog_iterator(variant_catalog_json_or_bed):
    """Takes the path of a JSON or BED file and returns an iterator over variant catalog records parsed from that file.

    Args:
        variant_catalog_json_or_bed (str): path to a JSON or BED file containing variant catalog records

    Yields:
        dict: a variant catalog record parsed from the file
    """

    if ".json" in variant_catalog_json_or_bed:
        with open_file(variant_catalog_json_or_bed, is_text_file=True) as f:
            for record in ijson.items(f, "item"):
                yield record
    else:
        with open_file(variant_catalog_json_or_bed, is_text_file=True) as input_variant_catalog:
            for line in input_variant_catalog:
                fields = line.strip().split("\t")
                unmodified_chrom = fields[0]
                chrom = unmodified_chrom.replace("chr", "")
                start_0based = int(fields[1])
                end_1based = int(fields[2])
                motif = fields[3].strip("()*+").upper()
                if not ACGT_REGEX.match(motif):
                    print(f"WARNING: skipping line with invalid motif: {line.strip()}")
                    continue

                record = {
                    "LocusId": f"{chrom}-{start_0based + 1}-{end_1based}-{motif}",
                    "ReferenceRegion": f"{unmodified_chrom}:{start_0based}-{end_1based}",
                    "LocusStructure": f"({motif})*",
                    "VariantType": "Repeat",
                }
                yield record
