import ijson
import re
from str_analysis.utils.file_utils import open_file
from str_analysis.utils.misc_utils import parse_interval
from tqdm import tqdm

ACGTN_REGEX = re.compile("^[ACGTN]+$", re.IGNORECASE)

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


def get_variant_catalog_iterator(variant_catalog_json_or_bed, show_progress_bar=False):
    """Takes the path of a JSON or BED file and returns an iterator over variant catalog records parsed from that file.

    Args:
        variant_catalog_json_or_bed (str): path to a JSON or BED file containing variant catalog records
        show_progress_bar (bool): whether to show a tqdm progress bar
    Yields:
        dict: a variant catalog record parsed from the file
    """

    if ".json" in variant_catalog_json_or_bed:
        with open_file(variant_catalog_json_or_bed, is_text_file=True) as f:
            if show_progress_bar:
                f = tqdm(f, unit=" records", unit_scale=True)

            for record in ijson.items(f, "item", use_float=True):
                yield record
    else:
        with open_file(variant_catalog_json_or_bed, is_text_file=True) as input_variant_catalog:
            if show_progress_bar:
                input_variant_catalog = tqdm(input_variant_catalog, unit=" records", unit_scale=True)

            for line in input_variant_catalog:
                fields = line.strip().split("\t")
                try:
                    unmodified_chrom = fields[0]
                    chrom = unmodified_chrom.replace("chr", "")
                    start_0based = int(fields[1])
                    end_1based = int(fields[2])

                    if "ID=" in fields[3] and "STRUC=" in fields[3]:
                        # This is a TRGT format line, skip it
                        info = dict(key_value.split("=") for key_value in fields[3].split(";"))
                        locus_id = info["ID"]
                        locus_structure = info["STRUC"].replace("n", "*")
                    else:
                        motif = fields[3].strip("()*+").upper()
                        locus_id = f"{chrom}-{start_0based}-{end_1based}-{motif}"
                        locus_structure = f"({motif})*"
                        if not ACGTN_REGEX.match(motif):
                            print(f"WARNING: skipping line with invalid motif: {line.strip()}")
                            continue

                    record = {
                        "LocusId": locus_id,
                        "ReferenceRegion": f"{unmodified_chrom}:{start_0based}-{end_1based}",
                        "LocusStructure": locus_structure,
                        "VariantType": "Repeat",
                    }
                    yield record
                except Exception as e:
                    raise ValueError(f"Unable to parse line {fields}. Error: {e}")