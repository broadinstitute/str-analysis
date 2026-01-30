import collections
import ijson
import itertools
import re
from tqdm import tqdm

from str_analysis.utils.file_utils import open_file
from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif

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
            iterator = ijson.items(f, "item", use_float=True)
            if show_progress_bar:
                iterator = tqdm(iterator, unit=" records", unit_scale=True)

            for record in iterator:
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
                        motif = fields[3].split(":")[0].strip("()*+").upper()
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


def compute_repeat_unit_id(canonical_repeat_unit):
    """Compute a unique identifier for a canonical repeat unit. Repeat units with the same id are treated as the same repeat unit.
    
    Args:
        canonical_repeat_unit (str): canonical repeat unit
    """
    
    if len(canonical_repeat_unit) <= 6:
        return canonical_repeat_unit
    else:
        return len(canonical_repeat_unit)


def group_overlapping_loci(
        catalog_record_iterator,
        only_group_loci_with_similar_motifs=False,
        min_overlap_size=-1,
        verbose=False,
):
    """Convert catalog to groups of overlapping tandem repeat loci.

    Args:
        catalog_record_iterator (iter): Iterator over variant catalog records represented as dictionaries. The records
            must already be sorted by genomic coordinates.
        only_group_loci_with_similar_motifs (bool): only group loci that have similar motifs - defined as same
            canonical motif for STRs and same motif length for VNTRs
        min_overlap_size (int): The min overlap size in base pairs. If it's negative, loci that are that many base pairs
            apart will be grouped.
        verbose (bool): if True, print verbose output

    Yield:
        list: a list of records that overlap each other
    """

    def get_chrom_from_record(record):
        if isinstance(record["ReferenceRegion"], list):
            raise ValueError("This method does not support locus definitions with a list of ReferenceRegions like " + str(record["ReferenceRegion"]))

        chrom = record["ReferenceRegion"].split(":")[0]
        chrom = chrom.replace("chr", "")
        return chrom

    # process records one chromosome at a time
    seen_chroms = set()
    for chrom, records_in_chrom in itertools.groupby(catalog_record_iterator, key=get_chrom_from_record):
        if chrom in seen_chroms:
            raise ValueError(f"Input catalog is not sorted: chromosome {chrom} appears non-consecutively. "
                             f"All records for a chromosome must be grouped together.")
        motif_id_to_record_group = collections.defaultdict(list)
        motif_id_to_record_group_end_1based = collections.defaultdict(int)
        previous_record_coordinates = None
        for record in records_in_chrom:
            if only_group_loci_with_similar_motifs:
                motifs = parse_motifs_from_locus_structure(record["LocusStructure"])
                if len(motifs) != 1:
                    raise ValueError(f"This method requires exactly one motif in the "
                                     f"LocusStructure, but got {len(motifs)} motifs in {record['LocusStructure']}")

                record_motif = motifs[0]
                current_repeat_unit_id = compute_repeat_unit_id(
                    compute_canonical_motif(record_motif, include_reverse_complement=True))
            else:
                current_repeat_unit_id = 0

            _, record_start_0based, record_end_1based = parse_interval(record["ReferenceRegion"])

            current_record_coordinates = (record_start_0based, record_end_1based)
            if previous_record_coordinates is not None and previous_record_coordinates > current_record_coordinates:
                raise ValueError(f"Input catalog is not sorted: "
                                 f"{chrom}:{previous_record_coordinates[0]}-{previous_record_coordinates[1]} "
                                 f"should be after "
                                 f"{chrom}:{current_record_coordinates[0]}-{current_record_coordinates[1]}")

            previous_record_coordinates = current_record_coordinates

            if len(motif_id_to_record_group[current_repeat_unit_id]) == 0:
                motif_id_to_record_group_end_1based[current_repeat_unit_id] = record_end_1based
                motif_id_to_record_group[current_repeat_unit_id].append(record)
                continue

            # check if the current allele overlaps with the current group
            if record_start_0based <= motif_id_to_record_group_end_1based[current_repeat_unit_id] - min_overlap_size:
                motif_id_to_record_group[current_repeat_unit_id].append(record)
                motif_id_to_record_group_end_1based[current_repeat_unit_id] = max(motif_id_to_record_group_end_1based[current_repeat_unit_id], record_end_1based)
                continue

            # add the current group
            record_group = motif_id_to_record_group[current_repeat_unit_id]
            yield record_group

            # start a new group
            motif_id_to_record_group[current_repeat_unit_id] = [record]
            motif_id_to_record_group_end_1based[current_repeat_unit_id] = record_end_1based

        if len(motif_id_to_record_group) > 0:
            for current_repeat_unit_id, record_group in motif_id_to_record_group.items():
                yield record_group

        seen_chroms.add(chrom)
