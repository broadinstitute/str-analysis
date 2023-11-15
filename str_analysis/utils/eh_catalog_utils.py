
def parse_motifs_from_locus_structure(locus_structure):
    """Takes an ExpansionHunter LocusStructure like "(CAG)*AGAC(GCC)+" and returns a list of motifs ["CAG", "GCC"]"""
    return [
        locus_structure_part.split("(")[-1] for locus_structure_part in locus_structure.split(")")[:-1]
    ]


def convert_json_records_to_bed_format_tuples(json_records):
    """Takes an iterator over json records in ExpansionHunter format, and coverts them to BED format tuples

    Yield:
        5-tuples: (chrom, start_0based, end_1based, motif, repeat_count)
    """
    for record in output_catalog_record_iter:
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
