"""Utility functions for TRGT analysis"""


def convert_trgt_locus_to_expansion_hunter_format(
    fasta_obj, chrom, start_0based, end_1based, locus_structure, locus_id_prefix=None, verbose=False):
    """Convert a TRGT locus definition to one or more ExpansionHunter catalog records"""

    unprocessed_reference_sequence = fasta_obj[chrom][start_0based:end_1based]
    unprocessed_locus_structure = locus_structure
    current_start_coord_0based = start_0based
    while "(" in unprocessed_locus_structure:  # continue while there are more motifs
        current_end_coord = current_start_coord_0based
        left_flank_end_i = unprocessed_locus_structure.find("(")
        if left_flank_end_i != -1:
            left_flank = unprocessed_locus_structure[:left_flank_end_i]
            if not unprocessed_reference_sequence.startswith(left_flank):
                print(f"WARNING: left flank {left_flank} doesn't match reference sequence "
                      f"{unprocessed_reference_sequence[:left_flank_end_i]} @ locus {locus_id}. Skipping...")
                break

            unprocessed_reference_sequence = unprocessed_reference_sequence[len(left_flank):]
            unprocessed_locus_structure = unprocessed_locus_structure[left_flank_end_i:]
            current_start_coord_0based += len(left_flank)
            current_end_coord += len(left_flank)
        motif_end_i = unprocessed_locus_structure.find(")n")
        if motif_end_i == -1:
            break

        motif = unprocessed_locus_structure[1:motif_end_i]
        unprocessed_locus_structure = unprocessed_locus_structure[motif_end_i + 2:]

        # count number of repeats in reference that exactly match the motif
        for i in range(0, len(unprocessed_reference_sequence), len(motif)):
            if not unprocessed_reference_sequence.startswith(motif):
                break
            unprocessed_reference_sequence = unprocessed_reference_sequence[len(motif):]
            current_end_coord += len(motif)

        if current_start_coord_0based != current_end_coord:
            current_locus_id = locus_id_prefix if locus_id_prefix else f"{chrom}-{current_start_coord_0based}-{current_end_coord}"
            current_locus_id = f"{current_locus_id}-{motif}"

            if verbose:
                print(f"Adding locus {current_locus_id} with motif {motif} and reference region "
                      f"{chrom}:{current_start_coord_0based}-{current_end_coord}")

            yield {
                "LocusId": current_locus_id,
                "ReferenceRegion": f"{chrom}:{current_start_coord_0based}-{current_end_coord}",
                "LocusStructure": f"({motif})*",
                "VariantType": "Repeat",
            }

            current_start_coord_0based = current_end_coord
