
def parse_motifs_from_locus_structure(locus_structure):
    """Takes an ExpansionHunter LocusStructure like "(CAG)*AGAC(GCC)+" and returns a list of motifs ["CAG", "GCC"]"""
    return [
        locus_structure_part.split("(")[-1] for locus_structure_part in locus_structure.split(")")[:-1]
    ]
