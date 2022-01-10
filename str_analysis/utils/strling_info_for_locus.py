

def parse_strling_info_for_locus(strling_genotype_df, locus_chrom, locus_start, locus_end, margin=700):
    """Extract info related to a specific locus from the STRling genotype table.
    See https://strling.readthedocs.io/en/latest/outputs.html for more details.

    Args:
        strling_genotype_df (pd.DataFrame): parsed STRling genotype table
        locus_chrom (str): locus chromosome
        locus_start (int): locus start coord
        locus_end (int): locus end coord
        margin (int): when looking for matching STRling results, include calls this many base pairs away from the locus.
            This should be approximately the fragment-length or slightly larger (700 is a reasonable value for
            Illumina short read data).
    """

    locus_chrom = locus_chrom.replace("chr", "")
    strling_results_for_locus = strling_genotype_df[
        (strling_genotype_df.chrom.str.replace("chr", "") == locus_chrom)
        & (strling_genotype_df.start_1based < locus_end + margin)
        & (strling_genotype_df.end_1based > locus_start - margin)
    ]

    """
    Example row:
    
    {'allele1_est': 0.0,
     'allele2_est': 22.09,
     'anchored_reads': 8,
     'chrom': 'chr4',
     'depth': 32.0,
     'expected_spanning_pairs': 57.83,
     'left': 39348477,
     'left_clips': 5,
     'repeatunit': 'AAGGG',
     'right': 39348477,
     'right_clips': 0,
     'spanning_pairs': 0,
     'spanning_pairs_pctl': 0.0,
     'spanning_reads': 0,
     'sum_str_counts': 265,                
     'unplaced_pairs': 0}
    """

    return [row.to_dict() for _, row in strling_results_for_locus.iterrows()]
