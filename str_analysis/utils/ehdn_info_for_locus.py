
def parse_ehdn_info_for_locus(ehdn_profile, locus_chrom, locus_start, locus_end, margin=700, motifs_of_interest=None):
    """Extract info related to a specific locus from an ExpansionHunterDenovo profile.

    Args:
        ehdn_profile (dict): dictionary representing the data from an EHdn str_profile.json file.
        locus_chrom (str): locus chromosome
        locus_start (int): locus start coord
        locus_end (int): locus end coord
        margin (int): when looking for anchored-IRR regions, include regions this many base pairs away from the locus.
            This should be approximately the fragment-length or slightly larger.
        motifs_of_interest (set): optionally, a set of motifs to include in the results even if EHdn found only
            paired-IRRs, and no anchored IRRs near the given locus.
    Returns:
        list of dictionaries where each dictionary represents an :
               {
                "region": "chr18:52204909-52204910",   # EHdn region containing anchored IRRs
                "repeat_unit": "CAG",
                "n_anchored_regions_for_this_repeat_unit": 3, # number of anchored regions for this locus
                "anchored_irr_count_for_this_repeat_unit_and_region": 5,  # number of IRRs found
                "total_anchored_irr_count_for_this_repeat_unit": 10,  # number of IRRs found across all regions that have this same repeat unit
                "paired_irr_count_for_this_repeat_unit": 5,  # number of paired IRRs
                "total_irr_count_for_this_repeat_unit_and_region":  7.5   #  anchored_irr_count_for_this_repeat_unit_and_region plus the total_anchored_irr_count_for_this_repeat_unit weighted by the fraction of anchored IRRs at the locus of interest vs. all other loci
                "sample_read_depth": 30.6, # overall sample coverage computed by EHdn
            }
    """
    locus_chrom = locus_chrom.replace("chr", "")

    sample_read_depth = ehdn_profile.pop("Depth")
    sample_read_length = ehdn_profile.pop("ReadLength")

    records = []
    for repeat_unit, irr_counts in ehdn_profile.items():  # contains keys: IrrPairCounts, RegionsWithIrrAnchors
        total_anchored_irr_count = irr_counts.get("AnchoredIrrCount", 0)
        irr_pair_count = irr_counts.get("IrrPairCount", 0)
        anchored_irr_regions = irr_counts.get("RegionsWithIrrAnchors", {})
        for region, read_count in anchored_irr_regions.items():
            chrom, start_and_end = region.split(":")
            chrom = chrom.replace("chr", "")
            start, end = map(int, start_and_end.split("-"))

            overlaps_locus = ((chrom == locus_chrom) and
                              (end >= locus_start - margin) and
                              (start <= locus_end + margin))
            if not overlaps_locus:
                continue

            records.append({
                "region": region,
                "repeat_unit": repeat_unit,
                "n_anchored_regions_for_this_repeat_unit": len(anchored_irr_regions),
                "anchored_irr_count_for_this_repeat_unit_and_region": read_count,
                "total_anchored_irr_count_for_this_repeat_unit": total_anchored_irr_count,
                "paired_irr_count_for_this_repeat_unit": irr_pair_count,
                "total_irr_count_for_this_repeat_unit_and_region":  read_count + irr_pair_count * read_count / float(total_anchored_irr_count),
                "sample_read_depth": sample_read_depth,
            })
            break
        else:
            # if there are irr_pairs for a known repeat unit, record them even though no achored reads are present
            if motifs_of_interest and repeat_unit in motifs_of_interest:
                records.append({
                    "region": None,
                    "repeat_unit": repeat_unit,
                    "n_anchored_regions_for_this_repeat_unit": 0,
                    "anchored_irr_count_for_this_repeat_unit_and_region": 0,
                    "total_anchored_irr_count_for_this_repeat_unit": total_anchored_irr_count,
                    "paired_irr_count_for_this_repeat_unit": irr_pair_count,
                    "total_irr_count_for_this_repeat_unit_and_region":  irr_pair_count,
                    "sample_read_depth": sample_read_depth,
                })

    return records, sample_read_depth, sample_read_length
