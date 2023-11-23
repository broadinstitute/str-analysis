"""NOTE: These unittests should be run manually since they take a relatively long time and require large data files"""

import simplejson as json
import os
from pprint import pformat
import pysam
import unittest

from add_adjacent_loci_to_expansion_hunter_catalog import process_input_record, get_interval_tree_for_chrom, \
    get_min_and_max_coords_for_chrom
from str_analysis.utils.get_adjacent_repeats import MAX_DISTANCE_BETWEEN_REPEATS, MAX_TOTAL_ADJACENT_REGION_SIZE, \
    MAX_OVERLAP_BETWEEN_ADJACENT_REPEATS


class GetAdjacentRepeatsTests(unittest.TestCase):

    def setUp(self):
        self.pysam_fasta_file = pysam.FastaFile(os.path.expanduser(f"~/hg38.fa"))
        self.source_of_adjacent_loci = os.path.expanduser(
            "~/code/str-truth-set/ref/other/repeat_specs_GRCh38_without_mismatches.including_homopolymers.sorted.at_least_9bp.bed.gz")
        self.variant_catalog_path = os.path.expanduser(
            "~/code/str-analysis/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json")

        self.variant_catalog = [
            {
                "LocusId": "ATXN7",
                "LocusStructure": "(GCA)*",
                "ReferenceRegion": "chr3:63912684-63912714",
                "ExpectedLocusStructure": "(GCA)*(GCC)*",
                "ExpectedReferenceRegion": [
                    "chr3:63912684-63912714",
                    "chr3:63912714-63912723"
                ],
            },
            {
                "LocusId": "ATXN8OS",
                "LocusStructure": "(CTG)*",
                "ReferenceRegion":  "chr13:70139383-70139428",
                "ExpectedLocusStructure": "(CTA)*(CTG)*",
                "ExpectedReferenceRegion": [
                    "chr13:70139353-70139383",
                    "chr13:70139383-70139428"
                ],
            },
            {
                "LocusId": "CNBP",
                "LocusStructure": "(CAGG)*",
                "ReferenceRegion": "chr3:129172576-129172656",
                "ExpectedLocusStructure": "(CAGG)*(CAGA)*(CA)*",
                "ExpectedReferenceRegion": [
                    "chr3:129172576-129172656",
                    "chr3:129172656-129172696",
                    "chr3:129172696-129172732"
                ],
            },
            {
                "LocusId": "FXN",
                "LocusStructure": "(GAA)*",
                "ReferenceRegion": "chr9:69037286-69037304",
                "ExpectedLocusStructure": "(A)*(GAA)*",
                "ExpectedReferenceRegion": [
                    "chr9:69037270-69037286",
                    "chr9:69037286-69037304"
                ],
            },
            {
                "LocusId": "NOP56",
                "LocusStructure": "(GGCCTG)*",
                "ReferenceRegion": "chr20:2652733-2652757",
                "ExpectedLocusStructure": "(GGCCTG)*(CGCCTG)*",
                "ExpectedReferenceRegion": [
                    "chr20:2652733-2652757",
                    "chr20:2652757-2652775"
                ],
            },
            {
                "LocusId": "HTT",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr4:3074876-3074933",
                "ExpectedLocusStructure": "(CAG)*CAACAGCCGCCA(CCG)*",
                "ExpectedReferenceRegion": [
                    "chr4:3074876-3074933",
                    "chr4:3074945-3074966"
                ],
            },
            {
                "LocusId": "PABPN1",
                "LocusStructure": "(GCG)*",
                "ReferenceRegion": "chr14:23321472-23321490",
                "ExpectedLocusStructure": "(GCG)*(GCA)*",
                "ExpectedReferenceRegion": [
                    "chr14:23321472-23321490",
                    "chr14:23321490-23321499",
                ],
            },
            {
                "LocusId": "FGF14",
                "LocusStructure": "(AAG)*",
                "ReferenceRegion": "chr13:102161574-102161724",
                "ExpectedLocusStructure": "(AAGA)*(AAG)*",
                "ExpectedReferenceRegion": [
                    "chr13:102161566-102161574",
                    "chr13:102161574-102161724",
                ],
            },
            {
                "LocusId": "NUTM2B-AS1",
                "LocusStructure": "(GGC)*",
                "ReferenceRegion": "chr10:79826383-79826404",
                "ExpectedLocusStructure": "(GGA)*AGCGGCGG(GGC)*",
                "ExpectedReferenceRegion": [
                    "chr10:79826366-79826375",
                    "chr10:79826383-79826404",
                ],
            },
            {
                "LocusId": "ATXN2",
                "LocusStructure": "(GCT)*",
                "ReferenceRegion": "chr12:111598949-111599018",
                "ExpectedLocusStructure": "(GCG)*(GCT)*",
                "ExpectedReferenceRegion": [
                    "chr12:111598943-111598949",
                    "chr12:111598949-111599018",
                ],
            },
            {
                "LocusId": "RILPL1",
                "LocusStructure": "(GGC)*",
                "ReferenceRegion": "chr12:123533720-123533750",
                "ExpectedLocusStructure": "(GGC)*AGCGGGGAGG(GC)*",
                "ExpectedReferenceRegion": [
                    "chr12:123533720-123533750",
                    "chr12:123533760-123533768",
                ],
            },
            {
                "LocusId": "TCF4",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr18:55586155-55586227",
                "ExpectedLocusStructure": "(GAG)*(CAG)*",
                "ExpectedReferenceRegion": [
                    "chr18:55586137-55586155",
                    "chr18:55586155-55586227",
                ],
            },
            {
                "LocusId": "TBX1",
                "LocusStructure": "(GCN)*",
                "ReferenceRegion": "chr22:19766762-19766807",
                "ExpectedLocusStructure": "(ACC)*CCGTGAGTCCA(GCC)*",
                "ExpectedReferenceRegion": [
                    "chr22:19766736-19766751",
                    "chr22:19766762-19766807",
                ],
            },
            {
                "LocusId": "SAMD12",
                "LocusStructure": "(TAAAA)*",
                "ReferenceRegion": "chr8:118366815-118366880",
                "ExpectedLocusStructure": "(TAAAA)*T(AATAA)*",
                "ExpectedReferenceRegion": [
                    "chr8:118366815-118366880",
                    "chr8:118366881-118366916",
                ],
            },
            {
                "LocusId": "AFF2",
                "LocusStructure": "(GCC)*",
                "ReferenceRegion": "chrX:148500631-148500691",
                "ExpectedLocusStructure": "(GCC)*CCGGCT(GCCGC)*",
                "ExpectedReferenceRegion": [
                    "chrX:148500631-148500691",
                    "chrX:148500697-148500712",
                ],
            },
            {
                "LocusId": "AR",
                "LocusStructure": "(GCA)*",
                "ReferenceRegion": "chrX:67545316-67545385",
                "ExpectedLocusStructure": "(GCT)*(GCA)*",
                "ExpectedReferenceRegion": [
                    "chrX:67545307-67545316",
                    "chrX:67545316-67545385",
                ],
            },
            {
                "LocusId": "DAB1",
                "LocusStructure": "(AAAAT)*",
                "ReferenceRegion": "chr1:57367043-57367118",
                "ExpectedLocusStructure": "(AAAAT)*(AAAT)*",
                "ExpectedReferenceRegion": [
                    "chr1:57367043-57367118",
                    "chr1:57367118-57367126",
                ],
            },
            {
                "LocusId": "FXN",
                "LocusStructure": "(GAA)*",
                "ReferenceRegion": "chr9:69037286-69037304",
                "ExpectedLocusStructure": "(A)*(GAA)*",
                "ExpectedReferenceRegion": [
                    "chr9:69037270-69037286",
                    "chr9:69037286-69037304",
                ],
            },
            {
                "LocusId": "HTT",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr4:3074876-3074933",
                "ExpectedLocusStructure": "(CAG)*CAACAGCCGCCA(CCG)*",
                "ExpectedReferenceRegion": [
                    "chr4:3074876-3074933",
                    "chr4:3074945-3074966",
                ],
            },
            {
                "LocusId": "NOP56",
                "LocusStructure": "(GGCCTG)*",
                "ReferenceRegion": "chr20:2652733-2652757",
                "ExpectedLocusStructure": "(GGCCTG)*(CGCCTG)*",
                "ExpectedReferenceRegion": [
                    "chr20:2652733-2652757",
                    "chr20:2652757-2652775",
                ],
            },
            {
                "LocusId": "FRA10AC1",
                "LocusStructure": "(CCG)*",
                "ReferenceRegion": "chr10:93702522-93702546",
                "ExpectedLocusStructure": "(CCA)*(CCG)*",
                "ExpectedReferenceRegion": [
                    "chr10:93702516-93702522",
                    "chr10:93702522-93702546",
                ],
            },
            {
                "LocusId": "TMEM185A",
                "LocusStructure": "(CGC)*",
                "ReferenceRegion": "chrX:149631735-149631780",
                "ExpectedLocusStructure": "(CGCCGT)*(CGC)*",
                "ExpectedReferenceRegion": [
                    "chrX:149631723-149631735",
                    "chrX:149631735-149631780",
                ],
            },
            {
                "LocusId": "AFF3",
                "LocusStructure": "(GCC)*",
                "ReferenceRegion": "chr2:100104798-100104822",
                "ExpectedLocusStructure": "(GCC)*GCGGTGCTCTG(CGCC)*",
                "ExpectedReferenceRegion": [
                    "chr2:100104798-100104822",
                    "chr2:100104833-100104849",
                ],
            }
        ]

        locus_ids_already_in_catalog = {record["LocusId"] for record in self.variant_catalog}
        with open(self.variant_catalog_path) as f:
            self.full_variant_catalog = json.load(f)
        for record in self.full_variant_catalog:
            if record["LocusId"] in locus_ids_already_in_catalog or isinstance(record["ReferenceRegion"], list):
                continue
            record["ExpectedLocusStructure"] = record["LocusStructure"]
            record["ExpectedReferenceRegion"] = record["ReferenceRegion"]
            self.variant_catalog.append(record)

        self.variant_catalog.sort(key=lambda record: (
            record["ReferenceRegion"][0].split(":")[0]
            if isinstance(record["ReferenceRegion"], list) else
            record["ReferenceRegion"].split(":")[0]))

    def test_known_disease_associated_loci(self):
        current_chrom = None
        current_interval_tree = None
        for record in self.variant_catalog:
            print("="*100)
            print(f"Testing {record['LocusId']}")
            chrom = record["ReferenceRegion"].split(":")[0]
            max_total_adjacent_region_size = 1000
            if current_chrom != chrom:
                current_chrom = chrom
                min_coord, max_coord = get_min_and_max_coords_for_chrom(
                    self.variant_catalog, chrom)
                min_coord -= max_total_adjacent_region_size
                max_coord += max_total_adjacent_region_size

                current_interval_tree = get_interval_tree_for_chrom(
                    self.source_of_adjacent_loci, chrom, min_coord, max_coord)

            output_record = process_input_record(
                record,
                self.pysam_fasta_file,
                current_interval_tree,
                max_distance_between_adjacent_repeats=15,
                max_total_adjacent_region_size=MAX_TOTAL_ADJACENT_REGION_SIZE,
                max_overlap_between_adjacent_repeats=MAX_OVERLAP_BETWEEN_ADJACENT_REPEATS,
            )
            print("Output record: ", pformat(output_record))

            self.assertEqual(output_record["LocusStructure"], record["ExpectedLocusStructure"])
            self.assertEqual(output_record["ReferenceRegion"], record["ExpectedReferenceRegion"])
            print("Confirmed that the output record has the expected ReferenceRegion and LocusStructure")
