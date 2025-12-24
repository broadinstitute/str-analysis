import unittest

from str_analysis.utils.eh_catalog_utils import convert_json_records_to_bed_format_tuples
from utils.eh_catalog_utils import group_overlapping_loci


class Tests(unittest.TestCase):

    def test_convert_json_records_to_bed_format_tuples(self):

        json_records = [
            {
                "LocusStructure": "(CAG)*AGAC(GCC)+",
                "ReferenceRegion": [
                    "chr1:100-110",
                    "chr1:115-121",
                ]
            },
            {
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr22:12104420-12104435"
            }
        ]

        expected = [
            ("chr1", 100, 110, "CAG", "3.3"),
            ("chr1", 115, 121, "GCC", "2.0"),
            ("chr22", 12104420, 12104435, "CAG", "5.0")
        ]

        actual = list(convert_json_records_to_bed_format_tuples(json_records))
        self.assertEqual(actual, expected)

        json_records = [
            {
                "LocusStructure": "(CAG)*AGAC(GCC)+",
                "ReferenceRegion": "chr1:100-110",
            }
        ]

        with self.assertRaises(ValueError):
            list(convert_json_records_to_bed_format_tuples(json_records))

    def test_group_overlapping_loci(self):
        json_records = [
            { "LocusStructure": "(CAG)*", "ReferenceRegion": "chrX:1000-2000"},
            { "LocusStructure": "(CCG)*", "ReferenceRegion": "chrX:1000-2000"},
            { "LocusStructure": "(CAG)*", "ReferenceRegion": "chrX:1000-2000"},
            { "LocusStructure": "(CCG)*", "ReferenceRegion": "chrX:1998-3000"},
            { "LocusStructure": "(CAG)*", "ReferenceRegion": "chrX:1999-3000"},
            { "LocusStructure": "(CAG)*", "ReferenceRegion": "chrX:2998-3000"},
        ]

        groups = group_overlapping_loci(
            json_records,
            only_group_loci_with_similar_motifs=True,
            min_overlap_size=2,
            verbose=True,
        )

        groups = list(groups)
        self.assertEqual(len(groups), 3)
        self.assertEqual(len(groups[0]), 2)
        self.assertEqual(len(groups[1]), 2)
        self.assertEqual(len(groups[2]), 2)

        self.assertEqual(groups[0][0]["LocusStructure"], "(CAG)*")
        self.assertEqual(groups[0][1]["LocusStructure"], "(CAG)*")

        self.assertEqual(groups[1][0]["LocusStructure"], "(CAG)*")
        self.assertEqual(groups[1][1]["LocusStructure"], "(CAG)*")

        self.assertEqual(groups[2][0]["LocusStructure"], "(CCG)*")
        self.assertEqual(groups[2][1]["LocusStructure"], "(CCG)*")

        from pprint import pprint
        pprint(list(groups))