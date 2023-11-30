import unittest

from str_analysis.utils.eh_catalog_utils import convert_json_records_to_bed_format_tuples


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
