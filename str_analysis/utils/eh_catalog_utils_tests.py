import gzip
import os
import tempfile
import unittest

import simplejson

from str_analysis.utils.eh_catalog_utils import convert_json_records_to_bed_format_tuples
from str_analysis.utils.eh_catalog_utils import group_overlapping_loci
from str_analysis.utils.eh_catalog_utils import JsonListFileWriter


class JsonListFileWriterTests(unittest.TestCase):

    def test_output_is_the_same_as_dumping_the_whole_list(self):
        records = [
            {"LocusId": "1-100-110-A", "ReferenceRegion": "1:100-110", "NumRepeatsInReference": 10,
             "ReferenceRepeatPurity": 0.95, "Nested": {"a": [1, 2.5, None]}, "Empty": [], "NaN": float("nan")},
            {"LocusId": "2-5-11-CAG", "ReferenceRegion": ["2:5-8", "2:8-11"], "Text": 'with "quotes" and\nnewline'},
            {"LocusId": "X-1-2-T"},
        ]
        with tempfile.TemporaryDirectory() as temp_dir:
            for records_subset in records, records[:1], []:
                for use_gzip in False, True:
                    output_path = os.path.join(
                        temp_dir, f"output_with_{len(records_subset)}_records.json" + (".gz" if use_gzip else ""))
                    writer = JsonListFileWriter(output_path, use_gzip=use_gzip)
                    self.assertFalse(os.path.exists(output_path))
                    for record in records_subset:
                        writer.write(record)
                    writer.close()

                    with (gzip.open if use_gzip else open)(output_path, "rt") as f:
                        self.assertEqual(f.read(), simplejson.dumps(records_subset, indent=4, ignore_nan=True))
                    self.assertFalse(os.path.exists(f"{output_path}.incomplete"))

                    if use_gzip:
                        # the gzip header stores the same file name as gzip.open(output_path) would: the output
                        # file's name without .gz, after the 10-byte fixed header
                        with open(output_path, "rb") as f:
                            header = f.read(200)
                        self.assertTrue(header[3] & 0x08)
                        stored_file_name = header[10:header.index(b"\x00", 10)]
                        self.assertEqual(stored_file_name.decode(), os.path.basename(output_path)[:-len(".gz")])


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