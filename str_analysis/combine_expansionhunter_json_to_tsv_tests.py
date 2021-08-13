import json
import unittest

from combine_expansionhunter_json_to_tsv import parse_read_count_tuples, convert_expansionhunter_json_to_tsv_columns


VARIANT_CATALOG_CONTENTS = json.loads("""[{
  "LocusId": "X-149631736-149631780-TMEM185A",
  "LocusStructure": "(CGC)*",
  "ReferenceRegion": "chrX:149631735-149631780",
  "VariantType": "Repeat",
  "IsOfficial": false,
  "PMID": "7874164"
},
{
  "LocusId": "ATXN7",
  "LocusStructure": "(GCA)*(GCC)+",
  "ReferenceRegion": [
     "chr3:63912684-63912714",
     "chr3:63912714-63912726"
  ],
  "VariantId": [
     "ATXN7",
     "ATXN7_GCC_ADJACENT"
  ],
  "VariantType": [
     "Repeat",
     "Repeat"
  ],
  "IsOfficial": true,
  "PMID": "29398703"
},
{
  "LocusId": "X-25013650-25013697-ARX",
  "LocusStructure": "(CGC)*",
  "ReferenceRegion": "chrX:25013649-25013697",
  "VariantType": "Repeat",
  "IsOfficial": false,
  "PMID": "29946432"
}]
""")

EHv4_JSON = json.loads("""{
  "LocusResults": {
    "1-146228800-146228820-NOTCH2NLC": {
      "AlleleCount": 2,
      "Coverage": 50.469442942130875,
      "FragmentLength": 364,
      "LocusId": "1-146228800-146228820-NOTCH2NLC",
      "ReadLength": 151,
      "Variants": {
        "1-146228800-146228820-NOTCH2NLC": {
          "CountsOfSpanningReads": "(7, 60), (19, 4)",
          "CountsOfFlankingReads": "(1, 1), (2, 1), (4, 5), (6, 4), (7, 1)",
          "CountsOfInrepeatReads": "()",
          "Genotype": "7/7",
          "GenotypeConfidenceInterval": "7-7/7-7",
          "ReferenceRegion": "chr1:146228799-146228820",
          "RepeatUnit": "CGC",
          "VariantId": "1-146228800-146228820-NOTCH2NLC",
          "VariantSubtype": "Repeat",
          "VariantType": "Repeat"
        }
      }
    },
    "ATXN7": {
      "AlleleCount": 2,
      "Coverage": 57.574364521362895,
      "FragmentLength": 357,
      "LocusId": "ATXN7",
      "ReadLength": 151,
      "Variants": {
        "ATXN7": {
          "CountsOfSpanningReads": "(10, 29)",
          "CountsOfFlankingReads": "(1, 1), (2, 3), (3, 1), (5, 1), (6, 1), (7, 3), (8, 2), (9, 4), (10, 4)",
          "CountsOfInrepeatReads": "()",
          "Genotype": "10/10",
          "GenotypeConfidenceInterval": "10-10/10-10",
          "ReferenceRegion": "chr3:63912684-63912714",
          "RepeatUnit": "GCA",
          "VariantId": "ATXN7",
          "VariantSubtype": "Repeat",
          "VariantType": "Repeat"
        },
        "ATXN7_GCC_ADJACENT": {
          "CountsOfFlankingReads": "(1, 2), (3, 4), (4, 1)",
          "CountsOfInrepeatReads": "()",
          "CountsOfSpanningReads": "(4, 37)",
          "Genotype": "4/4",
          "GenotypeConfidenceInterval": "4-4/4-4",
          "ReferenceRegion": "chr3:63912714-63912726",
          "RepeatUnit": "GCC",
          "VariantId": "ATXN7_GCC_ADJACENT",
          "VariantSubtype": "Repeat",
          "VariantType": "Repeat"
        }
      }
    },
    "X-25013650-25013697-ARX": {
      "AlleleCount": 1,
      "Coverage": 25.31638723634397,
      "FragmentLength": 343,
      "LocusId": "X-25013650-25013697-ARX",
      "ReadLength": 151,
      "Variants": {
        "X-25013650-25013697-ARX": {
          "CountsOfSpanningReads": "(14, 1), (16, 11)",
          "CountsOfFlankingReads": "(1, 2), (6, 2), (9, 1), (10, 2), (11, 2), (13, 1), (14, 2), (15, 2), (16, 5)",
          "CountsOfInrepeatReads": "()",
          "Genotype": "16",
          "GenotypeConfidenceInterval": "16-16",
          "ReferenceRegion": "chrX:25013649-25013697",
          "RepeatUnit": "CGC",
          "VariantId": "X-25013650-25013697-ARX",
          "VariantSubtype": "Repeat",
          "VariantType": "Repeat"
        }
      }
    }
  },
  "SampleParameters": {
    "SampleId": "TestSample",
    "Sex": "Male"
  }
}""")


def get_expected_columns(row, with_allele_records=False, with_json_file_path=False):
    return (
    [] if not with_json_file_path else [
        'Dirname',
        'Filename',
    ]) + [
        'Sex',
        'SampleId',
        'LocusId',
        'ReadLength',
        'Coverage',
        'AlleleCount',
        'InVariantCatalog',
    ] + ([
        'VariantCatalog_LocusStructure',
        'VariantCatalog_IsOfficial',
        'VariantCatalog_PMID',
        'VariantCatalog_NumOfftargetRegions',
    ] if row['InVariantCatalog'] else []) + [
        'VariantId',
        'VariantSubtype',
        'RepeatUnit',
        'RepeatUnitLength',
        'ReferenceRegion',
        'CountsOfSpanningReads',
        'CountsOfFlankingReads',
        'CountsOfInrepeatReads',
        'NumSpanningReads',
        'NumFlankingReads',
        'NumInrepeatReads',
        'NumReadsTotal',
        'NumAllelesSupportedBySpanningReads',
        'NumAllelesSupportedByFlankingReads',
        'NumAllelesSupportedByInrepeatReads',
        'NumAllelesSupportedTotal',
        'Genotype',
        'GenotypeConfidenceInterval',
    ] + ([
        'Allele Number: Allele 1',
        'Num Repeats: Allele 1',
        'Repeat Size (bp): Allele 1',
        'CI start: Allele 1',
        'CI end: Allele 1',
        'CI size: Allele 1',
        'NumSpanningReadsThatSupportGenotype: Allele 1',
        'NumFlankingReadsThatSupportGenotype: Allele 1',
        'NumInrepeatReadsThatSupportGenotype: Allele 1',
        'NumReadsTotalThatSupportGenotype: Allele 1',
        'FractionOfReadsThatSupportsGenotype: Allele 1',
    ] + ([] if row['AlleleCount'] == 1 else [
        'Allele Number: Allele 2',
        'Num Repeats: Allele 2',
        'Repeat Size (bp): Allele 2',
        'CI start: Allele 2',
        'CI end: Allele 2',
        'CI size: Allele 2',
        'NumSpanningReadsThatSupportGenotype: Allele 2',
        'NumFlankingReadsThatSupportGenotype: Allele 2',
        'NumInrepeatReadsThatSupportGenotype: Allele 2',
        'NumReadsTotalThatSupportGenotype: Allele 2',
        'FractionOfReadsThatSupportsGenotype: Allele 2',
    ]) if not with_allele_records else [
        'Allele Number',
        'Num Repeats',
        'Repeat Size (bp)',
        'CI start',
        'CI end',
        'CI size',
        'NumSpanningReadsThatSupportGenotype',
        'NumFlankingReadsThatSupportGenotype',
        'NumInrepeatReadsThatSupportGenotype',
        'NumReadsTotalThatSupportGenotype',
        'FractionOfReadsThatSupportsGenotype',
    ])


class Tests(unittest.TestCase):

    def test_parse_read_count_tuples(self):
        result = parse_read_count_tuples("(11, 4), (16, 2), (17, 13)")
        self.assertListEqual(result, [(11, 4), (16, 2), (17, 13)])

        result = parse_read_count_tuples("()")
        self.assertListEqual(result, [])

    def test_convert_expansionhunter_json_to_tsv_columns_for_variants(self):
        variant_rows = convert_expansionhunter_json_to_tsv_columns(
            EHv4_JSON,
            variant_catalog_contents=VARIANT_CATALOG_CONTENTS,
            json_file_path="",
            return_allele_records=False,
        )
        self.assertEqual(len(variant_rows), 4)
        for row in variant_rows:
            self.assertIn('InVariantCatalog', row)
            self.assertIn('AlleleCount', row)
            self.assertListEqual(list(row.keys()), get_expected_columns(row, with_allele_records=False, with_json_file_path=False))

            self.assertIn(row["LocusId"], {"1-146228800-146228820-NOTCH2NLC", "ATXN7", "X-25013650-25013697-ARX"})
            self.assertIn(row["VariantId"], {"1-146228800-146228820-NOTCH2NLC", "ATXN7", "ATXN7_GCC_ADJACENT", "X-25013650-25013697-ARX"})

            self.assertEqual(row["Sex"], "Male")
            self.assertEqual(row["SampleId"], "TestSample")

        self.assertEqual(variant_rows[0]["Genotype"], "7/7")
        self.assertEqual(variant_rows[0]["GenotypeConfidenceInterval"], "7-7/7-7")
        self.assertEqual(variant_rows[0]["CountsOfSpanningReads"], "(7, 60), (19, 4)")
        self.assertEqual(variant_rows[0]["CountsOfFlankingReads"], "(1, 1), (2, 1), (4, 5), (6, 4), (7, 1)")
        self.assertEqual(variant_rows[0]["CountsOfInrepeatReads"], "()")
        self.assertEqual(variant_rows[0]["NumSpanningReads"], 64)
        self.assertEqual(variant_rows[0]["NumFlankingReads"], 12)
        self.assertEqual(variant_rows[0]["NumInrepeatReads"], 0)
        self.assertEqual(variant_rows[0]["NumReadsTotal"], 76)
        self.assertEqual(variant_rows[0]["NumAllelesSupportedBySpanningReads"], 2)
        self.assertEqual(variant_rows[0]["NumAllelesSupportedByFlankingReads"], 5)
        self.assertEqual(variant_rows[0]["NumAllelesSupportedByInrepeatReads"], 0)

        for allele_num in 1, 2:
            self.assertAlmostEqual(variant_rows[0][f"NumSpanningReadsThatSupportGenotype: Allele {allele_num}"], 30.0)
            self.assertAlmostEqual(variant_rows[0][f"NumFlankingReadsThatSupportGenotype: Allele {allele_num}"], 0.5)
            self.assertEqual(variant_rows[0][f"NumInrepeatReadsThatSupportGenotype: Allele {allele_num}"], 0)
            self.assertAlmostEqual(variant_rows[0][f"NumReadsTotalThatSupportGenotype: Allele {allele_num}"], 30.5)
            self.assertAlmostEqual(variant_rows[0][f"FractionOfReadsThatSupportsGenotype: Allele {allele_num}"], 30.5/76)

    def test_convert_expansionhunter_json_to_tsv_columns_for_alleles(self):
        allele_rows = convert_expansionhunter_json_to_tsv_columns(
            EHv4_JSON,
            variant_catalog_contents=VARIANT_CATALOG_CONTENTS,
            json_file_path="/temp/file.json",
            return_allele_records=True,
        )
        self.assertEqual(len(allele_rows), 7)
        for row in allele_rows:
            self.assertIn('InVariantCatalog', row)
            self.assertIn('AlleleCount', row)
            self.assertListEqual(list(row.keys()), get_expected_columns(row, with_allele_records=True, with_json_file_path=True))

            self.assertIn(row["LocusId"], {"1-146228800-146228820-NOTCH2NLC", "ATXN7", "X-25013650-25013697-ARX"})
            self.assertIn(row["VariantId"], {"1-146228800-146228820-NOTCH2NLC", "ATXN7", "ATXN7_GCC_ADJACENT", "X-25013650-25013697-ARX"})

            self.assertEqual(row["Sex"], "Male")
            self.assertEqual(row["SampleId"], "TestSample")

        self.assertEqual(allele_rows[0]["Genotype"], "7/7")
        self.assertEqual(allele_rows[0]["GenotypeConfidenceInterval"], "7-7/7-7")
        self.assertEqual(allele_rows[0]["CountsOfSpanningReads"], "(7, 60), (19, 4)")
        self.assertEqual(allele_rows[0]["CountsOfFlankingReads"], "(1, 1), (2, 1), (4, 5), (6, 4), (7, 1)")
        self.assertEqual(allele_rows[0]["CountsOfInrepeatReads"], "()")
        self.assertEqual(allele_rows[0]["NumSpanningReads"], 64)
        self.assertEqual(allele_rows[0]["NumInrepeatReads"], 0)
        self.assertEqual(allele_rows[0]["NumAllelesSupportedBySpanningReads"], 2)
        self.assertEqual(allele_rows[0]["NumAllelesSupportedByInrepeatReads"], 0)

        self.assertAlmostEqual(allele_rows[0]["NumSpanningReadsThatSupportGenotype"], 30.0)
        self.assertAlmostEqual(allele_rows[0]["NumFlankingReadsThatSupportGenotype"], 0.5)
        self.assertEqual(allele_rows[0]["NumInrepeatReadsThatSupportGenotype"], 0)
        self.assertAlmostEqual(allele_rows[0]["NumReadsTotalThatSupportGenotype"], 30.5)
        self.assertAlmostEqual(allele_rows[0]["FractionOfReadsThatSupportsGenotype"], 30.5/76)

        self.assertEqual(allele_rows[-1]["Genotype"], "16")
        self.assertEqual(allele_rows[-1]["CountsOfSpanningReads"], "(14, 1), (16, 11)")
        self.assertEqual(allele_rows[-1]["CountsOfFlankingReads"], "(1, 2), (6, 2), (9, 1), (10, 2), (11, 2), (13, 1), (14, 2), (15, 2), (16, 5)")
        self.assertEqual(allele_rows[-1]["CountsOfInrepeatReads"], "()")
        self.assertEqual(allele_rows[-1]["NumSpanningReadsThatSupportGenotype"], 11)
        self.assertEqual(allele_rows[-1]["NumFlankingReadsThatSupportGenotype"], 5)
        self.assertEqual(allele_rows[-1]["NumInrepeatReadsThatSupportGenotype"], 0)
        self.assertEqual(allele_rows[-1]["NumReadsTotalThatSupportGenotype"], 16)
        self.assertAlmostEqual(allele_rows[-1]["FractionOfReadsThatSupportsGenotype"], 16/31.)
