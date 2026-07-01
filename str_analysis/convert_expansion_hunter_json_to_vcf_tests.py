import os
import simplejson as json
import tempfile
import unittest

from str_analysis.convert_expansion_hunter_json_to_vcf import parse_expansion_hunter_json

EH_JSON = {
    "SampleParameters": {"SampleId": "HG002.pcr_free", "Sex": "Male"},
    "LocusResults": {
        # het call: ref count = (591751-591733)/1 = 18; alleles 20 and 97
        "1-591733-591751-A": {
            "LocusId": "1-591733-591751-A",
            "Variants": {
                "1-591733-591751-A": {
                    "VariantId": "1-591733-591751-A",
                    "RepeatUnit": "A",
                    "ReferenceRegion": "chr1:591733-591751",
                    "Genotype": "20/97",
                    "GenotypeConfidenceInterval": "20-20/95-141",
                }
            }
        },
        # hom-ref call: ref count = (650423-650383)/2 = 20; both alleles 20 -> GT 0/0, no ALT
        "1-650383-650423-AC": {
            "LocusId": "1-650383-650423-AC",
            "Variants": {
                "1-650383-650423-AC": {
                    "VariantId": "1-650383-650423-AC",
                    "RepeatUnit": "AC",
                    "ReferenceRegion": "chr1:650383-650423",
                    "Genotype": "20/20",
                    "GenotypeConfidenceInterval": "20-20/20-20",
                }
            }
        },
        # hemizygous call (single allele, e.g. male chrX): must be expanded to a duplicated diploid genotype
        "X-1000-1010-A": {
            "LocusId": "X-1000-1010-A",
            "Variants": {
                "X-1000-1010-A": {
                    "VariantId": "X-1000-1010-A",
                    "RepeatUnit": "A",
                    "ReferenceRegion": "chrX:1000-1010",
                    "Genotype": "15",
                    "GenotypeConfidenceInterval": "14-16",
                }
            }
        },
    },
}


class Tests(unittest.TestCase):

    def setUp(self):
        self.json_file = tempfile.NamedTemporaryFile(mode="wt", suffix=".json", delete=False)
        json.dump(EH_JSON, self.json_file)
        self.json_file.flush()

    def test_parse(self):
        records = list(parse_expansion_hunter_json(self.json_file.name))
        self.assertEqual(len(records), 3)
        by_locus = {r[4]: r for r in records}

        # hemizygous call is duplicated into a diploid genotype (avoids EnsembleTR's diploid-only EH scorer crash)
        _, _, _, _, _, _, _, hemi_allele_counts, hemi_ci_strings = by_locus["X-1000-1010-A"]
        self.assertEqual(hemi_allele_counts, [15, 15])
        self.assertEqual(hemi_ci_strings, ["14-16", "14-16"])

        # het locus: POS = start0+1 = 591734, end = 591751, ref_count = 18, alleles [20, 97]
        sample_id, chrom, pos, end, locus_id, motif, ref_count, allele_counts, ci_strings = by_locus["1-591733-591751-A"]
        self.assertEqual(sample_id, "HG002.pcr_free")
        self.assertEqual((chrom, pos, end, motif, ref_count), ("chr1", 591734, 591751, "A", 18))
        self.assertEqual(allele_counts, [20, 97])
        self.assertEqual(ci_strings, ["20-20", "95-141"])

        # hom-ref locus: ref_count = 20, alleles [20, 20]
        _, _, _, _, _, motif2, ref_count2, allele_counts2, _ = by_locus["1-650383-650423-AC"]
        self.assertEqual((motif2, ref_count2, allele_counts2), ("AC", 20, [20, 20]))


if __name__ == "__main__":
    unittest.main()
