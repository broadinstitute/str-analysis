import argparse
import collections
import unittest

from str_analysis.extend_str_catalog_boundaries import build_output_records, catalog_record_key, \
    compute_output_path, extend_catalog_record


class FakeFasta:
    """The little bit of the pysam.FastaFile interface that extend_catalog_record actually uses."""

    def __init__(self, sequence):
        self.sequence = sequence
        self.references = ["chr1"]

    def get_reference_length(self, chrom):
        return len(self.sequence)

    def fetch(self, chrom, start, end):
        return self.sequence[max(0, start):end]


def record(locus_id, reference_region, motif="CAT"):
    return {"LocusId": locus_id, "ReferenceRegion": reference_region,
            "LocusStructure": f"({motif})*", "VariantType": "Repeat"}


def definitions(output_records):
    """The (region, motif) identity of each output record, in output order."""
    return [(r["ReferenceRegion"], r["LocusStructure"]) for r in output_records]


class TestBuildOutputRecords(unittest.TestCase):

    def setUp(self):
        self.counters = collections.Counter()

    def build(self, pairs, keep_originals):
        return build_output_records(pairs, keep_originals, self.counters)

    def test_unextended_locus_passes_through_in_both_modes(self):
        pairs = [(record("a", "chr1:100-110"), None)]
        for keep_originals in (False, True):
            output = self.build(pairs, keep_originals)
            self.assertEqual(definitions(output), [("chr1:100-110", "(CAT)*")])

    def test_extended_definition_replaces_the_original_by_default(self):
        pairs = [(record("a", "chr1:100-110"), record("a", "chr1:100-130"))]
        output = self.build(pairs, keep_originals=False)
        self.assertEqual(definitions(output), [("chr1:100-130", "(CAT)*")])

    def test_original_is_kept_alongside_the_extended_definition_when_asked(self):
        pairs = [(record("a", "chr1:100-110"), record("a", "chr1:100-130"))]
        output = self.build(pairs, keep_originals=True)
        self.assertEqual(definitions(output), [("chr1:100-110", "(CAT)*"), ("chr1:100-130", "(CAT)*")])

    def test_two_loci_extending_to_the_same_interval_add_only_one_definition(self):
        pairs = [
            (record("a", "chr1:100-110"), record("a", "chr1:100-130")),
            (record("b", "chr1:115-125"), record("b", "chr1:100-130")),
        ]
        output = self.build(pairs, keep_originals=True)
        self.assertEqual(definitions(output), [
            ("chr1:100-110", "(CAT)*"), ("chr1:100-130", "(CAT)*"), ("chr1:115-125", "(CAT)*")])
        self.assertEqual(sum(count for label, count in self.counters.items()
                             if "another locus extended to the same interval" in label), 1)

    def test_same_interval_but_a_different_motif_is_not_a_duplicate(self):
        pairs = [
            (record("a", "chr1:100-110", motif="CAT"), record("a", "chr1:100-130", motif="CAT")),
            (record("b", "chr1:115-125", motif="CAG"), record("b", "chr1:100-130", motif="CAG")),
        ]
        output = self.build(pairs, keep_originals=False)
        self.assertEqual(definitions(output),
                         [("chr1:100-130", "(CAT)*"), ("chr1:100-130", "(CAG)*")])

    def test_extended_definition_matching_a_locus_already_in_the_catalog_is_not_added(self):
        # the locus it collides with appears later in the file, so this only works if the whole
        # catalog is considered before deciding
        pairs = [
            (record("a", "chr1:100-110"), record("a", "chr1:100-130")),
            (record("b", "chr1:100-130"), None),
        ]
        output = self.build(pairs, keep_originals=True)
        self.assertEqual(definitions(output), [("chr1:100-110", "(CAT)*"), ("chr1:100-130", "(CAT)*")])
        self.assertEqual(sum(count for label, count in self.counters.items()
                             if "already in the catalog" in label), 1)

    def test_replace_mode_drops_an_extension_that_duplicates_a_locus_already_in_the_output(self):
        # "a" extends onto what "b" already is, so only "b"'s definition survives
        pairs = [
            (record("a", "chr1:100-110"), record("a", "chr1:100-130")),
            (record("b", "chr1:100-130"), None),
        ]
        output = self.build(pairs, keep_originals=False)
        self.assertEqual(definitions(output), [("chr1:100-130", "(CAT)*")])

    def test_two_loci_extending_onto_each_other_collapse_to_one_definition(self):
        pairs = [
            (record("a", "chr1:100-110"), record("a", "chr1:100-130")),
            (record("b", "chr1:115-125"), record("b", "chr1:100-130")),
        ]
        output = self.build(pairs, keep_originals=False)
        self.assertEqual(definitions(output), [("chr1:100-130", "(CAT)*")])
        self.assertEqual([r["LocusId"] for r in output], ["a"])

    def test_extended_definitions_are_renamed_after_the_interval_they_now_cover(self):
        pairs = [(record("some_gene", "chr1:100-110"), record("some_gene", "chr1:100-130"))]
        output = self.build(pairs, keep_originals=True)
        self.assertEqual([r["LocusId"] for r in output], ["some_gene", "1-100-130-CAT"])

    def test_replacing_a_locus_preserves_its_locus_id(self):
        # the original is gone, so its name is free, and a name like HTT is worth more than a
        # coordinate string
        pairs = [(record("HTT", "chr1:100-110"), record("HTT", "chr1:100-130"))]
        output = self.build(pairs, keep_originals=False)
        self.assertEqual([r["LocusId"] for r in output], ["HTT"])
        self.assertEqual(definitions(output), [("chr1:100-130", "(CAT)*")])

    def test_output_locus_ids_stay_unique_when_a_generated_name_is_already_taken(self):
        pairs = [
            (record("1-100-130-CAT", "chr1:200-210"), None),
            (record("a", "chr1:100-110"), record("a", "chr1:100-130")),
        ]
        output = self.build(pairs, keep_originals=True)
        locus_ids = [r["LocusId"] for r in output]
        self.assertEqual(len(set(locus_ids)), len(locus_ids))
        self.assertIn("1-100-130-CAT_extended", locus_ids)

    def test_key_ignores_the_locus_id(self):
        self.assertEqual(catalog_record_key(record("a", "chr1:100-110")),
                         catalog_record_key(record("b", "chr1:100-110")))


class TestExtendCatalogRecord(unittest.TestCase):

    def setUp(self):
        self.counters = collections.Counter()
        self.args = argparse.Namespace(max_repeats_in_gap=2, min_purity_of_new_sequence=0.9,
                                       verbose=False)
        # 7 clean CAT copies at 100-121, then 3 more in the right flank, and nothing repeat-like left
        self.fasta = FakeFasta("T" * 100 + "CAT" * 10 + "G" * 100)

    def test_a_locus_is_extended_over_the_repeat_that_continues_past_its_boundary(self):
        extended = extend_catalog_record(record("a", "chr1:100-121"), self.fasta, self.args,
                                         self.counters)
        self.assertEqual(extended["ReferenceRegion"], "chr1:100-130")

    def test_every_other_field_is_passed_through_unmodified(self):
        # ReferenceRegion is the only field this tool rewrites: optional fields, annotations included,
        # are the caller's to refresh, and inventing policy for them here would surprise more than help
        input_record = record("a", "chr1:100-121")
        input_record.update({"VariantId": "a", "MainReferenceRegion": "chr1:100-121",
                             "NumRepeatsInReference": 7, "ReferenceRepeatPurity": 1.0,
                             "Diseases": [{"Name": "Example", "NormalMax": 20}], "Gene": "SOME_GENE"})
        extended = extend_catalog_record(input_record, self.fasta, self.args, self.counters)
        self.assertEqual(extended["ReferenceRegion"], "chr1:100-130")
        for field, value in input_record.items():
            if field != "ReferenceRegion":
                self.assertEqual(extended[field], value, f"{field} should have been passed through")

    def test_a_locus_the_rule_declines_to_extend_returns_none(self):
        self.assertIsNone(extend_catalog_record(record("a", "chr1:10-20"), self.fasta, self.args,
                                                self.counters))

    def test_several_motifs_sharing_one_interval_is_a_hard_error(self):
        # how a TRGT-format BED row parses: one interval, several motifs. Where each motif starts
        # inside it is unknown, so the locus can neither be extended nor written back as one BED row.
        trgt_style = {"LocusId": "t", "ReferenceRegion": "chr1:100-121",
                      "LocusStructure": "(CAT)*(CAG)*", "VariantType": "Repeat"}
        with self.assertRaises(ValueError):
            extend_catalog_record(trgt_style, self.fasta, self.args, self.counters)


class TestOutputPath(unittest.TestCase):

    def test_output_path_is_always_a_bgzipped_bed(self):
        self.assertEqual(compute_output_path("catalog.bed"), "catalog.extended.bed.gz")
        self.assertEqual(compute_output_path("catalog.bed.gz"), "catalog.extended.bed.gz")
        self.assertEqual(compute_output_path("dir/catalog.bed.gz"), "catalog.extended.bed.gz")
        self.assertEqual(compute_output_path("my.catalog.v2.bed.bgz"), "my.catalog.v2.extended.bed.gz")


if __name__ == "__main__":
    unittest.main()
