import argparse
import collections
import unittest

from str_analysis.extend_str_catalog_boundaries import build_output_records, catalog_record_key, \
    compute_output_path, extend_catalog_record, rename_variant_ids


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

    def test_replace_mode_keeps_a_locus_whose_extension_would_duplicate_another(self):
        # "a" extends onto what "b" already is. The extended definition is dropped as a duplicate, and
        # "a" falls back to its own original definition rather than disappearing from the catalog.
        pairs = [
            (record("a", "chr1:100-110"), record("a", "chr1:100-130")),
            (record("b", "chr1:100-130"), None),
        ]
        output = self.build(pairs, keep_originals=False)
        self.assertEqual(definitions(output),
                         [("chr1:100-110", "(CAT)*"), ("chr1:100-130", "(CAT)*")])
        self.assertEqual([r["LocusId"] for r in output], ["a", "b"])

    def test_no_locus_is_ever_dropped_when_two_loci_extend_onto_each_other(self):
        # both extend to the same interval: one contributes the extended definition, the other falls
        # back to its original, so the catalog still describes two loci
        pairs = [
            (record("a", "chr1:100-110"), record("a", "chr1:100-130")),
            (record("b", "chr1:115-125"), record("b", "chr1:100-130")),
        ]
        output = self.build(pairs, keep_originals=False)
        self.assertEqual(len(output), 2)
        self.assertEqual([r["LocusId"] for r in output], ["a", "b"])
        self.assertEqual(definitions(output),
                         [("chr1:100-130", "(CAT)*"), ("chr1:115-125", "(CAT)*")])

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

    def test_variant_id_follows_the_locus_id_when_the_extended_definition_is_renamed(self):
        original = record("a", "chr1:100-110")
        original["VariantId"] = "a"
        extended = record("a", "chr1:100-130")
        extended["VariantId"] = "a"
        output = self.build([(original, extended)], keep_originals=True)
        variant_ids = [r["VariantId"] for r in output]
        self.assertEqual(len(set(variant_ids)), 2, f"VariantIds collide: {variant_ids}")
        self.assertEqual(variant_ids[1], output[1]["LocusId"])

    def test_multi_region_records_are_keyed_on_all_of_their_regions(self):
        multi = {"LocusId": "m", "ReferenceRegion": ["chr1:100-110", "chr1:110-120"],
                 "LocusStructure": "(CAT)*(CAG)*", "VariantType": "Repeat"}
        self.assertEqual(catalog_record_key(multi),
                         (("chr1:100-110", "chr1:110-120"), ("CAT", "CAG")))

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

    def test_main_reference_region_moves_with_the_boundary_it_names(self):
        input_record = record("a", "chr1:100-121")
        input_record["MainReferenceRegion"] = "chr1:100-121"
        extended = extend_catalog_record(input_record, self.fasta, self.args, self.counters)
        self.assertEqual(extended["MainReferenceRegion"], "chr1:100-130")


    def test_annotations_computed_from_the_old_interval_are_dropped(self):
        # listed out in full rather than derived from the constant, so that dropping a field from the
        # constant makes this test fail instead of quietly narrowing with it
        stale_annotations = ["NumRepeatsInReference", "ReferenceRepeatPurity", "NsInFlanks",
                             "LeftFlankMappability", "RightFlankMappability", "FlanksAndLocusMappability",
                             "HighestPurityMotif", "HighestPurityMotifPurity", "HighestPurityMotifQuality",
                             "TandemRepeatFinderMotif", "TandemRepeatFinderMotifQuality",
                             "KnownDiseaseAssociatedLocus"]
        input_record = record("a", "chr1:100-121")
        input_record.update({annotation: 1 for annotation in stale_annotations})
        input_record["Gene"] = "SOME_GENE"
        extended = extend_catalog_record(input_record, self.fasta, self.args, self.counters)
        for annotation in stale_annotations:
            self.assertNotIn(annotation, extended)
        # Gene does not depend on the exact boundary, so it survives
        self.assertEqual(extended["Gene"], "SOME_GENE")
        # and the input record itself is untouched, since its own interval did not move
        self.assertEqual(input_record["NumRepeatsInReference"], 1)

    def test_a_locus_the_rule_declines_to_extend_returns_none(self):
        self.assertIsNone(extend_catalog_record(record("a", "chr1:10-20"), self.fasta, self.args,
                                                self.counters))

    def test_several_motifs_sharing_one_region_are_passed_through_instead_of_crashing(self):
        # how a TRGT-format BED row parses: one spanning region, several motifs. Where each motif
        # starts inside the span is unknown, so the locus cannot be extended correctly.
        trgt_style = {"LocusId": "t", "ReferenceRegion": "chr1:100-121",
                      "LocusStructure": "(CAT)*(CAG)*", "VariantType": "Repeat"}
        self.assertIsNone(extend_catalog_record(trgt_style, self.fasta, self.args, self.counters))
        self.assertEqual(sum(count for label, count in self.counters.items()
                             if "do not correspond one-to-one" in label), 1)


class TestExtendMultiRegionRecord(unittest.TestCase):
    """A record with adjacent regions must move only its two OUTER boundaries."""

    def setUp(self):
        self.counters = collections.Counter()
        self.args = argparse.Namespace(max_repeats_in_gap=2, min_purity_of_new_sequence=0.9,
                                       verbose=False)
        # 20 clean CAT copies at 0-60, then nothing repeat-like, so both outer boundaries can move
        self.fasta = FakeFasta("CAT" * 20 + "G" * 50)
        self.record = {"LocusId": "m", "ReferenceRegion": ["chr1:30-39", "chr1:39-48"],
                       "LocusStructure": "(CAT)*(CAT)*", "VariantType": "Repeat"}

    def test_only_the_outer_boundaries_move(self):
        extended = extend_catalog_record(self.record, self.fasta, self.args, self.counters)
        # the interior boundary at 39 is untouched, while 30 moves left and 48 moves right
        self.assertEqual(extended["ReferenceRegion"], ["chr1:0-39", "chr1:39-60"])

    def test_main_reference_region_naming_the_first_region_follows_the_left_boundary(self):
        input_record = dict(self.record, MainReferenceRegion="chr1:30-39")
        extended = extend_catalog_record(input_record, self.fasta, self.args, self.counters)
        self.assertEqual(extended["MainReferenceRegion"], "chr1:0-39")

    def test_main_reference_region_naming_the_last_region_follows_the_right_boundary(self):
        input_record = dict(self.record, MainReferenceRegion="chr1:39-48")
        extended = extend_catalog_record(input_record, self.fasta, self.args, self.counters)
        self.assertEqual(extended["MainReferenceRegion"], "chr1:39-60")


class TestRenameVariantIds(unittest.TestCase):

    def test_a_plain_variant_id_is_replaced(self):
        self.assertEqual(rename_variant_ids("HTT", "HTT", "4-100-200-CAG"), "4-100-200-CAG")

    def test_a_list_keeps_its_shape_and_its_per_region_suffixes(self):
        self.assertEqual(rename_variant_ids(["HTT", "HTT_adjacent1"], "HTT", "NEW"),
                         ["NEW", "NEW_adjacent1"])

    def test_none_is_returned_when_nothing_mentions_the_old_id(self):
        # substitution cannot tell the two records apart here, so the caller drops the field rather
        # than leave the same VariantId on both
        self.assertIsNone(rename_variant_ids(["OTHER"], "HTT", "NEW"))
        self.assertIsNone(rename_variant_ids("OTHER", "HTT", "NEW"))

    def test_unrelated_entries_survive_when_another_entry_does_mention_the_old_id(self):
        self.assertEqual(rename_variant_ids(["HTT", "OTHER"], "HTT", "NEW"), ["NEW", "OTHER"])


class TestOutputPath(unittest.TestCase):

    def test_output_path_keeps_the_input_format_and_is_always_gzipped(self):
        self.assertEqual(compute_output_path("catalog.bed"), "catalog.extended.bed.gz")
        self.assertEqual(compute_output_path("catalog.bed.gz"), "catalog.extended.bed.gz")
        self.assertEqual(compute_output_path("dir/catalog.json"), "catalog.extended.json.gz")
        self.assertEqual(compute_output_path("catalog.json.gz"), "catalog.extended.json.gz")
        self.assertEqual(compute_output_path("my.catalog.v2.bed.bgz"), "my.catalog.v2.extended.bed.gz")


if __name__ == "__main__":
    unittest.main()
