import argparse
import collections
import unittest

from str_analysis.extend_str_catalog_boundaries import build_output_records, catalog_record_key, \
    compute_output_path, discard_extensions_outdone_by_the_input, extend_catalog_record, \
    select_longest_non_overlapping, sort_catalog_records


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
                             if "overlapping extended definition that was kept" in label), 1)

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
        self.assertEqual([r["LocusId"] for r in output], ["1-100-130-CAT"])

    def test_a_locus_whose_extension_is_dropped_is_put_back_when_it_would_otherwise_be_lost(self):
        # b's extension loses to a's, which is longer, but a's stops at 300 and b covered up to 350
        pairs = [
            (record("a", "chr1:100-110"), record("a", "chr1:100-300")),
            (record("b", "chr1:280-350"), record("b", "chr1:250-360")),
        ]
        output = self.build(pairs, keep_originals=False)
        self.assertEqual(definitions(output),
                         [("chr1:100-300", "(CAT)*"), ("chr1:280-350", "(CAT)*")])

    def test_a_locus_whose_extension_is_dropped_is_not_put_back_when_it_is_covered(self):
        pairs = [
            (record("a", "chr1:100-110"), record("a", "chr1:100-300")),
            (record("b", "chr1:150-160"), record("b", "chr1:120-280")),
        ]
        output = self.build(pairs, keep_originals=False)
        self.assertEqual(definitions(output), [("chr1:100-300", "(CAT)*")])

    def test_keeping_the_originals_never_puts_one_back_twice(self):
        pairs = [
            (record("a", "chr1:100-110"), record("a", "chr1:100-300")),
            (record("b", "chr1:280-350"), record("b", "chr1:250-360")),
        ]
        output = self.build(pairs, keep_originals=True)
        self.assertEqual(definitions(output), [("chr1:100-110", "(CAT)*"), ("chr1:100-300", "(CAT)*"),
                                               ("chr1:280-350", "(CAT)*")])

    def test_a_locus_is_not_put_back_when_some_other_kept_definition_already_is_it(self):
        # r's extension displaces p's, but s extended onto exactly what p covered, so putting p back
        # would write the same definition twice under two different LocusIds
        pairs = [
            (record("p", "chr1:200-230"), record("p", "chr1:200-330")),
            (record("s", "chr1:210-224"), record("s", "chr1:200-230")),
            (record("r", "chr1:310-360"), record("r", "chr1:300-500")),
        ]
        output = self.build(pairs, keep_originals=False)
        self.assertEqual(definitions(output),
                         [("chr1:200-230", "(CAT)*"), ("chr1:300-500", "(CAT)*")])

    def test_a_locus_is_not_put_back_when_two_kept_definitions_cover_it_between_them(self):
        # a and c are book-ended at 200, so together they describe every base b's original covered
        pairs = [
            (record("a", "chr1:80-90"), record("a", "chr1:50-200")),
            (record("b", "chr1:195-205"), record("b", "chr1:150-250")),
            (record("c", "chr1:340-350"), record("c", "chr1:200-350")),
        ]
        output = self.build(pairs, keep_originals=False)
        self.assertEqual(definitions(output), [("chr1:50-200", "(CAT)*"), ("chr1:200-350", "(CAT)*")])

    def test_a_gap_between_two_kept_definitions_still_puts_the_locus_back(self):
        # the same shape, except a and c leave 199-201 undescribed, which is where b's original sits
        pairs = [
            (record("a", "chr1:80-90"), record("a", "chr1:50-199")),
            (record("b", "chr1:195-205"), record("b", "chr1:150-250")),
            (record("c", "chr1:340-350"), record("c", "chr1:201-350")),
        ]
        output = self.build(pairs, keep_originals=False)
        self.assertEqual(definitions(output), [("chr1:50-199", "(CAT)*"), ("chr1:195-205", "(CAT)*"),
                                               ("chr1:201-350", "(CAT)*")])

    def test_a_rotation_of_the_same_motif_is_the_same_repeat_when_deduplicating(self):
        pairs = [
            (record("a", "chr1:100-110", motif="CAG"), record("a", "chr1:100-130", motif="CAG")),
            (record("b", "chr1:120-130", motif="GCA"), record("b", "chr1:99-140", motif="GCA")),
        ]
        output = self.build(pairs, keep_originals=False)
        self.assertEqual(definitions(output), [("chr1:99-140", "(GCA)*")])

    def test_extended_definitions_are_renamed_after_the_interval_they_now_cover(self):
        pairs = [(record("some_gene", "chr1:100-110"), record("some_gene", "chr1:100-130"))]
        output = self.build(pairs, keep_originals=True)
        self.assertEqual([r["LocusId"] for r in output], ["some_gene", "1-100-130-CAT"])

    def test_replacing_a_locus_still_renames_it_after_the_new_interval(self):
        # the original name is free once the original is gone, but it described the locus before it
        # grew, so the extended definition is named after what it now covers instead
        pairs = [(record("HTT", "chr1:100-110"), record("HTT", "chr1:100-130"))]
        output = self.build(pairs, keep_originals=False)
        self.assertEqual([r["LocusId"] for r in output], ["1-100-130-CAT"])
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
                                       min_reference_purity=0.9, verbose=False)
        # 7 clean CAT copies at 100-121, then 3 more in the right flank, and nothing repeat-like
        # around them. The padding is all G, a base CAT doesn't contain, so that neither boundary can
        # pick up a partial copy and the plain extension is what's being tested.
        self.fasta = FakeFasta("G" * 100 + "CAT" * 10 + "G" * 100)

    def test_a_locus_is_extended_over_the_repeat_that_continues_past_its_boundary(self):
        extended = extend_catalog_record(record("a", "chr1:100-121"), self.fasta, self.args,
                                         self.counters)
        self.assertEqual(extended["ReferenceRegion"], "chr1:100-130")

    def test_every_other_field_is_passed_through_unmodified(self):
        # ReferenceRegion is the only field extend_catalog_record rewrites: optional fields,
        # annotations included, are the caller's to refresh, and inventing policy for them here would
        # surprise more than help. The motif is never re-phased, so LocusStructure is passed through
        # like everything else, "+" quantifier and all.
        input_record = record("a", "chr1:100-121")
        input_record.update({"VariantId": "a", "MainReferenceRegion": "chr1:100-121",
                             "LocusStructure": "(CAT)+", "NumRepeatsInReference": 7,
                             "ReferenceRepeatPurity": 1.0,
                             "Diseases": [{"Name": "Example", "NormalMax": 20}], "Gene": "SOME_GENE"})
        extended = extend_catalog_record(input_record, self.fasta, self.args, self.counters)
        self.assertEqual(extended["ReferenceRegion"], "chr1:100-130")
        for field, value in input_record.items():
            if field != "ReferenceRegion":
                self.assertEqual(extended[field], value, f"{field} should have been passed through")

    def test_a_locus_the_rule_declines_to_extend_returns_none(self):
        self.assertIsNone(extend_catalog_record(record("a", "chr1:10-20"), self.fasta, self.args,
                                                self.counters))

    def test_an_extended_locus_is_polished_over_the_partial_copy_at_its_new_right_boundary(self):
        # 10 CAT copies at 101-131 with a leftover "CA" after them, so the extension itself stops at
        # 131 and the polish widens it over the "CA"
        fasta = FakeFasta("G" * 100 + "T" + "CAT" * 10 + "CA" + "G" * 100)
        extended = extend_catalog_record(record("a", "chr1:101-122"), fasta, self.args, self.counters)
        self.assertEqual(extended["ReferenceRegion"], "chr1:101-133")
        # the "T" at 100 continues the repeat too, but taking it would start the locus part-way
        # through the motif and make it a TCA locus, so the left boundary stops short of it
        self.assertEqual(extended["LocusStructure"], "(CAT)*")

    def test_the_left_boundary_moves_only_by_whole_copies_so_the_motif_keeps_its_frame(self):
        # EP400's shape: 10 CAT copies at 101-131, preceded by the "T" that ends the copy before
        # them. The locus grows left over the 3 whole copies it is short of and stops at 101.
        fasta = FakeFasta("G" * 100 + "T" + "CAT" * 10 + "G" * 100)
        extended = extend_catalog_record(record("a", "chr1:110-131"), fasta, self.args, self.counters)
        self.assertEqual(extended["ReferenceRegion"], "chr1:101-131")
        self.assertEqual(extended["LocusStructure"], "(CAT)*")

    def test_a_locus_the_rule_declines_to_extend_is_left_alone_even_next_to_a_partial_copy(self):
        # the "CA" just right of this locus does continue the repeat, but polishing is only for loci
        # the rule actually extended, and nothing here reaches far enough to be extended
        fasta = FakeFasta("G" * 100 + "CAT" * 3 + "CA" + "G" * 100)
        self.assertIsNone(extend_catalog_record(record("a", "chr1:100-109"), fasta, self.args,
                                                self.counters))

    def test_a_locus_that_is_a_poor_match_for_its_own_motif_is_left_alone(self):
        # three clean CAT copies sit immediately right of this locus, at 121, but the locus itself
        # matches its own motif at only 3 of its 21 bases, so the rule is not asked where its
        # boundaries belong
        fasta = FakeFasta("G" * 100 + "CATGGGGGGGGGGGGGGGGGG" + "CAT" * 3 + "G" * 100)
        self.assertIsNone(extend_catalog_record(record("a", "chr1:100-121"), fasta, self.args,
                                                self.counters))
        self.assertEqual(sum(count for label, count in self.counters.items()
                             if "not pure enough" in label), 1)

    def test_lowering_the_reference_purity_threshold_lets_that_same_locus_through(self):
        fasta = FakeFasta("G" * 100 + "CATGGGGGGGGGGGGGGGGGG" + "CAT" * 3 + "G" * 100)
        self.args.min_reference_purity = 0.1
        extended = extend_catalog_record(record("a", "chr1:100-121"), fasta, self.args,
                                         self.counters)
        self.assertEqual(extended["ReferenceRegion"], "chr1:100-130")

    def test_a_locus_with_a_few_interruptions_is_still_pure_enough(self):
        # 1 mismatched base out of 21 is 0.952, comfortably over the default threshold
        fasta = FakeFasta("G" * 100 + "CATCATCAGCATCATCATCAT" + "CAT" * 3 + "G" * 100)
        extended = extend_catalog_record(record("a", "chr1:100-121"), fasta, self.args,
                                         self.counters)
        self.assertEqual(extended["ReferenceRegion"], "chr1:100-130")

    def test_several_motifs_sharing_one_interval_is_a_hard_error(self):
        # how a TRGT-format BED row parses: one interval, several motifs. Where each motif starts
        # inside it is unknown, so the locus can neither be extended nor written back as one BED row.
        trgt_style = {"LocusId": "t", "ReferenceRegion": "chr1:100-121",
                      "LocusStructure": "(CAT)*(CAG)*", "VariantType": "Repeat"}
        with self.assertRaises(ValueError):
            extend_catalog_record(trgt_style, self.fasta, self.args, self.counters)


class TestDiscardingExtensionsOutdoneByTheInput(unittest.TestCase):

    def setUp(self):
        self.counters = collections.Counter()

    def surviving(self, *pairs):
        """The locus ids whose extension is still there after the input catalog has had its say."""
        return [record["LocusId"]
                for record, extended in discard_extensions_outdone_by_the_input(list(pairs),
                                                                               self.counters)
                if extended is not None]

    def test_a_longer_overlapping_input_locus_takes_the_extension_away(self):
        self.assertEqual(self.surviving(
            (record("a", "chr1:100-110"), record("a", "chr1:100-130")),
            (record("long", "chr1:120-200"), None)), [])
        self.assertEqual(sum(count for label, count in self.counters.items()
                             if "already covers that repeat" in label), 1)

    def test_a_locus_own_original_never_takes_its_extension_away(self):
        self.assertEqual(self.surviving((record("a", "chr1:100-110"),
                                         record("a", "chr1:100-130"))), ["a"])

    def test_an_input_locus_that_is_no_longer_than_the_extension_leaves_it_alone(self):
        # same length, so the extension still says as much about the repeat as the catalog did
        self.assertEqual(self.surviving(
            (record("a", "chr1:100-110"), record("a", "chr1:100-130")),
            (record("same_length", "chr1:120-150"), None)), ["a"])

    def test_a_longer_input_locus_that_does_not_overlap_leaves_it_alone(self):
        self.assertEqual(self.surviving(
            (record("a", "chr1:100-110"), record("a", "chr1:100-130")),
            (record("long", "chr1:130-300"), None)), ["a"])

    def test_a_longer_input_locus_with_a_different_repeat_leaves_it_alone(self):
        self.assertEqual(self.surviving(
            (record("a", "chr1:100-110", motif="CAT"), record("a", "chr1:100-130", motif="CAT")),
            (record("long", "chr1:120-200", motif="CAG"), None)), ["a"])

    def test_a_rotation_of_the_motif_is_the_same_repeat(self):
        self.assertEqual(self.surviving(
            (record("a", "chr1:100-110", motif="CAG"), record("a", "chr1:100-130", motif="CAG")),
            (record("long", "chr1:120-200", motif="GCA"), None)), [])

    def test_a_longer_input_locus_on_another_chromosome_leaves_it_alone(self):
        self.assertEqual(self.surviving(
            (record("a", "chr1:100-110"), record("a", "chr1:100-130")),
            (record("long", "chr2:120-200"), None)), ["a"])

    def test_an_input_locus_that_was_itself_extended_still_counts(self):
        # what the catalog already said is what matters, whatever became of that locus afterwards
        self.assertEqual(self.surviving(
            (record("a", "chr1:100-110"), record("a", "chr1:100-130")),
            (record("long", "chr1:120-200"), record("long", "chr1:120-260"))), ["long"])

    def test_one_input_locus_can_take_away_several_extensions(self):
        self.assertEqual(self.surviving(
            (record("a", "chr1:100-110"), record("a", "chr1:100-130")),
            (record("b", "chr1:150-160"), record("b", "chr1:140-180")),
            (record("long", "chr1:120-300"), None)), [])

    def test_the_scan_does_not_stop_at_the_first_extension_it_leaves_alone(self):
        # walking back over "huge", which this input locus is not long enough to take away, to reach
        # "target", which it is
        self.assertEqual(self.surviving(
            (record("target", "chr1:400-420"), record("target", "chr1:400-430")),
            (record("huge", "chr1:500-1100"), record("huge", "chr1:500-1200")),
            (record("long", "chr1:300-800"), None)), ["huge"])

    def test_the_scan_stops_once_nothing_further_back_can_reach(self):
        # "early" ends well before this input locus starts, so it is never in the running
        self.assertEqual(self.surviving(
            (record("early", "chr1:100-115"), record("early", "chr1:100-120")),
            (record("target", "chr1:400-420"), record("target", "chr1:400-430")),
            (record("long", "chr1:300-800"), None)), ["early"])


class TestSelectLongestNonOverlapping(unittest.TestCase):

    def dropped(self, *records):
        """The ids of the records that lose their place to a longer overlapping one."""
        kept, _ = select_longest_non_overlapping(list(records))
        return sorted(r["LocusId"] for index, r in enumerate(records) if index not in kept)

    def kept_intervals(self, *records):
        return dict(select_longest_non_overlapping(list(records))[1])

    def test_nothing_is_dropped_when_nothing_overlaps(self):
        self.assertEqual(self.dropped(record("a", "chr1:100-200"), record("b", "chr1:200-300")), [])

    def test_the_shorter_of_two_overlapping_definitions_is_dropped(self):
        self.assertEqual(self.dropped(record("a", "chr1:100-130"), record("b", "chr1:120-200")), ["a"])

    def test_a_chain_of_overlaps_only_collapses_where_it_has_to(self):
        # b is dropped for c, but a only reached c through b, so a keeps its place
        self.assertEqual(self.dropped(record("a", "chr1:100-200"), record("b", "chr1:150-260"),
                                      record("c", "chr1:250-500")), ["b"])

    def test_a_definition_straddling_two_kept_ones_is_dropped(self):
        # a short definition sitting across where two long, book-ended ones meet
        self.assertEqual(self.dropped(record("a", "chr1:199-201"), record("b", "chr1:100-200"),
                                      record("c", "chr1:200-300")), ["a"])

    def test_rotations_of_one_motif_are_the_same_repeat(self):
        self.assertEqual(self.dropped(record("a", "chr1:100-130", motif="CAG"),
                                      record("b", "chr1:105-140", motif="GCA")), ["a"])

    def test_a_different_canonical_motif_is_a_different_repeat(self):
        self.assertEqual(self.dropped(record("a", "chr1:100-130", motif="CAT"),
                                      record("b", "chr1:105-140", motif="CAG")), [])

    def test_the_same_interval_on_another_chromosome_does_not_overlap(self):
        self.assertEqual(self.dropped(record("a", "chr1:100-200"), record("b", "chr2:100-200")), [])

    def test_equally_long_definitions_are_resolved_by_coordinate_not_input_order(self):
        for records in ([record("a", "chr1:110-140"), record("b", "chr1:100-130")],
                        [record("b", "chr1:100-130"), record("a", "chr1:110-140")]):
            self.assertEqual(self.dropped(*records), ["a"])

    def test_identical_definitions_are_resolved_by_locus_id_not_input_order(self):
        for records in ([record("b", "chr1:100-130"), record("a", "chr1:100-130")],
                        [record("a", "chr1:100-130"), record("b", "chr1:100-130")]):
            self.assertEqual(self.dropped(*records), ["b"])

    def test_what_was_kept_is_reported_per_chromosome_and_canonical_motif(self):
        self.assertEqual(self.kept_intervals(record("a", "chr1:100-200", motif="CAG"),
                                             record("b", "chr1:300-400", motif="GCA"),
                                             record("c", "chr2:100-200", motif="CAT")),
                         {("chr1", "AGC"): ([100, 300], [200, 400]),
                          ("chr2", "ATC"): ([100], [200])})


class TestIupacMotifsAreNotExtended(unittest.TestCase):
    """A motif written with IUPAC codes matches more sequence than a concrete one, so the rule that
    reads the flank for more copies of the motif cannot be trusted on it.
    """

    def setUp(self):
        self.counters = collections.Counter()
        self.args = argparse.Namespace(
            max_repeats_in_gap=2, min_purity_of_new_sequence=0.9, min_reference_purity=0.9,
            verbose=False)

    def extend(self, motif, sequence):
        return extend_catalog_record(
            record("locus", "chr1:0-6", motif=motif), FakeFasta(sequence), self.args, self.counters)

    def test_a_gcn_motif_is_left_alone(self):
        # GCAGCA followed by more of the same: a concrete GCA motif extends here
        self.assertIsNone(self.extend("GCN", "GCAGCAGCAGCA"))
        self.assertEqual(self.counters["loci not considered because the motif contains IUPAC codes"], 1)

    def test_the_same_locus_with_a_concrete_motif_is_extended(self):
        extended = self.extend("GCA", "GCAGCAGCAGCA")
        self.assertIsNotNone(extended)
        self.assertEqual(self.counters["loci not considered because the motif contains IUPAC codes"], 0)

    def test_other_iupac_codes_are_caught_too(self):
        for motif in ("AARRG", "CWG", "NNN"):
            self.assertIsNone(self.extend(motif, "AAAAGAAAAGAAAAG"))


class TestOnlyExtendedDefinitions(unittest.TestCase):
    """--only-output-extended-loci turns the output into the set of definitions the rule adds to the
    input catalog, for use as a separate source of loci appended to that catalog.
    """

    def setUp(self):
        self.counters = collections.Counter()

    def build(self, pairs):
        return build_output_records(pairs, False, self.counters, only_extended_definitions=True)

    def test_only_the_extended_definitions_are_written(self):
        pairs = [(record("a", "chr1:100-110"), record("a", "chr1:100-130")),
                 (record("b", "chr1:500-510"), None)]
        output = self.build(pairs)
        self.assertEqual(definitions(output), [("chr1:100-130", "(CAT)*")])

    def test_an_unextended_catalog_produces_nothing(self):
        pairs = [(record("a", "chr1:100-110"), None), (record("b", "chr1:500-510"), None)]
        self.assertEqual(self.build(pairs), [])

    def test_a_definition_the_input_already_has_is_still_dropped(self):
        # nothing here may duplicate the catalog this output gets appended to
        pairs = [(record("a", "chr1:100-110"), record("a", "chr1:100-130")),
                 (record("b", "chr1:100-130"), None)]
        self.assertEqual(self.build(pairs), [])

    def test_a_displaced_locus_is_not_put_back(self):
        # in whole-catalog mode 'b' would come back so the catalog keeps covering 300-350; here the
        # input catalog still holds it, so putting it back would duplicate a locus that never left
        pairs = [(record("a", "chr1:100-110"), record("a", "chr1:100-300")),
                 (record("b", "chr1:320-350"), record("b", "chr1:150-350"))]
        output = self.build(pairs)
        self.assertEqual([r["LocusId"] for r in output], ["1-100-300-CAT"])

    def test_a_definition_matching_the_original_of_an_extended_locus_is_dropped(self):
        # 'p' is extended, so in whole-catalog mode its own definition would be gone and something
        # else growing onto it would be worth writing. Here 'p' stays in the input catalog, so 's'
        # growing into exactly 'p' has nothing to add. 'p' is dropped in favor of 'r' first, which is
        # what leaves 's' as the definition that would otherwise be written.
        pairs = [(record("p", "chr1:200-230"), record("p", "chr1:200-330")),
                 (record("s", "chr1:210-224"), record("s", "chr1:200-230")),
                 (record("r", "chr1:310-360"), record("r", "chr1:300-500"))]
        output = self.build(pairs)
        self.assertEqual(definitions(output), [("chr1:300-500", "(CAT)*")])


class TestSortCatalogRecords(unittest.TestCase):
    """The JSON output feeds group_overlapping_loci downstream, which raises on any pair of records
    whose coordinates decrease, so the catalog has to leave here in order.
    """

    def test_a_left_extension_is_placed_before_the_locus_it_grew_from(self):
        # what build_output_records emits with keep_originals: the original, then its extension,
        # which starts 24bp earlier, a whole number of motif copies
        records = [record("orig", "chr12:132062548-132062611"),
                   record("ext", "chr12:132062524-132062611")]
        self.assertEqual([r["ReferenceRegion"] for r in sort_catalog_records(records)],
                         ["chr12:132062524-132062611", "chr12:132062548-132062611"])

    def test_chromosomes_stay_contiguous_and_keep_the_order_they_arrived_in(self):
        records = [record("a", "chr10:500-510"), record("b", "chr2:100-110"),
                   record("c", "chr10:100-110")]
        self.assertEqual([r["ReferenceRegion"] for r in sort_catalog_records(records)],
                         ["chr10:100-110", "chr10:500-510", "chr2:100-110"])

    def test_records_at_the_same_interval_are_ordered_deterministically(self):
        records = [record("b", "chr1:100-110", motif="CAT"), record("a", "chr1:100-110", motif="CAT")]
        self.assertEqual([r["LocusId"] for r in sort_catalog_records(records)], ["a", "b"])
        # sorting the reversed input gives the same answer
        self.assertEqual([r["LocusId"] for r in sort_catalog_records(records[::-1])], ["a", "b"])


class TestOutputPath(unittest.TestCase):

    def test_a_bed_input_gives_a_bgzipped_bed(self):
        self.assertEqual(compute_output_path("catalog.bed"), "catalog.extended.bed.gz")
        self.assertEqual(compute_output_path("catalog.bed.gz"), "catalog.extended.bed.gz")
        self.assertEqual(compute_output_path("dir/catalog.bed.gz"), "catalog.extended.bed.gz")
        self.assertEqual(compute_output_path("my.catalog.v2.bed.bgz"), "my.catalog.v2.extended.bed.gz")

    def test_a_json_input_gives_a_gzipped_json(self):
        self.assertEqual(compute_output_path("catalog.json"), "catalog.extended.json.gz")
        self.assertEqual(compute_output_path("catalog.json.gz"), "catalog.extended.json.gz")
        self.assertEqual(compute_output_path("dir/catalog.json.gz"), "catalog.extended.json.gz")
        self.assertEqual(compute_output_path("my.catalog.v2.json.bgz"), "my.catalog.v2.extended.json.gz")


class TestExtendedLocusIds(unittest.TestCase):
    """Which extended definitions get renamed, now that a JSON output puts locus ids in the file."""

    def setUp(self):
        self.counters = collections.Counter()

    def build(self, pairs, keep_originals):
        return build_output_records(pairs, keep_originals, self.counters)

    def test_an_extended_definition_is_named_after_its_new_interval(self):
        for keep_originals in (False, True):
            pairs = [(record("1-100-110-CAT", "chr1:100-110"), record("1-100-110-CAT", "chr1:100-130"))]
            output = build_output_records(pairs, keep_originals, collections.Counter())
            self.assertEqual([r["LocusId"] for r in output][-1], "1-100-130-CAT")

    def test_a_name_is_replaced_by_coordinates_too(self):
        for keep_originals in (False, True):
            pairs = [(record("HTT", "chr1:100-110"), record("HTT", "chr1:100-130"))]
            output = build_output_records(pairs, keep_originals, collections.Counter())
            self.assertEqual([r["LocusId"] for r in output][-1], "1-100-130-CAT")

    def test_keeping_the_original_leaves_its_own_id_alone(self):
        pairs = [(record("HTT", "chr1:100-110"), record("HTT", "chr1:100-130"))]
        output = self.build(pairs, keep_originals=True)
        self.assertEqual([r["LocusId"] for r in output], ["HTT", "1-100-130-CAT"])

    def test_an_unextended_locus_keeps_its_id(self):
        pairs = [(record("1-100-110-CAT", "chr1:100-110"), None), (record("HTT", "chr1:200-210"), None)]
        output = self.build(pairs, keep_originals=False)
        self.assertEqual([r["LocusId"] for r in output], ["1-100-110-CAT", "HTT"])


if __name__ == "__main__":
    unittest.main()
