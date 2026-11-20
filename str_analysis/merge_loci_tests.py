#!/usr/bin/env python3

"""Comprehensive tests for merge_loci.py"""

import argparse
import collections
import os
import sys
import tempfile
import unittest
from unittest import mock
from intervaltree import IntervalTree, Interval

from str_analysis.merge_loci import (
    parse_args,
    parse_variant_catalog_paths_arg,
    check_for_sufficient_overlap_and_motif_match,
    check_whether_to_merge_adjacent_loci,
    compute_locus_definition_for_exact_match,
    add_variant_catalog_to_interval_trees,
    convert_interval_trees_to_output_records,
    REQUIRED_OUTPUT_FIELDS,
    SEPARATOR_FOR_MULTIPLE_SOURCES,
)


class TestConstants(unittest.TestCase):
    """Test module constants."""

    def test_required_output_fields(self):
        """Test REQUIRED_OUTPUT_FIELDS constant."""
        expected_fields = {"LocusId", "ReferenceRegion", "LocusStructure", "VariantType"}
        self.assertEqual(REQUIRED_OUTPUT_FIELDS, expected_fields)

    def test_separator_for_multiple_sources(self):
        """Test SEPARATOR_FOR_MULTIPLE_SOURCES constant."""
        self.assertEqual(SEPARATOR_FOR_MULTIPLE_SOURCES, " ||| ")


class TestParseArgs(unittest.TestCase):
    """Test parse_args function."""

    def setUp(self):
        """Create temporary catalog files."""
        self.temp_dir = tempfile.mkdtemp()
        self.catalog1 = os.path.join(self.temp_dir, "catalog1.json")
        self.catalog2 = os.path.join(self.temp_dir, "catalog2.bed")

        with open(self.catalog1, 'w') as f:
            f.write('[]')
        with open(self.catalog2, 'w') as f:
            f.write('chr1\t1000\t1100\tCAG\n')

    def test_parse_args_minimal(self):
        """Test parsing minimal required arguments."""
        sys.argv = ['test', self.catalog1, self.catalog2]
        args, paths = parse_args()

        self.assertEqual(len(paths), 2)
        self.assertEqual(args.overlap_fraction, 0.66)
        self.assertEqual(args.merge_type, 'union')
        self.assertEqual(args.overlapping_loci_action, 'keep-first')

    def test_parse_args_with_overlap_fraction(self):
        """Test parsing with custom overlap fraction."""
        sys.argv = ['test', '-f', '0.8', self.catalog1, self.catalog2]
        args, paths = parse_args()

        self.assertEqual(args.overlap_fraction, 0.8)

    def test_parse_args_with_jaccard_similarity(self):
        """Test parsing with Jaccard similarity threshold."""
        sys.argv = ['test', '--min-jaccard-similarity', '0.5', self.catalog1, self.catalog2]
        args, paths = parse_args()

        self.assertEqual(args.min_jaccard_similarity, 0.5)

    def test_parse_args_invalid_jaccard_similarity(self):
        """Test error with invalid Jaccard similarity."""
        sys.argv = ['test', '--min-jaccard-similarity', '1.5', self.catalog1, self.catalog2]
        with self.assertRaises(SystemExit):
            parse_args()

    def test_parse_args_merge_type_union(self):
        """Test parsing with merge type union."""
        sys.argv = ['test', '--merge-type', 'union', self.catalog1, self.catalog2]
        args, paths = parse_args()

        self.assertEqual(args.merge_type, 'union')

    def test_parse_args_merge_type_intersection(self):
        """Test parsing with merge type intersection."""
        sys.argv = ['test', '--merge-type', 'intersection', self.catalog1, self.catalog2]
        args, paths = parse_args()

        self.assertEqual(args.merge_type, 'intersection')

    def test_parse_args_overlapping_loci_actions(self):
        """Test parsing different overlapping loci actions."""
        for action in ['keep-first', 'keep-last', 'keep-both', 'keep-narrow', 'keep-wider', 'merge']:
            sys.argv = ['test', '--overlapping-loci-action', action, self.catalog1, self.catalog2]
            args, paths = parse_args()
            self.assertEqual(args.overlapping_loci_action, action)

    def test_parse_args_merge_adjacent_loci(self):
        """Test --merge-adjacent-loci-with-same-motif flag."""
        sys.argv = ['test', '-m', self.catalog1, self.catalog2]
        args, paths = parse_args()

        self.assertTrue(args.merge_adjacent_loci_with_same_motif)

    def test_parse_args_output_formats(self):
        """Test different output formats."""
        # JSON format
        sys.argv = ['test', '--output-format', 'JSON', self.catalog1, self.catalog2]
        args, paths = parse_args()
        self.assertEqual(args.output_format, 'JSON')

        # BED format
        sys.argv = ['test', '--output-format', 'BED', self.catalog1, self.catalog2]
        args, paths = parse_args()
        self.assertEqual(args.output_format, 'BED')

    def test_parse_args_add_found_in_fields_requires_labels(self):
        """Test that --add-found-in-fields requires catalog labels."""
        sys.argv = ['test', '--add-found-in-fields', self.catalog1, self.catalog2]
        with self.assertRaises(SystemExit):
            parse_args()

    def test_parse_args_conflicting_options(self):
        """Test conflicting option combinations."""
        # merge-adjacent with write-outer-join-table
        sys.argv = ['test', '-m', '--write-outer-join-table', self.catalog1, self.catalog2]
        with self.assertRaises(SystemExit):
            parse_args()

        # merge action with add-found-in-fields
        catalog1_labeled = f"cat1:{self.catalog1}"
        catalog2_labeled = f"cat2:{self.catalog2}"
        sys.argv = ['test', '--overlapping-loci-action', 'merge', '--add-found-in-fields',
                   catalog1_labeled, catalog2_labeled]
        with self.assertRaises(SystemExit):
            parse_args()

    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        shutil.rmtree(self.temp_dir)


class TestParseVariantCatalogPathsArg(unittest.TestCase):
    """Test parse_variant_catalog_paths_arg function."""

    def setUp(self):
        """Create temporary catalog files."""
        self.temp_dir = tempfile.mkdtemp()
        self.json_catalog = os.path.join(self.temp_dir, "catalog.json")
        self.bed_catalog = os.path.join(self.temp_dir, "catalog.bed")
        self.bed_gz_catalog = os.path.join(self.temp_dir, "catalog.bed.gz")

        with open(self.json_catalog, 'w') as f:
            f.write('[]')
        with open(self.bed_catalog, 'w') as f:
            f.write('chr1\t1000\t1100\tCAG\n')
        with open(self.bed_gz_catalog, 'w') as f:
            f.write('chr1\t1000\t1100\tCAG\n')

    def test_parse_paths_basic(self):
        """Test basic path parsing."""
        args = argparse.Namespace(
            variant_catalog_json_or_bed=[self.json_catalog, self.bed_catalog],
            add_found_in_fields=False
        )
        parser = argparse.ArgumentParser()

        paths = parse_variant_catalog_paths_arg(args, parser)

        self.assertEqual(len(paths), 2)
        self.assertEqual(paths[0][1], self.json_catalog)
        self.assertEqual(paths[0][2], 'JSON')
        self.assertEqual(paths[1][1], self.bed_catalog)
        self.assertEqual(paths[1][2], 'BED')

    def test_parse_paths_with_labels(self):
        """Test path parsing with catalog labels."""
        args = argparse.Namespace(
            variant_catalog_json_or_bed=[
                f"catalog1:{self.json_catalog}",
                f"catalog2:{self.bed_catalog}"
            ],
            add_found_in_fields=True
        )
        parser = argparse.ArgumentParser()

        paths = parse_variant_catalog_paths_arg(args, parser)

        self.assertEqual(len(paths), 2)
        self.assertEqual(paths[0][0], 'catalog1')
        self.assertEqual(paths[1][0], 'catalog2')

    def test_parse_paths_duplicate_catalog_name(self):
        """Test error with duplicate catalog names."""
        args = argparse.Namespace(
            variant_catalog_json_or_bed=[
                f"dup:{self.json_catalog}",
                f"dup:{self.bed_catalog}"
            ],
            add_found_in_fields=False
        )
        parser = argparse.ArgumentParser()

        with self.assertRaises(SystemExit):
            parse_variant_catalog_paths_arg(args, parser)

    def test_parse_paths_file_not_found(self):
        """Test error when file doesn't exist."""
        args = argparse.Namespace(
            variant_catalog_json_or_bed=['/nonexistent/file.json'],
            add_found_in_fields=False
        )
        parser = argparse.ArgumentParser()

        with self.assertRaises(SystemExit):
            parse_variant_catalog_paths_arg(args, parser)

    def test_parse_paths_separator_in_name(self):
        """Test error when catalog name contains separator."""
        args = argparse.Namespace(
            variant_catalog_json_or_bed=[
                f"cat |||  1:{self.json_catalog}"
            ],
            add_found_in_fields=False
        )
        parser = argparse.ArgumentParser()

        with self.assertRaises(SystemExit):
            parse_variant_catalog_paths_arg(args, parser)

    def test_parse_paths_unrecognized_extension(self):
        """Test error with unrecognized file extension."""
        bad_file = os.path.join(self.temp_dir, "catalog.txt")
        with open(bad_file, 'w') as f:
            f.write('test')

        args = argparse.Namespace(
            variant_catalog_json_or_bed=[bad_file],
            add_found_in_fields=False
        )
        parser = argparse.ArgumentParser()

        with self.assertRaises(SystemExit):
            parse_variant_catalog_paths_arg(args, parser)

    def test_parse_paths_various_extensions(self):
        """Test parsing paths with various valid extensions."""
        args = argparse.Namespace(
            variant_catalog_json_or_bed=[
                self.json_catalog,
                self.bed_catalog,
                self.bed_gz_catalog
            ],
            add_found_in_fields=False
        )
        parser = argparse.ArgumentParser()

        paths = parse_variant_catalog_paths_arg(args, parser)

        self.assertEqual(len(paths), 3)
        self.assertEqual(paths[0][2], 'JSON')
        self.assertEqual(paths[1][2], 'BED')
        self.assertEqual(paths[2][2], 'BED')

    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        shutil.rmtree(self.temp_dir)


class TestCheckForSufficientOverlapAndMotifMatch(unittest.TestCase):
    """Test check_for_sufficient_overlap_and_motif_match function."""

    def test_no_overlap(self):
        """Test when intervals don't overlap."""
        existing = Interval(1000, 1100, data={
            "LocusStructure": "(CAG)*",
            "ReferenceRegion": "chr1:1000-1100"
        })
        new = Interval(2000, 2100, data={
            "LocusStructure": "(CAG)*",
            "ReferenceRegion": "chr1:2000-2100"
        })

        result = check_for_sufficient_overlap_and_motif_match(existing, new)
        self.assertIsNone(result)

    def test_sufficient_overlap_same_motif(self):
        """Test sufficient overlap with same canonical motif."""
        existing = Interval(1000, 1100, data={
            "LocusStructure": "(CAG)*",
            "ReferenceRegion": "chr1:1000-1100"
        })
        new = Interval(1050, 1150, data={
            "LocusStructure": "(CAG)*",
            "ReferenceRegion": "chr1:1050-1150"
        })

        result = check_for_sufficient_overlap_and_motif_match(
            existing, new, min_overlap_fraction=0.5
        )
        self.assertEqual(result, existing)

    def test_insufficient_overlap(self):
        """Test insufficient overlap."""
        existing = Interval(1000, 1100, data={
            "LocusStructure": "(CAG)*",
            "ReferenceRegion": "chr1:1000-1100"
        })
        # Overlap of only 3bp = 1 motif length, less than required 2 motifs
        new = Interval(1097, 1150, data={
            "LocusStructure": "(CAG)*",
            "ReferenceRegion": "chr1:1097-1150"
        })

        result = check_for_sufficient_overlap_and_motif_match(
            existing, new, min_overlap_fraction=0.8
        )
        self.assertIsNone(result)

    def test_different_canonical_motifs(self):
        """Test overlap with different canonical motifs."""
        existing = Interval(1000, 1100, data={
            "LocusStructure": "(CAG)*",
            "ReferenceRegion": "chr1:1000-1100"
        })
        new = Interval(1050, 1150, data={
            "LocusStructure": "(AT)*",
            "ReferenceRegion": "chr1:1050-1150"
        })

        result = check_for_sufficient_overlap_and_motif_match(
            existing, new, min_overlap_fraction=0.5
        )
        self.assertIsNone(result)

    def test_overlap_by_two_motif_lengths(self):
        """Test overlap by at least 2 motif lengths."""
        existing = Interval(1000, 1100, data={
            "LocusStructure": "(CAG)*",
            "ReferenceRegion": "chr1:1000-1100"
        })
        # Overlap of 6bp = 2 * 3bp motif
        new = Interval(1094, 1150, data={
            "LocusStructure": "(CAG)*",
            "ReferenceRegion": "chr1:1094-1150"
        })

        result = check_for_sufficient_overlap_and_motif_match(
            existing, new, min_overlap_fraction=0.9  # High threshold
        )
        self.assertEqual(result, existing)

    def test_jaccard_similarity_filtering(self):
        """Test Jaccard similarity threshold."""
        existing = Interval(1000, 1100, data={
            "LocusStructure": "(CAG)*",
            "ReferenceRegion": "chr1:1000-1100"
        })
        new = Interval(1050, 1150, data={
            "LocusStructure": "(CAG)*",
            "ReferenceRegion": "chr1:1050-1150"
        })

        # Jaccard = 50 / 150 = 0.33, should fail with min=0.5
        result = check_for_sufficient_overlap_and_motif_match(
            existing, new, min_jaccard_similarity=0.5
        )
        self.assertIsNone(result)

        # Should pass with min=0.3
        result = check_for_sufficient_overlap_and_motif_match(
            existing, new, min_jaccard_similarity=0.3
        )
        self.assertEqual(result, existing)

    def test_vntr_motif_length_match(self):
        """Test VNTR motif length matching."""
        existing = Interval(1000, 1100, data={
            "LocusStructure": "(AAAAAAG)*",  # 7bp motif
            "ReferenceRegion": "chr1:1000-1100"
        })
        new = Interval(1050, 1150, data={
            "LocusStructure": "(CCCCCCT)*",  # 7bp motif, different sequence
            "ReferenceRegion": "chr1:1050-1150"
        })

        # Should fail without VNTR flag
        result = check_for_sufficient_overlap_and_motif_match(
            existing, new, min_overlap_fraction=0.5
        )
        self.assertIsNone(result)

        # Should pass with VNTR flag
        result = check_for_sufficient_overlap_and_motif_match(
            existing, new, min_overlap_fraction=0.5,
            motif_length_match_sufficient_for_VNTRs=True
        )
        self.assertEqual(result, existing)

    def test_adjacent_loci_same_structure(self):
        """Test loci with same LocusStructure (multiple reference regions)."""
        existing = Interval(1000, 1100, data={
            "LocusStructure": "(CAG)*ACGT(CTG)*",
            "ReferenceRegion": ["chr1:1000-1050", "chr1:1055-1100"]
        })
        new = Interval(1020, 1120, data={
            "LocusStructure": "(CAG)*ACGT(CTG)*",
            "ReferenceRegion": ["chr1:1020-1070", "chr1:1075-1120"]
        })

        result = check_for_sufficient_overlap_and_motif_match(
            existing, new, min_overlap_fraction=0.5
        )
        self.assertEqual(result, existing)

    def test_counter_tracking(self):
        """Test that counters are updated correctly."""
        counters = collections.Counter()
        existing = Interval(1000, 1100, data={
            "LocusStructure": "(CAG)*",
            "ReferenceRegion": "chr1:1000-1100"
        })
        new = Interval(1050, 1150, data={
            "LocusStructure": "(CAG)*",
            "ReferenceRegion": "chr1:1050-1150"
        })

        check_for_sufficient_overlap_and_motif_match(
            existing, new, counters=counters, min_overlap_fraction=0.5
        )

        self.assertGreater(len(counters), 0)
        self.assertIn('same canonical motif', str(counters))


class TestRequireExactMatchOfLocusDefinitions(unittest.TestCase):
    """Test the --require-exact-match-of-locus-definitions option."""

    @staticmethod
    def _interval(start, end, locus_structure, reference_region=None):
        return Interval(start, end, data={
            "LocusStructure": locus_structure,
            "ReferenceRegion": reference_region or f"chr1:{start}-{end}",
        })

    def _check(self, existing, new):
        return check_for_sufficient_overlap_and_motif_match(
            existing, new, require_exact_match_of_locus_definitions=True)

    def test_same_start_end_and_motif_match(self):
        existing = self._interval(1000, 1100, "(CAG)*")
        self.assertEqual(self._check(existing, self._interval(1000, 1100, "(CAG)*")), existing)

    def test_same_canonical_motif_after_cyclic_rotation_matches(self):
        existing = self._interval(1000, 1100, "(CAG)*")
        self.assertEqual(self._check(existing, self._interval(1000, 1100, "(AGC)*")), existing)

    def test_different_end_does_not_match(self):
        self.assertIsNone(self._check(self._interval(1000, 1100, "(CAG)*"), self._interval(1000, 1103, "(CAG)*")))

    def test_different_start_does_not_match(self):
        self.assertIsNone(self._check(self._interval(1000, 1100, "(CAG)*"), self._interval(1003, 1100, "(CAG)*")))

    def test_different_canonical_motif_does_not_match(self):
        self.assertIsNone(self._check(self._interval(1000, 1100, "(CAG)*"), self._interval(1000, 1100, "(CCG)*")))

    def test_overlap_that_passes_the_default_criteria_does_not_match(self):
        existing = self._interval(1000, 1100, "(CAG)*")
        new = self._interval(1000, 1099, "(CAG)*")
        self.assertEqual(check_for_sufficient_overlap_and_motif_match(existing, new), existing)
        self.assertIsNone(self._check(existing, new))

    def test_compound_locus_does_not_match_non_compound_locus_with_a_matching_component(self):
        compound = self._interval(1000, 1200, "(CAG)*(CCG)*", ["chr1:1000-1100", "chr1:1100-1200"])
        non_compound = self._interval(1000, 1100, "(CAG)*")
        self.assertIsNone(self._check(compound, non_compound))
        self.assertIsNone(self._check(non_compound, compound))

    def test_identical_compound_loci_match(self):
        regions = ["chr1:1000-1100", "chr1:1100-1200"]
        existing = self._interval(1000, 1200, "(CAG)*(CCG)*", regions)
        self.assertEqual(self._check(existing, self._interval(1000, 1200, "(CAG)*(CCG)*", regions)), existing)

    def test_compound_loci_with_different_component_boundaries_do_not_match(self):
        existing = self._interval(1000, 1200, "(CAG)*(CCG)*", ["chr1:1000-1100", "chr1:1100-1200"])
        new = self._interval(1000, 1200, "(CAG)*(CCG)*", ["chr1:1000-1103", "chr1:1103-1200"])
        self.assertIsNone(self._check(existing, new))

    def test_locus_definition_is_distinct_for_compound_and_non_compound_loci(self):
        self.assertNotEqual(
            compute_locus_definition_for_exact_match(
                {"ReferenceRegion": ["chr1:1000-1100", "chr1:1100-1200"], "LocusStructure": "(CAG)*(CCG)*"}),
            compute_locus_definition_for_exact_match({"ReferenceRegion": "chr1:1000-1100", "LocusStructure": "(CAG)*"}))

    def test_option_can_not_be_combined_with_min_jaccard_similarity(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            catalog = os.path.join(temp_dir, "catalog.json")
            with open(catalog, "w") as f:
                f.write("[]")
            sys.argv = ["test", "--require-exact-match-of-locus-definitions", "--min-jaccard-similarity", "0.5",
                        catalog, catalog]
            with self.assertRaises(SystemExit):
                with mock.patch("sys.stderr"):
                    parse_args()


class TestMergeActionAndAdjacentLociRegressions(unittest.TestCase):
    """Regression tests for the overlap merge action, adjacent-locus merging and catalog path parsing."""

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    def _write_catalog(self, name, regions_and_motifs):
        path = os.path.join(self.temp_dir.name, name)
        with open(path, "w") as f:
            f.write("".join(f"chr1\t{start}\t{end}\t{motif}\n" for start, end, motif in regions_and_motifs))
        return path

    def _add_catalogs(self, catalog_paths, **kwargs):
        interval_trees = collections.defaultdict(IntervalTree)
        outer_join_overlap_table = collections.defaultdict(dict)
        for path in catalog_paths:
            add_variant_catalog_to_interval_trees(
                os.path.basename(path), path, "BED", interval_trees, outer_join_overlap_table, **kwargs)
        return interval_trees, outer_join_overlap_table

    def test_merge_action_stores_the_merged_span_and_can_merge_a_third_overlapping_locus(self):
        catalogs = [
            self._write_catalog("a.bed", [(1000, 1100, "CAG")]),
            self._write_catalog("b.bed", [(1050, 1150, "CAG")]),
            self._write_catalog("c.bed", [(1120, 1200, "CAG")]),
        ]
        interval_trees, _ = self._add_catalogs(
            catalogs, overlapping_loci_action="merge", min_overlap_fraction=0.01,
            only_compare_loci_from_different_catalogs=True)
        output_records = list(convert_interval_trees_to_output_records(interval_trees))
        self.assertEqual([r["ReferenceRegion"] for r in output_records], ["chr1:1000-1200"])
        self.assertEqual({(i.begin, i.end) for i in interval_trees["1"] if not i.data.get("Discarded")}, {(1000, 1200)})

    def test_merge_action_with_intersection_keeps_the_merged_locus_present_in_every_catalog(self):
        catalogs = [self._write_catalog("a.bed", [(1000, 1100, "CAG")]),
                    self._write_catalog("b.bed", [(1050, 1150, "CAG")])]
        interval_trees, outer_join_overlap_table = self._add_catalogs(
            catalogs, overlapping_loci_action="merge", min_overlap_fraction=0.01)
        output_records = list(convert_interval_trees_to_output_records(
            interval_trees, outer_join_overlap_table, only_loci_present_in_n_catalogs=2))
        self.assertEqual([r["ReferenceRegion"] for r in output_records], ["chr1:1000-1150"])

    def test_merge_adjacent_loci_does_not_drop_or_resurrect_loci_around_a_discarded_locus(self):
        # the middle locus is a duplicate of the first, so keep-first discards it
        catalogs = [self._write_catalog("a.bed", [(1000, 1100, "CAG"), (2000, 2100, "CAG")]),
                    self._write_catalog("b.bed", [(1000, 1100, "CAG")])]
        interval_trees, _ = self._add_catalogs(catalogs, overlapping_loci_action="keep-first")
        output_records = list(convert_interval_trees_to_output_records(
            interval_trees, merge_adjacent_loci_with_same_motif=True))
        self.assertEqual([r["ReferenceRegion"] for r in output_records], ["chr1:1000-1100", "chr1:2000-2100"])

    def test_catalog_presence_is_recorded_when_the_only_overlapping_locus_was_discarded(self):
        catalogs = [self._write_catalog("a.bed", [(1000, 1100, "CAG")]),
                    self._write_catalog("b.bed", [(1000, 1100, "CAG")]),
                    self._write_catalog("c.bed", [(1000, 1100, "CAG")])]
        _, outer_join_overlap_table = self._add_catalogs(catalogs, overlapping_loci_action="keep-first")
        self.assertEqual(outer_join_overlap_table["1-1000-1100-AGC"],
                         {"a.bed": "Yes", "b.bed": "Yes", "c.bed": "Yes"})

    def test_merged_locus_source_is_split_into_catalog_names_in_the_presence_table(self):
        catalogs = [self._write_catalog("a.bed", [(1000, 1100, "CAG")]),
                    self._write_catalog("b.bed", [(1050, 1150, "CAG")]),
                    self._write_catalog("c.bed", [(1120, 1200, "CAG")])]
        interval_trees, outer_join_overlap_table = self._add_catalogs(
            catalogs, overlapping_loci_action="merge", min_overlap_fraction=0.01)
        all_catalog_names = {name for presence in outer_join_overlap_table.values() for name in presence}
        self.assertEqual(all_catalog_names, {"a.bed", "b.bed", "c.bed"})
        # the locus is in 3 catalogs, so requiring 4 catalogs must not output it
        self.assertEqual(list(convert_interval_trees_to_output_records(
            interval_trees, outer_join_overlap_table, only_loci_present_in_n_catalogs=4)), [])

    def test_zero_width_loci_are_skipped_instead_of_aborting_the_merge(self):
        catalogs = [self._write_catalog("a.bed", [(1000, 1000, "CAG"), (2000, 2100, "CAG")])]
        interval_trees, _ = self._add_catalogs(catalogs)
        self.assertEqual([(i.begin, i.end) for i in interval_trees["1"]], [(2000, 2100)])

    def test_http_catalog_paths_are_not_checked_with_file_exists(self):
        args = argparse.Namespace(variant_catalog_json_or_bed=["https://example.org/a.json"], add_found_in_fields=False)
        with mock.patch("str_analysis.merge_loci.file_exists", return_value=False):
            self.assertEqual(parse_variant_catalog_paths_arg(args, argparse.ArgumentParser()),
                             [("a.json", "https://example.org/a.json", "JSON")])

    def test_comma_separated_compound_bed_rows_are_skipped(self):
        catalogs = [self._write_catalog("a.bed", [(1000, 1100, "CAG,CCG"), (2000, 2100, "CAG")])]
        interval_trees, _ = self._add_catalogs(catalogs)
        self.assertEqual([(i.begin, i.end) for i in interval_trees["1"]], [(2000, 2100)])

    def test_uri_catalog_paths_are_not_split_at_the_scheme_colon(self):
        for path, expected_name, expected_path in [
                ("gs://bucket/dir/a.json", "a.json", "gs://bucket/dir/a.json"),
                ("catalog1:gs://bucket/dir/a.json", "catalog1", "gs://bucket/dir/a.json"),
                ("catalog1:/local/a.json", "catalog1", "/local/a.json")]:
            args = argparse.Namespace(variant_catalog_json_or_bed=[path], add_found_in_fields=False)
            with mock.patch("str_analysis.merge_loci.file_exists", return_value=True):
                self.assertEqual(parse_variant_catalog_paths_arg(args, argparse.ArgumentParser()),
                                 [(expected_name, expected_path, "JSON")])


class TestCheckWhetherToMergeAdjacentLoci(unittest.TestCase):
    """Test check_whether_to_merge_adjacent_loci function."""

    def test_adjacent_loci_same_motif(self):
        """Test merging adjacent loci with same motif."""
        previous = Interval(1000, 1100, data={
            "LocusStructure": "(CAG)*",
            "ReferenceRegion": "chr1:1000-1100",
            "Source": "test"
        })
        current = Interval(1101, 1200, data={  # 1bp gap
            "LocusStructure": "(CAG)*",
            "ReferenceRegion": "chr1:1101-1200",
            "Source": "test"
        })

        should_merge, merged_interval = check_whether_to_merge_adjacent_loci(previous, current)
        self.assertTrue(should_merge)
        self.assertIsNotNone(merged_interval)

    def test_non_adjacent_loci(self):
        """Test non-adjacent loci are not merged."""
        previous = Interval(1000, 1100, data={
            "LocusStructure": "(CAG)*",
            "ReferenceRegion": "chr1:1000-1100",
            "Source": "test"
        })
        current = Interval(1105, 1200, data={  # 5bp gap
            "LocusStructure": "(CAG)*",
            "ReferenceRegion": "chr1:1105-1200",
            "Source": "test"
        })

        should_merge, merged_interval = check_whether_to_merge_adjacent_loci(previous, current)
        self.assertFalse(should_merge)
        self.assertIsNone(merged_interval)

    def test_adjacent_loci_different_motif(self):
        """Test adjacent loci with different motifs are not merged."""
        previous = Interval(1000, 1100, data={
            "LocusStructure": "(CAG)*",
            "ReferenceRegion": "chr1:1000-1100",
            "Source": "test"
        })
        current = Interval(1101, 1200, data={
            "LocusStructure": "(AT)*",  # Different motif
            "ReferenceRegion": "chr1:1101-1200",
            "Source": "test"
        })

        should_merge, merged_interval = check_whether_to_merge_adjacent_loci(previous, current)
        self.assertFalse(should_merge)
        self.assertIsNone(merged_interval)

    def test_overlapping_loci(self):
        """Test overlapping loci are not merged."""
        previous = Interval(1000, 1100, data={
            "LocusStructure": "(CAG)*",
            "ReferenceRegion": "chr1:1000-1100",
            "Source": "test"
        })
        current = Interval(1050, 1150, data={
            "LocusStructure": "(CAG)*",
            "ReferenceRegion": "chr1:1050-1150",
            "Source": "test"
        })

        should_merge, merged_interval = check_whether_to_merge_adjacent_loci(previous, current)
        self.assertFalse(should_merge)
        self.assertIsNone(merged_interval)


class TestMainIntegration(unittest.TestCase):
    """Integration tests for main() function."""

    def setUp(self):
        """Create temporary files for testing."""
        self.temp_dir = tempfile.mkdtemp()

    def create_json_catalog(self, filename, records):
        """Helper to create JSON catalog file."""
        import simplejson
        path = os.path.join(self.temp_dir, filename)
        with open(path, 'w') as f:
            simplejson.dump(records, f)
        return path

    def create_bed_catalog(self, filename, records):
        """Helper to create BED catalog file."""
        path = os.path.join(self.temp_dir, filename)
        with open(path, 'w') as f:
            for chrom, start, end, motif in records:
                f.write(f"{chrom}\t{start}\t{end}\t{motif}\n")
        return path

    def test_main_merge_two_catalogs_union(self):
        """Test merging two catalogs with union."""
        catalog1 = self.create_json_catalog("cat1.json", [
            {
                "LocusId": "locus1",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr1:1000-1100",
                "VariantType": "Repeat"
            }
        ])

        catalog2 = self.create_json_catalog("cat2.json", [
            {
                "LocusId": "locus2",
                "LocusStructure": "(CTG)*",
                "ReferenceRegion": "chr1:2000-2100",
                "VariantType": "Repeat"
            }
        ])

        output_prefix = os.path.join(self.temp_dir, "merged")

        sys.argv = [
            'test',
            '--output-prefix', output_prefix,
            '--output-format', 'JSON',
            catalog1, catalog2
        ]

        from str_analysis.merge_loci import main
        main()

        # Check output exists (written as .json.gz by default)
        output_file = f"{output_prefix}.json.gz" if os.path.exists(f"{output_prefix}.json.gz") else f"{output_prefix}.json"
        self.assertTrue(os.path.exists(output_file))

        # Check both loci are in output
        import simplejson
        import gzip
        fopen = gzip.open if output_file.endswith('.gz') else open
        with fopen(output_file, 'rt') as f:
            output = simplejson.load(f)

        self.assertEqual(len(output), 2)

    def test_main_merge_overlapping_loci_keep_first(self):
        """Test merging with overlapping loci, keep first."""
        catalog1 = self.create_json_catalog("cat1.json", [
            {
                "LocusId": "locus1",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr1:1000-1100",
                "VariantType": "Repeat"
            }
        ])

        catalog2 = self.create_json_catalog("cat2.json", [
            {
                "LocusId": "locus2",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr1:1050-1150",  # Overlaps with locus1
                "VariantType": "Repeat"
            }
        ])

        output_prefix = os.path.join(self.temp_dir, "merged")

        sys.argv = [
            'test',
            '--output-prefix', output_prefix,
            '--output-format', 'JSON',
            '--overlapping-loci-action', 'keep-first',
            catalog1, catalog2
        ]

        from str_analysis.merge_loci import main
        main()

        import simplejson
        import gzip
        output_file = f"{output_prefix}.json.gz" if os.path.exists(f"{output_prefix}.json.gz") else f"{output_prefix}.json"
        fopen = gzip.open if output_file.endswith('.gz') else open
        with fopen(output_file, 'rt') as f:
            output = simplejson.load(f)

        # Should only have locus1 (first)
        self.assertEqual(len(output), 1)
        self.assertEqual(output[0]['LocusId'], 'locus1')

    def test_main_merge_overlapping_loci_keep_both(self):
        """Test merging with overlapping loci, keep both."""
        catalog1 = self.create_json_catalog("cat1.json", [
            {
                "LocusId": "locus1",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr1:1000-1100",
                "VariantType": "Repeat"
            }
        ])

        catalog2 = self.create_json_catalog("cat2.json", [
            {
                "LocusId": "locus2",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr1:1050-1150",
                "VariantType": "Repeat"
            }
        ])

        output_prefix = os.path.join(self.temp_dir, "merged")

        sys.argv = [
            'test',
            '--output-prefix', output_prefix,
            '--output-format', 'JSON',
            '--overlapping-loci-action', 'keep-both',
            catalog1, catalog2
        ]

        from str_analysis.merge_loci import main
        main()

        import simplejson
        import gzip
        output_file = f"{output_prefix}.json.gz" if os.path.exists(f"{output_prefix}.json.gz") else f"{output_prefix}.json"
        fopen = gzip.open if output_file.endswith('.gz') else open
        with fopen(output_file, 'rt') as f:
            output = simplejson.load(f)

        # Should have both loci
        self.assertEqual(len(output), 2)

    def test_main_merge_bed_format(self):
        """Test merging BED format catalogs."""
        catalog1 = self.create_bed_catalog("cat1.bed", [
            ("chr1", 1000, 1100, "CAG"),
            ("chr2", 2000, 2100, "CTG"),
        ])

        catalog2 = self.create_bed_catalog("cat2.bed", [
            ("chr3", 3000, 3100, "AT"),
        ])

        output_prefix = os.path.join(self.temp_dir, "merged")

        sys.argv = [
            'test',
            '--output-prefix', output_prefix,
            '--output-format', 'BED',
            catalog1, catalog2
        ]

        from str_analysis.merge_loci import main
        with mock.patch('os.system'):  # Mock bgzip/tabix calls
            main()

        # Check BED output exists
        self.assertTrue(os.path.exists(f"{output_prefix}.bed"))

    def test_main_merge_with_jaccard_filter(self):
        """Test merging with Jaccard similarity filter."""
        catalog1 = self.create_json_catalog("cat1.json", [
            {
                "LocusId": "locus1",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr1:1000-1200",  # 200bp
                "VariantType": "Repeat"
            }
        ])

        catalog2 = self.create_json_catalog("cat2.json", [
            {
                "LocusId": "locus2",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr1:1100-1300",  # 200bp, overlap=100bp
                "VariantType": "Repeat"
            }
        ])

        output_prefix = os.path.join(self.temp_dir, "merged")

        # Jaccard = 100/(300) = 0.33, should keep both with min=0.5
        sys.argv = [
            'test',
            '--output-prefix', output_prefix,
            '--output-format', 'JSON',
            '--min-jaccard-similarity', '0.5',
            '--overlapping-loci-action', 'keep-first',
            catalog1, catalog2
        ]

        from str_analysis.merge_loci import main
        main()

        import simplejson
        import gzip
        output_file = f"{output_prefix}.json.gz" if os.path.exists(f"{output_prefix}.json.gz") else f"{output_prefix}.json"
        fopen = gzip.open if output_file.endswith('.gz') else open
        with fopen(output_file, 'rt') as f:
            output = simplejson.load(f)

        # Should have both since Jaccard < 0.5
        self.assertEqual(len(output), 2)

    def test_main_merge_type_intersection(self):
        """Test merge type intersection."""
        catalog1 = self.create_json_catalog("cat1.json", [
            {
                "LocusId": "locus1",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr1:1000-1100",
                "VariantType": "Repeat"
            },
            {
                "LocusId": "locus2",
                "LocusStructure": "(CTG)*",
                "ReferenceRegion": "chr2:2000-2100",
                "VariantType": "Repeat"
            }
        ])

        catalog2 = self.create_json_catalog("cat2.json", [
            {
                "LocusId": "locus3",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr1:1050-1150",  # Overlaps locus1
                "VariantType": "Repeat"
            }
        ])

        output_prefix = os.path.join(self.temp_dir, "merged")

        sys.argv = [
            'test',
            '--output-prefix', output_prefix,
            '--output-format', 'JSON',
            '--merge-type', 'intersection',
            catalog1, catalog2
        ]

        from str_analysis.merge_loci import main
        main()

        import simplejson
        import gzip
        output_file = f"{output_prefix}.json.gz" if os.path.exists(f"{output_prefix}.json.gz") else f"{output_prefix}.json"
        fopen = gzip.open if output_file.endswith('.gz') else open
        with fopen(output_file, 'rt') as f:
            output = simplejson.load(f)

        # Should only have locus1 (present in both catalogs via overlap)
        self.assertEqual(len(output), 1)

    def test_main_merge_adjacent_loci(self):
        """Test --merge-adjacent-loci-with-same-motif."""
        catalog1 = self.create_json_catalog("cat1.json", [
            {
                "LocusId": "locus1",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr1:1000-1100",
                "VariantType": "Repeat",
                "Motif": "CAG"
            },
            {
                "LocusId": "locus2",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr1:1101-1200",  # Adjacent, 1bp gap
                "VariantType": "Repeat",
                "Motif": "CAG"
            }
        ])

        output_prefix = os.path.join(self.temp_dir, "merged")

        sys.argv = [
            'test',
            '--output-prefix', output_prefix,
            '--output-format', 'JSON',
            '-m',  # merge adjacent loci
            catalog1
        ]

        from str_analysis.merge_loci import main
        main()

        import simplejson
        import gzip
        output_file = f"{output_prefix}.json.gz" if os.path.exists(f"{output_prefix}.json.gz") else f"{output_prefix}.json"
        fopen = gzip.open if output_file.endswith('.gz') else open
        with fopen(output_file, 'rt') as f:
            output = simplejson.load(f)

        # Should merge to single locus
        self.assertEqual(len(output), 1)
        # Merged locus should span both
        self.assertIn('1000', output[0]['ReferenceRegion'])
        self.assertIn('1200', output[0]['ReferenceRegion'])

    def run_keep_first_merge_with_jaccard_threshold(self, catalog_paths, extra_args=()):
        """Runs main() with keep-first and a 0.66 Jaccard similarity threshold, and returns the output LocusIds."""
        import gzip
        import simplejson
        from str_analysis.merge_loci import main

        output_prefix = os.path.join(self.temp_dir, "merged")
        sys.argv = [
            'test',
            '--output-prefix', output_prefix,
            '--output-format', 'JSON',
            '--overlapping-loci-action', 'keep-first',
            '--overlap-fraction', '0.66',
            '--min-jaccard-similarity', '0.66',
            *extra_args,
            *catalog_paths,
        ]
        main()

        output_file = f"{output_prefix}.json.gz" if os.path.exists(f"{output_prefix}.json.gz") else f"{output_prefix}.json"
        fopen = gzip.open if output_file.endswith('.gz') else open
        with fopen(output_file, 'rt') as f:
            return [record['LocusId'] for record in simplejson.load(f)]

    def test_main_discarded_loci_dont_make_later_loci_redundant(self):
        """A record that only matches an already-discarded record is kept.

        locus2 (chr1:1000-1024) is discarded as redundant with locus1 (chr1:1000-1018, Jaccard 18/24). locus3
        (chr1:1000-1034) matches the discarded locus2 (Jaccard 24/34), but not locus1 (Jaccard 18/34), so nothing in
        the output describes it, and it's kept. This is the CHM13 catalog's chr1:1249185-1249219 (A)* locus.
        """
        catalog_paths = [
            self.create_json_catalog(f"cat{i}.json", [{
                "LocusId": f"locus{i}", "LocusStructure": "(A)*", "ReferenceRegion": region, "VariantType": "Repeat"}])
            for i, region in [(1, "chr1:1000-1018"), (2, "chr1:1000-1024"), (3, "chr1:1000-1034")]
        ]
        self.assertEqual(self.run_keep_first_merge_with_jaccard_threshold(catalog_paths), ["locus1", "locus3"])

    def test_main_discarded_loci_dont_make_later_loci_of_the_same_catalog_redundant(self):
        """Like the test above, but with all three records in one catalog."""
        catalog_path = self.create_json_catalog("cat.json", [
            {"LocusId": f"locus{i}", "LocusStructure": "(A)*", "ReferenceRegion": region, "VariantType": "Repeat"}
            for i, region in [(1, "chr1:1000-1018"), (2, "chr1:1000-1024"), (3, "chr1:1000-1034")]
        ])
        self.assertEqual(self.run_keep_first_merge_with_jaccard_threshold([catalog_path]), ["locus1", "locus3"])

    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        shutil.rmtree(self.temp_dir)


if __name__ == "__main__":
    unittest.main()
