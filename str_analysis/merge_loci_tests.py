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

    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        shutil.rmtree(self.temp_dir)


if __name__ == "__main__":
    unittest.main()
