#!/usr/bin/env python3

"""Comprehensive tests for annotate_and_filter_str_catalog.py"""

import argparse
import collections
import gzip
import io
import os
import sys
import tempfile
import unittest
from unittest import mock
from intervaltree import IntervalTree, Interval

from str_analysis.annotate_and_filter_str_catalog import (
    parse_args,
    parse_known_disease_associated_loci,
    get_overlap,
    output_tsv,
    count_Ns_in_flanks,
    VALID_GENE_REGIONS,
    KNOWN_DISEASE_ASSOCIATED_LOCI_COLUMNS,
)


class TestParseArgs(unittest.TestCase):
    """Test parse_args function."""

    def test_parse_args_minimal(self):
        """Test parsing minimal required arguments."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            f.write('[]')
            catalog_file = f.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write('>chr1\nACGT')
            ref_fasta = f.name

        try:
            sys.argv = ['test', '-R', ref_fasta, catalog_file]
            args, parser = parse_args()

            self.assertEqual(args.reference_fasta, ref_fasta)
            self.assertEqual(args.catalog_json_or_bed, catalog_file)
            self.assertIsNone(args.output_path)
            self.assertFalse(args.verbose)
        finally:
            os.unlink(catalog_file)
            os.unlink(ref_fasta)

    def test_parse_args_with_output_flags(self):
        """Test parsing with output format flags."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            f.write('[]')
            catalog_file = f.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write('>chr1\nACGT')
            ref_fasta = f.name

        try:
            sys.argv = ['test', '-R', ref_fasta, '--output-tsv', '--output-bed',
                       '--output-stats', catalog_file]
            args, parser = parse_args()

            self.assertTrue(args.output_tsv)
            self.assertTrue(args.output_bed)
            self.assertTrue(args.output_stats)
        finally:
            os.unlink(catalog_file)
            os.unlink(ref_fasta)

    def test_parse_args_with_filters(self):
        """Test parsing with filter arguments."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            f.write('[]')
            catalog_file = f.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write('>chr1\nACGT')
            ref_fasta = f.name

        try:
            sys.argv = ['test', '-R', ref_fasta, '--min-motif-size', '3',
                       '--max-motif-size', '6', '--min-interval-size-bp', '100',
                       catalog_file]
            args, parser = parse_args()

            self.assertEqual(args.min_motif_size, 3)
            self.assertEqual(args.max_motif_size, 6)
            self.assertEqual(args.min_interval_size_bp, 100)
        finally:
            os.unlink(catalog_file)
            os.unlink(ref_fasta)

    def test_parse_args_with_motif_filters(self):
        """Test parsing with motif-specific filters."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            f.write('[]')
            catalog_file = f.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write('>chr1\nACGT')
            ref_fasta = f.name

        try:
            sys.argv = ['test', '-R', ref_fasta, '-m', 'CAG', '-m', 'CTG',
                       '-ms', '3', '-xms', '1', catalog_file]
            args, parser = parse_args()

            self.assertEqual(args.motif, ['CAG', 'CTG'])
            self.assertEqual(args.motif_size, [3])
            self.assertEqual(args.exclude_motif_size, [1])
        finally:
            os.unlink(catalog_file)
            os.unlink(ref_fasta)

    def test_parse_args_locus_id_file(self):
        """Test parsing locus IDs from file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            f.write('[]')
            catalog_file = f.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write('>chr1\nACGT')
            ref_fasta = f.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write('locus1\nlocus2\nlocus3\n')
            locus_file = f.name

        try:
            sys.argv = ['test', '-R', ref_fasta, '--locus-id-file', locus_file, catalog_file]
            args, parser = parse_args()

            self.assertIn('locus1', args.locus_id)
            self.assertIn('locus2', args.locus_id)
            self.assertIn('locus3', args.locus_id)
            self.assertEqual(len(args.locus_id), 3)
        finally:
            os.unlink(catalog_file)
            os.unlink(ref_fasta)
            os.unlink(locus_file)

    def test_parse_args_missing_catalog_file_error(self):
        """Test error when catalog file doesn't exist."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write('>chr1\nACGT')
            ref_fasta = f.name

        try:
            sys.argv = ['test', '-R', ref_fasta, '/nonexistent/catalog.json']
            with self.assertRaises(SystemExit):
                parse_args()
        finally:
            os.unlink(ref_fasta)

    def test_parse_args_gene_region_validation(self):
        """Test that gene region filters require gene-models-source."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            f.write('[]')
            catalog_file = f.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write('>chr1\nACGT')
            ref_fasta = f.name

        try:
            # Test --region-type requires --gene-models-source
            sys.argv = ['test', '-R', ref_fasta, '-t', 'CDS', catalog_file]
            with self.assertRaises(SystemExit):
                parse_args()
        finally:
            os.unlink(catalog_file)
            os.unlink(ref_fasta)

    def test_valid_gene_regions_constant(self):
        """Test VALID_GENE_REGIONS constant."""
        expected_regions = {"CDS", "UTR", "5UTR", "3UTR", "promoter", "exon", "intron", "intergenic"}
        self.assertEqual(VALID_GENE_REGIONS, expected_regions)


class TestParseKnownDiseaseAssociatedLoci(unittest.TestCase):
    """Test parse_known_disease_associated_loci function."""

    def setUp(self):
        """Create temporary disease loci JSON file."""
        self.temp_dir = tempfile.mkdtemp()
        self.disease_json = os.path.join(self.temp_dir, "disease_loci.json")

        # Create a valid disease loci JSON
        disease_data = [
            {
                "LocusId": "HTT",
                "LocusStructure": "(CAG)*",
                "RepeatUnit": "CAG",
                "MainReferenceRegion": "chr4:3074876-3074933",
                "Gene": "HTT",
                "GeneRegion": "exon",
                "GeneId": "ENSG00000197386",
                "Diseases": ["Huntington disease"]
            },
            {
                "LocusId": "DMPK",
                "LocusStructure": "(CTG)*",
                "RepeatUnit": "CTG",
                "MainReferenceRegion": "chr19:45770204-45770264",
                "Gene": "DMPK",
                "GeneRegion": "3UTR",
                "GeneId": "ENSG00000104936",
                "Diseases": ["Myotonic dystrophy type 1"]
            },
            {
                "LocusId": "NotDisease",
                "LocusStructure": "(AT)*",
                "RepeatUnit": "AT",
                "MainReferenceRegion": "chr1:1000-1020",
                "Gene": "TEST",
                "GeneRegion": "intron",
                "GeneId": "ENSG00000000000",
                "Diseases": []  # No diseases - should be filtered out
            }
        ]

        import simplejson
        with open(self.disease_json, 'w') as f:
            simplejson.dump(disease_data, f)

    def test_parse_known_disease_loci(self):
        """Test parsing known disease loci from JSON."""
        args = argparse.Namespace(known_disease_associated_loci=self.disease_json)
        parser = argparse.ArgumentParser()

        interval_tree, canonical_motifs = parse_known_disease_associated_loci(args, parser)

        # Check that HTT and DMPK are in the interval tree, but NotDisease is not
        self.assertIn('4', interval_tree)  # chr4 normalized to '4'
        self.assertIn('19', interval_tree)

        # Check canonical motifs
        # CAG and CTG both have canonical motif 'AGC' (lexicographically smallest)
        self.assertIn('AGC', canonical_motifs)

    def test_parse_known_disease_loci_filters_no_diseases(self):
        """Test that loci with no diseases are filtered out."""
        args = argparse.Namespace(known_disease_associated_loci=self.disease_json)
        parser = argparse.ArgumentParser()

        interval_tree, canonical_motifs = parse_known_disease_associated_loci(args, parser)

        # Check that chr1 (NotDisease locus) is not in the tree
        self.assertNotIn('1', interval_tree)

    def test_parse_known_disease_loci_invalid_json(self):
        """Test error handling for invalid JSON."""
        invalid_json = os.path.join(self.temp_dir, "invalid.json")
        with open(invalid_json, 'w') as f:
            f.write("not valid json")

        args = argparse.Namespace(known_disease_associated_loci=invalid_json)
        parser = argparse.ArgumentParser()

        with self.assertRaises(SystemExit):
            parse_known_disease_associated_loci(args, parser)

    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        shutil.rmtree(self.temp_dir)


class TestGetOverlap(unittest.TestCase):
    """Test get_overlap function."""

    def setUp(self):
        """Create test interval trees."""
        self.interval_tree = collections.defaultdict(IntervalTree)

        # Add some test intervals
        self.interval_tree['1'].add(Interval(1000, 1100, data=argparse.Namespace(
            canonical_motif='CAG',
            motif='CAG',
            locus_id='HTT'
        )))

        self.interval_tree['2'].add(Interval(2000, 2100, data=argparse.Namespace(
            canonical_motif='CTG',
            motif='CTG',
            locus_id='DMPK'
        )))

    def test_get_overlap_exact_match(self):
        """Test exact match with same motif."""
        result = get_overlap(
            self.interval_tree, '1', 1000, 1100,
            canonical_motif='CAG', require_exact_locus_boundaries=True
        )
        self.assertEqual(result, 'HTT')

    def test_get_overlap_no_match_different_motif(self):
        """Test no match when motif differs."""
        result = get_overlap(
            self.interval_tree, '1', 1000, 1100,
            canonical_motif='CTG', require_exact_locus_boundaries=True
        )
        self.assertIsNone(result)

    def test_get_overlap_no_match_different_boundaries(self):
        """Test no match when boundaries differ."""
        result = get_overlap(
            self.interval_tree, '1', 1000, 1105,
            canonical_motif='CAG', require_exact_locus_boundaries=True
        )
        self.assertIsNone(result)

    def test_get_overlap_without_motif_check(self):
        """Test overlap check without motif matching."""
        result = get_overlap(
            self.interval_tree, '1', 1050, 1150,
            canonical_motif=None, require_exact_locus_boundaries=False
        )
        self.assertEqual(result, 'HTT')

    def test_get_overlap_no_overlap(self):
        """Test no overlap when intervals don't intersect."""
        result = get_overlap(
            self.interval_tree, '1', 2000, 2100,
            canonical_motif='CAG', require_exact_locus_boundaries=False
        )
        self.assertIsNone(result)

    def test_get_overlap_chrom_normalization(self):
        """Test that 'chr' prefix is removed."""
        result = get_overlap(
            self.interval_tree, 'chr1', 1000, 1100,
            canonical_motif='CAG', require_exact_locus_boundaries=True
        )
        self.assertEqual(result, 'HTT')

    def test_get_overlap_invalid_interval(self):
        """Test handling of invalid interval (end < start)."""
        with mock.patch('builtins.print') as mock_print:
            result = get_overlap(
                self.interval_tree, '1', 1100, 1000,
                canonical_motif='CAG', require_exact_locus_boundaries=False
            )
            # Should print warning and adjust end
            mock_print.assert_called_once()
            self.assertIn('WARNING', mock_print.call_args[0][0])


class TestOutputTsv(unittest.TestCase):
    """Test output_tsv function."""

    def setUp(self):
        """Create temporary output directory."""
        self.temp_dir = tempfile.mkdtemp()

    def test_output_tsv_single_reference_region(self):
        """Test TSV output with single reference regions."""
        output_path = os.path.join(self.temp_dir, "test.tsv")

        records = [
            {
                "LocusId": "locus1",
                "ReferenceRegion": "chr1:1000-1100",
                "Motif": "CAG",
                "NumRepeatsInReference": 33,
            },
            {
                "LocusId": "locus2",
                "ReferenceRegion": "chr2:2000-2100",
                "Motif": "CTG",
                "NumRepeatsInReference": 33,
            }
        ]

        output_tsv(output_path, records)

        # Check file was created
        self.assertTrue(os.path.exists(output_path))

        # Read and verify contents
        import pandas as pd
        df = pd.read_csv(output_path, sep='\t')
        self.assertEqual(len(df), 2)
        self.assertIn('LocusId', df.columns)
        self.assertIn('ReferenceRegion', df.columns)

    def test_output_tsv_multiple_reference_regions(self):
        """Test TSV output with multiple reference regions (adjacent repeats)."""
        output_path = os.path.join(self.temp_dir, "test_adjacent.tsv")

        records = [
            {
                "LocusId": "locus1",
                "ReferenceRegion": ["chr1:1000-1100", "chr1:1200-1300"],
                "NumRepeatsInReference": [33, 25],
                "VariantType": ["Repeat", "Repeat"],
                "ReferenceRepeatPurity": [0.95, 0.90],
            }
        ]

        output_tsv(output_path, records)

        # Read and verify contents
        import pandas as pd
        df = pd.read_csv(output_path, sep='\t')

        # Should expand to 2 rows
        self.assertEqual(len(df), 2)
        self.assertEqual(df.iloc[0]['ReferenceRegion'], 'chr1:1000-1100')
        self.assertEqual(df.iloc[1]['ReferenceRegion'], 'chr1:1200-1300')
        self.assertEqual(df.iloc[0]['NumRepeatsInReference'], 33)
        self.assertEqual(df.iloc[1]['NumRepeatsInReference'], 25)

    def test_output_tsv_gzipped(self):
        """Test TSV output to gzipped file."""
        output_path = os.path.join(self.temp_dir, "test.tsv.gz")

        records = [
            {
                "LocusId": "locus1",
                "ReferenceRegion": "chr1:1000-1100",
                "Motif": "CAG",
            }
        ]

        output_tsv(output_path, records)

        # Check gzipped file was created
        self.assertTrue(os.path.exists(output_path))

        # Read gzipped file
        import pandas as pd
        df = pd.read_csv(output_path, sep='\t')
        self.assertEqual(len(df), 1)

    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        shutil.rmtree(self.temp_dir)


class TestCountNsInFlanks(unittest.TestCase):
    """Test count_Ns_in_flanks function."""

    def setUp(self):
        """Create mock FASTA object."""
        self.mock_fasta = mock.Mock()

        # Mock chromosome size
        self.mock_fasta.lengths = [10000]
        self.mock_fasta.references = ['chr1']

    def test_count_Ns_in_flanks_no_Ns(self):
        """Test counting when there are no Ns."""
        # Mock get_reference_sequence_with_cache to return sequences without Ns
        with mock.patch('str_analysis.annotate_and_filter_str_catalog.get_reference_sequence_with_cache') as mock_get_seq, \
             mock.patch('str_analysis.annotate_and_filter_str_catalog.get_chromosome_size_with_cache', return_value=100000):

            mock_get_seq.side_effect = lambda fasta, chrom, start, end: 'A' * (end - start)

            count = count_Ns_in_flanks(self.mock_fasta, 'chr1', 5000, 6000, num_flanking_bases=1000)

            self.assertEqual(count, 0)

    def test_count_Ns_in_flanks_with_Ns(self):
        """Test counting when there are Ns in flanks."""
        with mock.patch('str_analysis.annotate_and_filter_str_catalog.get_reference_sequence_with_cache') as mock_get_seq, \
             mock.patch('str_analysis.annotate_and_filter_str_catalog.get_chromosome_size_with_cache', return_value=100000):

            def mock_sequence(fasta, chrom, start, end):
                # Left flank with 5 Ns
                if end == 5000:
                    return 'A' * 995 + 'NNNNN'
                # Right flank with 3 Ns
                else:
                    return 'NNN' + 'A' * 997

            mock_get_seq.side_effect = mock_sequence

            count = count_Ns_in_flanks(self.mock_fasta, 'chr1', 5000, 6000, num_flanking_bases=1000)

            self.assertEqual(count, 8)  # 5 + 3

    def test_count_Ns_in_flanks_at_chromosome_start(self):
        """Test counting Ns when repeat is near chromosome start."""
        with mock.patch('str_analysis.annotate_and_filter_str_catalog.get_reference_sequence_with_cache') as mock_get_seq, \
             mock.patch('str_analysis.annotate_and_filter_str_catalog.get_chromosome_size_with_cache', return_value=100000):

            mock_get_seq.return_value = 'A' * 1000

            # Repeat at position 500 - left flank should be truncated
            count = count_Ns_in_flanks(self.mock_fasta, 'chr1', 500, 600, num_flanking_bases=1000)

            # Check that left flank starts at 0, not negative
            calls = mock_get_seq.call_args_list
            self.assertEqual(calls[0][0][2], 0)  # start should be max(500-1000, 0) = 0

    def test_count_Ns_in_flanks_at_chromosome_end(self):
        """Test counting Ns when repeat is near chromosome end."""
        with mock.patch('str_analysis.annotate_and_filter_str_catalog.get_reference_sequence_with_cache') as mock_get_seq, \
             mock.patch('str_analysis.annotate_and_filter_str_catalog.get_chromosome_size_with_cache', return_value=10000):

            mock_get_seq.return_value = 'A' * 1000

            # Repeat at position 9500 - right flank should be truncated
            count = count_Ns_in_flanks(self.mock_fasta, 'chr1', 9400, 9500, num_flanking_bases=1000)

            # Check that right flank ends at chromosome size
            calls = mock_get_seq.call_args_list
            self.assertEqual(calls[1][0][3], 10000)  # end should be min(9500+1000, 10000) = 10000

    def test_count_Ns_in_flanks_custom_flank_size(self):
        """Test with custom flanking region size."""
        with mock.patch('str_analysis.annotate_and_filter_str_catalog.get_reference_sequence_with_cache') as mock_get_seq, \
             mock.patch('str_analysis.annotate_and_filter_str_catalog.get_chromosome_size_with_cache', return_value=100000):

            mock_get_seq.return_value = 'A' * 500

            count = count_Ns_in_flanks(self.mock_fasta, 'chr1', 5000, 6000, num_flanking_bases=500)

            # Check that flanks are 500bp
            calls = mock_get_seq.call_args_list
            left_flank_size = calls[0][0][3] - calls[0][0][2]
            right_flank_size = calls[1][0][3] - calls[1][0][2]

            self.assertEqual(left_flank_size, 500)
            self.assertEqual(right_flank_size, 500)


class TestConstants(unittest.TestCase):
    """Test module constants."""

    def test_known_disease_loci_columns(self):
        """Test KNOWN_DISEASE_ASSOCIATED_LOCI_COLUMNS constant."""
        expected_columns = [
            "LocusId", "LocusStructure", "RepeatUnit", "MainReferenceRegion",
            "Gene", "GeneRegion", "GeneId", "Diseases",
        ]
        self.assertEqual(KNOWN_DISEASE_ASSOCIATED_LOCI_COLUMNS, expected_columns)


class TestMainIntegration(unittest.TestCase):
    """Integration tests for main() function."""

    def setUp(self):
        """Create temporary files for testing."""
        self.temp_dir = tempfile.mkdtemp()

        # Create minimal reference FASTA
        self.ref_fasta = os.path.join(self.temp_dir, "ref.fa")
        with open(self.ref_fasta, 'w') as f:
            f.write(">chr1\n")
            f.write("CAGCAGCAGCAGCAGCAGCAGCAGCAGCAG" + "A"*1000 + "\n")
            f.write(">chr2\n")
            f.write("CTGCTGCTGCTGCTGCTGCTGCTGCTGCTG" + "A"*1000 + "\n")

        # Create FASTA index
        import pysam
        pysam.faidx(self.ref_fasta)

        # Create test catalog JSON
        self.catalog_json = os.path.join(self.temp_dir, "catalog.json")

    def create_catalog_file(self, records):
        """Helper to create catalog JSON file."""
        import simplejson
        with open(self.catalog_json, 'w') as f:
            simplejson.dump(records, f)

    @mock.patch('str_analysis.annotate_and_filter_str_catalog.download_local_copy')
    @mock.patch('str_analysis.annotate_and_filter_str_catalog.pyBigWig.open')
    def test_main_basic_processing(self, mock_bigwig, mock_download):
        """Test main() with basic catalog processing."""
        # Create simple catalog
        self.create_catalog_file([
            {
                "LocusId": "test1",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr1:0-30",
                "VariantType": "Repeat"
            }
        ])

        # Mock mappability bigwig
        mock_bw = mock.Mock()
        mock_bw.stats.return_value = (0.95,)
        mock_bigwig.return_value = mock_bw
        mock_download.return_value = "/tmp/mappability.bw"

        # Run main with minimal args
        sys.argv = [
            'test',
            '-R', self.ref_fasta,
            '--skip-gene-annotations',
            '--skip-disease-loci-annotations',
            '-o', os.path.join(self.temp_dir, 'output.json'),
            self.catalog_json
        ]

        from str_analysis.annotate_and_filter_str_catalog import main
        main()

        # Check output was created
        output_file = os.path.join(self.temp_dir, 'output.json')
        self.assertTrue(os.path.exists(output_file))

        # Verify output content
        import simplejson
        with open(output_file, 'r') as f:
            output_records = simplejson.load(f)

        self.assertEqual(len(output_records), 1)
        self.assertIn('NumRepeatsInReference', output_records[0])
        self.assertIn('ReferenceRepeatPurity', output_records[0])

    @mock.patch('str_analysis.annotate_and_filter_str_catalog.download_local_copy')
    @mock.patch('str_analysis.annotate_and_filter_str_catalog.pyBigWig.open')
    def test_main_motif_size_filter(self, mock_bigwig, mock_download):
        """Test main() with motif size filtering."""
        self.create_catalog_file([
            {
                "LocusId": "test1",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr1:0-30",
                "VariantType": "Repeat"
            },
            {
                "LocusId": "test2",
                "LocusStructure": "(AAAG)*",  # 4bp motif that won't simplify
                "ReferenceRegion": "chr1:100-200",
                "VariantType": "Repeat"
            }
        ])

        mock_bw = mock.Mock()
        mock_bw.stats.return_value = (0.95,)
        mock_bigwig.return_value = mock_bw
        mock_download.return_value = "/tmp/mappability.bw"

        sys.argv = [
            'test',
            '-R', self.ref_fasta,
            '--skip-gene-annotations',
            '--skip-disease-loci-annotations',
            '--min-motif-size', '3',
            '--max-motif-size', '3',
            '-o', os.path.join(self.temp_dir, 'output.json'),
            self.catalog_json
        ]

        from str_analysis.annotate_and_filter_str_catalog import main
        main()

        # Only 3bp motif should pass
        output_file = os.path.join(self.temp_dir, 'output.json')
        import simplejson
        with open(output_file, 'r') as f:
            output_records = simplejson.load(f)

        self.assertEqual(len(output_records), 1)
        self.assertEqual(output_records[0]['LocusId'], 'test1')

    @mock.patch('str_analysis.annotate_and_filter_str_catalog.download_local_copy')
    @mock.patch('str_analysis.annotate_and_filter_str_catalog.pyBigWig.open')
    def test_main_interval_size_filter(self, mock_bigwig, mock_download):
        """Test main() with interval size filtering."""
        self.create_catalog_file([
            {
                "LocusId": "short",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr1:0-30",
                "VariantType": "Repeat"
            },
            {
                "LocusId": "long",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr1:100-300",
                "VariantType": "Repeat"
            }
        ])

        mock_bw = mock.Mock()
        mock_bw.stats.return_value = (0.95,)
        mock_bigwig.return_value = mock_bw
        mock_download.return_value = "/tmp/mappability.bw"

        sys.argv = [
            'test',
            '-R', self.ref_fasta,
            '--skip-gene-annotations',
            '--skip-disease-loci-annotations',
            '--min-interval-size-bp', '100',
            '-o', os.path.join(self.temp_dir, 'output.json'),
            self.catalog_json
        ]

        from str_analysis.annotate_and_filter_str_catalog import main
        main()

        # Only long interval should pass
        output_file = os.path.join(self.temp_dir, 'output.json')
        import simplejson
        with open(output_file, 'r') as f:
            output_records = simplejson.load(f)

        self.assertEqual(len(output_records), 1)
        self.assertEqual(output_records[0]['LocusId'], 'long')

    @mock.patch('str_analysis.annotate_and_filter_str_catalog.download_local_copy')
    @mock.patch('str_analysis.annotate_and_filter_str_catalog.pyBigWig.open')
    def test_main_motif_simplification(self, mock_bigwig, mock_download):
        """Test main() simplifies complex motifs."""
        self.create_catalog_file([
            {
                "LocusId": "test1",
                "LocusStructure": "(CAGCAGCAG)*",  # Should simplify to (CAG)*
                "ReferenceRegion": "chr1:0-30",
                "VariantType": "Repeat"
            }
        ])

        mock_bw = mock.Mock()
        mock_bw.stats.return_value = (0.95,)
        mock_bigwig.return_value = mock_bw
        mock_download.return_value = "/tmp/mappability.bw"

        sys.argv = [
            'test',
            '-R', self.ref_fasta,
            '--skip-gene-annotations',
            '--skip-disease-loci-annotations',
            '-o', os.path.join(self.temp_dir, 'output.json'),
            self.catalog_json
        ]

        from str_analysis.annotate_and_filter_str_catalog import main
        main()

        output_file = os.path.join(self.temp_dir, 'output.json')
        import simplejson
        with open(output_file, 'r') as f:
            output_records = simplejson.load(f)

        # Motif should be simplified
        self.assertIn('(CAG)*', output_records[0]['LocusStructure'])

    @mock.patch('str_analysis.annotate_and_filter_str_catalog.download_local_copy')
    @mock.patch('str_analysis.annotate_and_filter_str_catalog.pyBigWig.open')
    def test_main_discard_invalid_bases_in_motif(self, mock_bigwig, mock_download):
        """Test main() filters loci with invalid bases in motif."""
        self.create_catalog_file([
            {
                "LocusId": "valid",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr1:0-30",
                "VariantType": "Repeat"
            },
            {
                "LocusId": "invalid",
                "LocusStructure": "(CAN)*",  # Contains N
                "ReferenceRegion": "chr1:100-200",
                "VariantType": "Repeat"
            }
        ])

        mock_bw = mock.Mock()
        mock_bw.stats.return_value = (0.95,)
        mock_bigwig.return_value = mock_bw
        mock_download.return_value = "/tmp/mappability.bw"

        sys.argv = [
            'test',
            '-R', self.ref_fasta,
            '--skip-gene-annotations',
            '--skip-disease-loci-annotations',
            '--discard-loci-with-non-ACGT-bases-in-motif',
            '-o', os.path.join(self.temp_dir, 'output.json'),
            self.catalog_json
        ]

        from str_analysis.annotate_and_filter_str_catalog import main
        main()

        output_file = os.path.join(self.temp_dir, 'output.json')
        import simplejson
        with open(output_file, 'r') as f:
            output_records = simplejson.load(f)

        # Only valid motif should pass
        self.assertEqual(len(output_records), 1)
        self.assertEqual(output_records[0]['LocusId'], 'valid')

    @mock.patch('str_analysis.annotate_and_filter_str_catalog.download_local_copy')
    @mock.patch('str_analysis.annotate_and_filter_str_catalog.pyBigWig.open')
    def test_main_set_locus_id(self, mock_bigwig, mock_download):
        """Test main() with --set-locus-id flag."""
        self.create_catalog_file([
            {
                "LocusId": "original_id",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr1:0-30",
                "VariantType": "Repeat"
            }
        ])

        mock_bw = mock.Mock()
        mock_bw.stats.return_value = (0.95,)
        mock_bigwig.return_value = mock_bw
        mock_download.return_value = "/tmp/mappability.bw"

        sys.argv = [
            'test',
            '-R', self.ref_fasta,
            '--skip-gene-annotations',
            '--skip-disease-loci-annotations',
            '--set-locus-id',
            '-o', os.path.join(self.temp_dir, 'output.json'),
            self.catalog_json
        ]

        from str_analysis.annotate_and_filter_str_catalog import main
        main()

        output_file = os.path.join(self.temp_dir, 'output.json')
        import simplejson
        with open(output_file, 'r') as f:
            output_records = simplejson.load(f)

        # LocusId should be reformatted as chrom-start-end-motif
        self.assertEqual(output_records[0]['LocusId'], '1-0-30-CAG')

    @mock.patch('str_analysis.annotate_and_filter_str_catalog.download_local_copy')
    @mock.patch('str_analysis.annotate_and_filter_str_catalog.pyBigWig.open')
    def test_main_output_bed_and_tsv(self, mock_bigwig, mock_download):
        """Test main() generates BED and TSV outputs."""
        self.create_catalog_file([
            {
                "LocusId": "test1",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr1:0-30",
                "VariantType": "Repeat"
            }
        ])

        mock_bw = mock.Mock()
        mock_bw.stats.return_value = (0.95,)
        mock_bigwig.return_value = mock_bw
        mock_download.return_value = "/tmp/mappability.bw"

        sys.argv = [
            'test',
            '-R', self.ref_fasta,
            '--skip-gene-annotations',
            '--skip-disease-loci-annotations',
            '--output-bed',
            '--output-tsv',
            '-o', os.path.join(self.temp_dir, 'output.json'),
            self.catalog_json
        ]

        from str_analysis.annotate_and_filter_str_catalog import main
        with mock.patch('os.system'):  # Mock bgzip/tabix calls
            main()

        # Check all outputs exist
        self.assertTrue(os.path.exists(os.path.join(self.temp_dir, 'output.json')))
        self.assertTrue(os.path.exists(os.path.join(self.temp_dir, 'output.tsv.gz')))
        self.assertTrue(os.path.exists(os.path.join(self.temp_dir, 'output.bed')))

    @mock.patch('str_analysis.annotate_and_filter_str_catalog.download_local_copy')
    @mock.patch('str_analysis.annotate_and_filter_str_catalog.pyBigWig.open')
    def test_main_locus_id_filter(self, mock_bigwig, mock_download):
        """Test main() filters by locus ID."""
        self.create_catalog_file([
            {
                "LocusId": "keep_this",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr1:0-30",
                "VariantType": "Repeat"
            },
            {
                "LocusId": "exclude_this",
                "LocusStructure": "(CTG)*",
                "ReferenceRegion": "chr2:0-30",
                "VariantType": "Repeat"
            }
        ])

        mock_bw = mock.Mock()
        mock_bw.stats.return_value = (0.95,)
        mock_bigwig.return_value = mock_bw
        mock_download.return_value = "/tmp/mappability.bw"

        sys.argv = [
            'test',
            '-R', self.ref_fasta,
            '--skip-gene-annotations',
            '--skip-disease-loci-annotations',
            '-l', 'keep_this',
            '-o', os.path.join(self.temp_dir, 'output.json'),
            self.catalog_json
        ]

        from str_analysis.annotate_and_filter_str_catalog import main
        main()

        output_file = os.path.join(self.temp_dir, 'output.json')
        import simplejson
        with open(output_file, 'r') as f:
            output_records = simplejson.load(f)

        self.assertEqual(len(output_records), 1)
        self.assertEqual(output_records[0]['LocusId'], 'keep_this')

    @mock.patch('str_analysis.annotate_and_filter_str_catalog.download_local_copy')
    @mock.patch('str_analysis.annotate_and_filter_str_catalog.pyBigWig.open')
    def test_main_exclude_chrom_filter(self, mock_bigwig, mock_download):
        """Test main() excludes specified chromosomes."""
        self.create_catalog_file([
            {
                "LocusId": "chr1_locus",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr1:0-30",
                "VariantType": "Repeat"
            },
            {
                "LocusId": "chr2_locus",
                "LocusStructure": "(CTG)*",
                "ReferenceRegion": "chr2:0-30",
                "VariantType": "Repeat"
            }
        ])

        mock_bw = mock.Mock()
        mock_bw.stats.return_value = (0.95,)
        mock_bigwig.return_value = mock_bw
        mock_download.return_value = "/tmp/mappability.bw"

        sys.argv = [
            'test',
            '-R', self.ref_fasta,
            '--skip-gene-annotations',
            '--skip-disease-loci-annotations',
            '-xc', 'chr2',
            '-o', os.path.join(self.temp_dir, 'output.json'),
            self.catalog_json
        ]

        from str_analysis.annotate_and_filter_str_catalog import main
        main()

        output_file = os.path.join(self.temp_dir, 'output.json')
        import simplejson
        with open(output_file, 'r') as f:
            output_records = simplejson.load(f)

        self.assertEqual(len(output_records), 1)
        self.assertEqual(output_records[0]['LocusId'], 'chr1_locus')

    @mock.patch('str_analysis.annotate_and_filter_str_catalog.download_local_copy')
    @mock.patch('str_analysis.annotate_and_filter_str_catalog.pyBigWig.open')
    def test_main_min_repeats_filter(self, mock_bigwig, mock_download):
        """Test main() filters by minimum number of repeats."""
        self.create_catalog_file([
            {
                "LocusId": "few_repeats",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr1:0-30",  # 10 repeats of CAG
                "VariantType": "Repeat"
            },
            {
                "LocusId": "many_repeats",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr1:100-300",  # 66 repeats of CAG
                "VariantType": "Repeat"
            }
        ])

        mock_bw = mock.Mock()
        mock_bw.stats.return_value = (0.95,)
        mock_bigwig.return_value = mock_bw
        mock_download.return_value = "/tmp/mappability.bw"

        sys.argv = [
            'test',
            '-R', self.ref_fasta,
            '--skip-gene-annotations',
            '--skip-disease-loci-annotations',
            '--min-repeats-in-reference', '20',
            '-o', os.path.join(self.temp_dir, 'output.json'),
            self.catalog_json
        ]

        from str_analysis.annotate_and_filter_str_catalog import main
        main()

        output_file = os.path.join(self.temp_dir, 'output.json')
        import simplejson
        with open(output_file, 'r') as f:
            output_records = simplejson.load(f)

        # Only locus with >= 20 repeats should pass
        self.assertEqual(len(output_records), 1)
        self.assertEqual(output_records[0]['LocusId'], 'many_repeats')

    @mock.patch('str_analysis.annotate_and_filter_str_catalog.download_local_copy')
    @mock.patch('str_analysis.annotate_and_filter_str_catalog.pyBigWig.open')
    def test_main_trim_loci(self, mock_bigwig, mock_download):
        """Test main() with --trim-loci flag."""
        self.create_catalog_file([
            {
                "LocusId": "test1",
                "LocusStructure": "(CAG)*",
                "ReferenceRegion": "chr1:0-31",  # 31bp is not divisible by 3
                "VariantType": "Repeat"
            }
        ])

        mock_bw = mock.Mock()
        mock_bw.stats.return_value = (0.95,)
        mock_bigwig.return_value = mock_bw
        mock_download.return_value = "/tmp/mappability.bw"

        sys.argv = [
            'test',
            '-R', self.ref_fasta,
            '--skip-gene-annotations',
            '--skip-disease-loci-annotations',
            '--trim-loci',
            '-o', os.path.join(self.temp_dir, 'output.json'),
            self.catalog_json
        ]

        from str_analysis.annotate_and_filter_str_catalog import main
        main()

        output_file = os.path.join(self.temp_dir, 'output.json')
        import simplejson
        with open(output_file, 'r') as f:
            output_records = simplejson.load(f)

        # ReferenceRegion should be trimmed to exact multiple of motif size
        self.assertEqual(output_records[0]['ReferenceRegion'], 'chr1:0-30')

    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        shutil.rmtree(self.temp_dir)


if __name__ == "__main__":
    unittest.main()
