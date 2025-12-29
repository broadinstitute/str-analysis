#!/usr/bin/env python3

"""Comprehensive tests for gtf_utils.py"""

import collections
import os
import tempfile
import unittest
from unittest import mock
from intervaltree import Interval, IntervalTree

from str_analysis.utils.gtf_utils import (
    parse_gtf_to_interval_trees,
    generate_gtf_records,
    compute_UTR_type,
    compute_genomic_region_of_interval,
    GENE_MODELS,
    PROMOTER_SIZE,
    GENCODE_INTERVAL_TREES,
    TRANSCRIPT_ID_TO_CDS_COORDS_MAP,
)


class TestConstants(unittest.TestCase):
    """Test module constants."""

    def test_gene_models_constant(self):
        """Test GENE_MODELS constant."""
        self.assertIn("gencode", GENE_MODELS)
        self.assertIn("mane", GENE_MODELS)
        self.assertIn("refseq", GENE_MODELS)

        # Check they're gs:// paths
        for model, path in GENE_MODELS.items():
            self.assertTrue(path.startswith("gs://"), f"{model} path should start with gs://")
            self.assertTrue(path.endswith(".gtf.gz"), f"{model} path should end with .gtf.gz")

    def test_promoter_size_constant(self):
        """Test PROMOTER_SIZE constant."""
        self.assertEqual(PROMOTER_SIZE, 1000)
        self.assertIsInstance(PROMOTER_SIZE, int)


class TestGenerateGtfRecords(unittest.TestCase):
    """Test generate_gtf_records function."""

    def setUp(self):
        """Create temporary GTF file."""
        self.temp_dir = tempfile.mkdtemp()
        self.gtf_file = os.path.join(self.temp_dir, "test.gtf")

    def create_gtf_file(self, content):
        """Helper to create GTF file with given content."""
        with open(self.gtf_file, 'w') as f:
            f.write(content)

    def test_generate_gtf_records_basic(self):
        """Test basic GTF record parsing."""
        gtf_content = """##gtf-version 2.2
chr1\tHAVANA\texon\t1000\t2000\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
"""
        self.create_gtf_file(gtf_content)

        records = list(generate_gtf_records(self.gtf_file))

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["chrom"], "chr1")
        self.assertEqual(records[0]["feature_type"], "exon")
        self.assertEqual(records[0]["start_1based"], 1000)
        self.assertEqual(records[0]["end_1based"], 2000)
        self.assertEqual(records[0]["strand"], "+")
        self.assertEqual(records[0]["gene_id"], "ENSG00000123")
        self.assertEqual(records[0]["transcript_id"], "ENST00000456")
        self.assertEqual(records[0]["gene_name"], "TEST1")

    def test_generate_gtf_records_skip_comments(self):
        """Test that comment lines are skipped."""
        gtf_content = """##gtf-version 2.2
# This is a comment
chr1\tHAVANA\texon\t1000\t2000\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
"""
        self.create_gtf_file(gtf_content)

        records = list(generate_gtf_records(self.gtf_file))
        self.assertEqual(len(records), 1)

    def test_generate_gtf_records_filter_feature_types(self):
        """Test that only relevant feature types are included."""
        gtf_content = """chr1\tHAVANA\texon\t1000\t2000\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
chr1\tHAVANA\tCDS\t1100\t1900\t.\t+\t0\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
chr1\tHAVANA\tUTR\t1000\t1099\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
chr1\tHAVANA\tgene\t1000\t5000\t.\t+\t.\tgene_id "ENSG00000123"; gene_name "TEST1";
chr1\tHAVANA\tstart_codon\t1100\t1102\t.\t+\t0\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
"""
        self.create_gtf_file(gtf_content)

        records = list(generate_gtf_records(self.gtf_file))

        # Should only have exon, CDS, and UTR (gene and start_codon are filtered out)
        self.assertEqual(len(records), 3)
        feature_types = [r["feature_type"] for r in records]
        self.assertIn("exon", feature_types)
        self.assertIn("CDS", feature_types)
        self.assertIn("UTR", feature_types)
        self.assertNotIn("gene", feature_types)
        self.assertNotIn("start_codon", feature_types)

    def test_generate_gtf_records_multiple_transcripts(self):
        """Test parsing GTF with multiple transcripts."""
        gtf_content = """chr1\tHAVANA\texon\t1000\t2000\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
chr1\tHAVANA\texon\t3000\t4000\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000457"; gene_name "TEST1";
chr2\tHAVANA\texon\t5000\t6000\t.\t-\t.\tgene_id "ENSG00000789"; transcript_id "ENST00000999"; gene_name "TEST2";
"""
        self.create_gtf_file(gtf_content)

        records = list(generate_gtf_records(self.gtf_file))

        self.assertEqual(len(records), 3)
        transcript_ids = [r["transcript_id"] for r in records]
        self.assertIn("ENST00000456", transcript_ids)
        self.assertIn("ENST00000457", transcript_ids)
        self.assertIn("ENST00000999", transcript_ids)

    def test_generate_gtf_records_negative_strand(self):
        """Test parsing records on negative strand."""
        gtf_content = """chr1\tHAVANA\texon\t1000\t2000\t.\t-\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
"""
        self.create_gtf_file(gtf_content)

        records = list(generate_gtf_records(self.gtf_file))

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["strand"], "-")

    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        shutil.rmtree(self.temp_dir)


class TestParseGtfToIntervalTrees(unittest.TestCase):
    """Test parse_gtf_to_interval_trees function."""

    def setUp(self):
        """Create temporary GTF file."""
        self.temp_dir = tempfile.mkdtemp()
        self.gtf_file = os.path.join(self.temp_dir, "test.gtf")

        # Clear any cached interval trees
        GENCODE_INTERVAL_TREES.clear()
        TRANSCRIPT_ID_TO_CDS_COORDS_MAP.clear()

    def create_gtf_file(self, content):
        """Helper to create GTF file."""
        with open(self.gtf_file, 'w') as f:
            f.write(content)

    def test_parse_gtf_to_interval_trees_basic(self):
        """Test basic GTF parsing into interval trees."""
        gtf_content = """chr1\tHAVANA\texon\t1000\t2000\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
chr2\tHAVANA\texon\t3000\t4000\t.\t+\t.\tgene_id "ENSG00000789"; transcript_id "ENST00000999"; gene_name "TEST2";
"""
        self.create_gtf_file(gtf_content)

        parse_gtf_to_interval_trees(self.gtf_file)

        # Check interval trees were created
        self.assertIn(self.gtf_file, GENCODE_INTERVAL_TREES)
        interval_trees = GENCODE_INTERVAL_TREES[self.gtf_file]

        # Check chromosome keys (without "chr" prefix)
        self.assertIn("1", interval_trees)
        self.assertIn("2", interval_trees)

        # Check intervals
        chr1_intervals = list(interval_trees["1"])
        self.assertEqual(len(chr1_intervals), 1)
        self.assertEqual(chr1_intervals[0].begin, 999)  # 0-based
        self.assertEqual(chr1_intervals[0].end, 2000)

    def test_parse_gtf_to_interval_trees_cds_caching(self):
        """Test that CDS coordinates are cached."""
        gtf_content = """chr1\tHAVANA\tCDS\t1100\t1900\t.\t+\t0\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
"""
        self.create_gtf_file(gtf_content)

        parse_gtf_to_interval_trees(self.gtf_file)

        # Check CDS coordinates were cached
        self.assertIn(self.gtf_file, TRANSCRIPT_ID_TO_CDS_COORDS_MAP)
        cds_map = TRANSCRIPT_ID_TO_CDS_COORDS_MAP[self.gtf_file]

        self.assertIn("ENST00000456", cds_map)
        cds_coords = cds_map["ENST00000456"]
        self.assertEqual(cds_coords[0], "chr1")  # chrom
        self.assertEqual(cds_coords[1], 1100)    # start_1based
        self.assertEqual(cds_coords[2], 1900)    # end_1based
        self.assertEqual(cds_coords[3], "+")     # strand

    def test_parse_gtf_to_interval_trees_with_limit(self):
        """Test parsing with record limit."""
        gtf_content = """chr1\tHAVANA\texon\t1000\t2000\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
chr1\tHAVANA\texon\t3000\t4000\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000457"; gene_name "TEST1";
chr1\tHAVANA\texon\t5000\t6000\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000458"; gene_name "TEST1";
"""
        self.create_gtf_file(gtf_content)

        parse_gtf_to_interval_trees(self.gtf_file, n=2)

        interval_trees = GENCODE_INTERVAL_TREES[self.gtf_file]
        chr1_intervals = list(interval_trees["1"])

        # Should only have 2 intervals (limited by n=2)
        self.assertEqual(len(chr1_intervals), 2)

    def test_parse_gtf_to_interval_trees_with_verbose(self):
        """Test parsing with verbose output."""
        gtf_content = """chr1\tHAVANA\texon\t1000\t1200\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
"""
        self.create_gtf_file(gtf_content)

        with mock.patch('builtins.print') as mock_print:
            parse_gtf_to_interval_trees(self.gtf_file, verbose=True)
            # Should print parsing and counter messages
            self.assertTrue(mock_print.called)
            call_args = [str(call) for call in mock_print.call_args_list]
            output = ' '.join(call_args)
            self.assertIn("Parsing", output)

    def test_generate_gtf_records_5utr_3utr_conversion(self):
        """Test that 5UTR and 3UTR are converted to UTR."""
        gtf_content = """chr1\tHAVANA\t5UTR\t1000\t1099\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
chr1\tHAVANA\t3UTR\t1900\t2000\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
"""
        self.create_gtf_file(gtf_content)
        records = list(generate_gtf_records(self.gtf_file))

        # Both should be converted to UTR
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0]["feature_type"], "UTR")
        self.assertEqual(records[1]["feature_type"], "UTR")

    def test_generate_gtf_records_unexpected_feature_warning(self):
        """Test warning for unexpected feature types."""
        gtf_content = """chr1\tHAVANA\tunknown_feature\t1000\t1100\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
chr1\tHAVANA\texon\t1000\t1100\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
"""
        self.create_gtf_file(gtf_content)

        with mock.patch('builtins.print') as mock_print:
            records = list(generate_gtf_records(self.gtf_file))
            # Should print warning for unknown_feature
            mock_print.assert_called()
            warning_call = str(mock_print.call_args)
            self.assertIn("WARNING", warning_call)
            self.assertIn("unknown_feature", warning_call)

        # Should only have the exon record
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["feature_type"], "exon")

    def test_generate_gtf_records_skip_invalid_coords(self):
        """Test that records with end <= start are skipped."""
        gtf_content = """chr1\tHAVANA\texon\t1100\t1000\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
chr1\tHAVANA\texon\t1000\t1100\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
"""
        self.create_gtf_file(gtf_content)
        records = list(generate_gtf_records(self.gtf_file))

        # Should skip the first record with end < start
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["start_1based"], 1000)
        self.assertEqual(records[0]["end_1based"], 1100)

    def test_generate_gtf_records_negative_strand_promoter(self):
        """Test promoter generation on negative strand."""
        gtf_content = """chr1\tHAVANA\ttranscript\t1000\t2000\t.\t-\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
"""
        self.create_gtf_file(gtf_content)
        records = list(generate_gtf_records(self.gtf_file))

        # Should generate both transcript and promoter
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0]["feature_type"], "transcript")
        self.assertEqual(records[1]["feature_type"], "promoter")
        # Promoter on negative strand is downstream of transcript end
        self.assertEqual(records[1]["start_1based"], 2001)

    def test_generate_gtf_records_invalid_strand(self):
        """Test error handling for invalid strand."""
        gtf_content = """chr1\tHAVANA\ttranscript\t1000\t2000\t.\tX\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
"""
        self.create_gtf_file(gtf_content)

        with self.assertRaises(ValueError) as cm:
            list(generate_gtf_records(self.gtf_file))
        self.assertIn("Unexpected strand value", str(cm.exception))

    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        shutil.rmtree(self.temp_dir)
        GENCODE_INTERVAL_TREES.clear()
        TRANSCRIPT_ID_TO_CDS_COORDS_MAP.clear()


class TestComputeUTRType(unittest.TestCase):
    """Test compute_UTR_type function."""

    def test_compute_utr_type_5prime_positive_strand(self):
        """Test 5' UTR on positive strand."""
        utr_record = {
            "feature_type": "UTR",
            "chrom": "chr1",
            "start_1based": 1000,
            "end_1based": 1099,
            "strand": "+",
            "transcript_id": "ENST00000456"
        }

        # CDS starts at 1100
        cds_coords = {
            "ENST00000456": ("chr1", 1100, 1900, "+")
        }

        result = compute_UTR_type(utr_record, cds_coords)
        self.assertEqual(result, "5' UTR")

    def test_compute_utr_type_3prime_positive_strand(self):
        """Test 3' UTR on positive strand."""
        utr_record = {
            "feature_type": "UTR",
            "chrom": "chr1",
            "start_1based": 1901,
            "end_1based": 2000,
            "strand": "+",
            "transcript_id": "ENST00000456"
        }

        # CDS ends at 1900
        cds_coords = {
            "ENST00000456": ("chr1", 1100, 1900, "+")
        }

        result = compute_UTR_type(utr_record, cds_coords)
        self.assertEqual(result, "3' UTR")

    def test_compute_utr_type_5prime_negative_strand(self):
        """Test 5' UTR on negative strand."""
        utr_record = {
            "feature_type": "UTR",
            "chrom": "chr1",
            "start_1based": 1901,
            "end_1based": 2000,
            "strand": "-",
            "transcript_id": "ENST00000456"
        }

        # CDS ends at 1900 (but on - strand, so 5' is downstream)
        cds_coords = {
            "ENST00000456": ("chr1", 1100, 1900, "-")
        }

        result = compute_UTR_type(utr_record, cds_coords)
        self.assertEqual(result, "5' UTR")

    def test_compute_utr_type_3prime_negative_strand(self):
        """Test 3' UTR on negative strand."""
        utr_record = {
            "feature_type": "UTR",
            "chrom": "chr1",
            "start_1based": 1000,
            "end_1based": 1099,
            "strand": "-",
            "transcript_id": "ENST00000456"
        }

        # CDS starts at 1100 (but on - strand, so 3' is upstream)
        cds_coords = {
            "ENST00000456": ("chr1", 1100, 1900, "-")
        }

        result = compute_UTR_type(utr_record, cds_coords)
        self.assertEqual(result, "3' UTR")

    def test_compute_utr_type_no_cds_found(self):
        """Test when CDS not found for transcript."""
        utr_record = {
            "feature_type": "UTR",
            "chrom": "chr1",
            "start_1based": 1000,
            "end_1based": 1099,
            "strand": "+",
            "transcript_id": "ENST00000999"
        }

        cds_coords = {}  # No CDS info

        with mock.patch('builtins.print') as mock_print:
            result = compute_UTR_type(utr_record, cds_coords)
            self.assertEqual(result, "UTR")
            mock_print.assert_called_once()
            self.assertIn("CDS not found", mock_print.call_args[0][0])

    def test_compute_utr_type_chrom_mismatch(self):
        """Test when UTR and CDS chromosomes don't match."""
        utr_record = {
            "feature_type": "UTR",
            "chrom": "chr2",
            "start_1based": 1000,
            "end_1based": 1099,
            "strand": "+",
            "transcript_id": "ENST00000456"
        }

        cds_coords = {
            "ENST00000456": ("chr1", 1100, 1900, "+")  # Different chrom
        }

        with mock.patch('builtins.print') as mock_print:
            result = compute_UTR_type(utr_record, cds_coords)
            self.assertEqual(result, "UTR")
            mock_print.assert_called_once()
            self.assertIn("chromosome", mock_print.call_args[0][0].lower())

    def test_compute_utr_type_strand_mismatch(self):
        """Test when UTR and CDS strands don't match."""
        utr_record = {
            "feature_type": "UTR",
            "chrom": "chr1",
            "start_1based": 1000,
            "end_1based": 1099,
            "strand": "+",
            "transcript_id": "ENST00000456"
        }

        cds_coords = {
            "ENST00000456": ("chr1", 1100, 1900, "-")  # Different strand
        }

        with mock.patch('builtins.print') as mock_print:
            result = compute_UTR_type(utr_record, cds_coords)
            self.assertEqual(result, "UTR")
            mock_print.assert_called_once()
            self.assertIn("strand", mock_print.call_args[0][0].lower())

    def test_compute_utr_type_invalid_feature(self):
        """Test error when feature type is not UTR."""
        utr_record = {
            "feature_type": "exon",  # Not UTR
            "chrom": "chr1",
            "start_1based": 1000,
            "end_1based": 1099,
            "strand": "+",
            "transcript_id": "ENST00000456"
        }

        cds_coords = {}

        # The actual code raises the feature_type string, which causes TypeError
        with self.assertRaises(TypeError):
            compute_UTR_type(utr_record, cds_coords)

    def test_compute_utr_type_ambiguous_utr_warning(self):
        """Test warning when UTR type cannot be determined."""
        utr_record = {
            "feature_type": "UTR",
            "chrom": "chr1",
            "start_1based": 1500,  # Overlaps CDS in an ambiguous way
            "end_1based": 1600,
            "strand": "+",
            "transcript_id": "ENST00000456"
        }

        # CDS is 1100-1900, so 1500-1600 overlaps CDS but in a way that's neither clearly 5' nor 3'
        cds_coords = {
            "ENST00000456": ("chr1", 1100, 1900, "+")
        }

        with mock.patch('builtins.print') as mock_print:
            result = compute_UTR_type(utr_record, cds_coords)
            # Should return generic "UTR" and print warning
            self.assertEqual(result, "UTR")
            mock_print.assert_called_once()
            warning = mock_print.call_args[0][0]
            self.assertIn("WARNING", warning)
            self.assertIn("Something wrong with UTR info", warning)


class TestComputeGenomicRegionOfInterval(unittest.TestCase):
    """Test compute_genomic_region_of_interval function."""

    def setUp(self):
        """Create temporary GTF file with test data."""
        self.temp_dir = tempfile.mkdtemp()
        self.gtf_file = os.path.join(self.temp_dir, "test.gtf")

        # Clear any cached interval trees
        GENCODE_INTERVAL_TREES.clear()
        TRANSCRIPT_ID_TO_CDS_COORDS_MAP.clear()

        # Create comprehensive test GTF
        gtf_content = """chr1\tHAVANA\ttranscript\t1000\t5000\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
chr1\tHAVANA\texon\t1000\t1200\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
chr1\tHAVANA\tCDS\t1100\t1900\t.\t+\t0\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
chr1\tHAVANA\texon\t1800\t2000\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
chr1\tHAVANA\tUTR\t1000\t1099\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
chr1\tHAVANA\tUTR\t1901\t2000\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
chr1\tHAVANA\tpromoter\t500\t999\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
"""
        with open(self.gtf_file, 'w') as f:
            f.write(gtf_content)

    def test_compute_genomic_region_cds(self):
        """Test interval in CDS region."""
        region, gene_name, gene_id, transcript_id = compute_genomic_region_of_interval(
            "chr1", 1150, 1200, genes_gtf_path=self.gtf_file
        )

        self.assertEqual(region, "CDS")
        self.assertEqual(gene_name, "TEST1")
        self.assertEqual(gene_id, "ENSG00000123")
        self.assertEqual(transcript_id, "ENST00000456")

    def test_compute_genomic_region_5prime_utr(self):
        """Test interval in 5' UTR."""
        region, gene_name, gene_id, transcript_id = compute_genomic_region_of_interval(
            "chr1", 1000, 1099, genes_gtf_path=self.gtf_file
        )

        self.assertEqual(region, "5' UTR")
        self.assertEqual(gene_name, "TEST1")

    def test_compute_genomic_region_3prime_utr(self):
        """Test interval in 3' UTR."""
        region, gene_name, gene_id, transcript_id = compute_genomic_region_of_interval(
            "chr1", 1901, 2000, genes_gtf_path=self.gtf_file
        )

        self.assertEqual(region, "3' UTR")
        self.assertEqual(gene_name, "TEST1")

    def test_compute_genomic_region_exon(self):
        """Test interval in exon (non-CDS)."""
        region, gene_name, gene_id, transcript_id = compute_genomic_region_of_interval(
            "chr1", 1800, 1850, genes_gtf_path=self.gtf_file
        )

        # Should be CDS since it overlaps CDS region
        self.assertEqual(region, "CDS")

    def test_compute_genomic_region_intron(self):
        """Test interval in intron."""
        # Region 2100-2200 is within transcript (1000-5000) but after all exons/CDS/UTR
        region, gene_name, gene_id, transcript_id = compute_genomic_region_of_interval(
            "chr1", 2100, 2200, genes_gtf_path=self.gtf_file
        )

        self.assertEqual(region, "intron")
        self.assertEqual(gene_name, "TEST1")

    def test_compute_genomic_region_promoter(self):
        """Test interval in promoter."""
        region, gene_name, gene_id, transcript_id = compute_genomic_region_of_interval(
            "chr1", 600, 800, genes_gtf_path=self.gtf_file
        )

        self.assertEqual(region, "promoter")
        self.assertEqual(gene_name, "TEST1")

    def test_compute_genomic_region_intergenic(self):
        """Test interval in intergenic region."""
        region, gene_name, gene_id, transcript_id = compute_genomic_region_of_interval(
            "chr1", 10000, 11000, genes_gtf_path=self.gtf_file
        )

        self.assertEqual(region, "intergenic")
        self.assertIsNone(gene_name)
        self.assertIsNone(gene_id)
        self.assertIsNone(transcript_id)

    def test_compute_genomic_region_chrom_normalization(self):
        """Test that chromosome names are normalized."""
        # Test with "chr" prefix
        region1, _, _, _ = compute_genomic_region_of_interval(
            "chr1", 1150, 1200, genes_gtf_path=self.gtf_file
        )

        # Test without "chr" prefix
        region2, _, _, _ = compute_genomic_region_of_interval(
            "1", 1150, 1200, genes_gtf_path=self.gtf_file
        )

        self.assertEqual(region1, region2)
        self.assertEqual(region1, "CDS")

    def test_compute_genomic_region_caching(self):
        """Test that interval trees are cached."""
        # First call
        compute_genomic_region_of_interval("chr1", 1150, 1200, genes_gtf_path=self.gtf_file)

        # Check cache
        self.assertIn(self.gtf_file, GENCODE_INTERVAL_TREES)

        # Second call should use cache
        region, _, _, _ = compute_genomic_region_of_interval(
            "chr1", 1150, 1200, genes_gtf_path=self.gtf_file
        )
        self.assertEqual(region, "CDS")

    def test_compute_genomic_region_default_gtf_path(self):
        """Test that default GTF path is used when genes_gtf_path=None."""
        # Mock the default GENE_MODELS to use our test GTF
        import str_analysis.utils.gtf_utils as gtf_utils
        original_gene_models = dict(gtf_utils.GENE_MODELS)
        original_value = gtf_utils.GENE_MODELS["gencode"]

        try:
            gtf_utils.GENE_MODELS["gencode"] = self.gtf_file

            # Call without specifying genes_gtf_path (should use default)
            region, gene_name, _, _ = compute_genomic_region_of_interval(
                "chr1", 1150, 1200, genes_gtf_path=None
            )

            self.assertEqual(region, "CDS")
            self.assertEqual(gene_name, "TEST1")
        finally:
            # Restore original
            gtf_utils.GENE_MODELS["gencode"] = original_value
            # Clear the cache that was created with test GTF
            if self.gtf_file in GENCODE_INTERVAL_TREES:
                del GENCODE_INTERVAL_TREES[self.gtf_file]
            if self.gtf_file in TRANSCRIPT_ID_TO_CDS_COORDS_MAP:
                del TRANSCRIPT_ID_TO_CDS_COORDS_MAP[self.gtf_file]

    def test_compute_genomic_region_priority_order(self):
        """Test that feature types are returned in priority order."""
        # Create GTF with overlapping features
        gtf_content = """chr1\tHAVANA\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
chr1\tHAVANA\texon\t1000\t2000\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
chr1\tHAVANA\tCDS\t1100\t1900\t.\t+\t0\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
chr1\tHAVANA\tUTR\t1000\t1099\t.\t+\t.\tgene_id "ENSG00000123"; transcript_id "ENST00000456"; gene_name "TEST1";
"""
        test_gtf = os.path.join(self.temp_dir, "priority.gtf")
        with open(test_gtf, 'w') as f:
            f.write(gtf_content)

        # Query interval that overlaps both CDS and UTR
        # CDS should take priority
        region, _, _, _ = compute_genomic_region_of_interval(
            "chr1", 1050, 1150, genes_gtf_path=test_gtf
        )

        # CDS has highest priority
        self.assertEqual(region, "CDS")

    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        shutil.rmtree(self.temp_dir)
        GENCODE_INTERVAL_TREES.clear()
        TRANSCRIPT_ID_TO_CDS_COORDS_MAP.clear()


if __name__ == "__main__":
    unittest.main()
