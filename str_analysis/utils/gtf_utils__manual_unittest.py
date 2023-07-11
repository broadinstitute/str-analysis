import unittest

from str_analysis.utils.file_utils import download_local_copy
from str_analysis.utils.gtf_utils import compute_genomic_region_of_interval, GENE_MODELS


class GTFUtilsTests(unittest.TestCase):

    def setUp(self):
        self.gencode_gtf_path = download_local_copy(GENE_MODELS["gencode"])

    def test_compute_genomic_region_of_interval(self):
        # intron
        result = compute_genomic_region_of_interval("chr1", 70_075_754, 70_075_973,
                                                    genes_gtf_path=self.gencode_gtf_path,
                                                    verbose=True)
        self.assertEqual(result, ('intron', 'LRRC7', 'ENSG00000033122', 'ENST00000310961'))

        # coding
        result = compute_genomic_region_of_interval("chr1", 70_038_218, 70_038_948,
                                                    genes_gtf_path=self.gencode_gtf_path,
                                                    verbose=True)
        self.assertEqual(result, ('CDS', 'LRRC7', 'ENSG00000033122', 'ENST00000310961'))

        # promoter
        result = compute_genomic_region_of_interval("chr1", 69_567_250, 69_567_319,
                                                    genes_gtf_path=self.gencode_gtf_path,
                                                    verbose=True)
        self.assertEqual(result, ('promoter', 'LRRC7', 'ENSG00000033122', 'ENST00000651989'))


        # 3' UTR
        result = compute_genomic_region_of_interval("chr1", 69_567_415, 69_568_261,
                                                    genes_gtf_path=self.gencode_gtf_path,
                                                    verbose=True)
        self.assertEqual(result, ("5' UTR", 'LRRC7', 'ENSG00000033122', 'ENST00000651989'))

        # 5' UTR
        result = compute_genomic_region_of_interval("chr1", 70_122_558, 70_127_178,
                                                    genes_gtf_path=self.gencode_gtf_path,
                                                    verbose=True)
        self.assertEqual(result, ("3' UTR", 'LRRC7', 'ENSG00000033122', 'ENST00000310961'))

    def test_compute_genomic_region_with_partial_overlap(self):
        # intron
        result = compute_genomic_region_of_interval("5", 70_049_769, 70_049_808,
                                                    genes_gtf_path=self.gencode_gtf_path,
                                                    verbose=True)
        self.assertEqual(result, ('intron', 'SMN2', 'ENSG00000205571', 'ENST00000380741'))

        # coding
        result = compute_genomic_region_of_interval("5", 70_049_743, 70_049_782,
                                                    genes_gtf_path=self.gencode_gtf_path,
                                                    verbose=True)
        self.assertEqual(result, ('CDS', 'SMN2', 'ENSG00000205571', 'ENST00000380741'))

        # promoter
        result = compute_genomic_region_of_interval("5", 70_049_480, 70_049_519,
                                                    genes_gtf_path=self.gencode_gtf_path,
                                                    verbose=True)
        self.assertEqual(result, ('promoter', 'SMN2', 'ENSG00000205571', 'ENST00000380741'))

        # 5' UTR
        result = compute_genomic_region_of_interval("5", 70_049_606, 70_049_663,
                                                    genes_gtf_path=self.gencode_gtf_path,
                                                    verbose=True)
        self.assertEqual(result, ("5' UTR", 'SMN2', 'ENSG00000205571', 'ENST00000506734'))

        # 3' UTR
        result = compute_genomic_region_of_interval("5", 70_077_528, 70_077_674,
                                                    genes_gtf_path=self.gencode_gtf_path,
                                                    verbose=True)
        self.assertEqual(result, ("3' UTR", 'SMN2', 'ENSG00000205571', 'ENST00000380742'))

        # intergenic
        result = compute_genomic_region_of_interval("5", 70_079_184, 70_080_487,
                                                    genes_gtf_path=self.gencode_gtf_path,
                                                    verbose=True)
        self.assertEqual(result, ("intergenic", None, None, None))

