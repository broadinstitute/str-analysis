import gzip
import os
import subprocess
import tempfile
import unittest
from unittest import mock

from str_analysis.extract_trid_metadata_from_TRGT_vcf import (
    DEFAULT_KNOWN_LOCI_CATALOG, load_known_locus_ids, main, no_sample_has_first_allele_called, parse_vcf_line)


def make_vcf_line(chrom, pos, trid, motif, genotypes, end=True):
    """Returns a TRGT-style VCF record line for an isolated repeat.

    Args:
        chrom (str): the CHROM column.
        pos (int): the POS column.
        trid (str): the TRID, which is also the LocusId for an isolated repeat.
        motif (str): the single motif in MOTIFS.
        genotypes (list): one GT string per sample.
        end (bool): whether INFO includes END, set to POS + 10.

    Returns:
        str: the tab-separated record, ending in a newline.
    """
    info = f"TRID={trid};MOTIFS={motif};STRUC=<TR1>" + (f";END={pos + 10}" if end else "")
    return "\t".join([chrom, str(pos), ".", "A", ".", ".", ".", info, "GT:AL"] +
                     [f"{gt}:10,10" for gt in genotypes]) + "\n"


def run_extractor(vcf_lines, extra_args=()):
    """Bgzip and index a VCF, run the extractor's main() on it, and return the rows it wrote.

    Args:
        vcf_lines (list): the VCF's lines, including its header lines.
        extra_args (list): additional command line args to pass to main().

    Returns:
        list: one list of column values per output row, without the header line.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        vcf_path = os.path.join(temp_dir, "test.vcf")
        with open(vcf_path, "wt") as f:
            f.writelines(vcf_lines)
        subprocess.run(["bgzip", vcf_path], check=True)
        subprocess.run(["tabix", "-p", "vcf", f"{vcf_path}.gz"], check=True)

        output_path = os.path.join(temp_dir, "trid_metadata.tsv.gz")
        argv = ["extract_trid_metadata_from_TRGT_vcf", "--input-vcf", f"{vcf_path}.gz",
                "--output-tsv", output_path, "--force"] + list(extra_args)
        with mock.patch("sys.argv", argv):
            main()

        with gzip.open(output_path, "rt") as f:
            next(f)
            return [line.rstrip("\n").split("\t") for line in f]


VCF_HEADER_LINES = [
    "##fileformat=VCFv4.2\n",
    "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "s1"]) + "\n",
]


def make_tmem185a_vcf_line(trid, genotype):
    """Returns a VCF record for the chrX:149631602-149631762 cluster from the AoU (Danzi) VCF."""
    return "\t".join(["chrX", "149631602", ".", "A", ".", ".", ".",
                      f"TRID={trid};END=149631762;MOTIFS=CGC,CGCCGT,TCC,GCT;STRUC=<VC272207>",
                      "GT:AL", f"{genotype}:10,10"]) + "\n"


TMEM185A_GENE_NAMED_TRID = ("X-149631602-149631617-TCC,X-149631685-149631694-GCT,TMEM185A_CGCCGT,"
                            "X-149631729-149631732-CGC")
TMEM185A_COORDINATE_TRID = ("X-149631602-149631617-TCC,X-149631685-149631694-GCT,"
                            "X-149631723-149631735-CGCCGT")


class Tests(unittest.TestCase):

    def test_record_is_skipped_only_when_no_sample_has_a_called_first_allele(self):
        for sample_fields in [["./.", "."], [".|.", "./0"], ["./0", ".|1"]]:
            self.assertTrue(no_sample_has_first_allele_called(["GT", "AL"], sample_fields), sample_fields)

        for sample_fields in [["0/.", "./."], ["./0", "1/1"], ["0"]]:
            self.assertFalse(no_sample_has_first_allele_called(["GT", "AL"], sample_fields), sample_fields)

    def test_gt_is_found_wherever_it_is_in_format(self):
        self.assertTrue(no_sample_has_first_allele_called(["AL", "GT"], ["10,12:./.", "10:."]))
        self.assertFalse(no_sample_has_first_allele_called(["AL", "GT"], ["10,12:./.", "10,12:0/1"]))

        # without a GT field there is nothing to decide on, so the record is kept
        self.assertFalse(no_sample_has_first_allele_called(["AL"], ["10,12"]))

    def test_parse_vcf_line_skips_records_trgt_lps_drops(self):
        self.assertEqual(parse_vcf_line(make_vcf_line("chr1", 1000, "1-1000-1010-AC", "AC", ["./0", "./."])), ())
        self.assertEqual(
            parse_vcf_line(make_vcf_line("chr1", 1000, "1-1000-1010-AC", "AC", ["0/.", "./."])),
            [("1-1000-1010-AC", "1-1000-1010-AC", "AC", "1:1000-1010", "")])

    def test_record_without_end_is_skipped(self):
        self.assertEqual(
            parse_vcf_line(make_vcf_line("chr1", 1000, "1-1000-1010-AC", "AC", ["0/1"], end=False)), ())
        self.assertEqual(
            parse_vcf_line(make_vcf_line("chr1", 1000, "1-1000-1010-AC", "AC", ["./."], end=False)), ())

    def test_known_locus_ids_map_gene_named_repeats_to_coordinates(self):
        known_locus_ids = load_known_locus_ids(DEFAULT_KNOWN_LOCI_CATALOG)

        self.assertEqual(known_locus_ids["TMEM185A_CGCCGT"], "X-149631723-149631735-CGCCGT")
        self.assertEqual(known_locus_ids["TMEM185A"], "X-149631735-149631780-CGC")
        self.assertEqual(known_locus_ids["HTT"], "4-3074876-3074933-CAG")
        self.assertEqual(known_locus_ids["HTT_CCG"], "4-3074939-3074966-CCG")
        self.assertEqual(known_locus_ids["EP400"], "12-132062548-132062611-CAG")

    def test_gene_named_ids_in_a_trid_get_metadata_rows(self):
        # record 2 from the AoU (Danzi) VCF, whose TRID names the TMEM185A repeat by gene
        trid = "X-149631602-149631617-TCC,X-149631685-149631694-GCT,TMEM185A_CGCCGT,X-149631729-149631732-CGC"
        line = "\t".join([
            "chrX", "149631602", ".", "A", ".", ".", ".",
            f"TRID={trid};END=149631762;MOTIFS=CGC,CGCCGT,TCC,GCT;STRUC=<VC272207>", "GT", "0/1"]) + "\n"

        rows = parse_vcf_line(line, load_known_locus_ids(DEFAULT_KNOWN_LOCI_CATALOG))

        self.assertEqual([(locus_id, motif) for _, locus_id, motif, _, _ in rows], [
            ("X-149631729-149631732-CGC", "CGC"),
            ("X-149631723-149631735-CGCCGT", "CGCCGT"),
            ("X-149631602-149631617-TCC", "TCC"),
            ("X-149631685-149631694-GCT", "GCT"),
        ])

        # without the map, the gene-named repeat gets no row
        self.assertNotIn("CGCCGT", [motif for _, _, motif, _, _ in parse_vcf_line(line)])

    def test_every_contig_in_the_index_is_extracted(self):
        # contig names without the "chr" prefix, plus an alt contig
        vcf_lines = VCF_HEADER_LINES + [
            make_vcf_line("1", 1000, "1-1000-1010-AC", "AC", ["0/1"]),
            make_vcf_line("1_KI270706v1_random", 2000, "1_KI270706v1_random-2000-2010-AG", "AG", ["0/0"]),
            make_vcf_line("X", 3000, "X-3000-3010-AAT", "AAT", ["1"]),
        ]

        self.assertEqual([row[1] for row in run_extractor(vcf_lines)], [
            "1-1000-1010-AC", "1_KI270706v1_random-2000-2010-AG", "X-3000-3010-AAT"])
        self.assertEqual(
            [row[1] for row in run_extractor(vcf_lines, ["--contig", "X", "--contig", "chrX"])],
            ["X-3000-3010-AAT"])

    def test_two_records_covering_the_same_repeats_both_get_rows(self):
        # the AoU (Danzi) VCF lists this cluster twice, once with a gene-named repeat id and once
        # with coordinate ids. Their TRIDs tell them apart, so both get rows here and the converter
        # decides what to write for the repeats they share.
        rows = run_extractor(VCF_HEADER_LINES + [
            make_tmem185a_vcf_line(TMEM185A_GENE_NAMED_TRID, "0/1"),
            make_tmem185a_vcf_line(TMEM185A_COORDINATE_TRID, "0/1"),
        ])

        self.assertEqual(sorted(row[1] for row in rows if row[0] == TMEM185A_GENE_NAMED_TRID), [
            "X-149631602-149631617-TCC", "X-149631685-149631694-GCT",
            "X-149631723-149631735-CGCCGT", "X-149631729-149631732-CGC"])
        self.assertEqual(sorted(row[1] for row in rows if row[0] == TMEM185A_COORDINATE_TRID), [
            "X-149631602-149631617-TCC", "X-149631685-149631694-GCT",
            "X-149631723-149631735-CGCCGT"])

    def test_the_same_record_twice_raises(self):
        with self.assertRaisesRegex(RuntimeError, "the VCF contains the same record twice"):
            run_extractor(VCF_HEADER_LINES + [
                make_tmem185a_vcf_line(TMEM185A_COORDINATE_TRID, "0/1"),
                make_tmem185a_vcf_line(TMEM185A_COORDINATE_TRID, "0/1"),
            ])


if __name__ == "__main__":
    unittest.main()
