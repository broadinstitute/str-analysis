import gzip
import os
import tempfile
import unittest

from str_analysis.extract_allele_sequences_from_vcf import (
    CALLED, NO_CALL, EXTRACTION_OK, EXTRACTION_REF_MISMATCH, EXTRACTION_MULTIALLELIC, OUTPUT_COLUMNS,
    extract_allele_sequences_from_vcf_record, get_locus_id, parse_info, parse_locus_interval, process_vcfs)

# TRGT: FORMAT is GT:AL:ALLR:SD:MC:MS:AP:AM, and the record carries a leading padding base, so POS is one before the
# repeat start given by TRID. Real example (HG002 30x PacBio):
#   chr1  591733  .  CAAAAAAAAAAAAAAAAAA  CAAA...(20 A) ,CAAA...(21 A)  TRID=1-591733-591751-A;END=591751  1/2:20,21:...
# Note there is no INFO/START: the repeat interval comes from TRID.
_TRGT_FORMAT = "GT:AL:ALLR:SD:MC:MS:AP:AM"
_TRGT_PAD = ":.:.:.:.:.:."  # the 6 trailing FORMAT fields after GT:AL

# ATaRVa: FORMAT is GT:AL:CN:LPM:AR:SD:DP:SN:SQ:MA:MR:DS:MV, and REF spans exactly [INFO START, INFO END) with no
# padding base (atarva/vcf_writer.py writes POS = locus_start + 1 and REF = ref.fetch(contig, locus_start, locus_end)).
_ATARVA_FORMAT = "GT:AL:CN:LPM:AR:SD:DP:SN:SQ:MA:MR:DS:MV"
_ATARVA_PAD = ":.:.:.:.:.:.:.:.:.:."  # the 10 trailing FORMAT fields after GT:AL:CN

# HipSTR: FORMAT starts GT:GB:Q:...; the locus id is in the ID column and INFO carries START (1-based) / END.
_HIPSTR_FORMAT = "GT:GB:Q:DP"


def trgt_record(trid, pos_1based, end, ref, alt, genotype, allele_sizes):
    return "\t".join([trid.split("-")[0], str(pos_1based), ".", ref, alt, ".", ".",
                      f"TRID={trid};END={end};MOTIFS=A;STRUC=(A)n",
                      _TRGT_FORMAT, f"{genotype}:{allele_sizes}" + _TRGT_PAD])


def atarva_record(chrom, start_0based, end, motif, ref, alt, genotype, filter_field="PASS", sample=None):
    return "\t".join([chrom, str(start_0based + 1), ".", ref, alt, "0", filter_field,
                      f"AC=1;AN=2;MOTIF={motif};START={start_0based};END={end};ID=.;REFCN=1",
                      _ATARVA_FORMAT, sample if sample is not None else f"{genotype}:.:." + _ATARVA_PAD])


def hipstr_record(locus_id, pos_1based, start_1based, end, ref, alt, genotype):
    return "\t".join([locus_id.split("-")[0], str(pos_1based), locus_id, ref, alt, ".", ".",
                      f"START={start_1based};END={end};PERIOD=4;DP=30",
                      _HIPSTR_FORMAT, f"{genotype}:0|0:1.00:30"])


# LongTR: a HipSTR fork, so it also echoes the catalog locus id in the ID column and can extend the genotyped region
# beyond the requested bed interval. Real example (HG002 30x PacBio):
#   chr1  591734  chr1-591733-591751-A  A(x18)  A(x20),A(x21)  START=591734;END=591751;MOTIF=A;PERIOD=1;BPDIFFS=2,3
_LONGTR_FORMAT = "GT:GB:Q:PQ:DP"


def longtr_record(locus_id, pos_1based, start_1based, end, ref, alt, genotype):
    return "\t".join([locus_id.split("-")[0], str(pos_1based), locus_id, ref, alt, ".", ".",
                      f"START={start_1based};END={end};MOTIF=A;PERIOD=1;DP=30",
                      _LONGTR_FORMAT, f"{genotype}:0|0:1.00:1.00:30"])


class Tests(unittest.TestCase):

    def test_parse_locus_interval(self):
        self.assertEqual(parse_locus_interval("1-591733-591751-A"), (591733, 591751))
        self.assertEqual(parse_locus_interval("chrX-100-140-AT"), (100, 140))
        # TRGT's own docstring example uses "_" separators
        self.assertEqual(parse_locus_interval("chr1_674823_674832"), (674823, 674832))
        # a locus id that isn't in "{chrom}-{start}-{end}-{motif}" form can't yield an interval
        self.assertIsNone(parse_locus_interval("some_custom_locus_name"))
        self.assertIsNone(parse_locus_interval("chr1-100"))

    def test_get_locus_id(self):
        info = parse_info("TRID=1-591733-591751-A;END=591751;MOTIFS=A")
        self.assertEqual(get_locus_id("TRGTv5", "chr1", ".", info), "1-591733-591751-A")
        self.assertEqual(get_locus_id("HipSTR", "chr1", "chr1-789568-789588-ATGGA", {"START": "789569"}),
                         "chr1-789568-789588-ATGGA")
        self.assertEqual(
            get_locus_id("ATaRVa", "chr1", ".", parse_info("MOTIF=CTT;START=15796;END=15849")),
            "chr1-15796-15849-CTT")
        # LongTR echoes the locus id from its regions bed into the ID column, the same way HipSTR does
        self.assertEqual(get_locus_id("LongTR", "chr1", "chr1-591733-591751-A", {"START": "591734"}),
                         "chr1-591733-591751-A")

    def test_longtr_alleles(self):
        # the first record of HG002 30x PacBio's LongTR vcf: REF spans the locus interval exactly, nothing to trim
        row = extract_allele_sequences_from_vcf_record(
            longtr_record("chr1-591733-591751-A", 591734, 591734, 591751, "A" * 18, f"{'A' * 20},{'A' * 21}", "1|2"),
            "LongTR")
        self.assertEqual(row["LocusId"], "1-591733-591751-A")
        self.assertEqual(row["AlleleSequence: Allele 1"], "A" * 20)
        self.assertEqual(row["AlleleSequence: Allele 2"], "A" * 21)
        self.assertEqual(row["AlleleStatus: Allele 1"], CALLED)
        self.assertEqual(row["ExtractionStatus"], EXTRACTION_OK)

    def test_longtr_extended_genotyping_region_is_trimmed(self):
        # LongTR inherits HipSTR's habit of extending the genotyped region past the requested bed interval (21 of the
        # 239 records in that same file do), so the flanks shared by REF and every ALT are trimmed to the locus.
        # INFO START is the start of the repeat, so it equals the locus id's start (100 0-based == 101 1-based) and
        # POS is what carries the 5bp extension. This extractor takes the interval from the locus id and never reads
        # INFO START, but convert_hipstr_vcf_to_expansion_hunter_json.py trims off INFO START/END, so a fixture that
        # disagreed with itself would stop the two implementations from being checkable against the same record.
        row = extract_allele_sequences_from_vcf_record(
            longtr_record("chr1-100-120-A", 96, 101, 120, "T" * 5 + "A" * 20 + "T" * 5, "T" * 5 + "A" * 18 + "T" * 5,
                          "0|1"),
            "LongTR")
        self.assertEqual(row["AlleleSequence: Allele 1"], "A" * 18)
        self.assertEqual(row["AlleleSequence: Allele 2"], "A" * 20)
        self.assertEqual(row["ExtractionStatus"], EXTRACTION_OK)

    def test_trgt_padding_base_is_trimmed(self):
        # TRGT emits a leading padding base: POS (591733) is one before the repeat start in TRID (591733 0-based, i.e.
        # 591734 1-based), so REF spans 1 extra base on the left and every allele must be trimmed by 1.
        # This also confirms the TRGT branch never needs an INFO/START field, which TRGT doesn't write.
        row = extract_allele_sequences_from_vcf_record(
            trgt_record("1-591733-591751-A", 591733, 591751, "C" + "A" * 18, "C" + "A" * 20 + "," + "C" + "A" * 21,
                        "1/2", "20,21"), "TRGTv5")
        self.assertEqual(row["LocusId"], "1-591733-591751-A")
        self.assertEqual(row["AlleleSequence: Allele 1"], "A" * 20)
        self.assertEqual(row["AlleleSequence: Allele 2"], "A" * 21)
        self.assertEqual(row["AlleleStatus: Allele 1"], CALLED)
        self.assertEqual(row["AlleleStatus: Allele 2"], CALLED)
        self.assertEqual(row["ExtractionStatus"], EXTRACTION_OK)

    def test_trgt_hom_and_size_sorting(self):
        # hom-alt (GT 1/1) -> both alleles are the single ALT
        row = extract_allele_sequences_from_vcf_record(
            trgt_record("1-651960-651977-A", 651960, 651977, "C" + "A" * 17, "C" + "A" * 16, "1/1", "16,16"), "TRGTv5")
        self.assertEqual(row["AlleleSequence: Allele 1"], "A" * 16)
        self.assertEqual(row["AlleleSequence: Allele 2"], "A" * 16)

        # hom-ref (GT 0/0, ALT ".") -> both alleles are the trimmed REF
        row = extract_allele_sequences_from_vcf_record(
            trgt_record("1-651960-651977-A", 651960, 651977, "C" + "A" * 17, ".", "0/0", "17,17"), "TRGTv5")
        self.assertEqual(row["AlleleSequence: Allele 1"], "A" * 17)
        self.assertEqual(row["AlleleSequence: Allele 2"], "A" * 17)
        self.assertEqual(row["AlleleStatus: Allele 1"], CALLED)

        # GT lists the long allele first; the output is size-sorted so Allele 1 is always the shorter one
        row = extract_allele_sequences_from_vcf_record(
            trgt_record("1-591733-591751-A", 591733, 591751, "C" + "A" * 18, "C" + "A" * 21 + "," + "C" + "A" * 20,
                        "1/2", "21,20"), "TRGTv5")
        self.assertEqual(row["AlleleSequence: Allele 1"], "A" * 20)
        self.assertEqual(row["AlleleSequence: Allele 2"], "A" * 21)

    def test_trgt_no_call(self):
        row = extract_allele_sequences_from_vcf_record(
            trgt_record("1-1030767-1030787-TG", 1030767, 1030787, "A" + "TG" * 10, ".", ".", "."), "TRGTv5")
        self.assertEqual(row["AlleleStatus: Allele 1"], NO_CALL)
        self.assertEqual(row["AlleleStatus: Allele 2"], NO_CALL)
        self.assertEqual(row["ExtractionStatus"], EXTRACTION_OK)

    def test_atarva_del_allele_differs_from_no_call(self):
        # <DEL> is a real zero-length allele: the sequence is empty but the status is "called", so downstream code
        # scores it against the truth sequence rather than dropping it.
        del_row = extract_allele_sequences_from_vcf_record(
            atarva_record("chr1", 9683125, 9683142, "A", "A" * 15 + "CC", "<DEL>", "0|1"), "ATaRVa")
        self.assertEqual(del_row["AlleleSequence: Allele 1"], "")
        self.assertEqual(del_row["AlleleSequence: Allele 2"], "A" * 15 + "CC")
        self.assertEqual(del_row["AlleleStatus: Allele 1"], CALLED)
        self.assertEqual(del_row["AlleleStatus: Allele 2"], CALLED)
        self.assertEqual(del_row["ExtractionStatus"], EXTRACTION_OK)

        # FILTER=LESS_READS: ATaRVa could not genotype the locus, so the empty sequence is NOT an allele
        no_call_row = extract_allele_sequences_from_vcf_record(
            atarva_record("chr1", 13089246, 13089273, "T", "T" * 27, ".", ".", filter_field="LESS_READS",
                          sample="." + ":." * 12), "ATaRVa")
        self.assertEqual(no_call_row["AlleleSequence: Allele 1"], "")
        self.assertEqual(no_call_row["AlleleStatus: Allele 1"], NO_CALL)
        self.assertEqual(no_call_row["AlleleStatus: Allele 2"], NO_CALL)

        # the two must be distinguishable, otherwise a failed locus would be scored as a zero-length allele
        self.assertNotEqual(
            [del_row[c] for c in OUTPUT_COLUMNS[1:]], [no_call_row[c] for c in OUTPUT_COLUMNS[1:]])

    def test_atarva_locus_id_chr_prefix_is_stripped(self):
        row = extract_allele_sequences_from_vcf_record(
            atarva_record("chr1", 15796, 15802, "CTT", "CTTCTT", "CTT", "0|1"), "ATaRVa")
        self.assertEqual(row["LocusId"], "1-15796-15802-CTT")

    def test_hipstr_padding_is_trimmed(self):
        # HipSTR can extend the genotyped region beyond the requested bed interval, so POS < INFO/START and the
        # record's REF runs past INFO/END. Both flanks must be trimmed off REF and every ALT.
        row = extract_allele_sequences_from_vcf_record(
            hipstr_record("chr1-1000-1016-AAAC", 999, 1001, 1016, "GG" + "AAAC" * 4 + "TT", "GG" + "AAAC" * 3 + "TT",
                          "0|1"), "HipSTR")
        self.assertEqual(row["LocusId"], "1-1000-1016-AAAC")
        self.assertEqual(row["AlleleSequence: Allele 1"], "AAAC" * 3)
        self.assertEqual(row["AlleleSequence: Allele 2"], "AAAC" * 4)

        # an unpadded record (POS - 1 == INFO/START - 1, REF ends at INFO/END) is passed through untrimmed
        row = extract_allele_sequences_from_vcf_record(
            hipstr_record("chr1-789568-789588-ATGGA", 789569, 789569, 789588, "ATGGA" * 4, "ATGGA" * 3, "0|1"),
            "HipSTR")
        self.assertEqual(row["AlleleSequence: Allele 1"], "ATGGA" * 3)
        self.assertEqual(row["AlleleSequence: Allele 2"], "ATGGA" * 4)

    def test_haploid_genotype_is_duplicated(self):
        # male chrX/chrY outside the PAR: ATaRVa writes a single GT index. The truth set stores HEMI loci as a
        # duplicated diploid genotype and add_tool_results_columns.py mirrors that, so do the same here.
        row = extract_allele_sequences_from_vcf_record(
            atarva_record("chrY", 2804694, 2804708, "T", "T" * 14, "T" * 13, "1", sample="1:13:13" + _ATARVA_PAD),
            "ATaRVa")
        self.assertEqual(row["AlleleSequence: Allele 1"], "T" * 13)
        self.assertEqual(row["AlleleSequence: Allele 2"], "T" * 13)
        self.assertEqual(row["AlleleStatus: Allele 1"], CALLED)
        self.assertEqual(row["AlleleStatus: Allele 2"], CALLED)

    def test_multiallelic(self):
        # ATaRVa's "multizygous" loci report more than 2 alleles, which can't be expressed as a diploid genotype
        row = extract_allele_sequences_from_vcf_record(
            atarva_record("chr3", 498, 540, "CAG", "CAG" * 14, "CAG" * 10 + "," + "CAG" * 12 + "," + "CAG" * 16,
                          "1/2/3", sample="1/2/3:.:." + _ATARVA_PAD), "ATaRVa")
        self.assertEqual(row["ExtractionStatus"], EXTRACTION_MULTIALLELIC)
        self.assertEqual(row["AlleleStatus: Allele 1"], NO_CALL)
        self.assertEqual(row["AlleleStatus: Allele 2"], NO_CALL)

        # more than 2 ALT alleles is also rejected, even with a diploid GT
        row = extract_allele_sequences_from_vcf_record(
            atarva_record("chr3", 498, 540, "CAG", "CAG" * 14, "CAG" * 10 + "," + "CAG" * 12 + "," + "CAG" * 16,
                          "1/2"), "ATaRVa")
        self.assertEqual(row["ExtractionStatus"], EXTRACTION_MULTIALLELIC)

    def test_ref_mismatch(self):
        # a REF that doesn't match the reference genome means the record and the truth locus don't describe the same
        # interval, so scoring its alleles would silently inflate every edit distance
        def get_reference_sequence(chrom, start_0based, end):
            return "A" * (end - start_0based)

        row = extract_allele_sequences_from_vcf_record(
            atarva_record("chr1", 100, 110, "AT", "GGGGGGGGGG", "GGGGGGGGGGGG", "0|1"), "ATaRVa",
            get_reference_sequence)
        self.assertEqual(row["ExtractionStatus"], EXTRACTION_REF_MISMATCH)
        self.assertEqual(row["AlleleStatus: Allele 1"], NO_CALL)

        # the same record passes when REF does match
        row = extract_allele_sequences_from_vcf_record(
            atarva_record("chr1", 100, 110, "A", "A" * 10, "A" * 12, "0|1"), "ATaRVa", get_reference_sequence)
        self.assertEqual(row["ExtractionStatus"], EXTRACTION_OK)

        # hg38 soft-masks tandem repeats in lower case, so the comparison must be case-insensitive
        row = extract_allele_sequences_from_vcf_record(
            atarva_record("chr1", 100, 110, "A", "a" * 10, "a" * 12, "0|1"), "ATaRVa",
            lambda chrom, start_0based, end: "A" * (end - start_0based))
        self.assertEqual(row["ExtractionStatus"], EXTRACTION_OK)
        self.assertEqual(row["AlleleSequence: Allele 2"], "A" * 12)

        # a record whose REF is narrower than the locus interval can't be reconciled by trimming either
        row = extract_allele_sequences_from_vcf_record(
            atarva_record("chr1", 100, 130, "A", "A" * 10, "A" * 12, "0|1"), "ATaRVa")
        self.assertEqual(row["ExtractionStatus"], EXTRACTION_REF_MISMATCH)

    def test_malformed_records_return_none(self):
        # a locus id that isn't "{chrom}-{start}-{end}-{motif}" gives no interval to trim to
        self.assertIsNone(extract_allele_sequences_from_vcf_record(
            trgt_record("custom_locus_name", 1000, 1010, "A" * 11, ".", "0/0", "10,10"), "TRGTv5"))

        # a GT index with no corresponding allele
        self.assertIsNone(extract_allele_sequences_from_vcf_record(
            atarva_record("chr1", 100, 110, "A", "A" * 10, "A" * 12, "1/2"), "ATaRVa"))

        # a truncated record with no sample column
        self.assertIsNone(extract_allele_sequences_from_vcf_record("chr1\t100\t.\tA\tT\t.\t.\t.", "ATaRVa"))

    def test_process_vcfs_writes_expected_table(self):
        vcf_lines = "\n".join([
            "##fileformat=VCFv4.2",
            "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "HG002"]),
            atarva_record("chr1", 100, 110, "A", "A" * 10, "A" * 12, "0|1"),
            atarva_record("chr1", 200, 210, "A", "A" * 10, "<DEL>", "0|1"),
            atarva_record("chr2", 300, 310, "A", "A" * 10, ".", ".", filter_field="LESS_READS",
                          sample="." + ":." * 12),
        ]) + "\n"

        with tempfile.TemporaryDirectory() as temp_dir:
            vcf_path = os.path.join(temp_dir, "test.vcf")
            with open(vcf_path, "wt") as f:
                f.write(vcf_lines)
            output_path = os.path.join(temp_dir, "test.allele_sequences.tsv.gz")
            counters = process_vcfs([vcf_path], "ATaRVa", output_path)

            with gzip.open(output_path, "rt") as f:
                header = f.readline().rstrip("\n").split("\t")
                rows = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in f]

        self.assertEqual(header, OUTPUT_COLUMNS)
        self.assertEqual(len(rows), 3)
        self.assertEqual([r["LocusId"] for r in rows],
                         ["1-100-110-A", "1-200-210-A", "2-300-310-A"])
        self.assertEqual(rows[0]["AlleleSequence: Allele 1"], "A" * 10)
        self.assertEqual(rows[0]["AlleleSequence: Allele 2"], "A" * 12)
        # the <DEL> allele sorts first (length 0) and stays "called"
        self.assertEqual(rows[1]["AlleleSequence: Allele 1"], "")
        self.assertEqual(rows[1]["AlleleStatus: Allele 1"], CALLED)
        self.assertEqual(rows[2]["AlleleStatus: Allele 1"], NO_CALL)

        self.assertEqual(counters["total"], 3)
        self.assertEqual(counters[EXTRACTION_OK], 3)
        self.assertEqual(counters[CALLED], 4)
        self.assertEqual(counters[NO_CALL], 2)


if __name__ == "__main__":
    unittest.main()
