import unittest

from str_analysis.extract_allele_sequences_from_vcf import CALLED, NO_CALL, EXTRACTION_OK, EXTRACTION_MULTIALLELIC
from str_analysis.extract_allele_sequences_from_expansion_hunter_json import extract_allele_sequences_from_locus_json


def locus(locus_id, variants):
    return {"LocusId": locus_id, "Variants": variants}


class Tests(unittest.TestCase):

    def test_het_call_with_consensus(self):
        row = extract_allele_sequences_from_locus_json(locus("chr1-591733-591751-A", {
            "1-591733-591751-A": {"Genotype": "18/20", "ConsensusSequences": ["A" * 18, "A" * 20]},
        }))
        self.assertEqual(row["LocusId"], "1-591733-591751-A")
        self.assertEqual(row["AlleleSequence: Allele 1"], "A" * 18)
        self.assertEqual(row["AlleleSequence: Allele 2"], "A" * 20)
        self.assertEqual(row["AlleleStatus: Allele 1"], CALLED)
        self.assertEqual(row["AlleleStatus: Allele 2"], CALLED)
        self.assertEqual(row["ExtractionStatus"], EXTRACTION_OK)

    def test_defensively_resorts_short_first(self):
        # ConsensusSequences reversed relative to Genotype order shouldn't happen in practice, but the extractor
        # re-sorts by length rather than trusting the input order
        row = extract_allele_sequences_from_locus_json(locus("chr1-1-10-A", {
            "1-1-10-A": {"Genotype": "5/8", "ConsensusSequences": ["A" * 8, "A" * 5]},
        }))
        self.assertEqual(row["AlleleSequence: Allele 1"], "A" * 5)
        self.assertEqual(row["AlleleSequence: Allele 2"], "A" * 8)

    def test_haploid_call_is_duplicated(self):
        # male chrX/chrY outside the PAR: a single Genotype value, no "/"
        row = extract_allele_sequences_from_locus_json(locus("chrX-1-10-A", {
            "1-1-10-A": {"Genotype": "9", "ConsensusSequences": ["A" * 9]},
        }))
        self.assertEqual(row["AlleleSequence: Allele 1"], "A" * 9)
        self.assertEqual(row["AlleleSequence: Allele 2"], "A" * 9)
        self.assertEqual(row["AlleleStatus: Allele 1"], CALLED)
        self.assertEqual(row["AlleleStatus: Allele 2"], CALLED)

    def test_homozygous_call_with_one_consensus_is_duplicated(self):
        # A homozygous diploid call reports "20/20" but only ONE consensus sequence for the two identical
        # alleles (both the full genotyper and the fast path do this). Requiring one sequence per genotype field
        # would drop every homozygous locus -- 6,417 of HG002's 60,394 consensus-bearing loci at 31x.
        row = extract_allele_sequences_from_locus_json(locus("chr1-711957-712094-TATATAT", {
            "1-711957-712094-TATATAT": {"Genotype": "20/20", "ConsensusSequences": ["TATATAT" * 20]},
        }))
        self.assertEqual(row["AlleleSequence: Allele 1"], "TATATAT" * 20)
        self.assertEqual(row["AlleleSequence: Allele 2"], "TATATAT" * 20)
        self.assertEqual(row["AlleleStatus: Allele 1"], CALLED)
        self.assertEqual(row["AlleleStatus: Allele 2"], CALLED)
        self.assertEqual(row["ExtractionStatus"], EXTRACTION_OK)

    def test_heterozygous_call_with_one_consensus_is_no_call(self):
        # a single consensus is only meaningful when the two alleles are identical; for a het call it's an
        # unexpected shape, so score nothing rather than duplicating one allele's sequence onto the other
        row = extract_allele_sequences_from_locus_json(locus("chr1-1-10-A", {
            "1-1-10-A": {"Genotype": "5/8", "ConsensusSequences": ["A" * 5]},
        }))
        self.assertEqual(row["AlleleStatus: Allele 1"], NO_CALL)
        self.assertEqual(row["AlleleStatus: Allele 2"], NO_CALL)

    def test_missing_consensus_sequences_is_no_call(self):
        # Genotype present (ExpansionHunter did call it) but no ConsensusSequences (e.g. consensus sequences were
        # disabled for this run, or read support was too poor to build one)
        row = extract_allele_sequences_from_locus_json(locus("chr1-1-10-A", {
            "1-1-10-A": {"Genotype": "5/8"},
        }))
        self.assertEqual(row["AlleleStatus: Allele 1"], NO_CALL)
        self.assertEqual(row["AlleleStatus: Allele 2"], NO_CALL)
        self.assertEqual(row["AlleleSequence: Allele 1"], "")
        self.assertEqual(row["AlleleSequence: Allele 2"], "")
        self.assertEqual(row["ExtractionStatus"], EXTRACTION_OK)

    def test_multiple_repeat_variants_is_multiallelic(self):
        row = extract_allele_sequences_from_locus_json(locus("chr1-1-10-A", {
            "1-1-10-A": {"Genotype": "5/8", "ConsensusSequences": ["A" * 5, "A" * 8]},
            "1-11-20-A": {"Genotype": "3/4", "ConsensusSequences": ["A" * 3, "A" * 4]},
        }))
        self.assertEqual(row["ExtractionStatus"], EXTRACTION_MULTIALLELIC)
        self.assertEqual(row["AlleleStatus: Allele 1"], NO_CALL)

    def test_locus_with_no_genotype_returns_none(self):
        # ExpansionHunter never emitted a call for this locus at all
        self.assertIsNone(extract_allele_sequences_from_locus_json(locus("chr1-1-10-A", {
            "1-1-10-A": {"RepeatUnit": "A"},
        })))

    def test_locus_with_no_locus_id_returns_none(self):
        self.assertIsNone(extract_allele_sequences_from_locus_json({"Variants": {}}))

    def test_more_than_two_genotype_alleles_is_multiallelic(self):
        row = extract_allele_sequences_from_locus_json(locus("chr1-1-10-A", {
            "1-1-10-A": {"Genotype": "5/8/9", "ConsensusSequences": ["A" * 5, "A" * 8, "A" * 9]},
        }))
        self.assertEqual(row["ExtractionStatus"], EXTRACTION_MULTIALLELIC)


if __name__ == "__main__":
    unittest.main()
