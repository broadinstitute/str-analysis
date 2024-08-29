import simplejson as json
import tempfile
import unittest

from str_analysis.convert_trgt_vcf_to_expansion_hunter_json import process_trgt_vcf

VCF_CONTENTS1 = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=TRID,Number=1,Type=String,Description="Tandem repeat ID">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=MOTIFS,Number=.,Type=String,Description="Motifs that the tandem repeat is composed of">
##INFO=<ID=STRUC,Number=1,Type=String,Description="Structure of the region">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AL,Number=.,Type=Integer,Description="Length of each allele">
##FORMAT=<ID=ALLR,Number=.,Type=String,Description="Length range per allele">
##FORMAT=<ID=SD,Number=.,Type=Integer,Description="Number of spanning reads supporting per allele">
##FORMAT=<ID=MC,Number=.,Type=String,Description="Motif counts per allele">
##FORMAT=<ID=MS,Number=.,Type=String,Description="Motif spans per allele">
##FORMAT=<ID=AP,Number=.,Type=Float,Description="Allele purity per allele">
##FORMAT=<ID=AM,Number=.,Type=Float,Description="Mean methylation level per allele">
##contig=<ID=chr1,length=248956422>
##trgtVersion=1.1.1-62f1f0e
##trgtCommand=trgt genotype --genome /io/local_copy/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta --reads /io/local_copy/fc-secure-ab18917b-873b-4082-8181-7653393999c9/GSS237629/LR_GS/pacbio/GSS237629.GCA_000001405_15.haplotagged.bam --repeats /io/local_copy/bw-proj/gnomad-bw/p1_iter4/2024_08_24__EH_and_TRGT_genome_wide_analysis_iter2/long_read_genome_wide_catalog.trgt.bed --karyotype XX --output-prefix GSS237629 --threads 16
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
chr1	71175	.	GAAAAAAAAAAA	GAAAAAAAAAAAAAA,GAAAAAAAAAAAAAAA	0	.	TRID=1-71175-71186-A;END=71186;MOTIFS=A;STRUC=(A)n	GT:AL:ALLR:SD:MC:MS:AP:AM	1/2:14,15:14-14,15-15:1,1:14,15:0(0-14),0(0-15):1.000000,1.000000:.,.
chr1	493790	.	TGGCAGGCAGGCAG	.	0	.	TRID=1-493790-493803-GGCA;END=493803;MOTIFS=GGCA;STRUC=(GGCA)n	GT:AL:ALLR:SD:MC:MS:AP:AM	.:.:.:.:.:.:.:.
chr1	624509	.	AAGGGAGGGAGGGAGGGAGGGAGGGA	AAGGAAGGAAGGAAGGGAGGGAGGGAGGGAGGGA	0	.	TRID=1-624509-624534-AGGG;END=624534;MOTIFS=AGGG;STRUC=(AGGG)n	GT:AL:ALLR:SD:MC:MS:AP:AM	0/1:25,33:25-25,33-33:2,6:6,8:0(0-25),0(0-33):0.960000,0.878788:.,.
chr1	29135	.	AGCGGCGGCG	AGCCGCGGCG	0	.	TRID=1-29135-29144-GCG;END=29144;MOTIFS=GCG;STRUC=(GCG)n	GT:AL:ALLR:SD:MC:MS:AP:AM	0/1:9,9:9-10,9-9:5,5:3,3:0(0-9),0(0-9):1.000000,0.888889:0.16,0.16
chr1	199861	.	CCCGCCGCCGC	.	0	.	TRID=1-199861-199871-CCG;END=199871;MOTIFS=CCG;STRUC=(CCG)n	GT:AL:ALLR:SD:MC:MS:AP:AM	0/0:10,10:10-10,10-10:2,6:4,4:0(0-10),0(0-10):0.833333,0.833333:0.37,0.37
chr1	788757	.	TAATGGAATGGAATGGAATGGAATGGAATGGAATGGAATGGAATGGAATGGA	TAATGGAATGGAATGGAATGGAATGGAATGGAATGGAATGGAATGGAATGGAATGGAATGGA	0	.	TRID=1-788757-788808-AATGG;END=788808;MOTIFS=AATGG;STRUC=(AATGG)n	GT:AL:ALLR:SD:MC:MS:AP:AM	0/1:51,61:51-61,61-61:3,2:10,12:0(0-51),0(0-61):0.980392,0.983607:.,.
chr1	670451	.	GTTATTTTATT	GTTATTTTATTTTATT	0	.	TRID=1-670451-670461-TTATT;END=670461;MOTIFS=TTATT;STRUC=(TTATT)n	GT:AL:ALLR:SD:MC:MS:AP:AM	0/1:10,15:10-10,10-15:10,3:2,3:0(0-10),0(0-15):1.000000,1.000000:.,."""

class Tests(unittest.TestCase):

    def setUp(self):
        # write vcf contents to temp file
        self.test_vcf1 = tempfile.NamedTemporaryFile(mode="wt", delete=False)
        self.test_vcf1.write(VCF_CONTENTS1)
        self.test_vcf1.flush()
        self.test_vcf1.seek(0)

    def test_parse_read_count_tuples(self):
        results = process_trgt_vcf(self.test_vcf1.name)
        self.assertEqual(len(results["LocusResults"]), 5)
        self.assertEqual(len(results["LocusResults"]["1-71175-71186-A"]["Variants"]), 1)
        self.assertEqual(len(results["LocusResults"]["1-624509-624534-AGGG"]["Variants"]), 1)
        self.assertEqual(len(results["LocusResults"]["1-29135-29144-GCG"]["Variants"]), 1)
        self.assertEqual(len(results["LocusResults"]["1-788757-788808-AATGG"]["Variants"]), 1)
        self.assertEqual(len(results["LocusResults"]["1-670451-670461-TTATT"]["Variants"]), 1)

        for locus_id, variant_info in results["LocusResults"].items():
            for variant_id, v in variant_info["Variants"].items():
                alleles = [v["Ref"]] + v["Alt"].split(",")
                genotype = v["GT"].split("/")
                motif = v["RepeatUnit"]
                purity = v["AP"]
                motif = v["RepeatUnit"]
                ref = v["Ref"]
                alelle1 = alleles[int(genotype[0])]
                alelle2 = alleles[int(genotype[1])]
                num_repeats = v["Genotype"].split("/")
                num_repeats1 = int(num_repeats[0])
                num_repeats2 = int(num_repeats[1])
                purity_scores = v["AP"].split(",")
                allele1_purity = float(purity_scores[0])
                allele2_purity = float(purity_scores[1])
                if allele1_purity > 0.99:
                    #print(alelle1[1:], motif*num_repeats1)
                    self.assertEqual(alelle1[1:], motif*num_repeats1)
                if allele2_purity > 0.99:
                    #print(alelle2[1:], motif*num_repeats2)
                    self.assertEqual(alelle2[1:], motif*num_repeats2)
