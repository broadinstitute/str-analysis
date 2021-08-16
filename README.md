# str-analysis
This package contains scripts and utilities related to analyzing short tandem repeats (STRs). 

---
## Scripts

### call_rfc1_canvas_alleles

RFC1 STR expansions have recently been linked to [CANVAS](https://www.omim.org/entry/614575) [ [Cortese 2019](https://pubmed.ncbi.nlm.nih.gov/30926972/) ].
This STR locus is unique in that it's both autosomal recessive and has a pathogenic repeat motif (AAGGG) 
that differs from the motif in the reference genome (AAAAG). Several other benign and pathogenic motifs have also been reported, such as  
AAAGG and ACAGG, as well as motifs of uncertain significances such as AAGAG and AGAGG [ [Akcimen 2019](https://pubmed.ncbi.nlm.nih.gov/31824583/) ]. 
It's not unusual for individuals to have one motif on one chromosome, and another motif on the other chromosome. 
Due to this multi-allelic nature, current STR genotyping tools like [ExpansionHunter](https://github.com/Illumina/ExpansionHunter) struggle to accurately 
genotype RFC1.  

The more recent tool [ExpansionHunter Denovo](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02017-z) 
is better at detecting RFC1 motif(s), but has low sensitivity, and is unable to distinguish unaffected carriers  
(who may have one benign AAAAG reference allele, and one pathogenic AAGGG allele) from 
affected individuals (who are homozygous for the AAGGG allele). In both cases, ExpansionHunter Denovo might detect the 
AAGGG allele, but there's no way to tell from the output whether a benign allele is also present. 
This script addresses some of these limitations. 

This script takes a whole genome (WGS) bam or cram file and outputs a .json file with information about the repeat 
motifs it detected at the RFC1/CANVAS STR locus. Currently, several repeat motifs are known to be benign 
(AAAAG, AAAGG) [ [Cortese 2019](https://pubmed.ncbi.nlm.nih.gov/30926972/) ], and several are pathogenic 
(AAGGG, ACAGG) [ [Scriba 2020](https://pubmed.ncbi.nlm.nih.gov/33103729/) ] when expanded beyond ~400 repeats. 
This script doesn't attempt to estimate the repeat size, but simply outputs a **call** field that indicates 
the kind of repeat motif(s) it detected at the locus (see detailed description below). 

In initial tests, we find that this approach is sufficient to distinguish affected from unaffected individuals. 
In a cohort of 4447 samples from individuals with different rare disease phenotypes as well as their unaffected
family members, this script identified 9 individuals as `call = PATHOGENIC MOTIF / PATHOGENIC MOTIF`
meaning that their read data supports only the pathogenic motifs, and no other motifs (ie. benign motifs, or 
motifs of uncertain significance). Of these 9 individuals, 2 are positive controls with previously-validated RFC1/CANVAS
pathogenic expansions, and 1 is an affected individual with a phenotype that is highly consistent with CANVAS (now 
proceeding to clinical validation). The other 6 cases are likely false-positives (or secondary findings).

Although the script works relatively well with our current data, we may need to adjust the script's internal thresholds 
if/when we test it on more positive controls or on simulated samples with different types of expansions.

These are all the output fields:

**sample_id**: *If this value is not specified as a command line arg, it is parsed from the input bam/cram file header or filename prefix.*    
**call**: *describes the alleles detected at the RFC1/CANVAS locus. Its format is analogous to a VCF genotype. Possible values are:*
* `PATHOGENIC MOTIF / PATHOGENIC MOTIF`: *only pathogenic allele(s) detected*
* `BENIGN MOTIF / BENIGN MOTIF`: *only benign allele(s) detected*
* `MOTIF OF UNCERTAIN SIGNIFICANCE / MOTIF OF UNCERTAIN SIGNIFICANCE`: *non-canonical allele(s) detected with unknown pathogenicity*
* `BENIGN MOTIF / PATHOGENIC MOTIF`: *heterozygous for a benign allele and a pathogenic allele, implying carrier status*
* `PATHOGENIC MOTIF / MOTIF OF UNCERTAIN SIGNIFICANCE`: *heterozygous for a pathogenic allele and a non-canonical allele(s) detected with unknown pathogenicity*
* `BENIGN MOTIF / MOTIF OF UNCERTAIN SIGNIFICANCE`: *heterozygous for a benign allele and a non-canonical allele(s) detected with unknown pathogenicity*
* `null`: *not enough evidence in the read data to support any of the above options*

**allele1_repeat_unit**: *the repeat unit that is supported by the most reads.*  
**allele1_read_count**: *the number of reads supporting allele1.*  
**allele1_n_occurrences**: *the total number of times allele1 occurs in the reads at the RFC1 locus.*  

**allele2_repeat_unit**: *the repeat unit that is supported by the next most reads, or null if all reads support allele1.*  
**allele2_read_count**: *see "allele1_read_count" description.*  
**allele2_n_occurrences**: *see "allele1_n_occurrences" description.*  

**left_flank_coverage**: *average read depth within a 2kb window immediately to the left of the RFC1 locus*  
**right_flank_coverage**: *average read depth within a 2kb window immediately to the right of the RFC1 locus*  
  
**found_n_reads_overlap_rfc1_locus**: *number of reads that overlap the AAAAG repeat in the reference genome 
at the RFC1 locus and have a MAPQ > 2*    
**found_repeats_in_n_reads**: *number of those reads that have a 5bp or 6bp repeat unit that covers > 70% of the overlapping read sequence*    
**found_repeats_in_fraction_of_reads**: `found_repeats_in_n_reads` / `found_n_reads_overlap_rfc1_locus`  

Also, this script optionally takes an ExpansionHunterDenovo profile for this sample and copies relevant fields to the
output. The ExpansionHunterDenovo profile isn't used in calculations. 

Example command line:

```
call_rfc1_canvas_alleles -e sample1.str_profile.json -g 38 sample1.cram
```

### combine_json_to_tsv

This script can combine the `call_rfc1_canvas_alleles` output json files for multiple samples into 
a single .tsv. The script takes the paths of the .json files as input, or, if none are provided, it searches for .json 
files in the current directory and subdirectories.

Example command line:
```
combine_json_to_tsv  sample1.rfc1_canvas_alleles.json  sample2.rfc1_canvas_alleles.json
```

## Installation

To install the scripts or utilities, run:

```
python3 -m pip install --upgrade str_analysis
```

Alternatively, you can use this docker image:

```
docker run -it weisburd/str-analysis:latest
```

