# str-analysis
This package contains scripts and utilities related to analyzing short tandem repeats (STRs). 

---
## Scripts

### call_rfc1_canvas_alleles

RFC1 STR expansions have recently been linked to [CANVAS](https://www.omim.org/entry/614575) [ [Cortese 2019](https://pubmed.ncbi.nlm.nih.gov/30926972/) ].
This STR locus is unique in that it's both autosomal recessive and has a pathogenic repeat motif (AAGGG) 
that differs from the motif in the reference genome (AAAAG). Several other benign and pathogenic motifs have also been reported, such as  
AAAGG and ACAGG, as well as motifs of uncertain significance such as AAGAG and AGAGG [ [Akcimen 2019](https://pubmed.ncbi.nlm.nih.gov/31824583/) ].
It's not unusual for individuals to have one motif on one chromosome, and another motif on the other chromosome. 
Due to this multi-allelic nature, current STR genotyping tools like [ExpansionHunter](https://github.com/Illumina/ExpansionHunter) struggle to accurately 
genotype RFC1.  

A more recent tool, [ExpansionHunter Denovo](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02017-z), 
is better at detecting multi-allelic expanded RFC1 motif(s), but is unable to distinguish unaffected carriers 
(such as those that have one benign AAAAG reference allele, and one pathogenic AAGGG allele) from 
affected individuals (homozygous for the AAGGG allele). In both cases, ExpansionHunter Denovo might detect the 
AAGGG allele, but would fail to detect the benign allele unless it is also highly expanded. 

This script addresses some of the above limitations. It takes a whole genome (WGS) bam or cram file and 
outputs a .json file with information about the repeat motifs it detects at the RFC1/CANVAS STR locus. 
The script doesn't attempt to estimate the repeat size, but simply outputs a **call** field that describes 
the type of repeat motif(s) it detected (more details below). 


In initial tests, we found this approach is sufficient to distinguish affected from unaffected individuals. 
In a cohort of 4447 samples from individuals with different rare disease phenotypes as well as their unaffected
family members, this script identified 9 individuals as `call = PATHOGENIC MOTIF / PATHOGENIC MOTIF`
meaning that their read data supports only the pathogenic motifs. Of these 9 individuals, 2 are positive controls with previously-validated RFC1/CANVAS
pathogenic expansions, and 1 is an affected individual with a phenotype that is highly consistent with CANVAS (now 
proceeding to clinical validation). The other 6 cases are likely false-positives (or secondary findings).
In the future, as we continue to test the script on additional positive controls as well as simulated data, we may 
adjust thresholds to optimize sensitivity/specificity.

Description of all fields in the output `*.rfc1_canvas_alleles.json`:

**sample_id**: *If this value is not specified as a command line arg, it is parsed from the input bam/cram file header or filename prefix.*    
**call**: *describes the alleles detected at the RFC1/CANVAS locus. Its format is analogous to a VCF genotype. Possible values are:*
* `PATHOGENIC MOTIF / PATHOGENIC MOTIF`: *only pathogenic allele(s) detected*
* `BENIGN MOTIF / BENIGN MOTIF`: *only benign allele(s) detected*
* `MOTIF OF UNCERTAIN SIGNIFICANCE / MOTIF OF UNCERTAIN SIGNIFICANCE`: *non-canonical allele(s) detected with unknown pathogenicity*
* `BENIGN MOTIF / PATHOGENIC MOTIF`: *heterozygous for a benign allele and a pathogenic allele, implying carrier status*
* `PATHOGENIC MOTIF / MOTIF OF UNCERTAIN SIGNIFICANCE`: *heterozygous for a pathogenic allele and a non-canonical allele(s) detected with unknown pathogenicity*
* `BENIGN MOTIF / MOTIF OF UNCERTAIN SIGNIFICANCE`: *heterozygous for a benign allele and a non-canonical allele(s) detected with unknown pathogenicity*
* `NO CALL`: *not enough evidence in the read data to support any of the above options*

**allele1_repeat_unit**: *the repeat unit that is supported by the most reads.*     
**allele1_read_count**: *the number of reads supporting allele1.*   
**allele1_normalized_read_count**: *same as allele1_read_count, but normalized by depth
of coverage in the flanking regions of the RFC1 locus*   
**allele1_n_occurrences**: *the total number of times allele1 occurs in the reads at the RFC1 locus.*     
**allele1_read_count_with_offtargets**: *the number of reads supporting allele1 within off-target regions 
    for this repeat unit. These are ~1kb regions where fully-repetitive (aka. IRR) reads may mismap to based on 
    experiments with simulated data.*   
**allele1_normalized_read_count_with_offtargets**: *same as allele1_read_count_with_offtargets, but normalized by depth 
    of coverage in the flanking regions of the RFC1 locus*

**allele2_repeat_unit**: *the repeat unit that is supported by the next most reads, or null if all reads support allele1.*    
**allele2_read_count**: *see "allele1_read_count" description.*  
**allele2_n_occurrences**: *see "allele1_n_occurrences" description.*  
...    
*NOTE:* allele2_* fields will only be generated if there is read support for more than 1 allele.    

**left_flank_coverage**: *average read depth within a 2kb window immediately to the left of the RFC1 locus*  
**right_flank_coverage**: *average read depth within a 2kb window immediately to the right of the RFC1 locus*  
  
**found_n_reads_overlap_rfc1_locus**: *number of reads that overlap the AAAAG repeat in the reference genome 
at the RFC1 locus and have a MAPQ > 2*    
**found_repeats_in_n_reads**: *number of reads that overlap the AAAAG repeat in the reference genome at the RFC1 locus, and have both MAPQ > 2 as well as some 5bp or 6bp repeat motif that covers > 70% of the overlapping read sequence (including any soft-clipped bases)*    
**found_repeats_in_fraction_of_reads**: `found_repeats_in_n_reads` / `found_n_reads_overlap_rfc1_locus`  

This script optionally takes an ExpansionHunterDenovo profile and copies relevant info to the output
`*.rfc1_canvas_alleles.json` file, adding more fields. 
The ExpansionHunterDenovo profile isn't used in calculating the fields listed above.  

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

