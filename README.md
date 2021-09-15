# str-analysis
This package contains scripts and utilities related to analyzing short tandem repeats (STRs). 

---
## Scripts

### call_non_ref_pathogenic_motifs

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
meaning that their read data supports only the pathogenic motifs. Of these 9 individuals, 2 are positive controls 
with previously-validated RFC1/CANVAS pathogenic expansions, and 1 is an affected individual with a phenotype that 
is highly consistent with CANVAS (now proceeding to clinical validation). The other 6 cases are likely false-positives 
(or secondary findings). In the future, as we continue to test the script on additional positive controls as well 
as simulated data, we may adjust thresholds to optimize sensitivity/specificity.

Description of the output `*.RFC1_motifs.json`:

**sample_id**: *If this value is not specified as a command line arg, it is parsed from the input bam/cram file header or filename prefix.*    
**call**: *describes the motifs detected at the RFC1/CANVAS locus. Its format is analogous to a VCF genotype. Possible values are:*
* `PATHOGENIC MOTIF / PATHOGENIC MOTIF`: *only pathogenic motif(s) detected*
* `BENIGN MOTIF / BENIGN MOTIF`: *only benign motif(s) detected*
* `MOTIF OF UNCERTAIN SIGNIFICANCE / MOTIF OF UNCERTAIN SIGNIFICANCE`: *non-canonical motif(s) detected with unknown pathogenicity*
* `BENIGN MOTIF / PATHOGENIC MOTIF`: *heterozygous for a benign motif and a pathogenic motif, implying carrier status*
* `PATHOGENIC MOTIF / MOTIF OF UNCERTAIN SIGNIFICANCE`: *heterozygous for a pathogenic motif and a non-canonical motif(s) detected with unknown pathogenicity*
* `BENIGN MOTIF / MOTIF OF UNCERTAIN SIGNIFICANCE`: *heterozygous for a benign motif and a non-canonical motif(s) detected with unknown pathogenicity*
* `NO CALL`: *not enough evidence in the read data to support any of the above options*

**motif1_repeat_unit**: *the repeat unit that is supported by the most reads.*     
**motif1_read_count**: *the number of reads supporting motif1.*   
**motif1_normalized_read_count**: *same as motif1_read_count, but normalized by depth
of coverage in the flanking regions of the RFC1 locus*   
**motif1_n_occurrences**: *the total number of times motif1 occurs in the reads at the RFC1 locus.*     
**motif1_read_count_with_offtargets**: *the number of reads supporting motif1 within off-target regions 
    for this repeat unit. These are ~1kb regions where fully-repetitive (aka. IRR) reads may mismap to based on 
    experiments with simulated data.*   
**motif1_normalized_read_count_with_offtargets**: *same as motif1_read_count_with_offtargets, but normalized by depth 
    of coverage in the flanking regions of the RFC1 locus*

**motif2_repeat_unit**: *the repeat unit that is supported by the next most reads, or null if all reads support motif1.*    
**motif2_read_count**: *see "motif1_read_count" description.*  
**motif2_n_occurrences**: *see "motif1_n_occurrences" description.*  
...    
*NOTE:* motif2_* fields will only be generated if there is read support for more than 1 motif.    

**left_flank_coverage**: *average read depth within a 2kb window immediately to the left of the RFC1 locus*  
**right_flank_coverage**: *average read depth within a 2kb window immediately to the right of the RFC1 locus*  
  
**found_n_reads_overlap_rfc1_locus**: *number of reads that overlap the AAAAG repeat in the reference genome 
at the RFC1 locus and have a MAPQ > 2*    
**found_repeats_in_n_reads**: *number of reads that overlap the AAAAG repeat in the reference genome at the RFC1 locus, and have both MAPQ > 2 as well as some 5bp or 6bp repeat motif that covers > 70% of the overlapping read sequence (including any soft-clipped bases)*    
**found_repeats_in_fraction_of_reads**: `found_repeats_in_n_reads` / `found_n_reads_overlap_rfc1_locus`  

This script optionally runs ExpansionHunterDenovo or takes an existing ExpansionHunterDenovo profile and 
copies relevant info to the output `*.RFC1_motifs.json` file. The ExpansionHunterDenovo profile isn't used in 
calculating the fields listed above, and this feature just makes it easier to collect information from multiple tools 
in a single output file. 

Similarly, the `--run-expansion-hunter` and `--run-reviewer` options tell the script to run ExpansionHunter and REViewer 
on the RFC1 locus. If these option(s) are specified, the script will first generate a custom variant catalog for the 
repeat motif(s) it detects, then run ExpansionHunter on those motifs, and then copy ExpansionHunter results to the 
output `*.RFC1_motifs.json` file.  

**9/15/2021** Besides RFC1, this tool now also supports calling 8 known pathogenic autosomal dominant loci with 
known non-ref pathogenic motifs. These loci are BEAN1, DAB1, MARCHF6, RAPGEF2, SAMD12, STARD7, TNRC6A, YEATS2. 
Use the `--locus` option to specify the names of one or more loci to call, or use the `--all-loci` option to generate 
calls for all 9 loci. 
See [[Depienne et. al](https://www.cell.com/ajhg/pdf/S0002-9297(21)00095-1.pdf)] for more information on these loci.

Example command lines:

```
# basic command 
call_non_ref_pathogenic_motifs -R hg38.fasta -g 38 sample1.cram --locus RFC1

# run ExpansionHunter and REViewer on all 9 loci with known non-ref pathogenic motifs
call_non_ref_pathogenic_motifs -R hg38.fasta --run-expansion-hunter --run-reviewer -g 38 sample1.cram --all-loci

# for 2 specific loci, run ExpansionHunter and REViewer + provide an existing ExpansionHunterDenovo profile
call_non_ref_pathogenic_motifs -R hg38.fasta --run-expansion-hunter --run-reviewer --ehdn-profile sample1.str_profile.json -g 38 sample1.cram --locus RFC1 --locus BEAN1
```

### combine_json_to_tsv

This script can combine the `call_non_ref_pathogenic_motifs` output json files for multiple samples into 
a single .tsv. The script takes the paths of the .json files as input, or, if none are provided, it searches for .json 
files in the current directory and subdirectories.

Example command line:
```
combine_json_to_tsv  sample1.RFC1_motifs.json  sample2.RFC1_motifs.json
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

