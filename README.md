# str-analysis
---
This package contains scripts and utilities related to analyzing short tandem repeats (STRs). 


## Installation

To install the scripts or utilities, run:

```
python3 -m pip install --upgrade str_analysis
```

Alternatively, you can use this docker image:

```
docker run -it weisburd/str-analysis:latest
```


---
## Scripts

### call_non_ref_pathogenic_motifs

RFC1 STR expansions have recently been linked to [CANVAS](https://www.omim.org/entry/614575) [ [Cortese 2019](https://pubmed.ncbi.nlm.nih.gov/30926972/) ].
This STR locus is unique in that it's both autosomal recessive and has a pathogenic repeat motif (AAGGG) 
that differs from the motif in the reference genome (AAAAG). Several other benign and pathogenic motifs have also been reported, such as  
AAAGG and ACAGG, as well as motifs of uncertain significance such as AAGAG and AGAGG [ [Akcimen 2019](https://pubmed.ncbi.nlm.nih.gov/31824583/) ].
It's not unusual for individuals to have one motif on one chromosome, and another motif on the other chromosome. 
Due to this multi-allelic nature, current STR genotyping tools like 
[ExpansionHunter](https://github.com/Illumina/ExpansionHunter) struggle to accurately genotype RFC1 since they require
the repeat unit to be specified apriori as input, and then only look for reads supporting that repeat unit.  

A more recent tool, [ExpansionHunter Denovo](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02017-z), 
is better at detecting multi-allelic expanded RFC1 motif(s), but is unable to distinguish unaffected carriers 
(such as those that have one benign AAAAG reference allele, and one pathogenic AAGGG allele) from 
affected individuals (homozygous for the AAGGG allele). In both cases, ExpansionHunter Denovo might detect the 
AAGGG allele, but would fail to detect the benign allele unless it is also highly expanded. 

This script addresses some of the above limitations. It takes a whole genome (WGS) bam or cram file and detects
which motif(s) are present at the RFC1/CANVAS STR locus. At a minimum, it then outputs a .json file with information on 
the detected motifs. Optionally it also runs ExpansionHunter using a custom variant catalog for each detected motif, 
and combines the results into a single diploid genotype and confidence interval. The custom variant catalog includes
off-target regions, which allow ExpansionHunter to accurately genotype expansions longer than fragment length. 
The script can also then run [REViewer](https://github.com/Illumina/REViewer) for each motif and combine the resulting 
images into a single diploid image where, if multiple motifs were detected, the short and long alleles are based on 
different motifs.  

Testing this approach on simulated data shows that it restores ExpansionHunter's accurately for non-reference motif 
expansions to the same level as for expansions with the reference (AAAAG) motif. The simulation results suggested that, 
to combine ExpansionHunter results for two different motifs, it makes sense to report the longest of the two long 
alleles, and then report the short allele from the other motif, so this is the rule this script uses.
See the Example Scenario section below for a concrete example.  

The script can also apply this approach to other known STR loci where the pathogenic motif differs from the reference 
such as DAB1, BEAN1, SAMD12, and others. Unlike RFC1/CANVAS, these other loci are autosomal dominant, which simplifies 
the task and makes it more ammenable to ExpansionHunterDenovo. To assist with this, the script also has a 
`--run-expansion-hunter-denovo` option to run ExpansionHunterDenovo on the input bam, or alternatively 
the `--ehdn-profile` to pass in an existing ExpansionHunterDenovo output file so that all relevant results from this 
script, ExpansionHunter, and ExpansionHunterDenovo can be combined into a single output file.  


**Example command lines:**

```
# basic command 
call_non_ref_pathogenic_motifs -R hg38.fasta -g 38 sample1.cram --locus RFC1

# run ExpansionHunter and REViewer on all 9 loci with known non-ref pathogenic motifs
call_non_ref_pathogenic_motifs -R hg38.fasta --run-expansion-hunter --run-reviewer -g 38 sample1.cram --all-loci

# for 2 specific loci, run ExpansionHunter and REViewer + provide an existing ExpansionHunterDenovo profile
call_non_ref_pathogenic_motifs -R hg38.fasta --run-expansion-hunter --run-reviewer --ehdn-profile sample1.str_profile.json -g 38 sample1.cram --locus RFC1 --locus BEAN1
```

** Command line arguments:**

```
positional arguments:
  bam_or_cram_path      bam or cram path

optional arguments:
  -h, --help            show this help message and exit
  -g {GRCh37,hg19,hg37,37,GRCh38,hg38,38}, --genome-version {GRCh37,hg19,hg37,37,GRCh38,hg38,38}
  -R REFERENCE_FASTA, --reference-fasta REFERENCE_FASTA
                        Reference fasta path. The reference fasta is sometimes
                        necessary for decoding cram files.
  -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        Output filename prefix
  -s SAMPLE_ID, --sample-id SAMPLE_ID
                        The sample id to put in the output json file. If not
                        specified, it will be retrieved from the bam/cram
                        header or filename prefix.
  --run-expansion-hunter-denovo
                        Optionally run ExpansionHunterDenovo and copy
                        information relevant to the locus from
                        ExpansionHunterDenovo results to the output json.
  --expansion-hunter-denovo-path EXPANSION_HUNTER_DENOVO_PATH
                        The path of the ExpansionHunterDenovo executable to
                        use when --expansion-hunter-denovo-path is specified.
  --expansion-hunter-denovo-profile EXPANSION_HUNTER_DENOVO_PROFILE
                        Optionally copy information relevant to the locus from
                        this ExpansionHunterDenovo profile to the output json.
                        This is instead of --run-expansion-hunter-denovo
  -r, --run-expansion-hunter
                        If this option is specified, this script will run
                        ExpansionHunter once for each of the motif(s) it
                        detects at the locus. ExpansionHunter doesn't
                        currently support genotyping multiallelic repeats such
                        as RFC1 where an individual may have 2 alleles with
                        motifs that differ from eachother (and from the
                        reference motif). Running ExpansionHunter separately
                        for each motif provides a work-around.
  --expansion-hunter-path EXPANSION_HUNTER_PATH
                        The path of the ExpansionHunter executable to use if
                        -r is specified. This must be ExpansionHunter version
                        3 or greater.
  --all-loci            Generate calls for all these loci: RFC1, BEAN1, DAB1,
                        MARCHF6, RAPGEF2, SAMD12, STARD7, TNRC6A, YEATS2
  -l {RFC1,BEAN1,DAB1,MARCHF6,RAPGEF2,SAMD12,STARD7,TNRC6A,YEATS2}, --locus {RFC1,BEAN1,DAB1,MARCHF6,RAPGEF2,SAMD12,STARD7,TNRC6A,YEATS2}
                        Generate calls for this specific locus. This argument
                        can be specified more than once to call multiple loci.
  --run-reviewer        Run the REViewer tool to visualize ExpansionHunter
                        output. --run-expansion-hunter must also be specified.
  --run-reviewer-for-pathogenic-calls
                        Run the REViewer tool to visualize ExpansionHunter
                        output only when this script calls a sample as having
                        PATHOGENIC MOTIF / PATHOGENIC MOTIF. --run-expansion-
                        hunter must also be specified.
  -v, --verbose         Print detailed log messages
```

**`*_motifs.json` output file:**

The `call_non_ref_pathogenic_motifs` script outputs a .json file with many fields summarizing what it found, as well as
ExpansionHunter, REViewer, and ExpansionHunterDenovo outputs when `--run-expansion-hunter`, `--run-reviewer`, and/or
`--run-expansion-hunter-denovo` args are used.
  
The key fields are:

**call**: (ex. `BENIGN MOTIF / PATHOGENIC MOTIF`)   
**motif1_repeat_unit**: (ex. `AAAGG`)  
**motif2_repeat_unit**: (ex. `AAGGG`)  
**expansion_hunter_call_repeat_unit**: (ex. `AAAGG / AAGGG`) The same information as above, but in genotype form.
**expansion_hunter_call_genotype**: (ex. `15/68`) The number of repeats of the motif(s) above.  
**expansion_hunter_call_CI**:  (ex. `15-15/55-87`) The confidence intervals for the genotype.   
**expansion_hunter_call_reviewer_svg**: (ex. `sample1.RFC1_AAGGG.expansion_hunter_reviewer.svg`) The REViewer read visualization image path.   


Output fields:

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
  
**motif2_repeat_unit**: *If more than one motif is detected, this field will contain the repeat unit that is supported by the next most reads*      
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
  
*NOTE: the expansion_hunter_* fields below will be added when `--run-expansion-hunter` is used.
    
**expansion_hunter_motif1_json_output_file** Path of the ExpansionHunter output json file when run on the motif1.    
**expansion_hunter_motif1_repeat_unit** Repeat unit passed to ExpansionHunter.  
**expansion_hunter_motif1_short_allele_genotype**  ExpansionHunter output genotype (number of repeats) of the short allele.   
**expansion_hunter_motif1_long_allele_genotype**  ExpansionHunter output genotype (number of repeats) of the long allele.  
**expansion_hunter_motif1_short_allele_CI_start**  ExpansionHunter output genotype confidence interval lower bound of the short allele.      
**expansion_hunter_motif1_short_allele_CI_end** ExpansionHunter output genotype confidence interval upper bound of the short allele.  
**expansion_hunter_motif1_long_allele_CI_start** ExpansionHunter output genotype confidence interval lower bound of the long allele.  
**expansion_hunter_motif1_long_allele_CI_end** ExpansionHunter output genotype confidence interval upper bound of the long allele.
**expansion_hunter_motif1_total_spanning_reads** ExpansionHunter output total number of spanning reads supporting the genotype for motif1.  
**expansion_hunter_motif1_total_flanking_reads** ExpansionHunter output total number of flanking reads supporting the genotype for motif1.  
**expansion_hunter_motif1_total_inrepeat_reads**  ExpansionHunter output total number of IRR reads supporting the genotype for motif1.  

**expansion_hunter_motif2_json_output_file** If a 2nd motif was detected at this locus, this will be the path of the
    ExpansionHunter output json file when run on the 2nd motif. All the fields listed above for motif1 will also be 
    present for motif2.
...

**expansion_hunter_call_repeat_unit** The repeat unit(s) used to run ExpansionHunter. If more than one was detected, this will have a format like `AAAAG / AAGGG`.  
**expansion_hunter_call_genotype** ExpansionHunter combined genotype based on the results of the ExpansionHunter run(s). 
**expansion_hunter_call_CI** ExpansionHunter output total number of spanning reads supporting the genotype for motif1.

*NOTE: the *_reviewer_svg fields below are only added when `--run-reviewer` is used.    
  
**expansion_hunter_motif1_reviewer_svg** Path of .svg image file generated by REViewer for this motif.  
**expansion_hunter_motif2_reviewer_svg** If a 2nd motif was detected at this locus, this will be the path of the .svg image file generated by REViewer for this other motif.  
**expansion_hunter_call_reviewer_svg** The final .svg image file. If more than one motif was detected, this will contain 
  a merged image with 1 short allele panel and 1 long allele panel selected from the motif1 and motif2 images above.  

*NOTE: the ehn_ fields below are only added when `--run-expansion-hunter-denovo` is used.

**ehdn_motif1_repeat_unit**  
**ehdn_motif1_anchored_irr_count**  
**ehdn_motif1_paired_irr_count**  
**ehdn_motif1_total_irr_count**    
**ehdn_sample_read_depth**  
**ehdn_motif1_n_anchored_regions**  


**Example Scenario:**

Lets say the script detects that sample1 contains two sets of reads at the RFC1 locus - some with the AAGGG motif and 
some with AAAAG. The script then runs ExpansionHunter for the AAGGG motif and then runs it again for the AAAAG motif. 
Lets say ExpansionHunter outputs 15/73 as the genotype for AAGGG and 15/22 as the genotype for AAAAG.
This script would then output:

**call**: (`BENIGN MOTIF / PATHOGENIC MOTIF`)  
**expansion_hunter_call_repeat_unit**: (`AAAAG / AAGGG`)  
**expansion_hunter_call_genotype**: (`15/73`) Selected from the two genotypes above.  
**expansion_hunter_call_CI**: (`15-15/55-87`) Selected from two sets of confidence intervals.  
**expansion_hunter_call_reviewer_svg** (`sample1.RFC1_AAGGG.expansion_hunter_reviewer.svg`) Merged from the two REViewer images for motifs 1 and 2.   


### combine_json_to_tsv

This script can combine the `call_non_ref_pathogenic_motifs` output json files for multiple samples into 
a single .tsv. The script takes the paths of the .json files as input, or, if none are provided, it searches for .json 
files in the current directory and subdirectories.

Example command line:
```
combine_json_to_tsv  sample1.RFC1_motifs.json  sample2.RFC1_motifs.json
```
