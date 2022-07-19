# str-analysis
---
This package contains scripts and utilities for analyzing short tandem repeats (STRs). 


## Installation

To install the scripts and utilities, run:

```
python3 -m pip install --upgrade str_analysis
```

Alternatively, you can use this docker image:

```
docker run -it weisburd/str-analysis:latest
```


---

### call_non_ref_pathogenic_motifs

This script takes a bam or cram file and determines which motifs are present
at known pathogenic STR loci (such as RFC1, BEAN1, DAB1, etc.) where several
motifs are known to segregate in the population. It then optionally runs
ExpansionHunterDenovo, ExpansionHunter, and/or STRling and gathers relevant
fields from their outputs that users can then compare or use for downstream analyses. It can also 
run REViewer to generate read visualization images based on the
ExpansionHunter outputs. Finally it generates a json file per locus that
contains all collected information as well as a "call" field indicating
whether pathogenic motifs were detected.

**Example command lines:**

```
# basic command 
call_non_ref_pathogenic_motifs -R hg38.fasta -g 38 sample1.cram --locus RFC1

# run ExpansionHunter and REViewer on all 9 loci with known non-ref pathogenic motifs
call_non_ref_pathogenic_motifs -R hg38.fasta --run-expansion-hunter --run-reviewer -g 38 sample1.cram --all-loci

# for 2 specific loci, run ExpansionHunter and REViewer + provide an existing ExpansionHunterDenovo profile
call_non_ref_pathogenic_motifs -R hg38.fasta --run-expansion-hunter --run-reviewer --ehdn-profile sample1.str_profile.json -g 38 sample1.cram --locus RFC1 --locus BEAN1
```

**Command-line args:**

```
positional arguments:
  bam_or_cram_path      bam or cram path.

optional arguments:
  -h, --help            show this help message and exit
  -g {GRCh37,hg19,hg37,37,GRCh38,hg38,38}, --genome-version {GRCh37,hg19,hg37,37,GRCh38,hg38,38}
                        Reference genome version
  -r REFERENCE_FASTA, --reference-fasta REFERENCE_FASTA
                        Reference fasta path.
  -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        Output filename prefix.
  -s SAMPLE_ID, --sample-id SAMPLE_ID
                        The sample id to put in the output json file. If not
                        specified, it will be retrieved from the bam/cram
                        header or filename prefix.
  --strling-genotype-table STRLING_GENOTYPE_TABLE
                        Optionally provide an existing STRling output file for
                        this sample. If specified, the script will skip
                        running STRling.
  --run-strling         Optionally run STRling and copy information relevant
                        to the locus from the STRling results to the json
                        output file.
  --strling-path STRLING_PATH
                        The path of the STRling executable to use.
  --strling-reference-index STRLING_REFERENCE_INDEX
                        Optionally provide the path of a pre-computed STRling
                        reference index file. If provided, it will save a step
                        and allow STRling to complete faster.
  --expansion-hunter-denovo-profile EXPANSION_HUNTER_DENOVO_PROFILE
                        Optionally copy information relevant to the locus from
                        this ExpansionHunterDenovo profile to the output json.
                        This is instead of --run-expansion-hunter-denovo.
  --run-expansion-hunter-denovo
                        Optionally run ExpansionHunterDenovo and copy
                        information relevant to the locus from
                        ExpansionHunterDenovo results to the output json.
  --expansion-hunter-denovo-path EXPANSION_HUNTER_DENOVO_PATH
                        The path of the ExpansionHunterDenovo executable to
                        use if --run-expansion-hunter-denovo is specified.
  --run-expansion-hunter
                        If this option is specified, this script will run
                        ExpansionHunter once for each of the motif(s) it
                        detects at the locus. ExpansionHunter doesn't
                        currently support genotyping multiallelic repeats such
                        as RFC1 where an individual may have 2 alleles with
                        motifs that differ from each other (and from the
                        reference motif). Running ExpansionHunter separately
                        for each motif provides a work-around.
  --expansion-hunter-path EXPANSION_HUNTER_PATH
                        The path of the ExpansionHunter executable to use if
                        --run-expansion-hunter is specified. This must be
                        ExpansionHunter version 3 or greater.
  --use-offtarget-regions
                        Optionally use off-target regions when counting reads
                        that support a motif, and when running
                        ExpansionHunter.
  --run-reviewer        Run the REViewer tool to visualize ExpansionHunter
                        output. --run-expansion-hunter must also be specified.
  --run-reviewer-for-pathogenic-calls
                        Run the REViewer tool to visualize ExpansionHunter
                        output only when this script calls a sample as having
                        PATHOGENIC MOTIF / PATHOGENIC MOTIF. --run-expansion-
                        hunter must also be specified.
  --all-loci            Generate calls for all these loci: RFC1, BEAN1, DAB1,
                        MARCHF6, RAPGEF2, SAMD12, STARD7, TNRC6A, YEATS2
  -l {RFC1,BEAN1,DAB1,MARCHF6,RAPGEF2,SAMD12,STARD7,TNRC6A,YEATS2}, --locus {RFC1,BEAN1,DAB1,MARCHF6,RAPGEF2,SAMD12,STARD7,TNRC6A,YEATS2}
                        Generate calls for this specific locus. This argument
                        can be specified more than once to call multiple loci.
  -v, --verbose         Print detailed log messages.
```


**Command line arguments:**

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
  --use-offtarget-regions  
                        Optionally use off-target regions when counting reads 
                        that support a motif, and when running ExpansionHunter.
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

**Summary of `*_motifs.json` output file:**

The `call_non_ref_pathogenic_motifs` script outputs a .json file with many fields summarizing what it found, as well as
ExpansionHunter, REViewer, ExpansionHunterDenovo, and STRling outputs when `--run-expansion-hunter`, `--run-reviewer`, and/or
`--run-expansion-hunter-denovo` or `--run-strling` args are used.
  
The key fields are:

**call**: (ex. `BENIGN MOTIF / PATHOGENIC MOTIF`)   
**motif1_repeat_unit**: (ex. `AAAGG`)  
**motif2_repeat_unit**: (ex. `AAGGG`)  
**expansion_hunter_call_repeat_unit**: (ex. `AAAGG / AAGGG`) The same information as above, but in genotype form.
**expansion_hunter_call_genotype**: (ex. `15/68`) The number of repeats of the motif(s) above.  
**expansion_hunter_call_CI**:  (ex. `15-15/55-87`) The confidence intervals for the genotype.   
**expansion_hunter_call_reviewer_svg**: (ex. `sample1.RFC1_AAGGG.expansion_hunter_reviewer.svg`) The REViewer read visualization image path.   


**Output fields:**

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
  
*NOTE:* the expansion_hunter_* fields below will be added when `--run-expansion-hunter` is used.
    
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

*NOTE:* the *_reviewer_svg fields below are only added when `--run-reviewer` is used.    
  
**expansion_hunter_motif1_reviewer_svg** Path of .svg image file generated by REViewer for this motif.  
**expansion_hunter_motif2_reviewer_svg** If a 2nd motif was detected at this locus, this will be the path of the .svg image file generated by REViewer for this other motif.  
**expansion_hunter_call_reviewer_svg** The final .svg image file. If more than one motif was detected, this will contain 
  a merged image with 1 short allele panel and 1 long allele panel selected from the motif1 and motif2 images above.  

*NOTE:* the ehn_ fields below are only added when `--run-expansion-hunter-denovo` is used.

**ehdn_motif1_repeat_unit**  
**ehdn_motif1_anchored_irr_count**  
**ehdn_motif1_paired_irr_count**  
**ehdn_motif1_total_irr_count**    
**ehdn_sample_read_depth**  
**ehdn_motif1_n_anchored_regions**  


**Example Output:**

Lets say the script detects that sample1 contains two sets of reads at the RFC1 locus - some with the AAGGG motif and 
some with AAAAG. The script then runs ExpansionHunter for the AAGGG motif and then runs it again for the AAAAG motif. 
Lets say ExpansionHunter outputs 15/73 as the genotype for AAGGG and 15/22 as the genotype for AAAAG.
This script would then output:

**call**: (`BENIGN MOTIF / PATHOGENIC MOTIF`)  
**expansion_hunter_call_repeat_unit**: (`AAAAG / AAGGG`)  
**expansion_hunter_call_genotype**: (`15/73`) Selected from the two genotypes above.  
**expansion_hunter_call_CI**: (`15-15/55-87`) Selected from two sets of confidence intervals.  
**expansion_hunter_call_reviewer_svg** (`sample1.RFC1_AAGGG.expansion_hunter_reviewer.svg`)    

