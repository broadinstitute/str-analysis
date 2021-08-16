# str-analysis
This package contains scripts and utilities related to analyzing short tandem repeats (STRs). 

---
## Scripts

### call_rfc1_canvas_alleles

This script takes a WGS bam or cram file and outputs a .json file containing details about alleles it 
detected at the RFC1/CANVAS STR locus. The main fields in the output dictionary are:

**sample_id**: *If this value is not specified as a command line arg, it is parsed from the input bam/cram file header or filename prefix.*    
**call**: *describes the alleles detected at the RFC1/CANVAS locus. Its format is analogous to a VCF genotype. Possible values are:*
* "PATHOGENIC MOTIF / PATHOGENIC MOTIF": *only pathogenic allele(s) detected*
* "BENIGN MOTIF / BENIGN MOTIF": *only benign allele(s) detected*
* "MOTIF OF UNCERTAIN SIGNIFICANCE / MOTIF OF UNCERTAIN SIGNIFICANCE": *non-canonical allele(s) detected with unknown pathogenicity*
* "BENIGN MOTIF / PATHOGENIC MOTIF": *heterozygous for a benign allele and a pathogenic allele, implying carrier status*
* "PATHOGENIC MOTIF / MOTIF OF UNCERTAIN SIGNIFICANCE": *heterozygous for a pathogenic allele and a non-canonical allele(s) detected with unknown pathogenicity*
* "BENIGN MOTIF / MOTIF OF UNCERTAIN SIGNIFICANCE": *heterozygous for a benign allele and a non-canonical allele(s) detected with unknown pathogenicity*
* null: *not enough evidence in the read data to support any of the above options*

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
**found_repeats_in_fraction_of_reads**: `found_repeats_in_n_reads` / `found_n_reads_overlap_rfc1_locus  

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

