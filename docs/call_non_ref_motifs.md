
### call_non_ref_motifs

This script takes a BAM or CRAM file and determines which motifs are present
at known disease-associated STR loci where non-reference 
motifs are seen in the population (RFC1, BEAN1, DAB1, EIF4A3, MARCHF6, RAPGEF2, SAMD12, STARD7, TNRC6A, YEATS2).  
It's main output is a TSV file with a "call" field that indicates whether known pathogenic motifs
were detected in a given sample. We have found this script to be useful for diagnosing pathogenic RFC1 expansions
in rare disease cohorts. Its utility for other loci remains to be seen.

Besides motif detection, the script can also optionally run ExpansionHunterDenovo, ExpansionHunter,
and/or STRling, collect their outputs, and include them in the single TSV output file.
It can also run REViewer to generate read visualization images for the ExpansionHunter calls.

**Example command lines:**

```
# basic command 
call_non_ref_motifs -r hg38.fasta -g 38 sample1.cram --known-loci --locus RFC1 

# run ExpansionHunter and REViewer on all 9 loci with known non-ref pathogenic motifs
call_non_ref_motifs -r hg38.fasta -g 38 --run-expansion-hunter --run-reviewer sample1.cram --known-loci

# for 2 specific loci, run ExpansionHunter and REViewer + provide an existing ExpansionHunterDenovo profile
call_non_ref_motifs -r hg38.fasta -g 38 --run-expansion-hunter --run-reviewer --ehdn-profile sample1.str_profile.json sample1.cram --known-loci --locus RFC1 --locus BEAN1

# run the script on other loci defined in a custom ExpansionHunter variant catalog
call_non_ref_motifs -r hg38.fasta -g 38 --variant-catalog custom_variant_catalog.json sample1.cram 
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

