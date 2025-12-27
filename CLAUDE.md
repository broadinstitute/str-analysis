## Overview

This repository contains tools for analyzing tandem repeat (TR) variation

## Principles

- Prefer correctness and reproducibility over speed.
- Avoid introducing new dependencies unless necessary.
- Favor clear, explicit code over clever abstractions.
- Maintain backward compatibility unless explicitly asked otherwise.

## Code Style

- Python 3.9+
- Avoid adding type hints
- Prefer pure functions
- Avoid global state
- Use pathlib, not os.path
- print() instead of logging
- Use docstring format from the Google style guide for python (https://google.github.io/styleguide/pyguide.html) 

## Preferred Tools

- pandas, numpy, scipy
- pysam for BAM/CRAM

## Clarifications

If something is ambiguous, ask one concise question before proceeding.

## Core Utility Modules (`str_analysis/utils/`)

- **`find_repeat_unit.py`**: Core algorithm to detect tandem repeats in sequences. Supports pure repeats and repeats with interruptions
- **`find_motif_utils.py`**: Functions for finding and analyzing repeat motifs
- **`canonical_repeat_unit.py`**: Converts repeat motifs to canonical form (handles rotations and reverse complements)
- **`eh_catalog_utils.py`**: Parse and manipulate ExpansionHunter variant catalogs (JSON/BED format)
- **`trf_runner.py`**: Interface to Tandem Repeats Finder (TRF) for repeat detection
- **`file_utils.py`**: File I/O utilities (handles .gz, .bgz compression transparently)
- **`fasta_utils.py`**: FASTA reference genome utilities using pyfaidx
- **`cram_bam_utils.py`**: BAM/CRAM file manipulation using pysam
- **`gtf_utils.py`**: GTF gene annotation parsing and overlap checking

## Main Command-Line Tools

Tools follow a consistent pattern:
1. Parse arguments with `argparse`
2. Use utilities from `str_analysis/utils/`
3. Process genomic data (VCF, BAM/CRAM, variant catalogs)
4. Output results (TSV tables, JSON, modified catalogs)

Major tool categories:

**Variant Filtering & Detection:**
- `filter_vcf_to_STR_variants.py` - Extract STR variants from VCF using brute-force k-mer search
- `filter_vcf_to_tandem_repeats.py` - Newer version of above algorithm
- `call_non_ref_motifs.py` - Detect non-reference motifs at disease loci (especially RFC1)

**Catalog Management:**
- `merge_loci.py` - Combine multiple STR catalogs, removing overlapping duplicates
- `annotate_and_filter_str_catalog.py` - Annotate loci with gene overlap, filter by criteria

## Terminology

- **TR** = tandem repeat.
- **Loci** are identified by chr:start-end (0-based, half-open).
- **Motif/Repeat Unit**: The sequence that repeats (e.g., "CAG")
- **Canonical Motif**: Normalized form considering rotations and reverse complement
- **Pure Repeat**: Uninterrupted tandem repeats
- **Interrupted Repeat**: Tandem repeats with insertions/deletions/substitutions
- **LocusStructure**: Pattern notation like `(CAG)*` or `ACGT(CAG)*GCTA(GCC)+`
- **Off-target Regions**: Genomic regions where reads with a given motif may spuriously map
- **Adjacent Loci**: Nearby STRs that should be genotyped together

## Test Organization

Tests use Python's `unittest` framework:
- Test files named `*_tests.py` or `*tests.py`
- Located alongside the modules they test in `str_analysis/` and `str_analysis/utils/`
- Some manual tests named `*__manual_unittest.py` (excluded from automated runs)
- Data files for tests in `str_analysis/data/tests/`
