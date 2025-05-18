# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

SIDEARM (Statistical Inference for Damage and Repair at Motifs) is a pipeline for analyzing DNA damage and repair around transcription factor binding sites. The code contains multiple Python scripts for processing, encoding, and analyzing genomic data including:

1. Processing transcription factor binding sites (TFBS)
2. Analyzing damage patterns (CPD-seq, BPDE)
3. Statistical analysis of damage/mutation profiles around binding motifs
4. Visualization of results with sequence logos and statistical metrics
5. Encoding genomic positions and motifs for efficient computation

## Key Components

- **Encoding Utilities**: Functions for efficiently encoding and decoding genomic positions using bit-packed integers
- **Intersection Analysis**: Core algorithms for identifying overlaps between damage sites and transcription factor binding sites
- **Statistical Framework**: Simulation methods to establish statistical significance with FDR correction
- **Visualization Tools**: Plot generation for damage profiles, Z-scores, p-values, and sequence logos

## Common Commands

### Running Intersection Analysis

To run the intersection analysis for a specific transcription factor:
```bash
python scripts/intersection_counting.py <TF_NAME> <RUN_ID>
```

### Running SLURM Batch Jobs

For large-scale analysis on SLURM-managed clusters:
```bash
sbatch scripts/slurm/normal_intersection.sh
```

### Running with Optimized Code (Newer Versions)

For optimized versions of the analysis:
```bash
python scripts/intersection_optimized.py <TF_NAME> <RUN_ID>
python scripts/redistribution_optimized.py <INPUT_FILE> <OUTPUT_FILE>
```

### Running with Generalized k-mer System

For analysis with configurable k-mer lengths:
```bash
# Extract k-mers from genomic intervals
python scripts/interval_extraction_generalized.py --genome <GENOME_FILE> --plus <PLUS_FILE> --minus <MINUS_FILE> --intervals <INTERVALS_FILE> --output-bin <OUTPUT_BIN> --output-index <OUTPUT_INDEX> --damage-type <DAMAGE_TYPE>

# Redistribute damages
python scripts/redistribution_generalized.py --id <START_ID> --runs <NUM_RUNS> --input <INPUT_BIN> --index <INDEX_FILE> --output <OUTPUT_DIR> --damage-type <DAMAGE_TYPE>

# Redistribute mutations
python scripts/redistribution_mut_generalized.py --input <INPUT_BIN> --index <INDEX_FILE> --output <OUTPUT_DIR> --runs <NUM_RUNS> --start_id <START_ID> --damage-type <DAMAGE_TYPE>
```

Where `<DAMAGE_TYPE>` can be 'CPD', 'BPDE', or '6MER' (or any custom type defined in kmer_config.py).

### Processing Mutation Data

For mutation data analysis:
```bash
python scripts/redistribution_mut_optimized.py <INPUT_FILE> <OUTPUT_FILE>
```

## Code Structure

### Core Modules

1. **Encoding/Decoding**:
   - Binary encoding of genomic positions for efficient storage
   - Functions like `decode()`, `decompress_mut()`, and `pack_tuple_nochrom()`
   - New generalized system in `kmer_config.py` for configurable k-mer lengths

2. **Matching & Overlap**:
   - Functions to match sequences to encoded representations
   - Functions to identify overlaps between damage sites and TF binding sites

3. **Statistical Analysis**:
   - Simulation framework for significance testing
   - Functions for calculating Z-scores and FDR-corrected p-values

4. **Visualization**:
   - Plot generation using Matplotlib and Logomaker
   - Visualization of damage profiles, sequence logos, and statistical metrics

### Important Files

- `scripts/program_utils.py`: Core shared utilities for file loading, sequence manipulation, and encoding
- `scripts/intersection_counting.py`: Main script for analyzing overlaps and significance
- `scripts/intersection_optimized.py`: Optimized version of the intersection counter
- `scripts/redistribution_optimized.py`: Optimized redistribution algorithms
- `scripts/kmer_config.py`: Configurable k-mer system for different damage types
- `scripts/interval_extraction_generalized.py`: Generalized k-mer extraction
- `scripts/redistribution_generalized.py`: Generalized damage redistribution
- `scripts/redistribution_mut_generalized.py`: Generalized mutation redistribution
- `scripts/22/encoding_utils.py` and `scripts/22/general_utils.py`: Newer utility functions
- `scripts/genome_utils/program_utils.py`: Utilities for genomic analysis

## Development Notes

- Many functions need a genome dictionary loaded from a FASTA file for sequence context
- Most analyses operate on processed BED files containing genomic positions and values
- File paths are often hard-coded and may need to be adjusted
- Analysis of both cellular and "naked" DNA damage is typically performed for comparison
- The code uses multiprocessing for parallel processing of simulations