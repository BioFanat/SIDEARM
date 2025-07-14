## Overview

SIDEARM is a computational framework for analyzing DNA damage patterns and their relationship to mutational signatures in genomic regions. The core methodology employs simulation-based inference to test whether observed mutational patterns at specific genomic features (e.g., transcription factor binding sites, enhancers) show enrichment or depletion compared to random distribution expectations.

## Core Concept: Simulation-Inference

The simulation-inference approach works by:

1. **Measuring actual damage/mutation patterns**: Count real damage events or mutations at specific genomic features
2. **Creating null distributions**: Simulate many random redistributions of the same total damage/mutation counts across all possible genomic locations
3. **Statistical comparison**: Compare actual patterns to simulated distributions to identify statistically significant enrichment or depletion

This approach isolates the contribution from sequence context biases (e.g., certain k-mers are more damage-prone) to directly model the contribution TF binding activity has on damage formation and mutagenesis.

## Key Components

- **Encoding Utilities**: Functions for efficiently encoding and decoding genomic positions using bit-packed integers
- **Intersection Analysis**: Core algorithms for identifying overlaps between damage sites and transcription factor binding sites
- **Statistical Framework**: Simulation methods to establish statistical significance with FDR correction
- **Visualization Tools**: Plot generation for damage profiles, Z-scores, p-values, and sequence logos

## Commands

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


This repository containing code for the SIDEARMv1 pipeline by Chi et. al., as submitted for review in the ISMB 2025 Proceedings. Version 2 of SIDEARM is currently in the works, soon to be released within the next two weeks!
