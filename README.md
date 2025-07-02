## Overview

SIDEARM is a computational framework for analyzing DNA damage patterns and their relationship to mutational signatures in genomic regions. The core methodology employs simulation-based inference to test whether observed mutational patterns at specific genomic features (e.g., transcription factor binding sites, enhancers) show enrichment or depletion compared to random distribution expectations.

## Core Concept: Simulation-Inference

The simulation-inference approach works by:

1. **Measuring actual damage/mutation patterns**: Count real damage events or mutations at specific genomic features
2. **Creating null distributions**: Simulate many random redistributions of the same total damage/mutation counts across all possible genomic locations
3. **Statistical comparison**: Compare actual patterns to simulated distributions to identify statistically significant enrichment or depletion

This approach isolates the contribution from sequence context biases (e.g., certain k-mers are more damage-prone) to directly model the contribution TF binding activity has on damage formation and mutagenesis.

This repository containing code for the SIDEARMv1 pipeline by Chi et. al., as submitted for review in the ISMB 2025 Proceedings. Version 2 of SIDEARM is currently in the works, soon to be released within the next two weeks!
