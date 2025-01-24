#!/bin/csh -f
#SBATCH --mem=12G
#SBATCH --output=script_out/remote_pulling.out # File to which stdout will be written
#SBATCH --error=script_out/remote_pulling.err

tabix https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0/hg38.archetype_motifs.v1.0.bed.gz -R ../../data/raw/a549_regions_sorted.bed > ../../data/tfbs/archetype_motifs_intersect_a549.bed
