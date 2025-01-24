#!/bin/csh -f
#SBATCH --mem=12G
#SBATCH --output=script_out/write_encodings_%a.out # File to which stdout will be written
#SBATCH --error=script_out/write_encodings_%a.err
#SBATCH --array=1-10


# Note the different syntax for command substitution and variable referencing in csh
python ../write_encodings.py $SLURM_ARRAY_TASK_ID

# bedtools intersect -a ../../data/raw/mutations_au_transitions_cor.bed -b ../../data/raw/atac_150bp_intervals_merged.bed > ../../data/raw/atac_mutations_transitions_correct.bed
# awk '{ if ($5 == "C" || $5 == "G") print }' ../../data/raw/atac_mutations_transitions_correct.bed > ../../data/raw/atac_mutations_transitions_C_only.bed
