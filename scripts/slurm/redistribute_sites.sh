#!/bin/bash
#SBATCH --mem=12G
#SBATCH --output=script_out/redistribute_sites_%a.out # File to which stdout will be written
#SBATCH --error=script_out/redistribute_sites_%a.err
#SBATCH --array=1-1

# bedtools intersect -a ../../data/raw/lung_mut_C_A.bed -b ../../data/raw/a549_regions_sorted.bed > ../../data/raw/lung_mut_C_A_accessible.bed
python ../redistribute_sites.py $SLURM_ARRAY_TASK_ID

# python ../redistribution_optimized.py \
#         --id $SLURM_ARRAY_TASK_ID \
#         --runs 1000 \
#         --input ../../data/encodings/atac_4mers_sorted_nb.bin \
#         --index ../../data/encodings/atac_nb_dmg_index.out \
#         --output /usr/xtmp/bc301/sim_uv_cpd_full/