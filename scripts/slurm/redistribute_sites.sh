#!/bin/bash
#SBATCH --mem=12G
#SBATCH --output=script_out/redistribute_sites_%a.out # File to which stdout will be written
#SBATCH --error=script_out/redistribute_sites_%a.err
#SBATCH --array=1-10

# bedtools intersect -a ../../data/raw/lung_mut_C_A.bed -b ../../data/raw/a549_regions_sorted.bed > ../../data/raw/lung_mut_C_A_accessible.bed
# python ../redistribute_sites.py $SLURM_ARRAY_TASK_ID

# DAMAGE
# python ../redistribution_optimized.py \
#         --id $SLURM_ARRAY_TASK_ID \
#         --runs 1000 \
#         --input ../../data/encodings/atac_4mers_sorted_nb.bin \
#         --index ../../data/encodings/atac_naked_dmg_index.out \
#         --output /usr/xtmp/bc301/sim_uv_cpd_naked/

# MUTATIONS
python ../redistribution_mut_optimized.py \
  --input ../../data/encodings/skin_mutations_packed_atac.bin \
  --index ../../data/mutagenic_skin_3mer_accessible.out \
  --output /usr/xtmp/bc301/sim_uv_mutagenic \
  --runs 1000 \
  --start_id $SLURM_ARRAY_TASK_ID


# CHECK
# python ../redistribution_mut_optimized.py \
#   --check /usr/xtmp/bc301/sim_uv_mut/mut_run_full_1.bin \
#   --index ../../data/mutation_skin_3mer_accessible.out
