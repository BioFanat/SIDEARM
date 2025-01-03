#!/bin/csh -f
#SBATCH --mem=24G
#SBATCH --output=script_out/naive_intersections_%a.out # File to which stdout will be written
#SBATCH --error=script_out/naive_intersections_%a.err  # File to which stderr will be written
#SBATCH --mincpus=40
#SBATCH --array=1-10

set NAME=`sed -n "${SLURM_ARRAY_TASK_ID}p" new_top_tfs.txt`
# python ../intersection_counting.py $NAME
python ../naive_sims.py $NAME