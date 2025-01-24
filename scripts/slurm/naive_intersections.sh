#!/bin/csh -f
#SBATCH --mem=24G
#SBATCH --output=script_out/runtime_experiment_%a.out # File to which stdout will be written
#SBATCH --error=script_out/runtime_experiment_%a.err  # File to which stderr will be written
#SBATCH --mincpus=40
#SBATCH --array=1-6

set NAME=`sed -n "${SLURM_ARRAY_TASK_ID}p" runtime.txt`
# python ../intersection_counting.py $NAME $SLURM_ARRAY_TASK_ID
python ../naive_sims.py $NAME $SLURM_ARRAY_TASK_ID