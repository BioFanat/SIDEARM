#!/bin/csh -f
#SBATCH --mem=48G
#SBATCH --output=script_out/all_tfs/atac_intersections_%a.out # File to which stdout will be written
#SBATCH --error=script_out/all_tfs/atac_intersections_%a.err  # File to which stderr will be written
#SBATCH --mincpus=40
#SBATCH --array=1-4

set NAME=`sed -n "${SLURM_ARRAY_TASK_ID}p" test_clusters.txt`
# python ../moods_processing.py ../../data/tfbs/fib_atac $NAME # preprocess as needed
python ../intersection_counting_naked.py $NAME $SLURM_ARRAY_TASK_ID 