#!/bin/csh -f
#SBATCH --mem=24G
#SBATCH --output=script_out/atac_intersections_%a.out # File to which stdout will be written
#SBATCH --error=script_out/atac_intersections_%a.err  # File to which stderr will be written
#SBATCH --mincpus=40
#SBATCH --array=1-10

set NAME=`sed -n "${SLURM_ARRAY_TASK_ID}p" top_tfs_uv.txt`
# python ../moods_processing.py ../../data/tfbs/a549 $NAME # preprocess as needed
python ../intersection_counting.py $NAME $SLURM_ARRAY_TASK_ID
