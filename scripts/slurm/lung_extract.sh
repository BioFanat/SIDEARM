#!/bin/csh -f
#SBATCH --mem=48G
#SBATCH --output=script_out/lung_extract.out # File to which stdout will be written
#SBATCH --error=script_out/lung_extract.err

#bedtools merge -i ../../data/raw/a549_regions_sorted.bed > ../../data/raw/a549_regions_merged.bed
#(awk '{mid=$2+1; print $1"\t"mid"\t"(mid+1)"\t"$6"\t1"}' ../../data/damages/A549_data/A549_BPDE_cells_1.bed; awk '{mid=$2+1; print $1"\t"mid"\t"(mid+1)"\t"$6"\t1"}' ../../data/damages/A549_data/A549_BPDE_cells_2.bed) | sort -k1,1 -k2,2n | awk '{key=$1"\t"$2"\t"$3"\t"$4; count[key]+=$5} END{for(k in count) print k"\t"count[k]}' > ../../data/damages/a549_bpde_1_2.bed
(awk '{mid=$2+1; print $1"\t"mid"\t"(mid+1)"\t"$6"\t1"}' ../../data/damages/A549_data/A549_BPDE_cells_1.bed; awk '{mid=$2+1; print $1"\t"mid"\t"(mid+1)"\t"$6"\t1"}' ../../data/damages/A549_data/A549_BPDE_cells_2.bed) | bedtools sort | bedtools groupby -g 1,2,3,4 -c 5 -o sum | bedtools intersect -a - -b ../../data/raw/a549_regions_merged.bed | awk '{if($4=="+") print > "../../data/damages/a549_bpde_1_2_plus.bed"; else print > "../../data/damages/a549_bpde_1_2_minus.bed"}'