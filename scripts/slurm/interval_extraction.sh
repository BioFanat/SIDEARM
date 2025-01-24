#!/bin/csh -f
#SBATCH --mem=12G
#SBATCH --output=script_out/intervals.out # File to which stdout will be written
#SBATCH --error=script_out/intervals.err


# bedtools sort -i ../../data/raw/atac_150bp_intervals.bed > ../../data/raw/atac_150bp_intervals_sorted.bed

# bedtools merge -i ../../data/raw/atac_150bp_intervals_sorted.bed > ../../data/raw/atac_150bp_intervals_merged.bed

# bedtools getfasta -fi ../../data/genome/hg38.fa -bed ../../data/raw/a549_regions_merged.bed -fo ../../data/raw/a549_regions_merged.fa.out
python ../interval_extraction.py
# python ../lung_extract.py /usr/xtmp/hiw4/for_bo/Lung_Tissue_data_ENCODE/ENCODE_intersected_DNase_seq_peaks_summits_3kb_intervals.bed > ../../data/raw/all_lung_150bp_intervals.bed
