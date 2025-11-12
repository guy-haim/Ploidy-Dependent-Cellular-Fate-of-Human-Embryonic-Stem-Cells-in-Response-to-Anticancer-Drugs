#!/bin/bash



trimmomatic SE /sci/labs/nissimb/guyhaim/icore-data/Xeno_transcriptome_alignment/fastq/10R_4c_p48_5_FU_2_5uM_48h_S9_L001_R1_001.fastq.gz /sci/labs/nissimb/guyhaim/icore-data/Xeno_transcriptome_alignment/trim/10R_4c_p48_5_FU_2_5uM_48h_S9_L001_R1_001/10R_4c_p48_5_FU_2_5uM_48h_S9_L001_R1_001_trimmed.fastq ILLUMINACLIP:/sci/labs/nissimb/guyhaim/icore-data/Xeno_transcriptome_alignment/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


#!/bin/bash

# Human:
STAR --runThreadN 8 --genomeDir /sci/labs/nissimb/guyhaim/icore-data/Xeno_transcriptome_alignment/genome_index --readFilesIn /sci/labs/nissimb/guyhaim/icore-data/Xeno_transcriptome_alignment/trim/10R_4c_p48_5_FU_2_5uM_48h_S9_L001_R1_001/10R_4c_p48_5_FU_2_5uM_48h_S9_L001_R1_001_trimmed.fastq --outSAMtype BAM SortedByCoordinate --twopassMode Basic --outSAMattributes NM --outFileNamePrefix /sci/labs/nissimb/guyhaim/icore-data/Xeno_transcriptome_alignment/star_out/10R_4c_p48_5_FU_2_5uM_48h_S9_L001_R1_001/10R_4c_p48_5_FU_2_5uM_48h_S9_L001_R1_001_

# Mouse:
#STAR# --runThreadN 8# --genomeDir /sci/labs/nissimb/guyhaim/icore-data/Xeno_transcriptome_alignment/genome_index_mouse# --readFilesIn /sci/labs/nissimb/guyhaim/icore-data/Xeno_transcriptome_alignment/trim/10R_4c_p48_5_FU_2_5uM_48h_S9_L001_R1_001/10R_4c_p48_5_FU_2_5uM_48h_S9_L001_R1_001_trimmed.fastq# --outSAMtype BAM SortedByCoordinate# --twopassMode Basic# --outSAMattributes NM# --outFileNamePrefix /sci/labs/nissimb/guyhaim/icore-data/Xeno_transcriptome_alignment/star_mouse_out/10R_4c_p48_5_FU_2_5uM_48h_S9_L001_R1_001/10R_4c_p48_5_FU_2_5uM_48h_S9_L001_R1_001_mouse_

human=/sci/labs/nissimb/guyhaim/icore-data/Xeno_transcriptome_alignment/star_out/10R_4c_p48_5_FU_2_5uM_48h_S9_L001_R1_001/10R_4c_p48_5_FU_2_5uM_48h_S9_L001_R1_001_Aligned.sortedByCoord.out.bam
mouse=/sci/labs/nissimb/guyhaim/icore-data/Xeno_transcriptome_alignment/star_mouse_out/10R_4c_p48_5_FU_2_5uM_48h_S9_L001_R1_001/10R_4c_p48_5_FU_2_5uM_48h_S9_L001_R1_001_mouse_Aligned.sortedByCoord.out.bam
destination=/sci/labs/nissimb/guyhaim/icore-data/Xeno_transcriptome_alignment/bam/10R_4c_p48_5_FU_2_5uM_48h_S9_L001_R1_001

Rscript --vanilla ./xenofilterr.R ${human} ${mouse} ${destination}

featureCounts -T 8 -a /sci/labs/nissimb/guyhaim/icore-data/Xeno_transcriptome_alignment/gencode.v38.primary_assembly.annotation.gtf -o /sci/labs/nissimb/guyhaim/icore-data/Xeno_transcriptome_alignment/feature_counts/10R_4c_p48_5_FU_2_5uM_48h_S9_L001_R1_001/10R_4c_p48_5_FU_2_5uM_48h_S9_L001_R1_001_counts.txt -s 2 /sci/labs/nissimb/guyhaim/icore-data/Xeno_transcriptome_alignment/star_out/10R_4c_p48_5_FU_2_5uM_48h_S9_L001_R1_001/10R_4c_p48_5_FU_2_5uM_48h_S9_L001_R1_001_Aligned.sortedByCoord.out.bam
