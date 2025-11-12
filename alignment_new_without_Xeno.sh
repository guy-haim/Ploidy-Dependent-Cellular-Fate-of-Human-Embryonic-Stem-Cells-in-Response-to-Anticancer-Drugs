#!/bin/bash
echo "please enter full path to the project dir:"
read PROJECT_DIR

echo "please enter full path to the fastq list file:"
read FASTQ_LIST

STAR_INDEX="/sci/labs/nissimb/guyhaim/icore-data/Xeno_transcriptome_alignment/genome_index"
STAR_MOUSE_INDEX="/sci/labs/nissimb/guyhaim/icore-data/Xeno_transcriptome_alignment/genome_index_mouse"
GTF="/sci/labs/nissimb/guyhaim/icore-data/Xeno_transcriptome_alignment/gencode.v38.primary_assembly.annotation.gtf"
TRUSEQ3="/sci/labs/nissimb/guyhaim/icore-data/Xeno_transcriptome_alignment/TruSeq3-SE.fa"
THREADS=8

module load trimmomatic R4 samtools STAR subread/2.0.1

mkdir ${PROJECT_DIR}/RunOutput
mkdir ${PROJECT_DIR}/fastq
mkdir ${PROJECT_DIR}/star_out
mkdir ${PROJECT_DIR}/bam
mkdir ${PROJECT_DIR}/star_mouse_out
mkdir ${PROJECT_DIR}/feature_counts
mkdir ${PROJECT_DIR}/trim

while read FASTQ
do
mkdir ${PROJECT_DIR}/RunOutput/${FASTQ}
mkdir ${PROJECT_DIR}/star_out/${FASTQ}
mkdir ${PROJECT_DIR}/bam/${FASTQ}
mkdir ${PROJECT_DIR}/star_mouse_out/${FASTQ}
mkdir ${PROJECT_DIR}/feature_counts/${FASTQ}
mkdir ${PROJECT_DIR}/trim/${FASTQ}

cat > xenofilterr.R << EOF

#!/usr/bin/env Rscript
library("XenofilteR")
args <- commandArgs(trailingOnly = TRUE)
sample.list <- data.frame(human = args[1], mouse=args[2])
bp.param <- SnowParam(workers = 4, type = "SOCK")
XenofilteR(sample.list, destination.folder = args[3], bp.param = bp.param)
EOF

cat > run_aligner.sh << EOF
#!/bin/bash



trimmomatic SE ${PROJECT_DIR}/fastq/${FASTQ}.fastq.gz ${PROJECT_DIR}/trim/${FASTQ}/${FASTQ}_trimmed.fastq ILLUMINACLIP:${TRUSEQ3}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


#!/bin/bash

# Human:
STAR\
 --runThreadN ${THREADS}\
 --genomeDir ${STAR_INDEX}\
 --readFilesIn ${PROJECT_DIR}/trim/${FASTQ}/${FASTQ}_trimmed.fastq\
 --outSAMtype BAM SortedByCoordinate\
 --twopassMode Basic\
 --outSAMattributes NM\
 --outFileNamePrefix ${PROJECT_DIR}/star_out/${FASTQ}/${FASTQ}_\


# Mouse:
#STAR\
# --runThreadN ${THREADS}\
# --genomeDir ${STAR_MOUSE_INDEX}\
# --readFilesIn ${PROJECT_DIR}/trim/${FASTQ}/${FASTQ}_trimmed.fastq\
# --outSAMtype BAM SortedByCoordinate\
# --twopassMode Basic\
# --outSAMattributes NM\
# --outFileNamePrefix ${PROJECT_DIR}/star_mouse_out/${FASTQ}/${FASTQ}_mouse_

human=${PROJECT_DIR}/star_out/${FASTQ}/${FASTQ}_Aligned.sortedByCoord.out.bam
mouse=${PROJECT_DIR}/star_mouse_out/${FASTQ}/${FASTQ}_mouse_Aligned.sortedByCoord.out.bam
destination=${PROJECT_DIR}/bam/${FASTQ}

Rscript --vanilla ./xenofilterr.R \${human} \${mouse} \${destination}

featureCounts\
 -T ${THREADS}\
 -a ${GTF}\
 -o ${PROJECT_DIR}/feature_counts/${FASTQ}/${FASTQ}_counts.txt\
 -s 2\
 ${PROJECT_DIR}/star_out/${FASTQ}/${FASTQ}_Aligned.sortedByCoord.out.bam
EOF

sbatch -c${THREADS} --mem=64G --time=4:0:0 -o ${PROJECT_DIR}/RunOutput/${FASTQ}/out ./run_aligner.sh
done < ${FASTQ_LIST}
