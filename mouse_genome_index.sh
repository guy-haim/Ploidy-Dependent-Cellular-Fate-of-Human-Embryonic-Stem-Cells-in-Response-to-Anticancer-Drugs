#!/bin/bash
#SBATCH --time=5:0:0
#SBATCH --mem=32G
#SBATCH -c8
#SBATCH -o index_out

# mkdir /vol/sci/bio/data/nissim.benvenisty/guyhaim/genome_index_mouse \

STAR\
 --runThreadN 8 \
 --runMode genomeGenerate \
 --genomeDir /vol/sci/bio/data/nissim.benvenisty/guyhaim/genome_index_mouse \
 --genomeFastaFiles /vol/sci/bio/data/nissim.benvenisty/guyhaim/GRCm38.primary_assembly.genome.fa \
 --sjdbGTFfile /vol/sci/bio/data/nissim.benvenisty/guyhaim/gencode.vM25.primary_assembly.annotation.gtf \
 --sjdbOverhang 100

