#!/bin/bash
SBATCH --time=5:0:0
SBATCH --mem=32G
SBATCH -c8
SBATCH -o index_out

mkdir /vol/sci/bio/data/nissim.benvenisty/guyhaim/genome_index \

STAR\
 --runThreadN 8 \
 --runMode genomeGenerate \
 --genomeDir /vol/sci/bio/data/nissim.benvenisty/guyhaim/genome_index \
 --genomeFastaFiles /vol/sci/bio/data/nissim.benvenisty/guyhaim/GRCh38.primary_assembly.genome.fa \
 --sjdbGTFfile /vol/sci/bio/data/nissim.benvenisty/guyhaim/gencode.v38.primary_assembly.annotation.gtf \
 --sjdbOverhang 100
