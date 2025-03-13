#!/bin/bash
# This script prepares all the directories 
ref=~/Project_chip_seq/ref_genome
raw_data=~/Project_chip_seq/raw_data  # Input file directory
OUTPUT_DIR=~/Project_chip_seq/Batch1      # Output file directory 
GENOME_INDEX=~/Project_chip_seq/ref_genome/yeast_index       # Bowtie2 index path
THREADS=8      

# Create output directories
mkdir -p $ref
mkdir -p $raw_data
mkdir -p $GENOME_INDEX
mkdir -p ${OUTPUT_DIR}/fastp
mkdir -p ${OUTPUT_DIR}/fastqc
mkdir -p ${OUTPUT_DIR}/aligned
mkdir -p ${OUTPUT_DIR}/dedup
mkdir -p ${OUTPUT_DIR}/peak_calling
mkdir -p ${OUTPUT_DIR}/coverage_tracks

# Download Reference Genome
cd $ref
echo "Downloading reference genome pombe...."
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/945/GCF_000002945.2_ASM294v3/GCF_000002945.2_ASM294v3_genomic.fna.gz
gunzip GCF_000002945.2_ASM294v3_genomic.fna.gz
mv GCF_000002945.2_ASM294v3_genomic.fna yeast.fa

# build bowtie2 index file (prefix)
echo "Building reference genome index...."
bowtie2-build --threads 8 $ref/yeast.fa yeast_index 

# Generate genome.size file
echo "Build genome.size file...."
samtools faidx $ref/yeast.fa
cut -f1,2 yeast.fa.fai > genome.size
sort -t'_' -k3,3n genome.size > genome_sorted.size