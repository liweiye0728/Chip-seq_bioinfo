#!/bin/bash
# Set working directory and parameters
ref=~/Project_chip_seq/ref_genome
raw_data=~/Project_chip_seq/raw_data  # Input file directory
OUTPUT_DIR=~/Project_chip_seq/Batch1      # Output file directory 
GENOME_INDEX=~/Project_chip_seq/ref_genome/yeast_index       # Bowtie2 index path
THREADS=8     

# Require py3.8
# Here run multiqc in python>=3.7 environment
cd ${OUTPUT_DIR}/fastqc

multiqc --title "FastQC summary" -d -dd 1 --module fastqc --filename multiqc_report.html -o ${OUTPUT_DIR}/multiqc $PWD


R1_FILES=(${raw_data}/*.R1.fastq.gz)
# Loop through each sample

for R1_FILE in "${R1_FILES[@]}"; do
    SAMPLE=$(basename "$R1_FILE" .R1.fastq.gz)
    
    # Construct the corresponding R2 file name
    R2_FILE="${raw_data}/${SAMPLE}.R2.fastq.gz"
    
    # Check if R2 file exists
    if [[ ! -f "$R2_FILE" ]]; then
        echo "Error: R2 file for ${SAMPLE} not found!"
        continue
    fi
    # 4. Use Picard MarkDuplicates to remove duplicates
    # Add ReadGroup information to the sorted BAM file if it is missing
    picard AddOrReplaceReadGroups \
           I="${OUTPUT_DIR}/aligned/${SAMPLE}.sorted.bam" \
           O="${OUTPUT_DIR}/aligned/${SAMPLE}.sorted_with_readgroup.bam" \
           RGID=${SAMPLE} \
           RGLB=lib1 \
           RGPL=illumina \
           RGPU=unit1 \
           RGSM=${SAMPLE}

    # Now run MarkDuplicates on the new BAM file
  
    picard MarkDuplicates \
           I="${OUTPUT_DIR}/aligned/${SAMPLE}.sorted_with_readgroup.bam" \
           O="${OUTPUT_DIR}/dedup/${SAMPLE}.dedup.bam" \
           M="${OUTPUT_DIR}/dedup/${SAMPLE}_metrics.txt" \
           REMOVE_DUPLICATES=true
    
    # Check if Picard was successful
    if [[ $? -ne 0 ]]; then
        echo "Error: Picard MarkDuplicates failed for ${SAMPLE}"
        continue
    fi
done