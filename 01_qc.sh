#!/bin/bash
# Environment: conda, chip_seq and chip_3.12(for picard)
# This script collects chip seq raw data and filter low quality Q<20, adapter 
# Set working directory and parameters
ref=~/Project_chip_seq/ref_genome
raw_data=~/Project_chip_seq/raw_data  # Input file directory
OUTPUT_DIR=~/Project_chip_seq/Batch1      # Output file directory 
GENOME_INDEX=~/Project_chip_seq/ref_genome/yeast_index       # Bowtie2 index path
THREADS=8     

cd $raw_data
fastqc -t 4 --outdir ${OUTPUT_DIR}/fastqc *.fastq.gz

# Get a list of all R1 files (assuming files end with .R1.fastq.gz)
R1_FILES=(${raw_data}/*.R1.fastq.gz)

# Loop through each sample
for R1_FILE in "${R1_FILES[@]}"; do
    # Extract sample name from the R1 file name
    # Example: Fxl1-447-in_L7_Q0061W0195.R1.fastq.gz -> Fxl1-447-in_L7_Q0061W0195
    SAMPLE=$(basename "$R1_FILE" .R1.fastq.gz)
    
    # Construct the corresponding R2 file name
    R2_FILE="${raw_data}/${SAMPLE}.R2.fastq.gz"
    
    # Check if R2 file exists
    if [[ ! -f "$R2_FILE" ]]; then
        echo "Error: R2 file for ${SAMPLE} not found!"
        continue
    fi
    
    echo "Processing sample: ${SAMPLE}"
    
    # 1. Use Fastp for quality trimming and adapter removal
    fastp -i "$R1_FILE" -I "$R2_FILE" \
          -o "${OUTPUT_DIR}/fastp/${SAMPLE}.trimmed.R1.fastq.gz" \
          -O "${OUTPUT_DIR}/fastp/${SAMPLE}.trimmed.R2.fastq.gz" \
          -q 30 -u 30 --length_required 35 -5 -3 \
          --detect_adapter_for_pe \
          --low_complexity_filter \
          -w ${THREADS} \
          -j "${OUTPUT_DIR}/fastp/${SAMPLE}_fastp.json" \
          -h "${OUTPUT_DIR}/fastp/${SAMPLE}_fastp.html" 
    
    # Check if Fastp was successful
    if [[ $? -ne 0 ]]; then
        echo "Error: Fastp failed for ${SAMPLE}"
        continue
    fi
    
    # # 2. Use Bowtie2 for alignment
    # bowtie2 -x ${GENOME_INDEX} \
    #         -1 "${OUTPUT_DIR}/fastp/${SAMPLE}.trimmed.R1.fastq.gz" \
    #         -2 "${OUTPUT_DIR}/fastp/${SAMPLE}.trimmed.R2.fastq.gz" \
    #         -p ${THREADS} \
    #         -S "${OUTPUT_DIR}/aligned/${SAMPLE}.sam"
    
    # # Check if Bowtie2 alignment was successful
    # if [[ $? -ne 0 ]]; then
    #     echo "Error: Bowtie2 alignment failed for ${SAMPLE}"
    #     continue
    # fi
    
    # # 3. Use SAMtools for conversion and sorting
    # samtools view -bS "${OUTPUT_DIR}/aligned/${SAMPLE}.sam" > "${OUTPUT_DIR}/aligned/${SAMPLE}.bam"
    # samtools sort "${OUTPUT_DIR}/aligned/${SAMPLE}.bam" -o "${OUTPUT_DIR}/aligned/${SAMPLE}.sorted.bam" -@ ${THREADS}
    
    # # Delete temporary SAM file to save space
    # rm "${OUTPUT_DIR}/aligned/${SAMPLE}.sam"
    
    # # Check if SAMtools processing was successful
    # if [[ $? -ne 0 ]]; then
    #     echo "Error: SAMtools processing failed for ${SAMPLE}"
    #     continue
    # fi
    
    # 4. Use Picard MarkDuplicates to remove duplicates
    # picard MarkDuplicates \
    #        I="${OUTPUT_DIR}/aligned/${SAMPLE}.sorted.bam" \
    #        O="${OUTPUT_DIR}/dedup/${SAMPLE}.dedup.bam" \
    #        M="${OUTPUT_DIR}/dedup/${SAMPLE}_metrics.txt" \
    #        REMOVE_DUPLICATES=true
    
    # # Check if Picard was successful
    # if [[ $? -ne 0 ]]; then
    #     echo "Error: Picard MarkDuplicates failed for ${SAMPLE}"
    #     continue
    # fi
    
    echo "Finished processing sample: ${SAMPLE}"
done

echo "All samples processed!"
