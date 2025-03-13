#!/bin/bash
# Goal: Generate coverage tracks, perform normalization, and convert formats to provide data for subsequent analysis and visualization
# Please ensure that samtools, MACS2, bedSort, and bedGraphToBigWig are in the PATH

# Set working directory and parameters
ref=~/Project_chip_seq/ref_genome # contain genome.size文件
raw_data=~/Project_chip_seq/raw_data  # Input file directory
OUTPUT_DIR=~/Project_chip_seq/Batch1      # Output file directory 
GENOME_INDEX=~/Project_chip_seq/ref_genome/yeast_index       # Bowtie2 index path
THREADS=8     



# Loop through each sample (based on deduplicated BAM files)
for dedup_bam in ${OUTPUT_DIR}/dedup/*.dedup.bam; do
    SAMPLE=$(basename "$dedup_bam" .dedup.bam)
    
    echo "Processing sample: ${SAMPLE}"
    
    # Calculate mapped reads count: Use samtools flagstat to extract the number of mapped reads from the deduplicated BAM file
    mappedReads=$(samtools flagstat "$dedup_bam" | grep " mapped (" | head -n1 | grep -Eo '^[0-9]+')
    
    if [[ -z "$mappedReads" ]]; then
       echo "Error: Unable to extract mapped reads count for ${SAMPLE}"
       continue
    fi
    echo "Mapped reads for ${SAMPLE}: ${mappedReads}"
    
    # Calculate normalization scale factor (based on million mapped reads)
    scale=$(perl -e "printf('%.3f', 1000000/$mappedReads)")
    echo "Normalization scale factor for ${SAMPLE}: $scale"
    
    # Generate coverage tracks (pileup), note that we use the deduplicated BAM file here
    macs2 pileup -f BAM --extsize 150 -i "$dedup_bam" -o "${OUTPUT_DIR}/coverage_tracks/${SAMPLE}_pileup.bdg"
    
    # Normalize using MACS2: Multiply the signal in the pileup file by the normalization factor
    echo "Normalizing ${SAMPLE}_pileup.bdg with factor $scale"
    macs2 bdgopt -i "${OUTPUT_DIR}/coverage_tracks/${SAMPLE}_pileup.bdg" -m multiply -p $scale -o "${OUTPUT_DIR}/coverage_tracks/temp_normalized.bdg"
    
    # Remove the first line of the normalized bedGraph file (if it contains extra headers or comments)
    sed -n '2,$p' "${OUTPUT_DIR}/coverage_tracks/temp_normalized.bdg" > "${OUTPUT_DIR}/coverage_tracks/${SAMPLE}_normalized.bdg"
    rm "${OUTPUT_DIR}/coverage_tracks/temp_normalized.bdg"
    
    # Convert the normalized bedGraph file to BigWig format for easier visualization in genome browsers
    bedSort "${OUTPUT_DIR}/coverage_tracks/${SAMPLE}_normalized.bdg" "${OUTPUT_DIR}/coverage_tracks/${SAMPLE}_normalized.bdg"
    bedGraphToBigWig "${OUTPUT_DIR}/coverage_tracks/${SAMPLE}_normalized.bdg" "$ref_genome/genome_sorted.size" "${OUTPUT_DIR}/coverage_tracks/${SAMPLE}_normalized.bw"
    
    echo "Normalization and BigWig conversion completed for ${SAMPLE}! Results saved in ${OUTPUT_DIR}/coverage_tracks/"
done
