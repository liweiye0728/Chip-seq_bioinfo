#!/bin/bash
set -e
# Goal: Perform peak calling and visualization
# Modification Notes:
# 1. Adjust input and output directories to be consistent with the normalization script (e.g., processed_data/dedup)
# 2. Set -g parameter to Schizosaccharomyces pombe genome size (1.26e7, 12.6 Mb)
# 3. Visualization prioritizes normalized BigWig files, with added BigWig generation step
# 4. Integrate MACS2 parameters optimized for 150bp paired-end reads and high duplication rate
# 5. Add BAM indexing and fix parameter conflicts (e.g., remove -p, correct --skipZeros)
# 6. Add logging for better tracking
# 7. Update sample-control matching to ignore trailing suffix differences

# Directory settings
OUTPUT_DIR=~/Project_chip_seq/Batch1/dedup      # Output file directory 
RESULT_DIR=~/Project_chip_seq/Batch1/peak_calling
COVERAGE_DIR=~/Project_chip_seq/Batch1/coverage_tracks
LOG_FILE="${RESULT_DIR}/peak_calling_$(date +%F_%H-%M-%S).log"

# # Create result directories (if they do not exist)
# mkdir -p "$RESULT_DIR" "$COVERAGE_DIR"

# Record start time
echo "Script started at $(date)" > "$LOG_FILE"

# Peak calling
echo "Starting peak calling..." | tee -a "$LOG_FILE"
for SAMPLE in "$OUTPUT_DIR"*.dedup.bam; do
    SAMPLE_NAME=$(basename "$SAMPLE" .dedup.bam)
    
    # 仅对含 "ChIp" 的文件进行峰调用，跳过对照文件本身
    if [[ "$SAMPLE_NAME" != *ChIp* ]]; then
        continue
    fi

    # 提取中间的数字 (如 3315、447、636 等)
    # 假设文件名格式类似：Fxl10-3315-ChIp_L7_Q0258W0196
    # 通过 sed 正则，匹配第一个 "-" 后的数字。
    # 注意：如有更多/更复杂的破折号，需要酌情调整。
    MIDDLE=$(echo "$SAMPLE_NAME" | sed -E 's/^[^-]+-([0-9]+).*/\1/')

    # 根据这段数字构造对照文件的匹配模式
    # 例如: "*-3315-in_*.dedup.bam"
    CONTROL_PATTERN="*-${MIDDLE}-in_*.dedup.bam"
    CONTROL_BAM=$(ls "$OUTPUT_DIR"/$CONTROL_PATTERN 2>/dev/null | head -n 1)

    if [[ -f "$CONTROL_BAM" ]]; then
        echo "Sample: $SAMPLE_NAME  -->  Control: $(basename "$CONTROL_BAM")"
        # 在此处调用 MACS2
        macs2 callpeak \
            -t "$SAMPLE" \
            -c "$CONTROL_BAM" \
            -f BAMPE -g 1.26e7 -q 0.01 \
            --nomodel --shift -75 --extsize 150 \
            --call-summits -B --SPMR \
            -n "${SAMPLE_NAME}_peaks" \
            --outdir "$RESULT_DIR"
    else
        echo "Warning: No control found for $SAMPLE_NAME (pattern: $CONTROL_PATTERN)"
    fi
done

# # Visualization - Coverage plots
# echo "Generating coverage plots..." | tee -a "$LOG_FILE"
# if [ -n "$(find "$COVERAGE_DIR" -name '*_normalized.bw' -print -quit)" ]; then
#     plotCoverage -b "$COVERAGE_DIR"/*_normalized.bw -o "${RESULT_DIR}/coverage_plots.pdf" \
#                  --skipZeros 2>&1 | tee -a "$LOG_FILE"
# else
#     echo "Normalized BigWig files not found in ${COVERAGE_DIR}, using dedup BAM files..." | tee -a "$LOG_FILE"
#     plotCoverage -b "$OUTPUT_DIR"/*.dedup.bam -o "${RESULT_DIR}/coverage_plots.pdf" \
#                  --skipZeros 2>&1 | tee -a "$LOG_FILE"
# fi


# 在您展示的代码中，computeMatrix和plotHeatmap生成的热图：

# **横轴**：代表每个峰的中心点(referencePoint center)前后各500bp区域(-b 500 -a 500)

# **纵轴**：代表每个找到的峰(narrowPeak文件中的峰)

# **参数解释**：
# - `reference-point`：以峰的中心点为参考
# - `--referencePoint center`：选择峰的中心作为参考点
# - `-b 500 -a 500`：分析中心点前后各500bp
# - `-R *_peaks.narrowPeak`：输入峰文件
# - `-S *_normalized.bw`：输入信号文件(BigWig)
# - `--skipZeros`：跳过没有信号的区域
# - `--sortRegions descend`：按信号强度降序排列
# - `--colorMap RdBu`：使用红蓝色图

# 热图显示了在每个峰附近区域的ChIP-seq信号分布模式。

# Visualization - Heatmap generation via computeMatrix and plotHeatmap
if [ -n "$(find "$RESULT_DIR" -name '*_peaks.narrowPeak' -print -quit)" ]; then
    echo "Computing matrix for heatmap..." | tee -a "$LOG_FILE"
    computeMatrix reference-point --referencePoint center \
                  -b 500 -a 500 \
                  -R "$RESULT_DIR"/*_peaks.narrowPeak \
                  -S "$COVERAGE_DIR"/*_normalized.bw \
                  --skipZeros -o "${RESULT_DIR}/matrix.gz" 2>&1 | tee -a "$LOG_FILE"
    
    echo "Generating heatmap..." | tee -a "$LOG_FILE"
    plotHeatmap -m "${RESULT_DIR}/matrix.gz" -out "${RESULT_DIR}/heatmap.pdf" \
                --sortRegions descend --colorMap RdBu 2>&1 | tee -a "$LOG_FILE"
else
    echo "No peak files found in ${RESULT_DIR}, skipping heatmap generation." | tee -a "$LOG_FILE"
fi

echo "Peak calling and visualization completed at $(date)! Results are saved in ${RESULT_DIR}" | tee -a "$LOG_FILE"
