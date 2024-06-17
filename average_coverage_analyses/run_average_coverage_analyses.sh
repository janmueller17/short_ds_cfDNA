#!/usr/bin/env bash
# Generate matrix, average line plots, and heatmaps for various genomic regions in the three styles used in the manuscript

# Non-standard tools required (not including respective dependencies): 
# deepTools - https://github.com/deeptools/deepTools

healthy_short_cfDNA_bigwig="/path/to/Healthy_S03_processed.bw"
healthy_regular_cfDNA_bigwig="/path/to/Healthy_S06_processed.bw"
healthy_cfMBD_bigwig="/path/to/Healthy_S25_processed.bw"
references_path="/path/to/references/"
out_path="/path/to/out/"
threads=4

# Function to run computeMatrix
run_compute_matrix() {
    local matrix_type="$1"
    shift
    local options=("$@")
    
    computeMatrix "${matrix_type}" "${options[@]}" || { echo "computeMatrix for ${matrix_type} failed"; exit 1; }
}

# Function to run plotProfile
run_plot_profile() {
    local matrix_file="$1"
    local output_file="$2"
    shift 2
    local options=("$@")
    
    plotProfile -m "${matrix_file}" -out "${output_file}" "${options[@]}" || { echo "plotProfile for ${matrix_file} failed"; exit 1; }
}

# Function to run plotHeatmap
run_plot_heatmap() {
    local matrix_file="$1"
    local output_file="$2"
    shift 2
    local options=("$@")
    
    plotHeatmap -m "${matrix_file}" -out "${output_file}" "${options[@]}" || { echo "plotHeatmap for ${matrix_file} failed"; exit 1; }
}

# TSS RefSeq
run_compute_matrix reference-point \
    --referencePoint center \
    -S "${healthy_short_cfDNA_bigwig}" "${healthy_regular_cfDNA_bigwig}" \
    --samplesLabel "Short cfDNA" "Regular cfDNA" \
    -R "${references_path}/refseq_hg19_TSS.bed" \
    --beforeRegionStartLength 1500 --afterRegionStartLength 1500 --skipZeros \
    -o "${out_path}/matrix_hg19-genes-TSS.mat.gz" \
    --numberOfProcessors "${threads}"

run_plot_profile \
    "${out_path}/matrix_hg19-genes-TSS.mat.gz" \
    "${out_path}/average_profile_plot_hg19-genes-TSS.pdf" \
    --perGroup --colors green purple --plotFileFormat pdf --plotType=lines --refPointLabel "TSS" \
    --plotWidth 6 --plotHeight 6 --yAxisLabel "CPM" --plotTitle "" --regionsLabel ""

# Full gene RefSeq
run_compute_matrix scale-regions \
    -S "${healthy_short_cfDNA_bigwig}" "${healthy_regular_cfDNA_bigwig}" \
    --samplesLabel "Short cfDNA" "Regular cfDNA" \
    -R "${references_path}/refseq_hg19_genes.bed" \
    --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
    -o "${out_path}/matrix_hg19-genes.mat.gz" \
    --numberOfProcessors "${threads}"

run_plot_profile \
    "${out_path}/matrix_hg19-genes.mat.gz" \
    "${out_path}/average_profile_plot_hg19-genes.pdf" \
    --perGroup --colors green purple --plotFileFormat pdf --plotType=lines --endLabel "TTS" --yAxisLabel "CPM" \
    --plotTitle "" --regionsLabel "" --labelRotation 90 --plotWidth 10 --plotHeight 7

# All ChIP TFBS from ENCODE3
run_compute_matrix reference-point \
    --referencePoint center \
    -S "${healthy_short_cfDNA_bigwig}" "${healthy_regular_cfDNA_bigwig}" \
    --samplesLabel "Short cfDNA" "Regular cfDNA" \
    -R "${references_path}/Txn_Factr_ChIP_E3.bed3" \
    --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --skipZeros \
    -o "${out_path}/matrix_all_ENCODE_TFBS.mat.gz" \
    --numberOfProcessors "${threads}" --missingDataAsZero

run_plot_profile \
    "${out_path}/matrix_all_ENCODE_TFBS.mat.gz" \
    "${out_path}/average_profile_plot_all_ENCODE_TFBS.pdf" \
    --colors green purple --plotWidth 5 --plotHeight 6 --labelRotation 90 \
    --yMin 0 --yMax 0.075 --samplesLabel "Short cfDNA" "Regular cfDNA"

# CpG islands with MBD-seq
run_compute_matrix scale-regions \
    -R "${references_path}/cpgIslandExt.hg19.chrnr.bed" \
    -S "${healthy_cfMBD_bigwig}" "${healthy_short_cfDNA_bigwig}" "${healthy_regular_cfDNA_bigwig}" \
    -b 1000 -a 1000 --regionBodyLength 500 \
    -o "${out_path}/matrix_CpGs.mat.gz" \
    --samplesLabel "cfMBD" "Short cfDNA" "Regular cfDNA" --numberOfProcessors "${threads}"

run_plot_heatmap \
    "${out_path}/matrix_CpGs.mat.gz" \
    "${out_path}/heatmap_plot_CpGs.pdf" \
    --clusterUsingSamples 1 --kmeans 2 --zMin 0 0 0 --zMax 2 0.15 0.05 \
    --yMin 0 0 0 --yMax 2 0.15 0.05 --endLabel "CpG end" --startLabel "CpG start" \
    --colorMap "inferno" --heatmapHeight 10 --whatToShow "heatmap and colorbar" --missingDataColor 0 --xAxisLabel "" --labelRotation 90

run_plot_profile \
    "${out_path}/matrix_CpGs.mat.gz" \
    "${out_path}/average_profile_plot_CpGs.pdf" \
    --clusterUsingSamples 1 --kmeans 2 --yMin 0 0 0 --yMax 2 0.15 0.05 --endLabel "CpG end" --startLabel "CpG start" \
    --labelRotation 90 --plotWidth 5 --plotHeight 7 --colors red blue
