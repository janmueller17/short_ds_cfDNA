#!/usr/bin/env bash
# Call transcription factor motif enrichments with AME from MEME suite

# Non-standard tools required (not including respective dependencies): 
# AME - https://meme-suite.org/meme/doc/ame.html

sample_path="/path/to/in"
tf_motif_reference="/path/to/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme"

# Function to run AME with appropriate checks
run_ame() {
    local out_path="$1"
    local fasta_file="$2"
    local control_file="$3"
    
    if [ -d "${out_path}" ]; then
        echo " "
        echo "Output directory '${out_path}' exists. Rename folder to prevent overwriting."
    else
        mkdir -p "${out_path}" || { echo "Failed to create output directory '${out_path}'"; exit 1; }
        ame --scoring avg --method fisher --control "${control_file}" -o "${out_path}" "${fasta_file}" "${tf_motif_reference}" || { echo "AME for ${fasta_file} with ${control_file} as control failed"; exit 1; }
    fi
}

# Enrichment in one healthy sample short cfDNA
run_ame "${sample_path}/healthy_single/" "${sample_path}/Healthy_S03_short_cfDNA_macs2_narrow_peaks.fasta" "--shuffle--"

# Differential enrichment for two conditions
# PDAC vs CRC
run_ame "${sample_path}/differential/CRC_vs_PDAC/CRC/" "${sample_path}/CRC_consensus_narrow_peaks.fasta" "${sample_path}/PDAC_consensus_narrow_peaks.fasta"
run_ame "${sample_path}/differential/CRC_vs_PDAC/PDAC/" "${sample_path}/PDAC_consensus_narrow_peaks.fasta" "${sample_path}/CRC_consensus_narrow_peaks.fasta"

# Sepsis vs Post-OP
run_ame "${sample_path}/differential/Sepsis_vs_PostOP/Sepsis/" "${sample_path}/Sepsis_consensus_narrow_peaks.fasta" "${sample_path}/PostOP_consensus_narrow_peaks.fasta"
run_ame "${sample_path}/differential/Sepsis_vs_PostOP/PostOP/" "${sample_path}/PostOP_consensus_narrow_peaks.fasta" "${sample_path}/Sepsis_consensus_narrow_peaks.fasta"
