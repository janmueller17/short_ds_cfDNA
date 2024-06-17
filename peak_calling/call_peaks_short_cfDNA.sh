#!/usr/bin/env bash
# Call peaks for short cfDNA Seq data with macs2

# Non-standard tools required (not including respective dependencies): 
# MACS2 - https://pypi.org/project/MACS2/

sample_path="/path/to/in/"
out_path_base="/path/to/out/"
sample_names=("sample_x" "sample_y" "sample_z")

for sample_name in "${sample_names[@]}"; do
  out_path="${out_path_base}/${sample_name}"
  if [ -d "${out_path}" ]; then
    echo " "
    echo "Output directory '${out_path}' exists. Rename folder to prevent overwriting."
    continue
  else
    echo "Processing ${sample_name} ..."
    mkdir -p "${out_path}" || { echo "Failed to create output directory '${out_path}'"; exit 1; }

    macs2 callpeak -t "${sample_path}/${sample_name}.bam" --outdir "${out_path}" 2> "${out_path}/${sample_name}_narrow_peakcall.log" -f BAM -g 2864785220 -n "${sample_name}_macs2_narrow" --nomodel --extsize 32 --call-summits --min-length 30 -q 0.05 || { echo "MACS2 narrow peak calling failed"; exit 1; }

    macs2 callpeak -t "${sample_path}/${sample_name}.bam" --outdir "${out_path}" 2> "${out_path}/${sample_name}_broad_peakcall.log" -f BAM -g 2864785220 -n "${sample_name}_macs2_broad" --nomodel --extsize 32 --broad --max-gap 100 --min-length 500 --broad-cutoff 0.05 || { echo "MACS2 broad peak calling failed"; exit 1; }
  fi
done
