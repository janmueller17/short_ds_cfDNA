#!/usr/bin/env bash
# Call nucleosome free regions for regular cfDNA-Seq data with nucleR and nucdyn

# Non-standard tools required (not including respective dependencies): 
# NucDyn - https://github.com/nucleosome-dynamics/nucleosome_dynamics

sample_path="/path/to/in/"
out_path_base="/path/to/out"
sample_names=("sample_x" "sample_y" "sample_z") # Add all your sample names here

for sample_name in "${sample_names[@]}"; do
  out_path="${out_path_base}/${sample_name}"
  if [ -d "${out_path}" ]; then
    echo " "
    echo "Output directory '${out_path}' exists. Rename folder to prevent overwriting."
    continue
  else
    echo "Processing ${sample_name} ..."
    mkdir -p "${out_path}" || { echo "Failed to create output directory '${out_path}'"; exit 1; }

    Rscript /path/to/readBAM.R --input "${sample_path}/${sample_name}.bam" --output "${out_path}/${sample_name}.RData" --type single || { echo "readBAM.R script failed"; exit 1; }
    Rscript /path/to/nucleR.R --input "${out_path}/${sample_name}.RData" --output "${out_path}/${sample_name}.nucleosome_calls.gff" --type single || { echo "nucleR.R script failed"; exit 1; }
    Rscript /path/to/NFR.R  --input "${out_path}/${sample_name}.nucleosome_calls.gff"  --output "${out_path}/${sample_name}.NFR.gff" || { echo "NFR.R script failed"; exit 1; }

    awk 'FS="\t"{print $1"\t"$4"\t"$5"\t"$8"\t"$6}' < "${out_path}/${sample_name}.nucleosome_calls.gff" > "${out_path}/${sample_name}.nucleosome_calls.bed" || { echo "AWK command for nucleosome_calls.bed failed"; exit 1; }
    awk 'FS="\t"{print $1"\t"$6"\t"$7}' < "${out_path}/${sample_name}.NFR.gff" > "${out_path}/${sample_name}.NFR.bed" || { echo "AWK command for NFR.bed failed"; exit 1; }
  fi
done
