#!/usr/bin/env bash
# Process and map regular cfDNA seq data

# Non-standard tools required (not including respective dependencies): 
# FastQC - https://github.com/s-andrews/FastQC
# bbduk - https://github.com/BioInfoTools/BBMap
# NextGenMap (ngm) - https://github.com/Cibiv/NextGenMap
# samtools - https://github.com/samtools/samtools

sample_path="/path/to/in/"
out_path_base="/path/to/out/"
threads=4
genome_reference="/path/to/genome_reference/genome.fasta"
seq_adapter_reference="/path/to/adapter_reference/adapter.fasta"
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

    # Initial QC of raw data
    fastqc -o "${out_path}" --noextract -t "${threads}" "${sample_path}/${sample_name}_R1.fastq.gz" || { echo "Initial QC failed"; exit 1; }
    echo "Initial QC finished"

    # First read processing: artifact filtering
    /path/to/bbduk.sh hdist=1 mink=12 threads=${threads} maxns=10 minlen=50 trimq=20 qtrim=t ktrim=r k=28 ref="${seq_adapter_reference}" in="${sample_path}/${sample_name}_R1.fastq.gz" out="${out_path}/${sample_name}_R1_processed.fastq.gz" || { echo "BBDuk failed"; exit 1; } # polyG trimming, size filter, adapter trimming
    fastqc -o "${out_path}" --noextract -t "${threads}" "${out_path}/${sample_name}_R1_processed.fastq.gz" || { echo "QC after processing failed"; exit 1; } # get QC stats after processing
    echo "First read processing finished"

    # Read mapping to genome and basic processing steps (remove unmapped, sort, index, and flag duplicates)
    ngm -t "${threads}" --no-progress -b -o "${out_path}/${sample_name}.bam" -r "${genome_reference}" -q "${out_path}/${sample_name}_R1_processed.fastq.gz" 2> "${out_path}/${sample_name}_ngm_mapping_report.txt" || { echo "NGM mapping failed"; exit 1; }
    samtools view -@ "${threads}" -b -f 4 "${out_path}/${sample_name}.bam" > "${out_path}/${sample_name}_unmapped.bam" || { echo "Samtools view failed"; exit 1; }
    samtools sort -@ "${threads}" -o "${out_path}/${sample_name}_sorted.bam" "${out_path}/${sample_name}.bam" || { echo "Samtools sort failed"; exit 1; }
    samtools index -@ "${threads}" "${out_path}/${sample_name}_sorted.bam" || { echo "Samtools index failed"; exit 1; }
    samtools markdup -@ "${threads}" "${out_path}/${sample_name}_sorted.bam" "${out_path}/${sample_name}_dup_flagged.bam" || { echo "Samtools markdup flag failed"; exit 1; } # flag duplicates 
    samtools flagstat -@ "${threads}" "${out_path}/${sample_name}_dup_flagged.bam" > "${out_path}/${sample_name}_stats.txt" || { echo "Samtools flagstat failed"; exit 1; }
    echo "Mapping finished"
  fi
done


