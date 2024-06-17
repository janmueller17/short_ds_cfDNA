#!/usr/bin/env bash
# Process and map short cfDNA seq data
# This version keeps all temp files (greppable with '_temp'), these can be removed if not needed.

# Non-standard tools required (not including respective dependencies): 
# FastQC - https://github.com/s-andrews/FastQC
# bbduk - https://github.com/BioInfoTools/BBMap
# prinseq-lite - https://github.com/uwb-linux/prinseq
# NextGenMap (ngm) - https://github.com/Cibiv/NextGenMap
# samtools - https://github.com/samtools/samtools
# deepTools - https://github.com/deeptools/deepTools

sample_path="/path/to/in/"
out_path_base="/path/to/out/"
threads=4
bbduk_path="/path/to/bbduk.sh"
prinseq_path="/path/to/prinseq-lite.pl"
genome_reference="/path/to/genome_reference/genome.fasta"
genome_reference_eff_size=2864785220 # hg19, see deepTools doc for exemplary values
seq_adapter_reference="/path/to/adapter_reference/adapter.fasta"
blacklist_bed="/path/to/ENCODE_consensusBlacklist.bed"
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

    # First read processing: artifact filtering (terminal A, polyG, etc.)
    "${bbduk_path}" trimpolyg=10 hdist=1 mink=12 threads=${threads} maxns=10 minlen=20 maxlength=60 trimq=20 qtrim=t ktrim=r k=28 ref="${seq_adapter_reference}" in="${sample_path}/${sample_name}_R1.fastq.gz" out="${out_path}/${sample_name}_R1_processed_temp1.fastq.gz" || { echo "BBDuk step 1 failed"; exit 1; } # polyG trimming, size filter, adapter trimming
    "${bbduk_path}" forcetrimright2=1 threads=${threads} ref="${seq_adapter_reference}" in="${out_path}/${sample_name}_R1_processed_temp1.fastq.gz" out="${out_path}/${sample_name}_R1_processed_temp2.fastq.gz" || { echo "BBDuk step 2 failed"; exit 1; } # 1 base right side trim to remove single A from library prep
    
    gunzip "${out_path}/${sample_name}_R1_processed_temp2.fastq.gz" || { echo "gunzip failed"; exit 1; }
    "${prinseq_path}" -lc_method dust -lc_threshold 7 -fastq "${out_path}/${sample_name}_R1_processed_temp2.fastq" -out_good "${out_path}/${sample_name}_R1_processed" -out_bad "${out_path}/${sample_name}_R1_processed_bad" || { echo "PRINSEQ-lite processing failed"; exit 1; } # low complexity filtering
    gzip "${out_path}/${sample_name}_R1_processed.fastq" || { echo "gzip failed"; exit 1; }
    
    fastqc -o "${out_path}" --noextract -t "${threads}" "${out_path}/${sample_name}_R1_processed.fastq.gz" || { echo "QC after processing failed"; exit 1; } # get QC stats before processing
    echo "First read processing finished"

    # Read mapping to genome and basic processing steps (remove unmapped, sort, and index)
    ngm -t "${threads}" --no-progress -b -o "${out_path}/${sample_name}.bam" -r "${genome_reference}" -q "${out_path}/${sample_name}_R1_processed.fastq.gz" 2> "${out_path}/${sample_name}_ngm_mapping_report.txt" || { echo "NGM mapping failed"; exit 1; }
    samtools view -@ "${threads}" -b -f 4 "${out_path}/${sample_name}.bam" > "${out_path}/${sample_name}_unmapped.bam" || { echo "Samtools view failed"; exit 1; }
    samtools sort -@ "${threads}" -o "${out_path}/${sample_name}_sorted.bam" "${out_path}/${sample_name}.bam" || { echo "Samtools sort failed"; exit 1; }
    samtools index -@ "${threads}" "${out_path}/${sample_name}_sorted.bam" || { echo "Samtools index failed"; exit 1; }
    samtools flagstat -@ "${threads}" "${out_path}/${sample_name}_sorted.bam" > "${out_path}/${sample_name}_stats_before_processing.txt" || { echo "Samtools flagstat before processing failed"; exit 1; }
    echo "Mapping finished"

    # Further processing
    samtools markdup -@ "${threads}" "${out_path}/${sample_name}_sorted.bam" "${out_path}/${sample_name}_dup_flagged.bam" || { echo "Samtools markdup flag failed"; exit 1; } # flag duplicates 
    samtools markdup -r -s -@ "${threads}" "${out_path}/${sample_name}_sorted.bam" "${out_path}/${sample_name}_dup_rem_temp.bam" || { echo "Samtools markdup remove failed"; exit 1; } # remove duplicates 
    samtools view -H "${out_path}/${sample_name}_dup_rem_temp.bam" | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - "${out_path}/${sample_name}_dup_rem_temp.bam" > "${out_path}/${sample_name}_chr_fix_temp.bam" || { echo "Reheadering failed"; exit 1; } # filter for autosomal and sex chromosomes only (chr 1 to 22 X and Y) and change chromosome names to UCSC style (1 -> chr1)
    samtools view "${out_path}/${sample_name}_chr_fix_temp.bam" -@ "${threads}" -o "${out_path}/${sample_name}_reads_in_bl_temp.bam" -U "${out_path}/${sample_name}_bl_rem_temp.bam" -L "${blacklist_bed}" || { echo "Blacklist removal failed"; exit 1; } # remove reads mapped to the ENCODE blacklist
    samtools view "${out_path}/${sample_name}_bl_rem_temp.bam" -@ "${threads}" -q 1 -o "${out_path}/${sample_name}_processed.bam" || { echo "Samtools view with mapQ filter failed"; exit 1; } # keep only reads with mapQ > 1
    samtools index -@ "${threads}" "${out_path}/${sample_name}_processed.bam" || { echo "Samtools index after processing failed"; exit 1; }
    samtools flagstat -@ "${threads}" "${out_path}/${sample_name}_processed.bam" > "${out_path}/${sample_name}_stats_after_processing.txt" || { echo "Samtools flagstat after processing failed"; exit 1; } # get QC stats after processing
    bamCoverage --bam "${out_path}/${sample_name}_processed.bam" -o "${out_path}/${sample_name}_processed.bw" --outFileFormat bigwig --normalizeUsing CPM --binSize 10 --numberOfProcessors "${threads}" --effectiveGenomeSize "${genome_reference_eff_size}" || { echo "bamCoverage failed"; exit 1; } # generate bigwig from processed reads
    echo "Processing finished"
  fi

done
