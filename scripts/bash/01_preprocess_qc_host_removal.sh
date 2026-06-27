#!/usr/bin/env bash
set -euo pipefail

source config/config.sh

mkdir -p \
  "$TRIMMED_DIR" "$FILTERED_DIR" "$HOST_REMOVED_DIR" \
  "$QC_DIR/pre_filter" "$QC_DIR/post_filter" "$QC_DIR/host_removal" "$LOG_DIR"

check_file() {
  local file="$1"
  if [[ ! -s "$file" ]]; then
    echo "ERROR: Missing or empty file: $file" >&2
    exit 1
  fi
}

check_tool() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: Tool not found in PATH: $1" >&2
    exit 1
  }
}

check_tool NanoPlot
check_tool chopper
check_tool minimap2
check_tool samtools
check_tool gzip

check_file "$SAMPLES_TSV"
check_file "$MOUSE_REF"

if [[ ! -f "${MOUSE_REF}.mmi" ]]; then
  echo "Indexing mouse reference..."
  minimap2 -d "${MOUSE_REF}.mmi" "$MOUSE_REF"
fi

tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r sample_id treatment timepoint replicate raw_fastq
do
  echo "Processing sample: $sample_id"
  check_file "$raw_fastq"

  trimmed_fastq="${TRIMMED_DIR}/${sample_id}.trimmed.fastq.gz"
  filtered_fastq="${FILTERED_DIR}/${sample_id}.filtered.fastq.gz"
  host_removed_fastq="${HOST_REMOVED_DIR}/${sample_id}.host_removed.fastq.gz"
  host_bam="${QC_DIR}/host_removal/${sample_id}.mouse_mapping.bam"
  host_stats="${QC_DIR}/host_removal/${sample_id}.host_removal_stats.txt"

  # Placeholder for adapter/barcode trimming. Use Dorado trim here if needed.
  if [[ ! -s "$trimmed_fastq" ]]; then
    echo "Phase 1: copying raw FASTQ as trimmed placeholder for $sample_id"
    cp "$raw_fastq" "$trimmed_fastq"
  fi

  if [[ ! -d "${QC_DIR}/pre_filter/${sample_id}_NanoPlot" ]]; then
    echo "Phase 2: initial NanoPlot QC for $sample_id"
    NanoPlot \
      --fastq "$trimmed_fastq" \
      --outdir "${QC_DIR}/pre_filter/${sample_id}_NanoPlot" \
      --threads "$THREADS"
  fi

  if [[ ! -s "$filtered_fastq" ]]; then
    echo "Phase 3: filtering reads for $sample_id"
    gzip -cd "$trimmed_fastq" | \
      chopper \
        --minlength "$MIN_LENGTH" \
        --quality "$MIN_Q" \
        --threads "$THREADS" | \
      gzip > "$filtered_fastq"
  fi
  check_file "$filtered_fastq"

  if [[ ! -s "$host_removed_fastq" ]]; then
    echo "Phase 4: removing mouse host reads for $sample_id"

    minimap2 \
      -t "$THREADS" \
      -ax map-ont \
      "${MOUSE_REF}.mmi" \
      "$filtered_fastq" | \
      samtools view -@ "$THREADS" -b -o "$host_bam"

    samtools flagstat "$host_bam" > "$host_stats"

    samtools view -@ "$THREADS" -b -f 4 "$host_bam" | \
      samtools fastq -@ "$THREADS" - | \
      gzip > "$host_removed_fastq"
  fi
  check_file "$host_removed_fastq"

  if [[ ! -d "${QC_DIR}/post_filter/${sample_id}_NanoPlot" ]]; then
    echo "Phase 5: post-filtering QC for $sample_id"
    NanoPlot \
      --fastq "$host_removed_fastq" \
      --outdir "${QC_DIR}/post_filter/${sample_id}_NanoPlot" \
      --threads "$THREADS"
  fi
done

echo "Phases 1-5 completed."
