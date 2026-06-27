#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "$PROJECT_ROOT"

source config/config.sh

mkdir -p \
  "$TRIMMED_DIR" "$FILTERED_DIR" "$HOST_REMOVED_DIR" \
  "$QC_DIR/pre_filter" "$QC_DIR/post_filter" "$QC_DIR/host_removal" "$LOG_DIR"

# Set TRIM_ADAPTERS=0 only when FASTQ files are already demultiplexed and adapter/barcode trimmed.
TRIM_ADAPTERS="${TRIM_ADAPTERS:-1}"
# Example: ADAPTER_TRIM_CMD='dorado trim --emit-fastq {input} | gzip -c > {output}'
# Example: ADAPTER_TRIM_CMD='porechop_abi -i {input} -o {output}'
ADAPTER_TRIM_CMD="${ADAPTER_TRIM_CMD:-}"

check_file() {
  local file="$1"
  if [[ ! -s "$file" ]]; then
    echo "ERROR: Missing or empty file: $file" >&2
    exit 1
  fi
}

check_dir_complete() {
  local dir="$1"
  if [[ ! -d "$dir" || ! -f "${dir}/.complete" ]]; then
    return 1
  fi
}

check_gzip() {
  local file="$1"
  check_file "$file"
  gzip -t "$file"
}

check_bam() {
  local file="$1"
  check_file "$file"
  samtools quickcheck "$file"
}

check_tool() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: Tool not found in PATH: $1" >&2
    exit 1
  }
}

resolve_path() {
  local path="$1"
  if [[ "$path" = /* ]]; then
    printf '%s\n' "$path"
  else
    printf '%s/%s\n' "$PROJECT_DIR" "$path"
  fi
}

cleanup_path() {
  local path="$1"
  if [[ -e "$path" ]]; then
    rm -rf "$path"
  fi
}

run_adapter_trim() {
  local input="$1"
  local output="$2"
  local tmp_output="${output}.tmp.$$"
  local cmd

  cleanup_path "$tmp_output"

  if [[ "$TRIM_ADAPTERS" == "0" ]]; then
    echo "Phase 1: FASTQ declared pre-trimmed; copying input for $sample_id"
    cp "$input" "$tmp_output"
  else
    if [[ -z "$ADAPTER_TRIM_CMD" ]]; then
      echo "ERROR: TRIM_ADAPTERS=1 but ADAPTER_TRIM_CMD is not set." >&2
      echo "Set ADAPTER_TRIM_CMD with {input} and {output} placeholders, or set TRIM_ADAPTERS=0 only for already-trimmed FASTQ." >&2
      exit 1
    fi

    echo "Phase 1: adapter/barcode trimming for $sample_id"
    cmd="${ADAPTER_TRIM_CMD//\{input\}/$input}"
    cmd="${cmd//\{output\}/$tmp_output}"
    bash -o pipefail -c "$cmd"
  fi

  check_gzip "$tmp_output"
  mv "$tmp_output" "$output"
}

run_nanoplot() {
  local fastq="$1"
  local outdir="$2"
  local tmp_outdir="${outdir}.tmp.$$"

  cleanup_path "$tmp_outdir"
  NanoPlot \
    --fastq "$fastq" \
    --outdir "$tmp_outdir" \
    --threads "$THREADS"
  touch "${tmp_outdir}/.complete"
  cleanup_path "$outdir"
  mv "$tmp_outdir" "$outdir"
}

validate_samples_tsv() {
  check_file "$SAMPLES_TSV"
  awk -F '\t' '
    NR == 1 {
      expected = "sample_id\ttreatment\ttimepoint\treplicate\traw_fastq"
      if ($0 != expected) {
        printf "ERROR: samples.tsv header must be: %s\n", expected > "/dev/stderr"
        exit 1
      }
      next
    }
    /^#/ || NF == 0 { next }
    NF != 5 {
      printf "ERROR: samples.tsv line %d has %d columns; expected 5\n", NR, NF > "/dev/stderr"
      exit 1
    }
    $1 == "" || $2 == "" || $3 == "" || $4 == "" || $5 == "" {
      printf "ERROR: samples.tsv line %d contains an empty required field\n", NR > "/dev/stderr"
      exit 1
    }
  ' "$SAMPLES_TSV"
}

check_tool NanoPlot
check_tool chopper
check_tool minimap2
check_tool samtools
check_tool gzip
check_tool awk

if [[ "$TRIM_ADAPTERS" != "0" ]]; then
  if [[ -z "$ADAPTER_TRIM_CMD" ]]; then
    echo "ERROR: adapter trimming is enabled but ADAPTER_TRIM_CMD is not set." >&2
    echo "Example: ADAPTER_TRIM_CMD='dorado trim --emit-fastq {input} | gzip -c > {output}'" >&2
    echo "Use TRIM_ADAPTERS=0 only if raw_fastq entries are already trimmed." >&2
    exit 1
  fi
fi

validate_samples_tsv
check_file "$MOUSE_REF"

if [[ ! -f "${MOUSE_REF}.mmi" ]]; then
  echo "Indexing mouse reference..."
  tmp_index="${MOUSE_REF}.mmi.tmp.$$"
  cleanup_path "$tmp_index"
  minimap2 -d "$tmp_index" "$MOUSE_REF"
  mv "$tmp_index" "${MOUSE_REF}.mmi"
fi

while IFS=$'\t' read -r sample_id treatment timepoint replicate raw_fastq
do
  [[ -z "${sample_id:-}" || "${sample_id:0:1}" == "#" ]] && continue

  echo "Processing sample: $sample_id"
  raw_fastq_path="$(resolve_path "$raw_fastq")"
  check_gzip "$raw_fastq_path"

  trimmed_fastq="${TRIMMED_DIR}/${sample_id}.trimmed.fastq.gz"
  filtered_fastq="${FILTERED_DIR}/${sample_id}.filtered.fastq.gz"
  host_removed_fastq="${HOST_REMOVED_DIR}/${sample_id}.host_removed.fastq.gz"
  host_bam="${QC_DIR}/host_removal/${sample_id}.mouse_mapping.bam"
  host_stats="${QC_DIR}/host_removal/${sample_id}.host_removal_stats.txt"

  if [[ ! -s "$trimmed_fastq" ]]; then
    run_adapter_trim "$raw_fastq_path" "$trimmed_fastq"
  fi
  check_gzip "$trimmed_fastq"

  pre_filter_dir="${QC_DIR}/pre_filter/${sample_id}_NanoPlot"
  if ! check_dir_complete "$pre_filter_dir"; then
    echo "Phase 2: initial NanoPlot QC for $sample_id"
    run_nanoplot "$trimmed_fastq" "$pre_filter_dir"
  fi

  if [[ ! -s "$filtered_fastq" ]]; then
    echo "Phase 3: filtering reads for $sample_id"
    tmp_filtered="${filtered_fastq}.tmp.$$"
    cleanup_path "$tmp_filtered"
    gzip -cd "$trimmed_fastq" | \
      chopper \
        --minlength "$MIN_LENGTH" \
        --quality "$MIN_Q" \
        --threads "$THREADS" | \
      gzip > "$tmp_filtered"
    check_gzip "$tmp_filtered"
    mv "$tmp_filtered" "$filtered_fastq"
  fi
  check_gzip "$filtered_fastq"

  if [[ ! -s "$host_removed_fastq" ]]; then
    echo "Phase 4: removing mouse host reads for $sample_id"
    tmp_host_bam="${host_bam}.tmp.$$"
    tmp_host_removed="${host_removed_fastq}.tmp.$$"
    tmp_host_stats="${host_stats}.tmp.$$"
    cleanup_path "$tmp_host_bam"
    cleanup_path "$tmp_host_removed"
    cleanup_path "$tmp_host_stats"

    minimap2 \
      -t "$THREADS" \
      -ax map-ont \
      "${MOUSE_REF}.mmi" \
      "$filtered_fastq" | \
      samtools view -@ "$THREADS" -b -o "$tmp_host_bam"

    check_bam "$tmp_host_bam"
    samtools flagstat "$tmp_host_bam" > "$tmp_host_stats"
    check_file "$tmp_host_stats"

    samtools view -@ "$THREADS" -b -f 4 "$tmp_host_bam" | \
      samtools fastq -@ "$THREADS" - | \
      gzip > "$tmp_host_removed"
    check_gzip "$tmp_host_removed"

    mv "$tmp_host_bam" "$host_bam"
    mv "$tmp_host_stats" "$host_stats"
    mv "$tmp_host_removed" "$host_removed_fastq"
  fi
  check_gzip "$host_removed_fastq"
  check_bam "$host_bam"
  check_file "$host_stats"

  post_filter_dir="${QC_DIR}/post_filter/${sample_id}_NanoPlot"
  if ! check_dir_complete "$post_filter_dir"; then
    echo "Phase 5: post-filtering QC for $sample_id"
    run_nanoplot "$host_removed_fastq" "$post_filter_dir"
  fi
done < <(tail -n +2 "$SAMPLES_TSV")

echo "Phases 1-5 completed."
