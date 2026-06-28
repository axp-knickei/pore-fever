#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "$PROJECT_ROOT"

source config/config.sh

mkdir -p \
  "$ASM_DIR" "$MAP_DIR" "$BIN_DIR" "$MAG_QC_DIR" \
  "$MAG_TAX_DIR" "$DREP_DIR" "$MAG_ABUND_DIR" "$LOG_DIR"

# Set FORCE_RERUN=1 to rebuild outputs that have completion markers.
FORCE_RERUN="${FORCE_RERUN:-0}"
# Optional file with one sample_id per line for pilot assemblies or selected high-depth samples.
ASSEMBLY_SAMPLE_LIST="${ASSEMBLY_SAMPLE_LIST:-}"
# SemiBin2 environment preset. Keep configurable because this is mouse gut, not human gut.
SEMIBIN_ENVIRONMENT="${SEMIBIN_ENVIRONMENT:-human_gut}"

check_file() {
  local file="$1"
  if [[ ! -s "$file" ]]; then
    echo "ERROR: Missing or empty file: $file" >&2
    exit 1
  fi
}

check_dir_complete() {
  local dir="$1"
  [[ -d "$dir" && -f "${dir}/.complete" ]]
}

check_gzip() {
  local file="$1"
  check_file "$file"
  gzip -t "$file"
}

check_bam() {
  local bam="$1"
  check_file "$bam"
  samtools quickcheck "$bam"
}

check_tool() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: Tool not found in PATH: $1" >&2
    exit 1
  }
}

cleanup_path() {
  local path="$1"
  if [[ -e "$path" ]]; then
    rm -rf "$path"
  fi
}

require_complete_or_rerun() {
  local path="$1"
  local label="$2"
  if [[ -e "$path" && "$FORCE_RERUN" != "1" ]]; then
    echo "ERROR: Found existing incomplete ${label}: $path" >&2
    echo "Remove it manually after inspection, or rerun with FORCE_RERUN=1 to rebuild." >&2
    exit 1
  fi
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

load_selected_samples() {
  SELECTED_SAMPLE_IDS=()
  if [[ -z "$ASSEMBLY_SAMPLE_LIST" ]]; then
    return 0
  fi

  check_file "$ASSEMBLY_SAMPLE_LIST"
  while IFS= read -r sample_id
  do
    [[ -z "$sample_id" || "${sample_id:0:1}" == "#" ]] && continue
    SELECTED_SAMPLE_IDS+=("$sample_id")
  done < "$ASSEMBLY_SAMPLE_LIST"

  if [[ "${#SELECTED_SAMPLE_IDS[@]}" -eq 0 ]]; then
    echo "ERROR: ASSEMBLY_SAMPLE_LIST did not contain any sample IDs: $ASSEMBLY_SAMPLE_LIST" >&2
    exit 1
  fi
}

sample_is_selected() {
  local sample_id="$1"
  local selected
  if [[ "${#SELECTED_SAMPLE_IDS[@]}" -eq 0 ]]; then
    return 0
  fi
  for selected in "${SELECTED_SAMPLE_IDS[@]}"
  do
    [[ "$sample_id" == "$selected" ]] && return 0
  done
  return 1
}

build_input_arrays() {
  SAMPLE_IDS=()
  HOST_FASTQS=()

  while IFS=$'\t' read -r sample_id treatment timepoint replicate raw_fastq
  do
    [[ -z "${sample_id:-}" || "${sample_id:0:1}" == "#" ]] && continue
    sample_is_selected "$sample_id" || continue

    fq="${HOST_REMOVED_DIR}/${sample_id}.host_removed.fastq.gz"
    check_gzip "$fq"
    SAMPLE_IDS+=("$sample_id")
    HOST_FASTQS+=("$fq")
  done < <(tail -n +2 "$SAMPLES_TSV")

  if [[ "${#HOST_FASTQS[@]}" -eq 0 ]]; then
    echo "ERROR: No host-removed FASTQ files were selected from samples.tsv." >&2
    exit 1
  fi
}

write_run_metadata() {
  local metadata_file="$1"
  {
    echo "Assembly and MAG pipeline run metadata"
    echo "Generated on: $(date)"
    echo
    echo "Project directory:"
    echo "$PROJECT_DIR"
    echo
    echo "Sample metadata:"
    echo "$SAMPLES_TSV"
    echo
    echo "Selected sample list:"
    echo "${ASSEMBLY_SAMPLE_LIST:-all samples from samples.tsv}"
    echo
    echo "Selected sample count:"
    echo "${#SAMPLE_IDS[@]}"
    echo
    echo "Selected samples:"
    printf '%s\n' "${SAMPLE_IDS[@]}"
    echo
    echo "Threads:"
    echo "$THREADS"
    echo
    echo "SemiBin2 environment:"
    echo "$SEMIBIN_ENVIRONMENT"
    echo
    echo "Tool versions:"
    echo "flye: $(flye --version 2>/dev/null || true)"
    echo "minimap2: $(minimap2 --version 2>/dev/null || true)"
    echo "samtools:"
    samtools --version | head -n 2 || true
    echo "seqkit: $(seqkit version 2>/dev/null || true)"
    echo "coverm: $(coverm --version 2>/dev/null || true)"
    echo "semibin2: $(semibin2 --version 2>/dev/null || true)"
    echo "checkm2: $(checkm2 --version 2>/dev/null || true)"
    echo "gtdbtk: $(gtdbtk --version 2>/dev/null || true)"
    echo "dRep: $(dRep --version 2>/dev/null || true)"
  } > "$metadata_file"
}

normalize_bin_catalog() {
  local source_dir="$1"
  local output_dir="$2"
  local tmp_dir="${output_dir}.tmp.$$"
  local genome base target count

  cleanup_path "$tmp_dir"
  mkdir -p "$tmp_dir"

  count=0
  for genome in "$source_dir"/*.fa "$source_dir"/*.fasta "$source_dir"/*.fna
  do
    [[ -s "$genome" ]] || continue
    base="$(basename "$genome")"
    base="${base%.*}"
    target="${tmp_dir}/${base}.fa"
    ln -s "$(realpath "$genome")" "$target"
    count=$((count + 1))
  done

  if [[ "$count" -eq 0 ]]; then
    echo "ERROR: No MAG bin FASTA files found in $source_dir with extensions .fa, .fasta, or .fna" >&2
    exit 1
  fi

  touch "${tmp_dir}/.complete"
  cleanup_path "$output_dir"
  mv "$tmp_dir" "$output_dir"
}

check_genomes_dir() {
  local dir="$1"
  local genomes=("$dir"/*.fa)
  if [[ ! -d "$dir" || "${#genomes[@]}" -eq 0 ]]; then
    echo "ERROR: Missing genome directory or no .fa files found: $dir" >&2
    exit 1
  fi
  local genome
  for genome in "${genomes[@]}"
  do
    check_file "$genome"
  done
}

check_tool flye
check_tool minimap2
check_tool samtools
check_tool seqkit
check_tool coverm
check_tool semibin2
check_tool checkm2
check_tool gtdbtk
check_tool dRep
check_tool gzip
check_tool awk
check_tool realpath

validate_samples_tsv
load_selected_samples
build_input_arrays

COASSEMBLY_FASTQ="${ASM_DIR}/all_samples.host_removed.fastq.gz"
COASSEMBLY_DONE="${COASSEMBLY_FASTQ}.complete"
ASSEMBLY_DIR="${ASM_DIR}/metaflye_coassembly"
ASSEMBLY_FASTA="${ASSEMBLY_DIR}/assembly.fasta"
ASSEMBLY_STATS="${ASM_DIR}/assembly_seqkit_stats.tsv"
RUN_METADATA="${ASM_DIR}/assembly_mag_pipeline_run_metadata.txt"
SEMIBIN_OUT="${BIN_DIR}/semibin2_output"
SEMIBIN_BINS="${SEMIBIN_OUT}/output_bins"
STANDARD_BINS="${BIN_DIR}/semibin2_bins_fa"
CHECKM2_OUT="${MAG_QC_DIR}/checkm2"
GTDB_OUT="${MAG_TAX_DIR}/gtdbtk"
DREP_OUT="${DREP_DIR}/drep_output"
DREP_GENOMES="${DREP_OUT}/dereplicated_genomes"
MAG_ABUNDANCE="${MAG_ABUND_DIR}/mag_abundance.tsv"

write_run_metadata "$RUN_METADATA"

if [[ "$FORCE_RERUN" == "1" || ! -s "$COASSEMBLY_FASTQ" || ! -f "$COASSEMBLY_DONE" ]]; then
  require_complete_or_rerun "$COASSEMBLY_FASTQ" "co-assembly FASTQ"
  echo "Creating co-assembly FASTQ from ${#HOST_FASTQS[@]} validated samples..."
  tmp_coassembly="${COASSEMBLY_FASTQ}.tmp.$$"
  cleanup_path "$tmp_coassembly"
  cat "${HOST_FASTQS[@]}" > "$tmp_coassembly"
  check_gzip "$tmp_coassembly"
  mv "$tmp_coassembly" "$COASSEMBLY_FASTQ"
  touch "$COASSEMBLY_DONE"
fi
check_gzip "$COASSEMBLY_FASTQ"

if [[ "$FORCE_RERUN" == "1" || ! -s "$ASSEMBLY_FASTA" || ! -f "${ASSEMBLY_DIR}/.complete" ]]; then
  require_complete_or_rerun "$ASSEMBLY_DIR" "metaFlye assembly directory"
  echo "Phase 9: running metaFlye co-assembly..."
  tmp_assembly_dir="${ASSEMBLY_DIR}.tmp.$$"
  cleanup_path "$tmp_assembly_dir"
  flye \
    --nano-raw "$COASSEMBLY_FASTQ" \
    --meta \
    --threads "$THREADS" \
    --out-dir "$tmp_assembly_dir"
  check_file "${tmp_assembly_dir}/assembly.fasta"
  touch "${tmp_assembly_dir}/.complete"
  cleanup_path "$ASSEMBLY_DIR"
  mv "$tmp_assembly_dir" "$ASSEMBLY_DIR"
fi
check_file "$ASSEMBLY_FASTA"

if [[ "$FORCE_RERUN" == "1" || ! -s "$ASSEMBLY_STATS" ]]; then
  echo "Phase 10: assembly quality summary"
  tmp_stats="${ASSEMBLY_STATS}.tmp.$$"
  cleanup_path "$tmp_stats"
  seqkit stats "$ASSEMBLY_FASTA" > "$tmp_stats"
  check_file "$tmp_stats"
  mv "$tmp_stats" "$ASSEMBLY_STATS"
fi
check_file "$ASSEMBLY_STATS"

MAP_BAMS=()
echo "Phase 11: mapping reads back to contigs and calculating coverage"
for idx in "${!SAMPLE_IDS[@]}"
do
  sample_id="${SAMPLE_IDS[$idx]}"
  fq="${HOST_FASTQS[$idx]}"
  bam="${MAP_DIR}/${sample_id}.to_contigs.bam"
  bai="${bam}.bai"

  if [[ "$FORCE_RERUN" == "1" || ! -s "$bam" || ! -s "$bai" ]]; then
    tmp_bam="${bam}.tmp.$$"
    tmp_bai="${tmp_bam}.bai"
    cleanup_path "$tmp_bam"
    cleanup_path "$tmp_bai"
    minimap2 \
      -t "$THREADS" \
      -ax map-ont \
      "$ASSEMBLY_FASTA" \
      "$fq" | \
      samtools sort -@ "$THREADS" -o "$tmp_bam"
    check_bam "$tmp_bam"
    samtools index "$tmp_bam"
    check_file "$tmp_bai"
    mv "$tmp_bam" "$bam"
    mv "$tmp_bai" "$bai"
  fi
  check_bam "$bam"
  check_file "$bai"
  MAP_BAMS+=("$bam")
done

if [[ "${#MAP_BAMS[@]}" -eq 0 ]]; then
  echo "ERROR: No read-to-contig BAM files were created or selected." >&2
  exit 1
fi

if [[ "$FORCE_RERUN" == "1" || ! -d "$SEMIBIN_OUT" || ! -f "${SEMIBIN_OUT}/.complete" ]]; then
  require_complete_or_rerun "$SEMIBIN_OUT" "SemiBin2 output directory"
  echo "Phase 12: running SemiBin2 binning"
  tmp_semibin="${SEMIBIN_OUT}.tmp.$$"
  cleanup_path "$tmp_semibin"
  semibin2 single_easy_bin \
    --input-fasta "$ASSEMBLY_FASTA" \
    --input-bam "${MAP_BAMS[@]}" \
    --environment "$SEMIBIN_ENVIRONMENT" \
    --output "$tmp_semibin" \
    --threads "$THREADS"
  if [[ ! -d "${tmp_semibin}/output_bins" ]]; then
    echo "ERROR: SemiBin2 completed but output_bins directory is missing: ${tmp_semibin}/output_bins" >&2
    exit 1
  fi
  touch "${tmp_semibin}/.complete"
  cleanup_path "$SEMIBIN_OUT"
  mv "$tmp_semibin" "$SEMIBIN_OUT"
fi

if [[ "$FORCE_RERUN" == "1" || ! -d "$STANDARD_BINS" || ! -f "${STANDARD_BINS}/.complete" ]]; then
  echo "Phase 13: standardizing MAG bin FASTA extensions for downstream tools"
  normalize_bin_catalog "$SEMIBIN_BINS" "$STANDARD_BINS"
fi
check_genomes_dir "$STANDARD_BINS"

echo "Phase 13: MAG refinement placeholder. Add DAS Tool, MAGScoT, or metaWRAP here if multiple binning outputs are available."

if [[ "$FORCE_RERUN" == "1" || ! -d "$CHECKM2_OUT" || ! -f "${CHECKM2_OUT}/.complete" ]]; then
  require_complete_or_rerun "$CHECKM2_OUT" "CheckM2 output directory"
  echo "Phase 14: running CheckM2"
  tmp_checkm2="${CHECKM2_OUT}.tmp.$$"
  cleanup_path "$tmp_checkm2"
  checkm2 predict \
    --threads "$THREADS" \
    --input "$STANDARD_BINS" \
    --output-directory "$tmp_checkm2"
  touch "${tmp_checkm2}/.complete"
  cleanup_path "$CHECKM2_OUT"
  mv "$tmp_checkm2" "$CHECKM2_OUT"
fi

if [[ "$FORCE_RERUN" == "1" || ! -d "$GTDB_OUT" || ! -f "${GTDB_OUT}/.complete" ]]; then
  require_complete_or_rerun "$GTDB_OUT" "GTDB-Tk output directory"
  echo "Phase 15: running GTDB-Tk"
  tmp_gtdb="${GTDB_OUT}.tmp.$$"
  cleanup_path "$tmp_gtdb"
  gtdbtk classify_wf \
    --genome_dir "$STANDARD_BINS" \
    --out_dir "$tmp_gtdb" \
    --extension fa \
    --cpus "$THREADS"
  touch "${tmp_gtdb}/.complete"
  cleanup_path "$GTDB_OUT"
  mv "$tmp_gtdb" "$GTDB_OUT"
fi

if [[ "$FORCE_RERUN" == "1" || ! -d "$DREP_OUT" || ! -f "${DREP_OUT}/.complete" ]]; then
  require_complete_or_rerun "$DREP_OUT" "dRep output directory"
  echo "Phase 16: running dRep dereplication"
  tmp_drep="${DREP_OUT}.tmp.$$"
  cleanup_path "$tmp_drep"
  genomes=("$STANDARD_BINS"/*.fa)
  dRep dereplicate \
    "$tmp_drep" \
    -g "${genomes[@]}" \
    -p "$THREADS" \
    -comp 50 \
    -con 10
  check_genomes_dir "${tmp_drep}/dereplicated_genomes"
  touch "${tmp_drep}/.complete"
  cleanup_path "$DREP_OUT"
  mv "$tmp_drep" "$DREP_OUT"
fi
check_genomes_dir "$DREP_GENOMES"

if [[ "$FORCE_RERUN" == "1" || ! -s "$MAG_ABUNDANCE" ]]; then
  echo "Phase 17: estimating MAG abundance with CoverM"
  tmp_mag_abundance="${MAG_ABUNDANCE}.tmp.$$"
  cleanup_path "$tmp_mag_abundance"
  coverm genome \
    --genome-fasta-directory "$DREP_GENOMES" \
    --reads "${HOST_FASTQS[@]}" \
    --min-read-percent-identity 90 \
    --min-read-aligned-percent 50 \
    --methods mean covered_fraction relative_abundance \
    --threads "$THREADS" \
    > "$tmp_mag_abundance"
  check_file "$tmp_mag_abundance"
  mv "$tmp_mag_abundance" "$MAG_ABUNDANCE"
fi
check_file "$MAG_ABUNDANCE"

cat <<'NOTE'
Phases 9-17 completed.
Downstream input for scripts/bash/04_functional_annotation.sh:
  results/mag_dereplication/drep_output/dereplicated_genomes/*.fa
NOTE
