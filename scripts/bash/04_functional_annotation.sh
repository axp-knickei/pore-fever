#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "$PROJECT_ROOT"

source config/config.sh

mkdir -p "$FUNC_DIR" "$GENE_CAT_DIR" "$LOG_DIR"

# Set FORCE_RERUN=1 to rebuild Bakta annotations, combined proteins, and MMseqs outputs.
FORCE_RERUN="${FORCE_RERUN:-0}"
# Optional Bakta database path. If unset, Bakta must be configured in the active environment.
BAKTA_DB="${BAKTA_DB:-}"
MMSEQS_MIN_SEQ_ID="${MMSEQS_MIN_SEQ_ID:-0.90}"
MMSEQS_COVERAGE="${MMSEQS_COVERAGE:-0.80}"

check_file() {
  local file="$1"
  if [[ ! -s "$file" ]]; then
    echo "ERROR: Missing or empty file: $file" >&2
    exit 1
  fi
}

check_fasta() {
  local file="$1"
  check_file "$file"
  if ! grep -q '^>' "$file"; then
    echo "ERROR: FASTA file has no header lines: $file" >&2
    exit 1
  fi
}

check_dir_complete() {
  local dir="$1"
  [[ -d "$dir" && -f "${dir}/.complete" ]]
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

check_bakta_database_hint() {
  if [[ -n "$BAKTA_DB" ]]; then
    if [[ ! -d "$BAKTA_DB" ]]; then
      echo "ERROR: BAKTA_DB is set but is not a directory: $BAKTA_DB" >&2
      exit 1
    fi
  else
    echo "WARNING: BAKTA_DB is not set. Bakta must find its database from the active environment/configuration." >&2
  fi
}

load_mag_genomes() {
  MAG_GENOMES=("$DREP_GENOMES"/*.fa)
  if [[ ! -d "$DREP_GENOMES" || "${#MAG_GENOMES[@]}" -eq 0 ]]; then
    echo "ERROR: Missing dereplicated MAG input files: ${DREP_GENOMES}/*.fa" >&2
    echo "Run scripts/bash/03_assembly_mag_pipeline.sh through dRep dereplication first." >&2
    exit 1
  fi

  local genome
  for genome in "${MAG_GENOMES[@]}"
  do
    check_fasta "$genome"
  done
}

find_bakta_protein_file() {
  local outdir="$1"
  local base="$2"
  local expected="${outdir}/${base}.faa"
  local proteins

  if [[ -s "$expected" ]]; then
    printf '%s\n' "$expected"
    return 0
  fi

  proteins=("$outdir"/*.faa)
  if [[ "${#proteins[@]}" -eq 1 && -s "${proteins[0]}" ]]; then
    printf '%s\n' "${proteins[0]}"
    return 0
  fi

  echo "ERROR: Could not find exactly one Bakta protein FASTA in $outdir" >&2
  return 1
}

write_run_metadata() {
  local metadata_file="$1"
  {
    echo "Functional annotation run metadata"
    echo "Generated on: $(date)"
    echo
    echo "Project directory:"
    echo "$PROJECT_DIR"
    echo
    echo "Input dereplicated MAG directory:"
    echo "$DREP_GENOMES"
    echo
    echo "Input MAG count:"
    echo "${#MAG_GENOMES[@]}"
    echo
    echo "Bakta annotation directory:"
    echo "$ANNOT_DIR"
    echo
    echo "Combined predicted proteins:"
    echo "$PROTEINS_ALL"
    echo
    echo "Nonredundant gene catalog:"
    echo "$REP_FASTA"
    echo
    echo "Bakta database:"
    echo "${BAKTA_DB:-environment/default configuration}"
    echo
    echo "MMseqs clustering thresholds:"
    echo "min_seq_id=${MMSEQS_MIN_SEQ_ID}"
    echo "coverage=${MMSEQS_COVERAGE}"
    echo
    echo "Tool versions:"
    echo "bakta: $(bakta --version 2>/dev/null || true)"
    echo "mmseqs: $(mmseqs version 2>/dev/null || true)"
  } > "$metadata_file"
}

write_placeholder_notes() {
  local pathway_note="${FUNC_DIR}/pathway_annotation_NOTE.txt"
  local abundance_note="${FUNC_DIR}/function_abundance_matrix_NOTE.txt"

  cat > "$pathway_note" <<NOTE
Phase 19 placeholder: pathway-level functional annotation has not been implemented in this scaffold.
Recommended next tools include DRAM, eggNOG-mapper, KOfamScan, HUMAnN, or another method selected for the study question and database availability.
Input candidates produced by this script:
- ${PROTEINS_ALL}
- ${REP_FASTA}
NOTE

  cat > "$abundance_note" <<NOTE
Phase 21 placeholder: function abundance quantification has not been implemented in this scaffold.
This script does not create results/functional_annotation/function_abundance_matrix.tsv.
Downstream R scripts will skip differential functional abundance until a real function abundance matrix is generated.
Recommended next step: map sample reads or genes to the nonredundant gene catalog and aggregate by gene family/pathway.
NOTE
}

check_tool bakta
check_tool mmseqs
check_tool grep
check_tool cat
check_tool awk

DREP_GENOMES="${DREP_DIR}/drep_output/dereplicated_genomes"
ANNOT_DIR="${FUNC_DIR}/bakta_mag_annotations"
PROTEINS_ALL="${GENE_CAT_DIR}/all_predicted_proteins.faa"
REP_FASTA="${GENE_CAT_DIR}/nonredundant_gene_catalog.faa"
RUN_METADATA="${FUNC_DIR}/functional_annotation_run_metadata.txt"

check_bakta_database_hint
load_mag_genomes
write_run_metadata "$RUN_METADATA"

mkdir -p "$ANNOT_DIR"

BAKTA_DB_ARGS=()
if [[ -n "$BAKTA_DB" ]]; then
  BAKTA_DB_ARGS=(--db "$BAKTA_DB")
fi

PROTEIN_FASTAS=()

echo "Phase 18: annotating ${#MAG_GENOMES[@]} MAGs with Bakta"
for genome in "${MAG_GENOMES[@]}"
do
  base="$(basename "$genome" .fa)"
  out="${ANNOT_DIR}/${base}"

  if [[ "$FORCE_RERUN" == "1" || ! -d "$out" || ! -f "${out}/.complete" ]]; then
    require_complete_or_rerun "$out" "Bakta annotation directory"
    tmp_out="${out}.tmp.$$"
    cleanup_path "$tmp_out"
    bakta \
      --threads "$THREADS" \
      --output "$tmp_out" \
      --prefix "$base" \
      "${BAKTA_DB_ARGS[@]}" \
      "$genome"
    protein_file="$(find_bakta_protein_file "$tmp_out" "$base")"
    check_fasta "$protein_file"
    touch "${tmp_out}/.complete"
    cleanup_path "$out"
    mv "$tmp_out" "$out"
  fi

  if ! check_dir_complete "$out"; then
    echo "ERROR: Bakta annotation directory is incomplete: $out" >&2
    exit 1
  fi

  protein_file="$(find_bakta_protein_file "$out" "$base")"
  check_fasta "$protein_file"
  PROTEIN_FASTAS+=("$protein_file")
done

if [[ "${#PROTEIN_FASTAS[@]}" -eq 0 ]]; then
  echo "ERROR: No Bakta protein FASTA files were found from completed annotations." >&2
  exit 1
fi

if [[ "$FORCE_RERUN" == "1" || ! -s "$PROTEINS_ALL" ]]; then
  echo "Collecting predicted proteins from completed Bakta annotations..."
  tmp_proteins="${PROTEINS_ALL}.tmp.$$"
  cleanup_path "$tmp_proteins"
  cat "${PROTEIN_FASTAS[@]}" > "$tmp_proteins"
  check_fasta "$tmp_proteins"
  mv "$tmp_proteins" "$PROTEINS_ALL"
fi
check_fasta "$PROTEINS_ALL"

write_placeholder_notes

echo "Phase 19: pathway-level annotation placeholder note written to ${FUNC_DIR}/pathway_annotation_NOTE.txt"

echo "Phase 20: building nonredundant gene catalog with MMseqs2"
if [[ "$FORCE_RERUN" == "1" || ! -s "$REP_FASTA" ]]; then
  tmp_run_dir="${GENE_CAT_DIR}/mmseqs_gene_catalog.tmp.$$"
  tmp_rep="${REP_FASTA}.tmp.$$"
  cleanup_path "$tmp_run_dir"
  cleanup_path "$tmp_rep"
  mkdir -p "$tmp_run_dir"

  cluster_prefix="${tmp_run_dir}/mmseqs_gene_catalog"
  mmseqs easy-linclust \
    "$PROTEINS_ALL" \
    "$cluster_prefix" \
    "${tmp_run_dir}/tmp" \
    --threads "$THREADS" \
    --min-seq-id "$MMSEQS_MIN_SEQ_ID" \
    -c "$MMSEQS_COVERAGE"

  check_fasta "${cluster_prefix}_rep_seq.fasta"
  cp "${cluster_prefix}_rep_seq.fasta" "$tmp_rep"
  check_fasta "$tmp_rep"
  mv "$tmp_rep" "$REP_FASTA"
  cleanup_path "$tmp_run_dir"
fi
check_fasta "$REP_FASTA"

echo "Phase 21: function quantification placeholder note written to ${FUNC_DIR}/function_abundance_matrix_NOTE.txt"
echo "Phases 18-21 completed."
