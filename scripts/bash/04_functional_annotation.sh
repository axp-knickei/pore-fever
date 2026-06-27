#!/usr/bin/env bash
set -euo pipefail

source config/config.sh

mkdir -p "$FUNC_DIR" "$GENE_CAT_DIR" "$LOG_DIR"

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

check_tool bakta
check_tool mmseqs

DREP_GENOMES="${DREP_DIR}/drep_output/dereplicated_genomes"
ANNOT_DIR="${FUNC_DIR}/bakta_mag_annotations"
PROTEINS_ALL="${GENE_CAT_DIR}/all_predicted_proteins.faa"
CLUSTER_DB="${GENE_CAT_DIR}/mmseqs_gene_catalog"
REP_FASTA="${GENE_CAT_DIR}/nonredundant_gene_catalog.faa"

if [[ ! -d "$DREP_GENOMES" ]]; then
  echo "ERROR: Missing dereplicated MAG directory: $DREP_GENOMES" >&2
  exit 1
fi

mkdir -p "$ANNOT_DIR"

echo "Phase 18: annotating MAGs with Bakta"
for genome in "${DREP_GENOMES}"/*.fa
do
  check_file "$genome"

  base=$(basename "$genome" .fa)
  out="${ANNOT_DIR}/${base}"

  if [[ ! -d "$out" ]]; then
    bakta \
      --threads "$THREADS" \
      --output "$out" \
      --prefix "$base" \
      "$genome"
  fi
done

if [[ ! -s "$PROTEINS_ALL" ]]; then
  echo "Collecting predicted proteins..."
  cat "${ANNOT_DIR}"/*/*.faa > "$PROTEINS_ALL"
fi
check_file "$PROTEINS_ALL"

echo "Phase 19: run DRAM, eggNOG-mapper, HUMAnN, or KOfamScan for pathway annotation."

echo "Phase 20: building nonredundant gene catalog with MMseqs2"
if [[ ! -s "$REP_FASTA" ]]; then
  mmseqs easy-linclust \
    "$PROTEINS_ALL" \
    "$CLUSTER_DB" \
    "${GENE_CAT_DIR}/tmp" \
    --threads "$THREADS" \
    --min-seq-id 0.90 \
    -c 0.80

  cp "${CLUSTER_DB}_rep_seq.fasta" "$REP_FASTA"
fi
check_file "$REP_FASTA"

echo "Phase 21: quantify genes/functions by mapping reads to the nonredundant catalog."
echo "Phases 18-21 completed."
