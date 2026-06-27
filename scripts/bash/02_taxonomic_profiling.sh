#!/usr/bin/env bash
set -euo pipefail

source config/config.sh

mkdir -p "$TAX_DIR" "$LOG_DIR"

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

check_tool sylph

check_file "$SAMPLES_TSV"
check_file "$SYLPH_DB"

MANIFEST="${TAX_DIR}/host_removed_fastq_manifest.tsv"
PROFILE_OUT="${TAX_DIR}/sylph_species_profile.tsv"

echo -e "sample_id\tfastq" > "$MANIFEST"

tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r sample_id treatment timepoint replicate raw_fastq
do
  fq="${HOST_REMOVED_DIR}/${sample_id}.host_removed.fastq.gz"
  check_file "$fq"
  echo -e "${sample_id}\t${fq}" >> "$MANIFEST"
done

echo "Phase 6: running sylph species profiling"
sylph profile \
  "$SYLPH_DB" \
  "${HOST_REMOVED_DIR}"/*.host_removed.fastq.gz \
  -t "$THREADS" \
  > "$PROFILE_OUT"

check_file "$PROFILE_OUT"

echo "Phase 7: read-level confirmation is optional after inspecting Phase 6."
echo "Phase 8: prepare clean abundance matrices with scripts/python/03_prepare_taxonomy_table.py"
echo "Phases 6-8 initial taxonomic profiling completed."
