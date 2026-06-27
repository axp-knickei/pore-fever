#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "$PROJECT_ROOT"

source config/config.sh

mkdir -p "$TAX_DIR" "$LOG_DIR"

# Set FORCE_RERUN=1 to rebuild the Sylph profile even when a complete output exists.
FORCE_RERUN="${FORCE_RERUN:-0}"

check_file() {
  local file="$1"
  if [[ ! -s "$file" ]]; then
    echo "ERROR: Missing or empty file: $file" >&2
    exit 1
  fi
}

check_gzip() {
  local file="$1"
  check_file "$file"
  gzip -t "$file"
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

check_sylph_db() {
  if [[ ! -s "$SYLPH_DB" ]]; then
    echo "ERROR: Sylph database is missing or empty: $SYLPH_DB" >&2
    echo "Edit SYLPH_DB in config/config.sh or place the database at the configured path." >&2
    exit 1
  fi
}

write_run_metadata() {
  local metadata_file="$1"
  local profile_file="$2"
  local manifest_file="$3"
  local sample_count="$4"
  shift 4
  local sylph_cmd=("$@")

  {
    echo "Sylph taxonomic profiling run metadata"
    echo "Generated on: $(date)"
    echo
    echo "Project directory:"
    echo "$PROJECT_DIR"
    echo
    echo "Sample metadata:"
    echo "$SAMPLES_TSV"
    echo
    echo "Sample count:"
    echo "$sample_count"
    echo
    echo "Input manifest:"
    echo "$manifest_file"
    echo
    echo "Sylph database:"
    echo "$SYLPH_DB"
    echo
    echo "Profile output:"
    echo "$profile_file"
    echo
    echo "Threads:"
    echo "$THREADS"
    echo
    echo "Sylph version:"
    sylph --version || true
    echo
    echo "Command:"
    printf '%q ' "${sylph_cmd[@]}"
    echo
    echo
    echo "Phase 7 read-level confirmation note:"
    echo "Optional follow-up after Phase 6. Candidate approaches include targeted read extraction plus minimap2 alignment to reference genomes, Kraken2/Centrifuge-style read classification, or BLAST/DIAMOND checks for selected taxa."
  } > "$metadata_file"
}

check_tool sylph
check_tool gzip
check_tool awk

validate_samples_tsv
check_sylph_db

MANIFEST="${TAX_DIR}/host_removed_fastq_manifest.tsv"
PROFILE_OUT="${TAX_DIR}/sylph_species_profile.tsv"
PROFILE_DONE="${PROFILE_OUT}.complete"
RUN_METADATA="${TAX_DIR}/sylph_run_metadata.txt"

FASTQ_FILES=()
TMP_MANIFEST="${MANIFEST}.tmp.$$"
cleanup_path "$TMP_MANIFEST"
echo -e "sample_id\tfastq" > "$TMP_MANIFEST"

while IFS=$'\t' read -r sample_id treatment timepoint replicate raw_fastq
do
  [[ -z "${sample_id:-}" || "${sample_id:0:1}" == "#" ]] && continue

  fq="${HOST_REMOVED_DIR}/${sample_id}.host_removed.fastq.gz"
  check_gzip "$fq"
  FASTQ_FILES+=("$fq")
  echo -e "${sample_id}\t${fq}" >> "$TMP_MANIFEST"
done < <(tail -n +2 "$SAMPLES_TSV")

if [[ "${#FASTQ_FILES[@]}" -eq 0 ]]; then
  echo "ERROR: No host-removed FASTQ files were found from samples.tsv." >&2
  exit 1
fi

mv "$TMP_MANIFEST" "$MANIFEST"

if [[ "$FORCE_RERUN" != "1" && -s "$PROFILE_OUT" && -f "$PROFILE_DONE" ]]; then
  echo "Phase 6: existing completed Sylph profile found; skipping rerun."
  echo "Use FORCE_RERUN=1 to rebuild: $PROFILE_OUT"
else
  echo "Phase 6: running sylph species profiling on ${#FASTQ_FILES[@]} samples"
  TMP_PROFILE="${PROFILE_OUT}.tmp.$$"
  TMP_METADATA="${RUN_METADATA}.tmp.$$"
  cleanup_path "$TMP_PROFILE"
  cleanup_path "$TMP_METADATA"
  cleanup_path "$PROFILE_DONE"

  SYLPH_CMD=(sylph profile "$SYLPH_DB" "${FASTQ_FILES[@]}" -t "$THREADS")
  "${SYLPH_CMD[@]}" > "$TMP_PROFILE"
  check_file "$TMP_PROFILE"

  write_run_metadata "$TMP_METADATA" "$PROFILE_OUT" "$MANIFEST" "${#FASTQ_FILES[@]}" "${SYLPH_CMD[@]}"
  mv "$TMP_PROFILE" "$PROFILE_OUT"
  mv "$TMP_METADATA" "$RUN_METADATA"
  touch "$PROFILE_DONE"
fi

check_file "$PROFILE_OUT"
check_file "$RUN_METADATA"

cat <<'NOTE'
Phase 7: read-level confirmation is optional after inspecting Phase 6.
Suggested follow-up for selected taxa: targeted read extraction plus minimap2, Kraken2/Centrifuge-style read classification, or BLAST/DIAMOND checks depending on the question.
Phase 8: prepare clean abundance matrices with scripts/python/03_prepare_taxonomy_table.py
Phases 6-8 initial taxonomic profiling completed.
NOTE
