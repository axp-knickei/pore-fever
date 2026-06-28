#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONFIG_FILE="${PROJECT_DIR}/config/config.sh"

if [[ ! -s "$CONFIG_FILE" ]]; then
  echo "ERROR: Missing or empty config file: $CONFIG_FILE" >&2
  exit 1
fi

cd "$PROJECT_DIR"
# shellcheck source=/dev/null
source "$CONFIG_FILE"
PROJECT_DIR="$(cd "$PROJECT_DIR" && pwd)"

REPORT="${PROJECT_DIR}/results/reproducibility_report.txt"
MANIFEST="${PROJECT_DIR}/results/reproducibility_manifest.tsv"
mkdir -p "${PROJECT_DIR}/results"

timestamp_utc="$(date -u '+%Y-%m-%dT%H:%M:%SZ')"

echo -e "section\titem\tstatus\tpath\tdetails" > "$MANIFEST"

append_manifest() {
  local section="$1"
  local item="$2"
  local status="$3"
  local path="$4"
  local details="$5"
  details="${details//$'\t'/ }"
  details="${details//$'\n'/; }"
  printf '%s\t%s\t%s\t%s\t%s\n' "$section" "$item" "$status" "$path" "$details" >> "$MANIFEST"
}

file_status() {
  local path="$1"
  if [[ -f "$path" && -s "$path" ]]; then
    echo "PRESENT"
  elif [[ -f "$path" ]]; then
    echo "EMPTY"
  elif [[ -d "$path" && -f "${path}/.complete" ]]; then
    echo "COMPLETE_DIR"
  elif [[ -d "$path" ]]; then
    echo "DIR_PRESENT"
  else
    echo "MISSING"
  fi
}

path_details() {
  local path="$1"
  if [[ -f "$path" ]]; then
    printf 'size_bytes=%s' "$(wc -c < "$path")"
  elif [[ -d "$path" ]]; then
    local count
    count="$(find "$path" -mindepth 1 -maxdepth 1 2>/dev/null | wc -l | tr -d ' ')"
    printf 'entries=%s' "$count"
  else
    printf 'not_found'
  fi
}

record_path() {
  local section="$1"
  local item="$2"
  local path="$3"
  local required="$4"
  local status
  local details
  status="$(file_status "$path")"
  details="$(path_details "$path"); required=${required}"
  append_manifest "$section" "$item" "$status" "$path" "$details"
  printf '%-46s %s\n' "${section} | ${item}:" "$status"
  printf '  path: %s\n' "$path"
  printf '  details: %s\n' "$details"
}

record_glob_count() {
  local section="$1"
  local item="$2"
  local pattern="$3"
  local required="$4"
  local count
  count="$(find "$(dirname "$pattern")" -maxdepth 1 -name "$(basename "$pattern")" 2>/dev/null | wc -l | tr -d ' ')"
  local status="PRESENT"
  if [[ "$count" == "0" ]]; then
    status="MISSING"
  fi
  append_manifest "$section" "$item" "$status" "$pattern" "matches=${count}; required=${required}"
  printf '%-46s %s\n' "${section} | ${item}:" "$status"
  printf '  path: %s\n' "$pattern"
  printf '  details: matches=%s; required=%s\n' "$count" "$required"
}

record_tool() {
  local label="$1"
  local executable="$2"
  shift 2
  printf '%s\n' "$label"
  if ! command -v "$executable" >/dev/null 2>&1; then
    echo "  status: MISSING"
    append_manifest "tool" "$label" "MISSING" "$executable" "command_not_found"
    echo
    return 0
  fi

  local exe_path
  exe_path="$(command -v "$executable")"
  echo "  status: FOUND"
  echo "  path: $exe_path"
  local output
  output="$({ "$executable" "$@"; } 2>&1 || true)"
  if [[ -z "$output" ]]; then
    output="version command produced no output"
  fi
  echo "$output" | sed -n '1,4p' | sed 's/^/  /'
  append_manifest "tool" "$label" "FOUND" "$exe_path" "$(echo "$output" | sed -n '1p')"
  echo
}

metadata_summary() {
  if [[ ! -s "$SAMPLES_TSV" ]]; then
    echo "Metadata summary: MISSING"
    append_manifest "input" "sample metadata validation" "MISSING" "$SAMPLES_TSV" "cannot_validate"
    return 0
  fi

  local header
  header="$(head -n 1 "$SAMPLES_TSV")"
  local required_cols=(sample_id treatment timepoint replicate raw_fastq)
  local missing=()
  local col
  for col in "${required_cols[@]}"; do
    if ! printf '%s\n' "$header" | tr '\t' '\n' | grep -Fxq "$col"; then
      missing+=("$col")
    fi
  done

  local sample_count
  sample_count="$(awk 'NR > 1 && $0 !~ /^#/ && NF > 0 {n++} END {print n+0}' "$SAMPLES_TSV")"
  echo "Metadata summary:"
  echo "  path: $SAMPLES_TSV"
  echo "  samples: $sample_count"
  if [[ "${#missing[@]}" -eq 0 ]]; then
    echo "  required columns: OK"
    append_manifest "input" "sample metadata validation" "PRESENT" "$SAMPLES_TSV" "samples=${sample_count}; required_columns=OK"
  else
    echo "  required columns: MISSING ${missing[*]}"
    append_manifest "input" "sample metadata validation" "INCOMPLETE" "$SAMPLES_TSV" "samples=${sample_count}; missing_columns=${missing[*]}"
  fi
}

git_summary() {
  echo "Git state:"
  if ! command -v git >/dev/null 2>&1 || ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "  git repository: unavailable"
    append_manifest "git" "repository" "MISSING" "$PROJECT_DIR" "git_not_available_or_not_a_repo"
    return 0
  fi

  local branch commit remote dirty recent
  branch="$(git branch --show-current 2>/dev/null || true)"
  commit="$(git rev-parse HEAD 2>/dev/null || true)"
  remote="$(git remote get-url origin 2>/dev/null || true)"
  if [[ -n "$(git status --porcelain 2>/dev/null || true)" ]]; then
    dirty="DIRTY"
  else
    dirty="CLEAN"
  fi
  echo "  branch: ${branch:-unknown}"
  echo "  commit: ${commit:-unknown}"
  echo "  remote: ${remote:-unknown}"
  echo "  working tree: $dirty"
  append_manifest "git" "current commit" "$dirty" "$PROJECT_DIR" "branch=${branch:-unknown}; commit=${commit:-unknown}; remote=${remote:-unknown}"
  echo "  recent commits:"
  recent="$(git log --oneline -5 2>/dev/null || true)"
  if [[ -n "$recent" ]]; then
    echo "$recent" | sed 's/^/    /'
  else
    echo "    unavailable"
  fi
}

write_metadata_snippet() {
  local title="$1"
  local path="$2"
  echo "$title"
  if [[ -s "$path" ]]; then
    sed -n '1,80p' "$path" | sed 's/^/  /'
  else
    echo "  MISSING: $path"
  fi
  echo
}

{
  echo "ONT Mouse Gut Metagenomics Reproducibility Report"
  echo "Generated UTC: $timestamp_utc"
  echo

  echo "Project and configuration"
  echo "Project directory: $PROJECT_DIR"
  echo "Config file: $CONFIG_FILE"
  echo "Sample metadata: $SAMPLES_TSV"
  echo "Mouse reference: $MOUSE_REF"
  echo "Sylph database: $SYLPH_DB"
  echo "THREADS=$THREADS"
  echo "MIN_LENGTH=$MIN_LENGTH"
  echo "MIN_Q=$MIN_Q"
  echo

  metadata_summary
  echo
  git_summary
  echo

  echo "Tool versions and availability"
  record_tool "minimap2" "minimap2" --version
  record_tool "samtools" "samtools" --version
  record_tool "NanoPlot" "NanoPlot" --version
  record_tool "chopper" "chopper" --version
  record_tool "sylph" "sylph" --version
  record_tool "flye" "flye" --version
  record_tool "coverm" "coverm" --version
  record_tool "SemiBin2" "SemiBin2" --version
  record_tool "checkm2" "checkm2" --version
  record_tool "gtdbtk" "gtdbtk" --version
  record_tool "dRep" "dRep" --version
  record_tool "bakta" "bakta" --version
  record_tool "seqkit" "seqkit" version
  record_tool "R" "R" --version
  record_tool "Rscript" "Rscript" --version
  record_tool "python3" "python3" --version
  record_tool "python" "python" --version
  record_tool "conda" "conda" --version
  record_tool "mamba" "mamba" --version

  echo "Input and database files"
  record_path "input" "config.sh" "$CONFIG_FILE" "yes"
  record_path "input" "samples.tsv" "$SAMPLES_TSV" "yes"
  record_path "input" "mouse reference FASTA" "$MOUSE_REF" "yes for host removal"
  record_path "input" "mouse reference minimap2 index" "${MOUSE_REF}.mmi" "recommended"
  record_path "input" "Sylph database" "$SYLPH_DB" "yes for taxonomic profiling"
  echo

  echo "Pipeline artifact checklist"
  record_path "phase_1_5" "host-removed FASTQ directory" "$HOST_REMOVED_DIR" "yes after preprocessing"
  record_path "phase_1_5" "QC directory" "$QC_DIR" "recommended"
  record_path "phase_1_5" "host-removed FASTQ QC summary" "${QC_DIR}/host_removed_fastq_qc_summary.tsv" "recommended after 05_summarize_fastq_qc.py"
  record_path "phase_6" "Sylph profile" "${TAX_DIR}/sylph_species_profile.tsv" "yes for taxonomy table"
  record_path "phase_6" "Sylph completion marker" "${TAX_DIR}/sylph_species_profile.tsv.complete" "recommended"
  record_path "phase_6" "Sylph FASTQ manifest" "${TAX_DIR}/host_removed_fastq_manifest.tsv" "recommended"
  record_path "phase_6" "Sylph run metadata" "${TAX_DIR}/sylph_run_metadata.txt" "recommended"
  record_path "phase_8" "taxa abundance matrix" "${TAX_DIR}/taxa_abundance_matrix.tsv" "yes for statistics/reporting"
  record_path "phase_8" "taxonomy preparation summary" "${TAX_DIR}/taxonomy_table_preparation_summary.txt" "recommended"
  record_path "phase_9_17" "assembly run metadata" "${ASM_DIR}/assembly_mag_pipeline_run_metadata.txt" "recommended"
  record_path "phase_9" "coassembly FASTQ" "${ASM_DIR}/all_samples.host_removed.fastq.gz" "recommended"
  record_path "phase_10" "metaFlye assembly FASTA" "${ASM_DIR}/metaflye_coassembly/assembly.fasta" "recommended"
  record_path "phase_10" "assembly seqkit stats" "${ASM_DIR}/assembly_seqkit_stats.tsv" "recommended"
  record_path "phase_12" "SemiBin2 output directory" "${BIN_DIR}/semibin2_output" "recommended"
  record_path "phase_13" "standardized bins directory" "${BIN_DIR}/semibin2_bins_fa" "recommended"
  record_path "phase_14" "CheckM2 output directory" "${MAG_QC_DIR}/checkm2" "recommended"
  record_path "phase_15" "GTDB-Tk output directory" "${MAG_TAX_DIR}/gtdbtk" "recommended"
  record_path "phase_16" "dRep output directory" "${DREP_DIR}/drep_output" "recommended"
  record_glob_count "phase_16" "dereplicated MAG FASTAs" "${DREP_DIR}/drep_output/dereplicated_genomes/*.fa" "recommended for functional annotation"
  record_path "phase_18_21" "functional run metadata" "${FUNC_DIR}/functional_annotation_run_metadata.txt" "recommended"
  record_path "phase_18_21" "function abundance placeholder note" "${FUNC_DIR}/function_abundance_matrix_NOTE.txt" "expected until function quantification exists"
  record_path "phase_22_24" "statistics run summary" "${STAT_DIR}/statistics_run_summary.txt" "recommended"
  record_path "phase_22" "alpha diversity table" "${STAT_DIR}/alpha_diversity_taxa.tsv" "recommended"
  record_path "phase_22" "PCoA coordinates" "${STAT_DIR}/beta_diversity_pcoa_bray_taxa_coordinates.tsv" "recommended"
  record_path "phase_23" "PERMANOVA report" "${STAT_DIR}/permanova_taxa_bray_treatment_timepoint_interaction.txt" "recommended"
  record_path "phase_25" "final Markdown report" "${PROJECT_DIR}/results/report/final_pipeline_report.md" "recommended"
  record_path "phase_25" "final report output manifest" "${PROJECT_DIR}/results/report/report_output_manifest.tsv" "recommended"
  record_path "phase_25" "taxonomic stacked barplot PNG" "${FIG_DIR}/taxonomic_stacked_barplot_top20.png" "recommended"
  echo

  echo "Upstream run metadata snippets"
  write_metadata_snippet "Sylph run metadata" "${TAX_DIR}/sylph_run_metadata.txt"
  write_metadata_snippet "Taxonomy table preparation summary" "${TAX_DIR}/taxonomy_table_preparation_summary.txt"
  write_metadata_snippet "Assembly/MAG run metadata" "${ASM_DIR}/assembly_mag_pipeline_run_metadata.txt"
  write_metadata_snippet "Functional annotation run metadata" "${FUNC_DIR}/functional_annotation_run_metadata.txt"
  write_metadata_snippet "Statistics run summary" "${STAT_DIR}/statistics_run_summary.txt"
  write_metadata_snippet "Final visualization/report manifest" "${PROJECT_DIR}/results/report/report_output_manifest.tsv"

  echo "Machine-readable manifest: $MANIFEST"
} > "$REPORT"

echo "Reproducibility report written to: $REPORT"
echo "Reproducibility manifest written to: $MANIFEST"
