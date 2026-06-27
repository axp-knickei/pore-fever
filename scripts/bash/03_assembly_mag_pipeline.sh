#!/usr/bin/env bash
set -euo pipefail

source config/config.sh

mkdir -p \
  "$ASM_DIR" "$MAP_DIR" "$BIN_DIR" "$MAG_QC_DIR" \
  "$MAG_TAX_DIR" "$DREP_DIR" "$MAG_ABUND_DIR" "$LOG_DIR"

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

check_tool flye
check_tool minimap2
check_tool samtools
check_tool seqkit
check_tool coverm
check_tool semibin2
check_tool checkm2
check_tool gtdbtk
check_tool dRep

COASSEMBLY_FASTQ="${ASM_DIR}/all_samples.host_removed.fastq.gz"
ASSEMBLY_DIR="${ASM_DIR}/metaflye_coassembly"
ASSEMBLY_FASTA="${ASSEMBLY_DIR}/assembly.fasta"

if [[ ! -s "$COASSEMBLY_FASTQ" ]]; then
  echo "Creating co-assembly FASTQ..."
  cat "${HOST_REMOVED_DIR}"/*.host_removed.fastq.gz > "$COASSEMBLY_FASTQ"
fi
check_file "$COASSEMBLY_FASTQ"

if [[ ! -s "$ASSEMBLY_FASTA" ]]; then
  echo "Phase 9: running metaFlye co-assembly..."
  flye \
    --nano-raw "$COASSEMBLY_FASTQ" \
    --meta \
    --threads "$THREADS" \
    --out-dir "$ASSEMBLY_DIR"
fi
check_file "$ASSEMBLY_FASTA"

echo "Phase 10: assembly quality summary"
seqkit stats "$ASSEMBLY_FASTA" > "${ASM_DIR}/assembly_seqkit_stats.tsv"

echo "Phase 11: mapping reads back to contigs and calculating coverage"
tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r sample_id treatment timepoint replicate raw_fastq
do
  fq="${HOST_REMOVED_DIR}/${sample_id}.host_removed.fastq.gz"
  bam="${MAP_DIR}/${sample_id}.to_contigs.bam"

  check_file "$fq"

  if [[ ! -s "$bam" ]]; then
    minimap2 \
      -t "$THREADS" \
      -ax map-ont \
      "$ASSEMBLY_FASTA" \
      "$fq" | \
      samtools sort -@ "$THREADS" -o "$bam"

    samtools index "$bam"
  fi
done

SEMIBIN_OUT="${BIN_DIR}/semibin2_output"

if [[ ! -d "$SEMIBIN_OUT" ]]; then
  echo "Phase 12: running SemiBin2 binning"
  semibin2 single_easy_bin \
    --input-fasta "$ASSEMBLY_FASTA" \
    --input-bam "${MAP_DIR}"/*.to_contigs.bam \
    --environment human_gut \
    --output "$SEMIBIN_OUT" \
    --threads "$THREADS"
fi

echo "Phase 13: add MAG refinement here if multiple binning outputs are available."

REFINED_BINS="${SEMIBIN_OUT}/output_bins"
if [[ ! -d "$REFINED_BINS" ]]; then
  echo "ERROR: Missing refined bin directory: $REFINED_BINS" >&2
  exit 1
fi

CHECKM2_OUT="${MAG_QC_DIR}/checkm2"
if [[ ! -d "$CHECKM2_OUT" ]]; then
  echo "Phase 14: running CheckM2"
  checkm2 predict \
    --threads "$THREADS" \
    --input "$REFINED_BINS" \
    --output-directory "$CHECKM2_OUT"
fi

GTDB_OUT="${MAG_TAX_DIR}/gtdbtk"
if [[ ! -d "$GTDB_OUT" ]]; then
  echo "Phase 15: running GTDB-Tk"
  gtdbtk classify_wf \
    --genome_dir "$REFINED_BINS" \
    --out_dir "$GTDB_OUT" \
    --extension fa \
    --cpus "$THREADS"
fi

if [[ ! -d "${DREP_DIR}/drep_output" ]]; then
  echo "Phase 16: running dRep dereplication"
  dRep dereplicate \
    "${DREP_DIR}/drep_output" \
    -g "${REFINED_BINS}"/*.fa \
    -p "$THREADS" \
    -comp 50 \
    -con 10
fi

DREP_GENOMES="${DREP_DIR}/drep_output/dereplicated_genomes"

if [[ ! -s "${MAG_ABUND_DIR}/mag_abundance.tsv" ]]; then
  echo "Phase 17: estimating MAG abundance with CoverM"
  coverm genome \
    --genome-fasta-directory "$DREP_GENOMES" \
    --reads "${HOST_REMOVED_DIR}"/*.host_removed.fastq.gz \
    --min-read-percent-identity 90 \
    --min-read-aligned-percent 50 \
    --methods mean covered_fraction relative_abundance \
    --threads "$THREADS" \
    > "${MAG_ABUND_DIR}/mag_abundance.tsv"
fi

echo "Phases 9-17 completed."
