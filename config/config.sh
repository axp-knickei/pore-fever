#!/usr/bin/env bash

# Project root. Run scripts from the repository root.
PROJECT_DIR="$(pwd)"

# Input metadata
SAMPLES_TSV="${PROJECT_DIR}/config/samples.tsv"

# Mouse reference genome
MOUSE_REF="${PROJECT_DIR}/config/mouse_reference/GRCm39.fa"

# Database paths. Replace these with local database locations.
SYLPH_DB="${PROJECT_DIR}/databases/sylph/gtdb-rs214-c200-dbv1.syldb"

# Threads
THREADS="${THREADS:-16}"

# Minimum filtering thresholds
MIN_LENGTH="${MIN_LENGTH:-1000}"
MIN_Q="${MIN_Q:-10}"

# Directories
RAW_DIR="${PROJECT_DIR}/data/raw_fastq"
TRIMMED_DIR="${PROJECT_DIR}/data/trimmed_fastq"
FILTERED_DIR="${PROJECT_DIR}/data/filtered_fastq"
HOST_REMOVED_DIR="${PROJECT_DIR}/data/host_removed_fastq"

QC_DIR="${PROJECT_DIR}/results/qc"
TAX_DIR="${PROJECT_DIR}/results/taxonomy"
ASM_DIR="${PROJECT_DIR}/results/assembly"
MAP_DIR="${PROJECT_DIR}/results/mapping"
BIN_DIR="${PROJECT_DIR}/results/bins"
MAG_QC_DIR="${PROJECT_DIR}/results/mag_qc"
MAG_TAX_DIR="${PROJECT_DIR}/results/mag_taxonomy"
DREP_DIR="${PROJECT_DIR}/results/mag_dereplication"
MAG_ABUND_DIR="${PROJECT_DIR}/results/mag_abundance"
FUNC_DIR="${PROJECT_DIR}/results/functional_annotation"
GENE_CAT_DIR="${PROJECT_DIR}/results/gene_catalog"
STAT_DIR="${PROJECT_DIR}/results/statistics"
FIG_DIR="${PROJECT_DIR}/results/figures"
LOG_DIR="${PROJECT_DIR}/logs"
