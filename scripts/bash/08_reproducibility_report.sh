#!/usr/bin/env bash
set -euo pipefail

source config/config.sh

REPORT="${PROJECT_DIR}/results/reproducibility_report.txt"
mkdir -p "${PROJECT_DIR}/results"

{
  echo "ONT Mouse Gut Metagenomics Reproducibility Report"
  echo "Generated on: $(date)"
  echo
  echo "Project directory:"
  echo "$PROJECT_DIR"
  echo
  echo "Sample metadata:"
  echo "$SAMPLES_TSV"
  echo
  echo "Mouse reference:"
  echo "$MOUSE_REF"
  echo
  echo "Filtering parameters:"
  echo "MIN_LENGTH=$MIN_LENGTH"
  echo "MIN_Q=$MIN_Q"
  echo
  echo "Tool versions:"
  echo
  echo "minimap2:"
  minimap2 --version || true
  echo
  echo "samtools:"
  samtools --version | head -n 2 || true
  echo
  echo "NanoPlot:"
  NanoPlot --version || true
  echo
  echo "chopper:"
  chopper --version || true
  echo
  echo "flye:"
  flye --version || true
  echo
  echo "coverm:"
  coverm --version || true
  echo
  echo "checkm2:"
  checkm2 --version || true
  echo
  echo "gtdbtk:"
  gtdbtk --version || true
  echo
  echo "R:"
  R --version | head -n 1 || true
  echo
  echo "Python:"
  python --version || true
} > "$REPORT"

echo "Reproducibility report written to: $REPORT"
