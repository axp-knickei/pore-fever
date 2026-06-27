# ONT Mouse Gut Metagenomics

Script scaffold for a 75-sample Oxford Nanopore mouse gut metagenomics study:

- 5 treatment groups: `control`, `placebo`, `treatment1`, `treatment2`, `treatment3`
- 5 time points: `T1` to `T5`
- 3 biological replicates per treatment/time point
- Input: basecalled ONT FASTQ files
- Main outputs: cleaned reads, QC summaries, taxonomic tables, MAG catalog, functional tables, statistics, and figures

## Layout

```text
config/                         sample metadata and shared settings
data/raw_fastq/                  input FASTQ files
data/trimmed_fastq/              adapter-trimmed reads
data/filtered_fastq/             length/quality-filtered reads
data/host_removed_fastq/         mouse-depleted reads
results/qc/                     read QC and host-removal summaries
results/taxonomy/               taxonomic profiles and abundance matrices
results/assembly/               co-assembly outputs
results/mapping/                read-to-contig mappings
results/bins/                   MAG binning outputs
results/mag_qc/                 MAG quality reports
results/mag_taxonomy/           MAG taxonomic assignments
results/mag_dereplication/      dereplicated MAG catalog
results/mag_abundance/          MAG abundance tables
results/functional_annotation/  MAG/gene/function annotations
results/gene_catalog/           nonredundant gene catalog
results/statistics/             diversity and differential analysis outputs
results/figures/                publication-style figures
scripts/bash/                   heavy workflow steps
scripts/python/                 validation and summary helpers
scripts/r/                      statistical analysis and visualization
logs/                           run logs
envs/                           environment notes or future conda files
```

## Before Running

1. Put raw FASTQ files under `data/raw_fastq/`.
2. Edit `config/samples.tsv` so every `raw_fastq` path matches your actual file names.
3. Put the mouse reference genome at `config/mouse_reference/GRCm39.fa`, or update `MOUSE_REF` in `config/config.sh`.
4. Update database paths in `config/config.sh`, especially `SYLPH_DB`.
5. Install the required bioinformatics tools in your HPC/conda environment.

## Recommended Order

```bash
bash scripts/bash/01_preprocess_qc_host_removal.sh
python scripts/python/05_summarize_fastq_qc.py
bash scripts/bash/02_taxonomic_profiling.sh
python scripts/python/03_prepare_taxonomy_table.py
bash scripts/bash/03_assembly_mag_pipeline.sh
bash scripts/bash/04_functional_annotation.sh
Rscript scripts/r/06_statistics_diversity_differential.R
Rscript scripts/r/07_visualization_reporting.R
bash scripts/bash/08_reproducibility_report.sh
```

For 75 ONT fecal samples, start with the read-based route first:

```text
FASTQ -> QC -> host removal -> taxonomic profiling -> diversity -> differential taxa
```

Then proceed to functional profiling and MAG reconstruction once read quality, depth, and group effects look reasonable.
