#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
})

project_dir <- getwd()

metadata_file <- file.path(project_dir, "config", "samples.tsv")
taxa_file <- file.path(project_dir, "results", "taxonomy", "taxa_abundance_matrix.tsv")
function_file <- file.path(project_dir, "results", "functional_annotation", "function_abundance_matrix.tsv")

outdir <- file.path(project_dir, "results", "statistics")
figdir <- file.path(project_dir, "results", "figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

stop_if_missing <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop(paste("Missing or empty file:", path))
  }
}

stop_if_missing(metadata_file)
metadata <- read_tsv(metadata_file, show_col_types = FALSE)

required_cols <- c("sample_id", "treatment", "timepoint", "replicate")
missing_cols <- setdiff(required_cols, colnames(metadata))
if (length(missing_cols) > 0) {
  stop(paste("Metadata missing columns:", paste(missing_cols, collapse = ", ")))
}

if (file.exists(taxa_file)) {
  taxa <- read_tsv(taxa_file, show_col_types = FALSE)

  taxon_col <- colnames(taxa)[1]
  taxa_mat <- taxa %>%
    column_to_rownames(taxon_col) %>%
    as.matrix()
  taxa_mat <- t(taxa_mat)

  common_samples <- intersect(rownames(taxa_mat), metadata$sample_id)
  if (length(common_samples) < 5) {
    stop("Too few matching samples between taxa table and metadata.")
  }

  taxa_mat <- taxa_mat[common_samples, , drop = FALSE]
  metadata2 <- metadata %>%
    filter(sample_id %in% common_samples) %>%
    arrange(match(sample_id, common_samples))

  shannon <- diversity(taxa_mat, index = "shannon")
  richness <- specnumber(taxa_mat)

  alpha_df <- metadata2 %>%
    mutate(
      shannon = shannon,
      richness = richness
    )

  write_tsv(alpha_df, file.path(outdir, "alpha_diversity_taxa.tsv"))

  p_alpha <- ggplot(alpha_df, aes(x = treatment, y = shannon)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.7) +
    facet_wrap(~ timepoint) +
    theme_bw() +
    labs(
      title = "Taxonomic alpha diversity",
      y = "Shannon diversity",
      x = "Treatment"
    )

  ggsave(
    file.path(figdir, "alpha_diversity_shannon_by_treatment_timepoint.png"),
    p_alpha,
    width = 10,
    height = 6,
    dpi = 300
  )

  bray <- vegdist(taxa_mat, method = "bray")
  permanova <- adonis2(bray ~ treatment + timepoint, data = metadata2)
  capture.output(
    permanova,
    file = file.path(outdir, "permanova_taxa_bray_treatment_timepoint.txt")
  )

  pcoa <- cmdscale(bray, k = 2, eig = TRUE)
  pcoa_df <- metadata2 %>%
    mutate(
      PCoA1 = pcoa$points[, 1],
      PCoA2 = pcoa$points[, 2]
    )

  p_pcoa <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = treatment, shape = timepoint)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_bw() +
    labs(
      title = "Beta diversity PCoA based on Bray-Curtis distance",
      x = "PCoA1",
      y = "PCoA2"
    )

  ggsave(
    file.path(figdir, "beta_diversity_pcoa_bray_taxa.png"),
    p_pcoa,
    width = 8,
    height = 6,
    dpi = 300
  )

  message("Phase 23: taxa table exists. Ready for ANCOM-BC2, ALDEx2, or MaAsLin covariate modeling.")
} else {
  message("Taxa abundance file not found. Skipping Phase 22 taxonomic diversity.")
}

if (file.exists(function_file)) {
  message("Phase 24: function abundance table exists. Ready for differential functional abundance.")
} else {
  message("Function abundance table not found. Skipping Phase 24.")
}

message("Statistics script completed.")
