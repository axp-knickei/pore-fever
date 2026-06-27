#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

project_dir <- getwd()
metadata_file <- file.path(project_dir, "config", "samples.tsv")
taxa_file <- file.path(project_dir, "results", "taxonomy", "taxa_abundance_matrix.tsv")
figdir <- file.path(project_dir, "results", "figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

stop_if_missing <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop(paste("Missing or empty file:", path))
  }
}

stop_if_missing(metadata_file)
metadata <- read_tsv(metadata_file, show_col_types = FALSE)

if (file.exists(taxa_file)) {
  taxa <- read_tsv(taxa_file, show_col_types = FALSE)
  taxon_col <- colnames(taxa)[1]

  taxa_long <- taxa %>%
    pivot_longer(
      cols = -all_of(taxon_col),
      names_to = "sample_id",
      values_to = "abundance"
    ) %>%
    left_join(metadata, by = "sample_id")

  top_taxa <- taxa_long %>%
    group_by(.data[[taxon_col]]) %>%
    summarise(total = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total)) %>%
    slice_head(n = 20) %>%
    pull(.data[[taxon_col]])

  taxa_plot <- taxa_long %>%
    mutate(
      taxon_plot = if_else(.data[[taxon_col]] %in% top_taxa, .data[[taxon_col]], "Other")
    ) %>%
    group_by(sample_id, treatment, timepoint, taxon_plot) %>%
    summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

  p <- ggplot(taxa_plot, aes(x = sample_id, y = abundance, fill = taxon_plot)) +
    geom_col() +
    facet_grid(timepoint ~ treatment, scales = "free_x", space = "free_x") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 6)
    ) +
    labs(
      title = "Top taxa across mouse gut ONT metagenomes",
      x = "Sample",
      y = "Relative abundance",
      fill = "Taxon"
    )

  ggsave(
    file.path(figdir, "taxonomic_stacked_barplot_top20.png"),
    p,
    width = 14,
    height = 8,
    dpi = 300
  )
} else {
  message("Taxa abundance file not found. Skipping taxonomic barplot.")
}

message("Visualization script completed.")
