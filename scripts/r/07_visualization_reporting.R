#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

args_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(args_file) > 0) sub("^--file=", "", args_file[1]) else NA_character_
if (is.na(script_path)) {
  project_dir <- getwd()
} else {
  project_dir <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
}

metadata_file <- file.path(project_dir, "config", "samples.tsv")
taxa_file <- file.path(project_dir, "results", "taxonomy", "taxa_abundance_matrix.tsv")
statistics_dir <- file.path(project_dir, "results", "statistics")
alpha_file <- file.path(statistics_dir, "alpha_diversity_taxa.tsv")
pcoa_file <- file.path(statistics_dir, "beta_diversity_pcoa_bray_taxa_coordinates.tsv")
permanova_file <- file.path(statistics_dir, "permanova_taxa_bray_treatment_timepoint_interaction.txt")
betadisper_file <- file.path(statistics_dir, "betadisper_taxa_bray_treatment_timepoint.txt")
betadisper_skipped_file <- file.path(statistics_dir, "betadisper_taxa_bray_treatment_timepoint_SKIPPED.txt")
statistics_summary_file <- file.path(statistics_dir, "statistics_run_summary.txt")

figdir <- file.path(project_dir, "results", "figures")
report_dir <- file.path(project_dir, "results", "report")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

report_lines <- character()
manifest <- tibble(output_file = character(), output_type = character(), description = character())

add_report <- function(...) {
  text <- paste(..., collapse = "")
  report_lines <<- c(report_lines, text)
  message(text)
}

record_output <- function(path, output_type, description) {
  manifest <<- bind_rows(
    manifest,
    tibble(
      output_file = normalizePath(path, mustWork = FALSE),
      output_type = output_type,
      description = description
    )
  )
}

write_note <- function(path, lines, description) {
  writeLines(lines, con = path)
  record_output(path, "note", description)
}

stop_if_missing <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop(paste("Missing or empty file:", path), call. = FALSE)
  }
}

validate_metadata <- function(metadata) {
  required_cols <- c("sample_id", "treatment", "timepoint", "replicate")
  missing_cols <- setdiff(required_cols, colnames(metadata))
  if (length(missing_cols) > 0) {
    stop(paste("Metadata missing columns:", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  for (col in required_cols) {
    if (any(is.na(metadata[[col]]) | trimws(as.character(metadata[[col]])) == "")) {
      stop(paste("Metadata contains empty values in column:", col), call. = FALSE)
    }
  }

  duplicated_samples <- metadata$sample_id[duplicated(metadata$sample_id)]
  if (length(duplicated_samples) > 0) {
    stop(paste("Duplicate sample_id values:", paste(unique(duplicated_samples), collapse = ", ")), call. = FALSE)
  }

  metadata2 <- metadata %>%
    mutate(
      sample_id = as.character(sample_id),
      treatment = factor(treatment, levels = c("control", "placebo", "treatment1", "treatment2", "treatment3")),
      timepoint = factor(timepoint, levels = c("T1", "T2", "T3", "T4", "T5")),
      replicate = factor(replicate)
    )

  if (any(is.na(metadata2$treatment))) {
    stop("Metadata contains treatment values outside expected levels: control, placebo, treatment1, treatment2, treatment3", call. = FALSE)
  }
  if (any(is.na(metadata2$timepoint))) {
    stop("Metadata contains timepoint values outside expected levels: T1, T2, T3, T4, T5", call. = FALSE)
  }

  metadata2
}

infer_abundance_scale <- function(sample_sums) {
  median_sum <- median(sample_sums)
  if (median_sum <= 1.5) {
    "proportion-like"
  } else if (median_sum <= 150) {
    "percent-like"
  } else {
    "count-or-depth-like"
  }
}

save_plot <- function(plot, filename_stem, width, height, description) {
  png_path <- file.path(figdir, paste0(filename_stem, ".png"))
  pdf_path <- file.path(figdir, paste0(filename_stem, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = 300)
  ggsave(pdf_path, plot, width = width, height = height)
  record_output(png_path, "figure_png", description)
  record_output(pdf_path, "figure_pdf", description)
}

prepare_taxa_long <- function(taxa_file, metadata) {
  stop_if_missing(taxa_file)
  taxa <- read_tsv(taxa_file, show_col_types = FALSE)
  if (ncol(taxa) < 2) {
    stop("Taxa abundance matrix must contain one taxon column plus sample columns.", call. = FALSE)
  }

  taxon_col <- colnames(taxa)[1]
  if (any(is.na(taxa[[taxon_col]]) | trimws(as.character(taxa[[taxon_col]])) == "")) {
    stop("Taxa abundance matrix contains empty taxon names.", call. = FALSE)
  }
  if (any(duplicated(taxa[[taxon_col]]))) {
    stop("Taxa abundance matrix contains duplicated taxon names.", call. = FALSE)
  }

  sample_cols <- setdiff(colnames(taxa), taxon_col)
  unmatched_cols <- setdiff(sample_cols, metadata$sample_id)
  if (length(unmatched_cols) > 0) {
    stop(
      paste("Taxa abundance matrix contains sample columns not found in metadata:", paste(unmatched_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  missing_samples <- setdiff(metadata$sample_id, sample_cols)
  if (length(missing_samples) > 0) {
    add_report("Warning: metadata samples without taxonomy columns: ", paste(missing_samples, collapse = ", "))
  }

  numeric_taxa <- taxa %>%
    mutate(across(all_of(sample_cols), ~ suppressWarnings(as.numeric(.x))))
  abundance_mat <- numeric_taxa %>%
    select(all_of(sample_cols)) %>%
    as.matrix()

  if (anyNA(abundance_mat)) {
    stop("Taxa abundance matrix contains NA/non-numeric abundance values.", call. = FALSE)
  }
  if (any(abundance_mat < 0)) {
    stop("Taxa abundance matrix contains negative abundance values.", call. = FALSE)
  }

  long <- numeric_taxa %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "sample_id",
      values_to = "abundance"
    ) %>%
    rename(taxon = all_of(taxon_col)) %>%
    mutate(taxon = as.character(taxon)) %>%
    left_join(metadata, by = "sample_id")

  if (any(is.na(long$treatment)) || any(is.na(long$timepoint))) {
    stop("Taxa abundance rows failed to join cleanly with metadata.", call. = FALSE)
  }

  sample_sums <- long %>%
    group_by(sample_id) %>%
    summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop")
  if (any(sample_sums$total_abundance <= 0)) {
    bad_samples <- sample_sums$sample_id[sample_sums$total_abundance <= 0]
    stop(paste("Taxa abundance matrix contains samples with zero total abundance:", paste(bad_samples, collapse = ", ")), call. = FALSE)
  }

  abundance_scale <- infer_abundance_scale(sample_sums$total_abundance)
  long <- long %>%
    left_join(sample_sums, by = "sample_id") %>%
    mutate(
      relative_abundance = abundance / total_abundance,
      sample_id = factor(sample_id, levels = metadata$sample_id),
      abundance_scale = abundance_scale
    )

  diagnostics <- sample_sums %>%
    left_join(metadata, by = "sample_id") %>%
    mutate(
      nonzero_taxa = map_int(as.character(sample_id), ~ sum(long$sample_id == .x & long$abundance > 0)),
      abundance_scale = abundance_scale
    )

  list(
    long = long,
    diagnostics = diagnostics,
    abundance_scale = abundance_scale,
    taxon_count = n_distinct(long$taxon),
    sample_count = n_distinct(long$sample_id)
  )
}

plot_taxa_barplot <- function(taxa_long, top_n = 20) {
  top_taxa <- taxa_long %>%
    group_by(taxon) %>%
    summarise(mean_relative_abundance = mean(relative_abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_relative_abundance), taxon) %>%
    slice_head(n = top_n) %>%
    pull(taxon)

  taxa_plot <- taxa_long %>%
    mutate(taxon_plot = if_else(taxon %in% top_taxa, taxon, "Other")) %>%
    group_by(sample_id, treatment, timepoint, taxon_plot) %>%
    summarise(relative_abundance = sum(relative_abundance, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      taxon_plot = fct_reorder(taxon_plot, relative_abundance, .fun = sum, .desc = TRUE),
      sample_id = fct_drop(sample_id)
    )

  ggplot(taxa_plot, aes(x = sample_id, y = relative_abundance, fill = taxon_plot)) +
    geom_col(width = 0.9) +
    facet_grid(timepoint ~ treatment, scales = "free_x", space = "free_x") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom",
      legend.key.size = grid::unit(0.35, "cm")
    ) +
    guides(fill = guide_legend(ncol = 2)) +
    labs(
      title = "Top taxa across mouse gut ONT metagenomes",
      subtitle = paste("Top", top_n, "taxa by mean relative abundance; remaining taxa grouped as Other"),
      x = "Sample",
      y = "Relative abundance within sample",
      fill = "Taxon"
    )
}

plot_alpha_if_available <- function(alpha_file) {
  if (!file.exists(alpha_file) || file.info(alpha_file)$size == 0) {
    write_note(
      file.path(report_dir, "alpha_diversity_reporting_SKIPPED.txt"),
      c(
        "Alpha diversity reporting plot skipped because the statistics output was not found.",
        paste("Expected file:", alpha_file),
        "Generate it with: Rscript scripts/r/06_statistics_diversity_differential.R"
      ),
      "Skipped alpha-diversity reporting note"
    )
    return(FALSE)
  }

  alpha <- read_tsv(alpha_file, show_col_types = FALSE)
  required <- c("sample_id", "treatment", "timepoint", "shannon", "richness")
  missing <- setdiff(required, colnames(alpha))
  if (length(missing) > 0) {
    stop(paste("Alpha diversity table missing columns:", paste(missing, collapse = ", ")), call. = FALSE)
  }

  alpha <- alpha %>%
    mutate(
      treatment = factor(treatment, levels = c("control", "placebo", "treatment1", "treatment2", "treatment3")),
      timepoint = factor(timepoint, levels = c("T1", "T2", "T3", "T4", "T5")),
      shannon = as.numeric(shannon),
      richness = as.numeric(richness)
    )
  if (anyNA(alpha$shannon) || anyNA(alpha$richness)) {
    stop("Alpha diversity table contains non-numeric shannon or richness values.", call. = FALSE)
  }

  alpha_long <- alpha %>%
    select(sample_id, treatment, timepoint, shannon, richness) %>%
    pivot_longer(cols = c(shannon, richness), names_to = "metric", values_to = "value") %>%
    mutate(metric = recode(metric, shannon = "Shannon diversity", richness = "Observed taxa"))

  p <- ggplot(alpha_long, aes(x = treatment, y = value, color = treatment)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.7, size = 1.8) +
    facet_grid(metric ~ timepoint, scales = "free_y") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = "Alpha diversity summary", x = "Treatment", y = "Diversity value")

  save_plot(p, "report_alpha_diversity_taxa", 12, 7, "Final report alpha-diversity figure")
  TRUE
}

plot_pcoa_if_available <- function(pcoa_file) {
  if (!file.exists(pcoa_file) || file.info(pcoa_file)$size == 0) {
    write_note(
      file.path(report_dir, "pcoa_reporting_SKIPPED.txt"),
      c(
        "PCoA reporting plot skipped because the statistics output was not found.",
        paste("Expected file:", pcoa_file),
        "Generate it with: Rscript scripts/r/06_statistics_diversity_differential.R"
      ),
      "Skipped PCoA reporting note"
    )
    return(FALSE)
  }

  pcoa <- read_tsv(pcoa_file, show_col_types = FALSE)
  required <- c("sample_id", "treatment", "timepoint", "PCoA1", "PCoA2")
  missing <- setdiff(required, colnames(pcoa))
  if (length(missing) > 0) {
    stop(paste("PCoA table missing columns:", paste(missing, collapse = ", ")), call. = FALSE)
  }

  pcoa <- pcoa %>%
    mutate(
      treatment = factor(treatment, levels = c("control", "placebo", "treatment1", "treatment2", "treatment3")),
      timepoint = factor(timepoint, levels = c("T1", "T2", "T3", "T4", "T5")),
      PCoA1 = as.numeric(PCoA1),
      PCoA2 = as.numeric(PCoA2)
    )
  if (anyNA(pcoa$PCoA1) || anyNA(pcoa$PCoA2)) {
    stop("PCoA table contains non-numeric coordinates.", call. = FALSE)
  }

  p <- ggplot(pcoa, aes(x = PCoA1, y = PCoA2, color = treatment, shape = timepoint)) +
    geom_point(size = 3, alpha = 0.85) +
    theme_bw() +
    labs(
      title = "Beta diversity PCoA summary",
      subtitle = "Bray-Curtis distances from the statistics step",
      x = "PCoA1",
      y = "PCoA2"
    )

  save_plot(p, "report_beta_diversity_pcoa_bray_taxa", 8, 6, "Final report beta-diversity PCoA figure")
  TRUE
}

copy_statistics_notes <- function() {
  stats_inputs <- c(
    permanova_file,
    betadisper_file,
    betadisper_skipped_file,
    statistics_summary_file
  )
  detected <- stats_inputs[file.exists(stats_inputs) & file.info(stats_inputs)$size > 0]
  if (length(detected) == 0) {
    write_note(
      file.path(report_dir, "statistics_outputs_SKIPPED.txt"),
      c(
        "No statistics text outputs were found for final reporting.",
        "Generate them with: Rscript scripts/r/06_statistics_diversity_differential.R"
      ),
      "Skipped statistics reporting note"
    )
    return(character())
  }
  detected
}

stop_if_missing(metadata_file)
metadata <- read_tsv(metadata_file, show_col_types = FALSE, col_types = cols(.default = col_character()))
metadata <- validate_metadata(metadata)

add_report("Final reporting project directory: ", project_dir)
add_report("Metadata samples: ", nrow(metadata))

if (file.exists(taxa_file) && file.info(taxa_file)$size > 0) {
  taxa_prepared <- prepare_taxa_long(taxa_file, metadata)
  taxa_long <- taxa_prepared$long
  diagnostics_file <- file.path(report_dir, "report_taxa_abundance_diagnostics.tsv")
  write_tsv(taxa_prepared$diagnostics, diagnostics_file)
  record_output(diagnostics_file, "table", "Taxa abundance diagnostics used by final reporting")

  plot_data_file <- file.path(report_dir, "report_taxa_relative_abundance_long.tsv")
  write_tsv(taxa_long, plot_data_file)
  record_output(plot_data_file, "table", "Long relative-abundance table used for final taxonomic barplot")

  p_taxa <- plot_taxa_barplot(taxa_long, top_n = 20)
  save_plot(p_taxa, "taxonomic_stacked_barplot_top20", 14, 8, "Top-20 taxonomic stacked barplot")

  add_report("Taxa reporting samples: ", taxa_prepared$sample_count)
  add_report("Taxa reporting features: ", taxa_prepared$taxon_count)
  add_report("Detected taxa abundance input scale: ", taxa_prepared$abundance_scale)
  add_report("Taxonomic plot normalized abundances within each sample before plotting.")
} else {
  write_note(
    file.path(report_dir, "taxa_visualization_SKIPPED.txt"),
    c(
      "Taxonomic visualization skipped because the taxa abundance matrix was not found or was empty.",
      paste("Expected file:", taxa_file),
      "Generate it with: python scripts/python/03_prepare_taxonomy_table.py"
    ),
    "Skipped taxonomic visualization note"
  )
  add_report("Taxa abundance matrix not found; taxonomic final figure skipped.")
}

alpha_created <- plot_alpha_if_available(alpha_file)
pcoa_created <- plot_pcoa_if_available(pcoa_file)
statistics_text_files <- copy_statistics_notes()

if (length(statistics_text_files) > 0) {
  add_report("Statistics text outputs detected for report: ", length(statistics_text_files))
  for (path in statistics_text_files) {
    record_output(path, "statistics_text", paste("Statistics output detected:", basename(path)))
  }
}
add_report("Alpha diversity report figure created: ", alpha_created)
add_report("PCoA report figure created: ", pcoa_created)

report_file <- file.path(report_dir, "final_pipeline_report.md")
report_body <- c(
  "# ONT Mouse Gut Metagenomics Final Report",
  "",
  "## Inputs",
  paste("- Metadata:", metadata_file),
  paste("- Taxa abundance matrix:", taxa_file),
  paste("- Statistics directory:", statistics_dir),
  "",
  "## Run Summary",
  paste("- Project directory:", project_dir),
  paste("- Metadata samples:", nrow(metadata)),
  report_lines,
  "",
  "## Main Figures",
  "- `results/figures/taxonomic_stacked_barplot_top20.png` and `.pdf` when taxonomy input is available.",
  "- `results/figures/report_alpha_diversity_taxa.png` and `.pdf` when statistics output is available.",
  "- `results/figures/report_beta_diversity_pcoa_bray_taxa.png` and `.pdf` when statistics output is available.",
  "",
  "## Interpretation Notes",
  "- Taxonomic barplots are normalized within sample for visualization, even when upstream abundance values are count/depth-like.",
  "- PERMANOVA, dispersion, and differential-abundance interpretation should be taken from the statistics outputs generated by `06_statistics_diversity_differential.R`.",
  "- Missing optional statistics inputs are recorded as skipped notes in this report directory."
)
writeLines(report_body, con = report_file)
record_output(report_file, "report", "Final Markdown report")

manifest_file <- file.path(report_dir, "report_output_manifest.tsv")
record_output(manifest_file, "manifest", "Final reporting output manifest")
write_tsv(manifest, manifest_file)
message("Visualization/reporting completed. Report written to: ", report_file)
message("Output manifest written to: ", manifest_file)
