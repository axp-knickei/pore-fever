#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
})

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (is.na(script_path)) {
  project_dir <- getwd()
} else {
  project_dir <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
}

metadata_file <- file.path(project_dir, "config", "samples.tsv")
taxa_file <- file.path(project_dir, "results", "taxonomy", "taxa_abundance_matrix.tsv")
function_file <- file.path(project_dir, "results", "functional_annotation", "function_abundance_matrix.tsv")

outdir <- file.path(project_dir, "results", "statistics")
figdir <- file.path(project_dir, "results", "figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

summary_lines <- character()
add_summary <- function(...) {
  text <- paste(..., collapse = "")
  summary_lines <<- c(summary_lines, text)
  message(text)
}

write_note <- function(path, lines) {
  writeLines(lines, con = path)
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

prepare_abundance_matrix <- function(table_file, metadata, label) {
  table <- read_tsv(table_file, show_col_types = FALSE)
  if (ncol(table) < 2) {
    stop(paste(label, "table must contain one feature column plus sample columns."), call. = FALSE)
  }

  feature_col <- colnames(table)[1]
  sample_cols <- setdiff(colnames(table), feature_col)
  unmatched_cols <- setdiff(sample_cols, metadata$sample_id)
  if (length(unmatched_cols) > 0) {
    stop(
      paste(label, "table contains sample columns not found in metadata:", paste(unmatched_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  numeric_table <- table %>% mutate(across(all_of(sample_cols), ~ suppressWarnings(as.numeric(.x))))
  mat <- numeric_table %>%
    column_to_rownames(feature_col) %>%
    as.matrix()

  if (anyNA(mat)) {
    stop(paste(label, "matrix contains NA/non-numeric values after parsing."), call. = FALSE)
  }
  if (any(mat < 0)) {
    stop(paste(label, "matrix contains negative abundance values."), call. = FALSE)
  }
  if (any(duplicated(rownames(mat)))) {
    stop(paste(label, "matrix contains duplicated feature names."), call. = FALSE)
  }

  mat <- t(mat)
  common_samples <- intersect(rownames(mat), metadata$sample_id)
  if (length(common_samples) < 5) {
    stop(paste("Too few matching samples between", label, "table and metadata."), call. = FALSE)
  }

  mat <- mat[common_samples, , drop = FALSE]
  metadata2 <- metadata %>%
    filter(sample_id %in% common_samples) %>%
    arrange(match(sample_id, common_samples))

  if (!identical(rownames(mat), metadata2$sample_id)) {
    stop(paste(label, "sample ordering mismatch after metadata alignment."), call. = FALSE)
  }

  feature_sums <- colSums(mat)
  keep_features <- feature_sums > 0
  dropped_features <- sum(!keep_features)
  mat <- mat[, keep_features, drop = FALSE]

  sample_sums <- rowSums(mat)
  diagnostics <- tibble(
    sample_id = rownames(mat),
    total_abundance = sample_sums,
    nonzero_features = rowSums(mat > 0)
  ) %>%
    left_join(metadata2, by = "sample_id")

  list(
    matrix = mat,
    metadata = metadata2,
    feature_col = feature_col,
    dropped_zero_features = dropped_features,
    sample_diagnostics = diagnostics
  )
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

run_diversity <- function(mat, metadata2) {
  shannon <- diversity(mat, index = "shannon")
  richness <- specnumber(mat)

  alpha_df <- metadata2 %>%
    mutate(
      shannon = shannon,
      richness = richness
    )
  write_tsv(alpha_df, file.path(outdir, "alpha_diversity_taxa.tsv"))

  p_alpha <- ggplot(alpha_df, aes(x = treatment, y = shannon, color = treatment)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.7) +
    facet_wrap(~ timepoint) +
    theme_bw() +
    theme(legend.position = "none") +
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

  bray <- vegdist(mat, method = "bray")

  permanova_formula <- bray ~ treatment * timepoint
  repeated_measure_note <- "No mouse_id column found; PERMANOVA was run without strata/blocking."
  strata_arg <- NULL
  if ("mouse_id" %in% colnames(metadata2)) {
    if (any(is.na(metadata2$mouse_id) | trimws(as.character(metadata2$mouse_id)) == "")) {
      stop("mouse_id column exists but contains empty values; cannot use it for repeated-measure blocking.", call. = FALSE)
    }
    repeated_measure_note <- "mouse_id column found; PERMANOVA used strata = mouse_id for repeated-measure blocking."
    strata_arg <- metadata2$mouse_id
  }

  if (is.null(strata_arg)) {
    permanova <- adonis2(permanova_formula, data = metadata2, permutations = 999)
  } else {
    permanova <- adonis2(permanova_formula, data = metadata2, permutations = 999, strata = strata_arg)
  }
  capture.output(
    list(
      model = "bray ~ treatment * timepoint",
      repeated_measure_note = repeated_measure_note,
      permanova = permanova
    ),
    file = file.path(outdir, "permanova_taxa_bray_treatment_timepoint_interaction.txt")
  )

  dispersion_group <- interaction(metadata2$treatment, metadata2$timepoint, drop = TRUE)
  dispersion_counts <- table(dispersion_group)
  if (length(dispersion_counts) >= 2 && all(dispersion_counts >= 2)) {
    beta_disp <- betadisper(bray, dispersion_group)
    beta_disp_perm <- permutest(beta_disp, permutations = 999)
    capture.output(
      list(
        note = "PERMANOVA can be confounded by unequal dispersion. Inspect this test before interpreting PERMANOVA.",
        group = "interaction(treatment, timepoint)",
        group_sizes = dispersion_counts,
        betadisper = beta_disp,
        permutest = beta_disp_perm
      ),
      file = file.path(outdir, "betadisper_taxa_bray_treatment_timepoint.txt")
    )
  } else {
    write_note(
      file.path(outdir, "betadisper_taxa_bray_treatment_timepoint_SKIPPED.txt"),
      c(
        "betadisper skipped because at least two treatment:timepoint groups with at least two samples each are required.",
        "Observed group sizes:",
        capture.output(print(dispersion_counts))
      )
    )
  }

  pcoa <- cmdscale(bray, k = 2, eig = TRUE)
  pcoa_df <- metadata2 %>%
    mutate(
      PCoA1 = pcoa$points[, 1],
      PCoA2 = pcoa$points[, 2]
    )
  write_tsv(pcoa_df, file.path(outdir, "beta_diversity_pcoa_bray_taxa_coordinates.tsv"))

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

  list(
    alpha_df = alpha_df,
    pcoa_df = pcoa_df,
    repeated_measure_note = repeated_measure_note
  )
}

write_differential_notes <- function(mat, metadata2) {
  note_file <- file.path(outdir, "differential_taxa_abundance_NOTE.txt")
  lines <- c(
    "Phase 23 compositional differential abundance note",
    "This scaffold does not force a single differential abundance method because method choice depends on table scale, sparsity, repeated-measure design, and installed packages.",
    "Recommended primary methods: ANCOM-BC2/ANCOMBC, ALDEx2, or MaAsLin2/MaAsLin3.",
    "Recommended formula: abundance ~ treatment * timepoint, with mouse_id as a random/blocking term if longitudinal repeated sampling exists.",
    "Always interpret differential abundance together with compositional assumptions, prevalence filtering, multiple testing correction, and dispersion diagnostics."
  )
  write_note(note_file, lines)

  if (requireNamespace("ALDEx2", quietly = TRUE)) {
    aldex_note <- file.path(outdir, "differential_taxa_aldex2_NOTE.txt")
    write_note(
      aldex_note,
      c(
        "ALDEx2 is installed, but automatic multi-factor treatment*timepoint execution is not implemented in this scaffold.",
        "Use the taxa abundance matrix with an explicit contrast design after selecting the biological comparison."
      )
    )
  }

  if (requireNamespace("Maaslin2", quietly = TRUE)) {
    maaslin_note <- file.path(outdir, "differential_taxa_maaslin2_NOTE.txt")
    write_note(
      maaslin_note,
      c(
        "MaAsLin2 is installed. Suggested fixed effects: treatment,timepoint,treatment:timepoint.",
        "If mouse_id exists, include it as a random effect after confirming the repeated-measure design."
      )
    )
  }

  if (requireNamespace("ANCOMBC", quietly = TRUE)) {
    ancom_note <- file.path(outdir, "differential_taxa_ancombc_NOTE.txt")
    write_note(
      ancom_note,
      c(
        "ANCOMBC is installed. Suggested formula: treatment * timepoint.",
        "Confirm whether your installed ANCOMBC version supports the exact repeated-measure/random-effect model you need before running production analysis."
      )
    )
  }
}

write_function_notes <- function(function_file) {
  if (file.exists(function_file)) {
    note <- file.path(outdir, "differential_functional_abundance_NOTE.txt")
    write_note(
      note,
      c(
        "Phase 24 functional abundance input exists, but functional differential abundance is not implemented in this scaffold yet.",
        paste("Input detected:", function_file),
        "Recommended next step: validate table format and apply a compositional method with treatment, timepoint, and repeated-measure terms as appropriate."
      )
    )
    add_summary("Function abundance file detected; wrote Phase 24 implementation note.")
  } else {
    note <- file.path(outdir, "functional_abundance_skipped_NOTE.txt")
    write_note(
      note,
      c(
        "Phase 24 skipped: function abundance table not found.",
        paste("Expected file:", function_file),
        "Current functional annotation scaffold writes function_abundance_matrix_NOTE.txt instead of a real matrix until quantification is implemented."
      )
    )
    add_summary("Function abundance file not found; skipped Phase 24 and wrote note.")
  }
}

stop_if_missing(metadata_file)
metadata <- read_tsv(metadata_file, show_col_types = FALSE, col_types = cols(.default = col_character()))
metadata <- validate_metadata(metadata)

add_summary("Statistics script project directory: ", project_dir)
add_summary("Metadata samples: ", nrow(metadata))

if (file.exists(taxa_file) && file.info(taxa_file)$size > 0) {
  taxa_prepared <- prepare_abundance_matrix(taxa_file, metadata, "Taxa abundance")
  taxa_mat <- taxa_prepared$matrix
  metadata2 <- taxa_prepared$metadata
  sample_sums <- rowSums(taxa_mat)
  abundance_scale <- infer_abundance_scale(sample_sums)

  write_tsv(taxa_prepared$sample_diagnostics, file.path(outdir, "taxa_abundance_matrix_diagnostics.tsv"))
  add_summary("Taxa matrix matched samples: ", nrow(taxa_mat))
  add_summary("Taxa matrix retained features: ", ncol(taxa_mat))
  add_summary("Dropped zero-sum taxa: ", taxa_prepared$dropped_zero_features)
  add_summary("Inferred taxa abundance scale: ", abundance_scale)

  run_result <- run_diversity(taxa_mat, metadata2)
  add_summary(run_result$repeated_measure_note)
  write_differential_notes(taxa_mat, metadata2)
  add_summary("Phase 22 diversity outputs and Phase 23 differential-abundance notes written.")
} else {
  note <- file.path(outdir, "taxa_abundance_skipped_NOTE.txt")
  write_note(
    note,
    c(
      "Phase 22 skipped: taxa abundance matrix not found or empty.",
      paste("Expected file:", taxa_file),
      "Generate this file with: python scripts/python/03_prepare_taxonomy_table.py"
    )
  )
  add_summary("Taxa abundance file not found; skipped Phase 22 and wrote note.")
}

write_function_notes(function_file)

summary_file <- file.path(outdir, "statistics_run_summary.txt")
write_note(summary_file, summary_lines)
message("Statistics script completed. Summary written to: ", summary_file)
