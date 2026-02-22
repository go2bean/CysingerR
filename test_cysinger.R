# test_cysinger.R - Verification script for Cysinger package
# Run with: Rscript test_cysinger.R

pkg_dir <- "/Users/tllee/Antigravity/MSD/Cysinger"
r_files <- list.files(file.path(pkg_dir, "R"), pattern = "\\.R$", full.names = TRUE)

# Load dependencies
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(viridis)
  library(scales)
})

# Source package files in dependency order
source_order <- c("codon_table.R", "cys_frequency.R", "fastq_filter.R",
                  "demultiplex.R", "quality_filter.R", "fasta_convert.R",
                  "alignment.R", "combine_reads.R", "substitution_table.R",
                  "classify_mutations.R", "data.R",
                  "plot_cys_heatmap.R", "plot_cys_barplot.R",
                  "plot_mutation_landscape.R", "plot_qc.R",
                  "plot_epitope_map.R")

for (f in source_order) {
  fp <- file.path(pkg_dir, "R", f)
  if (file.exists(fp)) {
    cat("Sourcing:", f, "\n")
    source(fp)
  } else {
    cat("WARNING: Missing file:", f, "\n")
  }
}
cat("\n=== All R files sourced successfully ===\n\n")

# Test 1: Load example data
cat("Test 1: Loading example data...\n")
cys_data <- read.delim(
  file.path(pkg_dir, "inst", "extdata", "Freq_CYS_mutants_res_level.tsv"),
  header = FALSE, sep = "\t", stringsAsFactors = FALSE, skip = 1
)
n_bins <- ncol(cys_data) - 3
names(cys_data) <- c("position", "wt_residue", "mut_residue",
                     paste0("bin_", 1:n_bins))
cat("  Loaded:", nrow(cys_data), "rows x", ncol(cys_data), "columns\n")
cat("  Positions:", min(cys_data$position), "to", max(cys_data$position), "\n\n")

# Test 2: Normalization
cat("Test 2: Normalizing data...\n")
cys_norm <- normalize_cys_frequency(cys_data, input_bin = 1, sorted_bins = 2:4)
cat("  Enrichment columns added:", sum(grepl("enrichment", names(cys_norm))), "\n")
cat("  Enrichment range bin_4:",
    sprintf("%.3f - %.3f", min(cys_norm$enrichment_bin_4),
            max(cys_norm$enrichment_bin_4)), "\n\n")

# Test 3: Translation
cat("Test 3: CcdB sequences...\n")
cat("  WT nucleotide length:", nchar(ccdB_wt_nucleotide()), "\n")
cat("  WT protein:", ccdB_wt_protein(), "\n")
cat("  Reverse complement test:", reverse_complement_seq("ATGC"), "\n\n")

# Test 4: Codon table
cat("Test 4: Codon table...\n")
ct <- get_codon_table()
cat("  ATG ->", ct["ATG"], "\n")
cat("  TGT ->", ct["TGT"], "\n")
cat("  TGC ->", ct["TGC"], "\n\n")

# Test 5: Generate figures
cat("Test 5: Generating figures...\n")
fig_dir <- file.path(pkg_dir, "inst", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

tryCatch({
  # Heatmap
  p1 <- plot_cys_heatmap(cys_data,
    bin_labels = c("Input", "Unsorted", "Low Sort", "High Sort"))
  ggsave(file.path(fig_dir, "heatmap.png"), p1, width = 14, height = 4, dpi = 200)
  cat("  [OK] heatmap.png\n")
}, error = function(e) cat("  [FAIL] heatmap:", conditionMessage(e), "\n"))

tryCatch({
  # Bar plot
  p2 <- plot_cys_enrichment(cys_data,
    bin_labels = c("Input", "Unsorted", "Low Sort", "High Sort"))
  ggsave(file.path(fig_dir, "enrichment_barplot.png"), p2, width = 14, height = 6, dpi = 200)
  cat("  [OK] enrichment_barplot.png\n")
}, error = function(e) cat("  [FAIL] barplot:", conditionMessage(e), "\n"))

tryCatch({
  # Mutation landscape
  p3 <- plot_mutation_landscape(cys_data, bin_to_plot = 4, bin_label = "High Sort")
  ggsave(file.path(fig_dir, "mutation_landscape.png"), p3, width = 14, height = 5, dpi = 200)
  cat("  [OK] mutation_landscape.png\n")
}, error = function(e) cat("  [FAIL] landscape:", conditionMessage(e), "\n"))

tryCatch({
  # Epitope profile
  p4 <- plot_epitope_profile(cys_norm,
    bin_labels = c("Unsorted", "Low Sort", "High Sort"),
    smooth = TRUE, span = 0.2,
    highlight_regions = list(
      list(start = 20, end = 30, label = "Loop 1"),
      list(start = 55, end = 65, label = "Active Site"),
      list(start = 90, end = 101, label = "C-terminal")
    ))
  ggsave(file.path(fig_dir, "epitope_profile.png"), p4, width = 14, height = 6, dpi = 200)
  cat("  [OK] epitope_profile.png\n")
}, error = function(e) cat("  [FAIL] epitope:", conditionMessage(e), "\n"))

tryCatch({
  # Surface accessibility
  p5 <- plot_surface_accessibility(cys_norm, enrichment_bin = 4)
  ggsave(file.path(fig_dir, "surface_accessibility.png"), p5, width = 8, height = 8, dpi = 200)
  cat("  [OK] surface_accessibility.png\n")
}, error = function(e) cat("  [FAIL] surface:", conditionMessage(e), "\n"))

tryCatch({
  # QC dashboard
  p6 <- plot_qc_summary(cys_data,
    bin_labels = c("Input", "Unsorted", "Low Sort", "High Sort"))
  ggsave(file.path(fig_dir, "qc_dashboard.png"), p6, width = 14, height = 10, dpi = 200)
  cat("  [OK] qc_dashboard.png\n")
}, error = function(e) cat("  [FAIL] qc:", conditionMessage(e), "\n"))

tryCatch({
  # Donut chart
  mock_class <- data.frame(
    class = c("wt", "single", "double", "triple", "multi"),
    count = c(45000, 28000, 8500, 2100, 1400),
    fraction = c(0.529, 0.329, 0.100, 0.025, 0.016)
  )
  p7 <- plot_mutation_class_donut(mock_class)
  ggsave(file.path(fig_dir, "mutation_donut.png"), p7, width = 7, height = 7, dpi = 200)
  cat("  [OK] mutation_donut.png\n")
}, error = function(e) cat("  [FAIL] donut:", conditionMessage(e), "\n"))

cat("\n=== VERIFICATION COMPLETE ===\n")
cat("Figures saved to:", fig_dir, "\n")
cat("Files:", paste(list.files(fig_dir), collapse = ", "), "\n")
