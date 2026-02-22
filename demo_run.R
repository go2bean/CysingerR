# ============================================================
# Cysinger Package — Example Demo Run
# ============================================================

pkg_dir <- "/Users/tllee/Antigravity/MSD/Cysinger"

# Load dependencies
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

# Source all package files
r_files <- c("codon_table.R", "cys_frequency.R", "fastq_filter.R",
             "demultiplex.R", "quality_filter.R", "fasta_convert.R",
             "alignment.R", "combine_reads.R", "substitution_table.R",
             "classify_mutations.R", "data.R", "Cysinger-package.R",
             "plot_cys_heatmap.R", "plot_cys_barplot.R",
             "plot_mutation_landscape.R", "plot_qc.R",
             "plot_epitope_map.R")
for (f in r_files) source(file.path(pkg_dir, "R", f), local = FALSE)

cat("=" |> strrep(60), "\n")
cat("  CYSINGER PACKAGE DEMO\n")
cat("=" |> strrep(60), "\n\n")

# ---- 1. Codon Table & WT Sequences ----
cat(">>> 1. Codon Table & WT Sequences\n")
cat("-" |> strrep(40), "\n")
ct <- get_codon_table()
cat("Codon table has", length(ct), "entries\n")
cat("Sample codons:\n")
cat("  ATG ->", ct["ATG"], "(Met/Start)\n")
cat("  TGT ->", ct["TGT"], "(Cys)\n")
cat("  TGC ->", ct["TGC"], "(Cys)\n")
cat("  TAA ->", ct["TAA"], "(Stop)\n\n")

cat("CcdB WT nucleotide sequence (306 nt):\n")
wt_nt <- ccdB_wt_nucleotide()
cat("  Length:", nchar(wt_nt), "nt\n")
cat("  First 60:", substr(wt_nt, 1, 60), "...\n")
cat("  Last  30: ...", substr(wt_nt, nchar(wt_nt)-29, nchar(wt_nt)), "\n\n")

cat("CcdB WT protein:\n")
wt_prot <- ccdB_wt_protein()
cat(" ", wt_prot, "\n")
cat("  Length:", nchar(wt_prot), "aa (including stop)\n\n")

# ---- 2. Translation ----
cat(">>> 2. Sequence Translation\n")
cat("-" |> strrep(40), "\n")
test_seqs <- c("ATGTGTAAATGA", "ATGGCTTACTAA", "TTTGGG")
for (s in test_seqs) {
  cat(sprintf("  %s -> %s\n", s, translate_sequence(s)))
}
cat("\n")

# ---- 3. Reverse Complement ----
cat(">>> 3. Reverse Complement\n")
cat("-" |> strrep(40), "\n")
rc_tests <- c("ATGC", "AAACCCGGG", "TATATATA")
for (s in rc_tests) {
  cat(sprintf("  %s -> %s\n", s, reverse_complement_seq(s)))
}
cat("\n")

# ---- 4. Load Example Data ----
cat(">>> 4. Loading CcdB Example Data\n")
cat("-" |> strrep(40), "\n")
cys_data <- read.delim(
  file.path(pkg_dir, "inst", "extdata", "Freq_CYS_mutants_res_level.tsv"),
  header = FALSE, sep = "\t", stringsAsFactors = FALSE, skip = 1
)
n_bins <- ncol(cys_data) - 3
names(cys_data) <- c("position", "wt_residue", "mut_residue",
                     paste0("bin_", 1:n_bins))

cat("  Dimensions:", nrow(cys_data), "positions x", ncol(cys_data), "columns\n")
cat("  Columns:", paste(names(cys_data), collapse=", "), "\n")
cat("  Position range:", min(cys_data$position), "-", max(cys_data$position), "\n\n")

cat("  First 10 rows:\n")
print(head(cys_data, 10))
cat("\n")

cat("  Summary statistics per bin:\n")
for (b in 1:n_bins) {
  bc <- paste0("bin_", b)
  cat(sprintf("    Bin %d: total=%s, mean=%.1f, median=%.0f, max=%d, zeros=%d\n",
              b, format(sum(cys_data[[bc]]), big.mark=","),
              mean(cys_data[[bc]]), median(cys_data[[bc]]),
              max(cys_data[[bc]]),
              sum(cys_data[[bc]] == 0)))
}
cat("\n")

# ---- 5. Normalization ----
cat(">>> 5. Normalizing Cysteine Frequencies\n")
cat("-" |> strrep(40), "\n")
cys_norm <- normalize_cys_frequency(cys_data, input_bin = 1, sorted_bins = 2:4)
cat("  Output columns:", paste(names(cys_norm), collapse=", "), "\n\n")

cat("  Enrichment ratio summary (Bin 4 / Bin 1):\n")
cat(sprintf("    Min:    %.3f\n", min(cys_norm$enrichment_bin_4)))
cat(sprintf("    Median: %.3f\n", median(cys_norm$enrichment_bin_4)))
cat(sprintf("    Mean:   %.3f\n", mean(cys_norm$enrichment_bin_4)))
cat(sprintf("    Max:    %.3f\n", max(cys_norm$enrichment_bin_4)))
cat("\n")

cat("  Top 10 enriched positions (Bin 4):\n")
top10 <- cys_norm[order(-cys_norm$enrichment_bin_4), ][1:10, ]
print(top10[, c("position", "wt_residue", "bin_1", "bin_4",
                "enrichment_bin_4", "log2_enrichment_bin_4")])
cat("\n")

cat("  Bottom 5 depleted positions (Bin 4):\n")
bot5 <- cys_norm[order(cys_norm$enrichment_bin_4), ][1:5, ]
print(bot5[, c("position", "wt_residue", "bin_1", "bin_4",
               "enrichment_bin_4", "log2_enrichment_bin_4")])
cat("\n")

# ---- 6. Mutation Classification (mock demo) ----
cat(">>> 6. Mutation Classification (Mock Data)\n")
cat("-" |> strrep(40), "\n")

# Build mock data matching build_substitution_table() output
set.seed(42)
n_reads <- 200
mock_sub <- data.frame(
  read_id = paste0("read_", 1:n_reads),
  bin = sample(1:4, n_reads, replace = TRUE),
  is_wt = c(rep(TRUE, 80), rep(FALSE, 120)),
  residue_positions = c(
    rep("wt", 80),                               # 80 WT reads
    sample(as.character(1:101), 50, replace=TRUE), # 50 single mutants
    paste(sample(1:101, 30, replace=TRUE),         # 30 double mutants
          sample(1:101, 30, replace=TRUE)),
    paste(sample(1:101, 20, replace=TRUE),         # 20 triple mutants
          sample(1:101, 20, replace=TRUE),
          sample(1:101, 20, replace=TRUE)),
    paste(sample(1:101, 20, replace=TRUE),         # 20 multi mutants
          sample(1:101, 20, replace=TRUE),
          sample(1:101, 20, replace=TRUE),
          sample(1:101, 20, replace=TRUE))
  ),
  stringsAsFactors = FALSE
)
classified <- classify_mutations(mock_sub)
cat("  Input: ", n_reads, "reads\n")
cat("  Classification summary:\n")
print(classified$summary)
cat("\n")

# ---- 7. Generate All Figures ----
cat(">>> 7. Generating Publication Figures\n")
cat("-" |> strrep(40), "\n")
fig_dir <- file.path(pkg_dir, "inst", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

bin_labels <- c("Input", "Unsorted", "Low Sort", "High Sort")

# Figure 1: Heatmap
p1 <- plot_cys_heatmap(cys_data, bin_labels = bin_labels,
                        title = "CcdB Cysteine Substitution Frequency")
ggsave(file.path(fig_dir, "demo_heatmap.png"), p1, width = 14, height = 4, dpi = 200)
cat("  [OK] demo_heatmap.png (14x4 in)\n")

# Figure 2: Enrichment bar plot
p2 <- plot_cys_enrichment(cys_data, bin_labels = bin_labels)
ggsave(file.path(fig_dir, "demo_barplot.png"), p2, width = 14, height = 6, dpi = 200)
cat("  [OK] demo_barplot.png (14x6 in)\n")

# Figure 3: Mutation landscape
p3 <- plot_mutation_landscape(cys_data, bin_to_plot = 4, bin_label = "High Sort")
ggsave(file.path(fig_dir, "demo_landscape.png"), p3, width = 14, height = 5, dpi = 200)
cat("  [OK] demo_landscape.png (14x5 in)\n")

# Figure 4: Epitope profile
p4 <- plot_epitope_profile(cys_norm,
  bin_labels = c("Unsorted", "Low Sort", "High Sort"),
  smooth = TRUE, span = 0.2,
  highlight_regions = list(
    list(start = 20, end = 30, label = "Loop 1"),
    list(start = 55, end = 65, label = "Active Site"),
    list(start = 90, end = 101, label = "C-terminal")
  ))
ggsave(file.path(fig_dir, "demo_epitope.png"), p4, width = 14, height = 6, dpi = 200)
cat("  [OK] demo_epitope.png (14x6 in)\n")

# Figure 5: Surface radar
p5 <- plot_surface_accessibility(cys_norm, enrichment_bin = 4)
ggsave(file.path(fig_dir, "demo_radar.png"), p5, width = 8, height = 8, dpi = 200)
cat("  [OK] demo_radar.png (8x8 in)\n")

# Figure 6: QC dashboard
p6 <- plot_qc_summary(cys_data, bin_labels = bin_labels)
ggsave(file.path(fig_dir, "demo_qc.png"), p6, width = 14, height = 10, dpi = 200)
cat("  [OK] demo_qc.png (14x10 in)\n")

# Figure 7: Donut chart
p7 <- plot_mutation_class_donut(classified$summary)
ggsave(file.path(fig_dir, "demo_donut.png"), p7, width = 7, height = 7, dpi = 200)
cat("  [OK] demo_donut.png (7x7 in)\n")

cat("\n")
cat("=" |> strrep(60), "\n")
cat("  DEMO COMPLETE — All 7 figures saved to:\n")
cat("  ", fig_dir, "\n")
cat("=" |> strrep(60), "\n")
