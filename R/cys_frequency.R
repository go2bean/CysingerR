#' Calculate Cysteine Mutation Frequency per Position
#'
#' @description From single-mutant reads, counts cysteine (TGT/TGC) mutations
#' at each residue position per bin at the codon level. Replaces
#' \code{get_CYS_freq.pl}.
#'
#' @param single_muts Data frame of single-mutant reads from
#'   \code{\link{classify_mutations}}.
#' @param wt_nt Character string. Wild-type nucleotide sequence.
#' @param n_bins Integer. Number of bins (MIDs).
#' @param cys_codons Character vector of cysteine codons (default
#'   \code{c("TGT", "TGC")}).
#'
#' @return A data frame with columns: \code{position}, \code{wt_residue},
#'   \code{cys_codon}, \code{mut_residue} ("C"), and one column per bin
#'   (\code{bin_1}, \code{bin_2}, ...) with read counts.
#' @export
calculate_cys_frequency <- function(single_muts, wt_nt, n_bins = 4L,
                                     cys_codons = c("TGT", "TGC")) {
  ct <- get_codon_table()
  protein_length <- nchar(wt_nt) %/% 3

  # Translate WT to get residue identities
  wt_protein <- translate_sequence(wt_nt)
  wt_aa <- strsplit(wt_protein, "")[[1]]

  # Parse single mutant data
  freq <- list()
  for (i in seq_len(nrow(single_muts))) {
    bin <- single_muts$bin[i]
    res_pos_str <- single_muts$residue_positions[i]
    mut_codon_str <- single_muts$mut_codons[i]

    if (res_pos_str == "wt" || is.na(res_pos_str)) next

    res_pos <- as.integer(strsplit(res_pos_str, "\\s+")[[1]])
    mut_codons <- strsplit(mut_codon_str, "\\s+")[[1]]

    # For single mutants, should be only one unique residue position
    if (length(res_pos) >= 1 && length(mut_codons) >= 1) {
      # Use only complete codons (3 chars)
      valid <- nchar(mut_codons) == 3
      if (any(valid)) {
        # Take the first valid codon for the residue position
        for (idx in which(valid)) {
          codon <- mut_codons[idx]
          key <- paste(res_pos[idx], bin, codon, sep = "_")
          freq[[key]] <- (freq[[key]] %||% 0L) + 1L
        }
      }
    }
  }

  # Build output table
  rows <- list()
  for (pos in seq_len(protein_length)) {
    wt_res <- if (pos <= length(wt_aa)) wt_aa[pos] else "X"
    for (cys_cod in cys_codons) {
      row_data <- list(
        position = pos,
        wt_residue = wt_res,
        cys_codon = cys_cod,
        mut_residue = "C"
      )
      for (b in seq_len(n_bins)) {
        key <- paste(pos, b, cys_cod, sep = "_")
        row_data[[paste0("bin_", b)]] <- freq[[key]] %||% 0L
      }
      rows[[length(rows) + 1]] <- as.data.frame(row_data, stringsAsFactors = FALSE)
    }
  }

  do.call(rbind, rows)
}

# Helper for NULL coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Aggregate Cysteine Frequency at Residue Level
#'
#' @description Sums cysteine mutation counts across synonymous codons
#' (TGT + TGC) to get residue-level frequency per position and bin. Replaces
#' \code{get_CYS_freq_reslevel.pl}.
#'
#' @param cys_freq Data frame from \code{\link{calculate_cys_frequency}}.
#' @param n_bins Integer. Number of bins.
#'
#' @return A data frame with columns: \code{position}, \code{wt_residue},
#'   \code{mut_residue} ("C"), and bin columns with summed counts.
#' @export
aggregate_cys_residue_level <- function(cys_freq, n_bins = 4L) {
  bin_cols <- paste0("bin_", seq_len(n_bins))

  # Group by position and sum across cysteine codons
  positions <- unique(cys_freq$position)
  results <- vector("list", length(positions))

  for (idx in seq_along(positions)) {
    pos <- positions[idx]
    subset <- cys_freq[cys_freq$position == pos, ]
    row_data <- list(
      position = pos,
      wt_residue = subset$wt_residue[1],
      mut_residue = "C"
    )
    for (bc in bin_cols) {
      row_data[[bc]] <- sum(subset[[bc]])
    }
    results[[idx]] <- as.data.frame(row_data, stringsAsFactors = FALSE)
  }

  do.call(rbind, results)
}


#' Normalize Cysteine Frequency Data
#'
#' @description Normalizes raw read counts to enrichment ratios. The
#' normalization follows Najar et al. (2017): divide each bin count by the
#' input library count, then normalize by the sum across positions.
#'
#' @param cys_res_level Data frame from \code{\link{aggregate_cys_residue_level}}.
#' @param input_bin Integer or character. Which bin column represents the input
#'   library (default 1, meaning \code{bin_1}).
#' @param sorted_bins Integer vector. Which bin columns represent sorted bins
#'   (default \code{2:4}).
#' @param pseudocount Numeric. Added to all counts before normalization to
#'   avoid division by zero (default 1).
#'
#' @return Data frame with additional columns: \code{enrichment_bin_X} for each
#'   sorted bin, and \code{log2_enrichment_bin_X}.
#' @export
normalize_cys_frequency <- function(cys_res_level,
                                     input_bin = 1L,
                                     sorted_bins = 2:4,
                                     pseudocount = 1) {
  input_col <- paste0("bin_", input_bin)

  # Add pseudocount
  input_counts <- cys_res_level[[input_col]] + pseudocount
  input_total <- sum(input_counts)
  input_freq <- input_counts / input_total

  for (sb in sorted_bins) {
    sort_col <- paste0("bin_", sb)
    sort_counts <- cys_res_level[[sort_col]] + pseudocount
    sort_total <- sum(sort_counts)
    sort_freq <- sort_counts / sort_total

    enrich_col <- paste0("enrichment_bin_", sb)
    log2_col <- paste0("log2_enrichment_bin_", sb)

    cys_res_level[[enrich_col]] <- sort_freq / input_freq
    cys_res_level[[log2_col]] <- log2(sort_freq / input_freq)
  }

  cys_res_level
}
