#' Demultiplex Reads by MID and Primer Combinations
#'
#' @description Assigns filtered reads to bins based on their MID × primer
#' combination. Replaces \code{bin_all_fwd.pl} and \code{bin_all_rev.pl}.
#'
#' @param reads_df Data frame from \code{\link{filter_fastq_reads}} (must have
#'   columns \code{mid_tag} and \code{primer}).
#' @param mid_tags Character vector of MID tag sequences (order defines bin
#'   numbering).
#' @param primers Character vector of primer sequences (order defines sub-bin
#'   lettering: a, b, c, ...).
#' @param mid_length Integer. Length of MID tag (default 6).
#' @param primer_length Integer. Length of primer (default 18).
#'
#' @return The input data frame with an additional \code{bin_id} column
#'   (e.g., "1a", "1b", "2a", "2b", ...) and \code{bin_numeric} (integer
#'   bin number).
#' @export
demultiplex_reads <- function(reads_df,
                               mid_tags,
                               primers,
                               mid_length = 6L,
                               primer_length = 18L) {
  # Only work with classified reads
  classified <- reads_df$classified
  reads_df$bin_id <- NA_character_
  reads_df$bin_numeric <- NA_integer_

  if (sum(classified) == 0) return(reads_df)

  # Build lookup table: MID+primer -> bin_id
  alpha <- letters[seq_along(primers)]
  bin_lookup <- list()
  for (m in seq_along(mid_tags)) {
    for (p in seq_along(primers)) {
      key <- paste0(mid_tags[m], primers[p])
      bin_lookup[[key]] <- list(
        bin_id = paste0(m, alpha[p]),
        bin_numeric = m
      )
    }
  }

  # Assign bins
  for (i in which(classified)) {
    seq_i <- reads_df$sequence[i]
    mid_i <- substr(seq_i, 1, mid_length)
    primer_i <- substr(seq_i, mid_length + 1, mid_length + primer_length)
    key_i <- paste0(mid_i, primer_i)

    if (!is.null(bin_lookup[[key_i]])) {
      reads_df$bin_id[i] <- bin_lookup[[key_i]]$bin_id
      reads_df$bin_numeric[i] <- bin_lookup[[key_i]]$bin_numeric
    }
  }

  reads_df
}


#' Combine Forward and Reverse Bin Files
#'
#' @description Merges forward and reverse read data frames for each bin.
#' Replaces \code{import_bins.pl}.
#'
#' @param fwd_reads Data frame of demultiplexed forward reads.
#' @param rev_reads Data frame of demultiplexed reverse reads.
#'
#' @return Combined data frame with all reads, sorted by bin.
#' @export
combine_bin_directions <- function(fwd_reads, rev_reads) {
  combined <- rbind(fwd_reads, rev_reads)
  combined <- combined[order(combined$bin_id), ]
  rownames(combined) <- NULL
  combined
}
