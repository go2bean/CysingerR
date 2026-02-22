#' Quality Filter Reads at Q20 Threshold
#'
#' @description Filters reads by Phred quality score, trimming at the first
#' base below the quality threshold. Reads shorter than a minimum length
#' after trimming are removed. Replaces \code{filter_scores_Q20.pl}.
#'
#' @param reads_df Data frame with at least \code{sequence} and \code{quality}
#'   columns.
#' @param min_quality Integer. Minimum Phred quality score (default 20).
#' @param min_length Integer. Minimum read length after trimming (default 75).
#'
#' @return A list with:
#'   \describe{
#'     \item{reads}{Data frame of quality-filtered reads with trimmed sequences.}
#'     \item{stats}{Data frame with QC statistics: total reads, reads with low
#'       quality bases, trimmed but used reads, short reads discarded, very
#'       good reads (no low quality bases).}
#'   }
#' @export
quality_filter_reads <- function(reads_df,
                                  min_quality = 20L,
                                  min_length = 75L) {
  # ASCII offset for Phred+33 encoding
  offset <- 33L

  # Characters with Phred < min_quality
  low_qual_chars <- rawToChar(as.raw(seq(offset, offset + min_quality - 1)),
                              multiple = TRUE)
  low_qual_pattern <- paste0("[", paste0("\\", low_qual_chars, collapse = ""),
                             "]")

  total <- nrow(reads_df)
  has_low_qual <- 0L
  used_trimmed <- 0L
  short_trimmed <- 0L
  very_good <- 0L
  used_good <- 0L
  short_good <- 0L
  keep <- logical(total)
  trimmed_seq <- character(total)
  trimmed_qual <- character(total)


  for (i in seq_len(total)) {
    qual <- reads_df$quality[i]
    seq_i <- reads_df$sequence[i]

    # Find first low-quality position
    m <- regexpr(low_qual_pattern, qual)
    first_low <- as.integer(m)

    if (first_low > 0) {
      # Has low quality bases
      has_low_qual <- has_low_qual + 1L
      trim_len <- first_low - 1L

      if (trim_len >= min_length) {
        used_trimmed <- used_trimmed + 1L
        trimmed_seq[i] <- substr(seq_i, 1, trim_len)
        trimmed_qual[i] <- substr(qual, 1, trim_len)
        keep[i] <- TRUE
      } else {
        short_trimmed <- short_trimmed + 1L
      }
    } else {
      # All good quality
      very_good <- very_good + 1L
      if (nchar(qual) >= min_length) {
        used_good <- used_good + 1L
        trimmed_seq[i] <- seq_i
        trimmed_qual[i] <- qual
        keep[i] <- TRUE
      } else {
        short_good <- short_good + 1L
      }
    }
  }

  filtered <- reads_df[keep, ]
  filtered$sequence <- trimmed_seq[keep]
  filtered$quality <- trimmed_qual[keep]
  rownames(filtered) <- NULL

  stats <- data.frame(
    metric = c("Total reads", "Reads with low quality bases",
               "Used reads (trimmed)", "Short reads (trimmed, discarded)",
               "Very good reads (no low quality)", "Used very good reads",
               "Short good reads (discarded)", "Total passed"),
    count = c(total, has_low_qual, used_trimmed, short_trimmed,
              very_good, used_good, short_good,
              used_trimmed + used_good),
    stringsAsFactors = FALSE
  )

  list(reads = filtered, stats = stats)
}
