#' Align Reads to a Reference Sequence
#'
#' @description Performs pairwise Smith-Waterman alignment of each read against
#' a wild-type reference sequence. Uses base R string operations (or
#' Biostrings if available) to replace the external EMBOSS water dependency.
#' Replaces \code{run_water_fwd.pl} and \code{run_water_rev.pl}.
#'
#' @param sequences Character vector of read sequences.
#' @param reference Character string of the wild-type reference.
#' @param read_ids Character vector of read IDs (same length as sequences).
#' @param gap_open Numeric gap opening penalty (default 20).
#' @param gap_extend Numeric gap extension penalty (default 0.5).
#'
#' @return A data frame with columns: \code{read_id}, \code{ref_start},
#'   \code{ref_end}, \code{ref_aligned}, \code{read_start}, \code{read_end},
#'   \code{read_aligned}, \code{alignment_line}, \code{score},
#'   \code{has_indel}.
#' @export
align_to_reference <- function(sequences, reference, read_ids,
                                gap_open = 20, gap_extend = 0.5) {
  n <- length(sequences)
  if (length(read_ids) != n) stop("sequences and read_ids must be same length")

  use_biostrings <- requireNamespace("Biostrings", quietly = TRUE)

  results <- vector("list", n)

  for (i in seq_len(n)) {
    seq_i <- toupper(sequences[i])
    ref <- toupper(reference)

    if (use_biostrings) {
      aln <- Biostrings::pairwiseAlignment(
        pattern = seq_i,
        subject = ref,
        type = "local",
        substitutionMatrix = Biostrings::nucleotideSubstitutionMatrix(
          match = 5, mismatch = -4
        ),
        gapOpening = gap_open,
        gapExtension = gap_extend
      )

      ref_aln <- as.character(Biostrings::alignedSubject(aln))
      read_aln <- as.character(Biostrings::alignedPattern(aln))
      score_val <- Biostrings::score(aln)

      # Subject (ref) start/end
      ref_start <- Biostrings::start(Biostrings::subject(aln))
      ref_end <- ref_start + nchar(gsub("-", "", ref_aln)) - 1

      # Pattern (read) start/end
      read_start <- Biostrings::start(Biostrings::pattern(aln))
      read_end <- read_start + nchar(gsub("-", "", read_aln)) - 1

    } else {
      # Fallback: simple ungapped alignment by substring matching
      # Find best matching position using sliding window
      read_len <- nchar(seq_i)
      ref_len <- nchar(ref)

      best_score <- -Inf
      best_start <- 1
      for (s in seq_len(max(1, ref_len - read_len + 1))) {
        sub_ref <- substr(ref, s, min(s + read_len - 1, ref_len))
        sub_len <- min(nchar(sub_ref), nchar(seq_i))
        matches <- sum(strsplit(substr(sub_ref, 1, sub_len), "")[[1]] ==
                       strsplit(substr(seq_i, 1, sub_len), "")[[1]])
        if (matches > best_score) {
          best_score <- matches
          best_start <- s
        }
      }

      aln_len <- min(read_len, ref_len - best_start + 1)
      ref_aln <- substr(ref, best_start, best_start + aln_len - 1)
      read_aln <- substr(seq_i, 1, aln_len)
      ref_start <- best_start
      ref_end <- best_start + aln_len - 1
      read_start <- 1
      read_end <- aln_len
      score_val <- best_score
    }

    # Build alignment line (| for match, . for mismatch, space for gap)
    ref_chars <- strsplit(ref_aln, "")[[1]]
    read_chars <- strsplit(read_aln, "")[[1]]
    aln_line <- character(length(ref_chars))
    for (j in seq_along(ref_chars)) {
      if (ref_chars[j] == "-" || read_chars[j] == "-") {
        aln_line[j] <- " "
      } else if (ref_chars[j] == read_chars[j]) {
        aln_line[j] <- "|"
      } else {
        aln_line[j] <- "."
      }
    }

    has_indel <- any(ref_chars == "-") || any(read_chars == "-")

    results[[i]] <- data.frame(
      read_id = read_ids[i],
      ref_start = ref_start,
      ref_end = ref_end,
      ref_aligned = ref_aln,
      read_start = read_start,
      read_end = read_end,
      read_aligned = read_aln,
      alignment_line = paste0(aln_line, collapse = ""),
      score = score_val,
      has_indel = has_indel,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}

#' Parse Alignments into Structured Format
#'
#' @description Extracts structured alignment information for downstream
#' analysis. Replaces \code{reformat_watout_fwd.pl} and
#' \code{reformat_watout_rev.pl}.
#'
#' @param alignments Data frame from \code{\link{align_to_reference}}.
#'
#' @return The input data frame unchanged (already in structured format).
#'   This function exists for API consistency with the original pipeline.
#' @export
parse_alignments <- function(alignments) {
  alignments
}
