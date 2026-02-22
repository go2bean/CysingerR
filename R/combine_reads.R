#' Combine Paired Forward and Reverse Reads
#'
#' @description Merges forward and reverse reads by read ID, handling
#' overlapping, non-overlapping, and subset cases. Replaces
#' \code{combine_all.pl}.
#'
#' @param fwd_alignments Data frame of forward read alignments (from
#'   \code{\link{align_to_reference}}).
#' @param rev_alignments Data frame of reverse read alignments (from
#'   \code{\link{align_to_reference}}).
#' @param wt_length Integer. Length of the wild-type nucleotide sequence.
#' @param wt_seq Character string. Wild-type nucleotide sequence.
#'
#' @return A list with:
#'   \describe{
#'     \item{combined}{Data frame of combined reads with columns: \code{read_id},
#'       \code{ref_start}, \code{ref_end}, \code{wt_aligned}, \code{read_aligned},
#'       \code{alignment_line}, \code{category} (overlap/no_overlap/fwd_only/rev_only).}
#'     \item{stats}{Data frame with summary counts.}
#'   }
#' @export
combine_paired_reads <- function(fwd_alignments, rev_alignments,
                                  wt_length, wt_seq) {
  # Remove reads with indels
  fwd_clean <- fwd_alignments[!fwd_alignments$has_indel, ]
  rev_clean <- rev_alignments[!rev_alignments$has_indel, ]

  fwd_indel <- sum(fwd_alignments$has_indel)
  rev_indel <- sum(rev_alignments$has_indel)

  fwd_ids <- fwd_clean$read_id
  rev_ids <- rev_clean$read_id
  both_ids <- intersect(fwd_ids, rev_ids)
  fwd_only_ids <- setdiff(fwd_ids, rev_ids)
  rev_only_ids <- setdiff(rev_ids, fwd_ids)

  results <- list()
  n_overlap <- 0L
  n_same_overlap <- 0L
  n_diff_overlap <- 0L
  n_no_overlap <- 0L
  n_subset <- 0L
  n_short <- 0L

  # Paired reads
  for (rid in both_ids) {
    f <- fwd_clean[fwd_clean$read_id == rid, , drop = FALSE][1, ]
    r <- rev_clean[rev_clean$read_id == rid, , drop = FALSE][1, ]

    f_start <- f$ref_start
    f_end <- f$ref_end
    r_start <- r$ref_start
    r_end <- r$ref_end

    # Both must span full range
    if (f_start == 1 && r_end == wt_length) {
      if (r_start <= f_end && r_end >= f_end) {
        # Overlap
        n_overlap <- n_overlap + 1L
        olap_len <- f_end - r_start + 1

        fwd_olap <- substr(f$read_aligned, r_start, f_end)
        rev_olap <- substr(r$read_aligned, 1, olap_len)

        if (identical(fwd_olap, rev_olap)) {
          n_same_overlap <- n_same_overlap + 1L
          # Combine: fwd + non-overlapping part of rev
          rev_new <- substr(r$read_aligned, olap_len + 1, nchar(r$read_aligned))
          wt_new <- substr(r$ref_aligned, olap_len + 1, nchar(r$ref_aligned))
          aln_new <- substr(r$alignment_line, olap_len + 1, nchar(r$alignment_line))

          results[[length(results) + 1]] <- data.frame(
            read_id = rid,
            ref_start = 1, ref_end = wt_length,
            wt_aligned = paste0(f$ref_aligned, wt_new),
            read_aligned = paste0(f$read_aligned, rev_new),
            alignment_line = paste0(f$alignment_line, aln_new),
            category = "overlap",
            stringsAsFactors = FALSE
          )
        } else {
          n_diff_overlap <- n_diff_overlap + 1L
        }
      } else if (r_start > f_end) {
        # No overlap — fill gap with WT
        n_no_overlap <- n_no_overlap + 1L
        gap <- substr(wt_seq, f_end + 1, r_start - 1)
        gap_aln <- paste0(rep("|", nchar(gap)), collapse = "")

        results[[length(results) + 1]] <- data.frame(
          read_id = rid,
          ref_start = 1, ref_end = wt_length,
          wt_aligned = paste0(f$ref_aligned, gap, r$ref_aligned),
          read_aligned = paste0(f$read_aligned, gap, r$read_aligned),
          alignment_line = paste0(f$alignment_line, gap_aln, r$alignment_line),
          category = "no_overlap",
          stringsAsFactors = FALSE
        )
      } else {
        n_subset <- n_subset + 1L
      }
    } else {
      # One or both reads don't span full range
      if (f_start == 1) {
        # Use forward read
        results[[length(results) + 1]] <- data.frame(
          read_id = rid,
          ref_start = f_start, ref_end = f_end,
          wt_aligned = f$ref_aligned,
          read_aligned = f$read_aligned,
          alignment_line = f$alignment_line,
          category = "fwd_partial",
          stringsAsFactors = FALSE
        )
      } else {
        n_short <- n_short + 1L
      }
    }
  }

  # Forward-only reads
  n_fwd_only_complete <- 0L
  n_fwd_only_partial <- 0L
  for (rid in fwd_only_ids) {
    f <- fwd_clean[fwd_clean$read_id == rid, , drop = FALSE][1, ]
    if (f$ref_start == 1) {
      if (f$ref_end == wt_length) n_fwd_only_complete <- n_fwd_only_complete + 1L
      else n_fwd_only_partial <- n_fwd_only_partial + 1L

      results[[length(results) + 1]] <- data.frame(
        read_id = rid,
        ref_start = f$ref_start, ref_end = f$ref_end,
        wt_aligned = f$ref_aligned,
        read_aligned = f$read_aligned,
        alignment_line = f$alignment_line,
        category = "fwd_only",
        stringsAsFactors = FALSE
      )
    } else {
      n_short <- n_short + 1L
    }
  }

  # Reverse-only reads
  n_rev_only_complete <- 0L
  n_rev_only_partial <- 0L
  for (rid in rev_only_ids) {
    r <- rev_clean[rev_clean$read_id == rid, , drop = FALSE][1, ]
    if (r$ref_end == wt_length) {
      if (r$ref_start == 1) n_rev_only_complete <- n_rev_only_complete + 1L
      else n_rev_only_partial <- n_rev_only_partial + 1L

      results[[length(results) + 1]] <- data.frame(
        read_id = rid,
        ref_start = r$ref_start, ref_end = r$ref_end,
        wt_aligned = r$ref_aligned,
        read_aligned = r$read_aligned,
        alignment_line = r$alignment_line,
        category = "rev_only",
        stringsAsFactors = FALSE
      )
    } else {
      n_short <- n_short + 1L
    }
  }

  combined <- if (length(results) > 0) do.call(rbind, results) else
    data.frame(read_id = character(), ref_start = integer(),
               ref_end = integer(), wt_aligned = character(),
               read_aligned = character(), alignment_line = character(),
               category = character(), stringsAsFactors = FALSE)
  rownames(combined) <- NULL

  stats <- data.frame(
    metric = c("Total paired", "With overlap", "Same overlap",
               "Different overlap", "No overlap", "Subset",
               "Short reads", "Fwd-only complete", "Fwd-only partial",
               "Rev-only complete", "Rev-only partial",
               "Fwd indels removed", "Rev indels removed"),
    count = c(length(both_ids), n_overlap, n_same_overlap,
              n_diff_overlap, n_no_overlap, n_subset,
              n_short, n_fwd_only_complete, n_fwd_only_partial,
              n_rev_only_complete, n_rev_only_partial,
              fwd_indel, rev_indel),
    stringsAsFactors = FALSE
  )

  list(combined = combined, stats = stats)
}
