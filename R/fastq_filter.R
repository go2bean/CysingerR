#' Filter FASTQ Reads by MID Tags and Primers
#'
#' @description Parses paired-end FASTQ files, identifies reads matching known
#' MID (Multiplex Identifier) tags and primer sequences, strips padding/spacer
#' nucleotides, and returns classified reads. Replaces the Perl scripts
#' \code{initial_filter_fwd.pl} and \code{initial_filter_rev.pl}.
#'
#' @param fastq_file Path to a FASTQ file.
#' @param mid_tags Character vector of MID tag sequences (6-mer barcodes).
#' @param primers Character vector of primer sequences.
#' @param sequencer_prefix Character string for the sequencer header prefix
#'   (default \code{"@M"}).
#' @param max_padding Integer. Maximum allowed random padding bases before the
#'   MID tag (default 4).
#' @param direction Character. Either \code{"forward"} or \code{"reverse"}.
#'
#' @return A data frame with columns: \code{read_id}, \code{sequence},
#'   \code{quality}, \code{mid_tag}, \code{primer}, \code{direction},
#'   \code{classified} (logical).
#' @export
#' @examples
#' # See vignette for full examples
filter_fastq_reads <- function(fastq_file,
                                mid_tags,
                                primers,
                                sequencer_prefix = "@M",
                                max_padding = 4L,
                                direction = c("forward", "reverse")) {
  direction <- match.arg(direction)

  lines <- readLines(fastq_file, warn = FALSE)
  n_lines <- length(lines)
  if (n_lines %% 4 != 0) {
    warning("FASTQ file does not have a number of lines divisible by 4. ",
            "Truncating to complete records.")
    n_lines <- (n_lines %/% 4) * 4
    lines <- lines[seq_len(n_lines)]
  }

  n_reads <- n_lines / 4
  # Pre-build regex pattern
  mid_pattern <- paste0(mid_tags, collapse = "|")
  pri_pattern <- paste0(primers, collapse = "|")
  full_pattern <- paste0("^([ATGC]{0,", max_padding, "})(", mid_pattern,
                         ")(G{0,3})(", pri_pattern, ")(.*)")

  results <- vector("list", n_reads)

  for (i in seq_len(n_reads)) {
    idx <- (i - 1) * 4
    header <- lines[idx + 1]
    seq_line <- lines[idx + 2]
    # plus_line <- lines[idx + 3]  # "+"
    qual_line <- lines[idx + 4]

    # Extract read ID from header
    read_id <- sub("^@", "", strsplit(header, " ")[[1]][1])

    # Try to match MID + primer pattern
    m <- regexec(full_pattern, seq_line)
    match_data <- regmatches(seq_line, m)

    if (length(match_data[[1]]) > 0) {
      # Matched: groups are: full, padding, mid, spacer_G, primer, insert
      padding <- match_data[[1]][2]
      mid <- match_data[[1]][3]
      spacer <- match_data[[1]][4]
      primer <- match_data[[1]][5]
      insert <- match_data[[1]][6]

      # Trimmed sequence = MID + primer + insert
      trimmed_seq <- paste0(mid, primer, insert)

      # Trim quality string correspondingly
      pad_len <- nchar(padding)
      spacer_len <- nchar(spacer)
      # Remove padding from start, spacer from quality
      qual_trimmed <- qual_line
      if (spacer_len > 0) {
        spacer_start <- nchar(padding) + nchar(mid) + 1  # 1-indexed start of spacer in original
        qual_trimmed <- paste0(
          substr(qual_line, pad_len + 1, spacer_start - 1),
          substr(qual_line, spacer_start + spacer_len, nchar(qual_line))
        )
      } else {
        qual_trimmed <- substr(qual_line, pad_len + 1, nchar(qual_line))
      }

      results[[i]] <- list(
        read_id = read_id,
        sequence = trimmed_seq,
        quality = qual_trimmed,
        mid_tag = mid,
        primer = primer,
        direction = direction,
        classified = TRUE
      )
    } else {
      # Unclassified read
      results[[i]] <- list(
        read_id = read_id,
        sequence = gsub("\n", "", seq_line),
        quality = qual_line,
        mid_tag = NA_character_,
        primer = NA_character_,
        direction = direction,
        classified = FALSE
      )
    }
  }

  do.call(rbind.data.frame, c(results, stringsAsFactors = FALSE))
}
