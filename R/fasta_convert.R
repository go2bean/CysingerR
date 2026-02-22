#' Convert Reads to FASTA Format
#'
#' @description Converts filtered reads into FASTA format strings. For
#' reverse reads, also computes the reverse complement. Replaces
#' \code{Convert2Fasta.pl} and \code{revComp_Convert2Fasta.pl}.
#'
#' @param reads_df Data frame with \code{read_id}, \code{sequence}, and
#'   \code{direction} columns.
#' @param mid_length Integer. Length of MID tag to strip from start of
#'   sequence (default 6).
#'
#' @return Data frame with columns: \code{read_id}, \code{fasta_header},
#'   \code{fasta_sequence} (reverse-complemented if direction is "reverse").
#' @export
reads_to_fasta <- function(reads_df, mid_length = 6L) {
  out <- data.frame(
    read_id = reads_df$read_id,
    fasta_header = paste0(">", reads_df$read_id),
    fasta_sequence = substr(reads_df$sequence, mid_length + 1,
                            nchar(reads_df$sequence)),
    direction = reads_df$direction,
    stringsAsFactors = FALSE
  )

  # Reverse complement for reverse reads
  rev_idx <- which(out$direction == "reverse")
  if (length(rev_idx) > 0) {
    out$fasta_sequence[rev_idx] <- vapply(
      out$fasta_sequence[rev_idx],
      reverse_complement_seq,
      character(1),
      USE.NAMES = FALSE
    )
  }

  out
}

#' Reverse Complement a DNA Sequence
#'
#' @description Computes the reverse complement of a DNA sequence string.
#'
#' @param seq Character string of DNA.
#' @return Character string of the reverse complement.
#' @export
#' @examples
#' reverse_complement_seq("ATGC")  # "GCAT"
reverse_complement_seq <- function(seq) {
  seq <- toupper(seq)
  chars <- strsplit(seq, "")[[1]]
  comp <- chartr("ACGT", "TGCA", chars)
  paste0(rev(comp), collapse = "")
}

#' Write FASTA Data Frame to File
#'
#' @description Writes a FASTA data frame to a file.
#'
#' @param fasta_df Data frame with \code{fasta_header} and
#'   \code{fasta_sequence} columns.
#' @param file Path to the output file.
#'
#' @return Invisible NULL. Called for side effect.
write_fasta <- function(fasta_df, file) {
  lines <- paste0(fasta_df$fasta_header, "\n", fasta_df$fasta_sequence)
  writeLines(lines, file)
  invisible(NULL)
}
