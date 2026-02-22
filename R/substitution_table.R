#' Build Substitution Table from Combined Reads
#'
#' @description Compares each combined read to the wild-type reference and
#' identifies nucleotide-level substitution positions, along with the WT and
#' mutant codons at each position. Replaces \code{get_substitution_table.pl}.
#'
#' @param combined_df Data frame from \code{\link{combine_paired_reads}} (the
#'   \code{combined} element) with columns: \code{read_id}, \code{ref_start},
#'   \code{wt_aligned}, \code{read_aligned}, \code{alignment_line}.
#' @param bin_id Character or integer identifying the current bin.
#'
#' @return A data frame with columns: \code{read_id}, \code{bin},
#'   \code{read_range}, \code{num_substitutions}, \code{positions},
#'   \code{residue_positions}, \code{wt_codons}, \code{mut_codons},
#'   \code{is_wt}.
#' @export
build_substitution_table <- function(combined_df, bin_id = 1) {
  n <- nrow(combined_df)
  results <- vector("list", n)

  for (i in seq_len(n)) {
    wt_seq <- combined_df$wt_aligned[i]
    read_seq <- combined_df$read_aligned[i]
    aln_line <- combined_df$alignment_line[i]
    read_id <- combined_df$read_id[i]
    ref_start <- combined_df$ref_start[i]

    wt_chars <- strsplit(wt_seq, "")[[1]]
    read_chars <- strsplit(read_seq, "")[[1]]
    aln_chars <- strsplit(aln_line, "")[[1]]

    # Find substitution positions (. in alignment line = mismatch)
    sub_positions <- which(aln_chars == ".")
    # Also check for direct mismatches
    for (j in seq_along(wt_chars)) {
      if (!is.na(wt_chars[j]) && !is.na(read_chars[j]) &&
          wt_chars[j] != read_chars[j] &&
          wt_chars[j] != "-" && read_chars[j] != "-") {
        sub_positions <- c(sub_positions, j)
      }
    }
    sub_positions <- sort(unique(sub_positions))

    num_subs <- length(sub_positions)
    is_wt <- (num_subs == 0)

    # Check for indels (spaces in alignment line)
    has_indel <- any(aln_chars == " ")

    if (num_subs > 0 && !has_indel) {
      # Compute nucleotide positions (1-based in WT)
      nt_positions <- sub_positions + ref_start - 1
      # Compute residue (codon) numbers
      res_positions <- ((nt_positions - 1) %/% 3) + 1

      # Get codons for each substitution position
      wt_codons <- character(num_subs)
      mut_codons <- character(num_subs)
      for (k in seq_len(num_subs)) {
        pos <- nt_positions[k]
        codon_start_nt <- ((pos - 1) %/% 3) * 3 + 1  # 1-based codon start in WT
        # Map back to alignment coordinates
        aln_pos <- codon_start_nt - ref_start + 1
        if (aln_pos >= 1 && aln_pos + 2 <= length(wt_chars)) {
          wt_codons[k] <- paste0(wt_chars[aln_pos:(aln_pos + 2)], collapse = "")
          mut_codons[k] <- paste0(read_chars[aln_pos:(aln_pos + 2)], collapse = "")
        } else {
          wt_codons[k] <- NA_character_
          mut_codons[k] <- NA_character_
        }
      }

      results[[i]] <- data.frame(
        read_id = read_id,
        bin = bin_id,
        read_range = paste0(ref_start, "_", ref_start + length(wt_chars) - 1),
        num_substitutions = num_subs,
        positions = paste(nt_positions, collapse = " "),
        residue_positions = paste(res_positions, collapse = " "),
        wt_codons = paste(wt_codons, collapse = " "),
        mut_codons = paste(mut_codons, collapse = " "),
        is_wt = FALSE,
        stringsAsFactors = FALSE
      )
    } else {
      results[[i]] <- data.frame(
        read_id = read_id,
        bin = bin_id,
        read_range = paste0(ref_start, "_", ref_start + length(wt_chars) - 1),
        num_substitutions = 0L,
        positions = "wt",
        residue_positions = "wt",
        wt_codons = "wt",
        mut_codons = "wt",
        is_wt = is_wt,
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, results)
}
